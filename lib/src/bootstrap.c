/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License along
 *   with this program; if not, write to the Free Software Foundation, Inc.,
 *   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "libgretl.h"
#include "bootstrap.h"

#define BDEBUG 0

enum {
    BOOT_CI          = 1 << 0,  /* compute confidence interval */
    BOOT_PVAL        = 1 << 1,  /* compute p-value */
    BOOT_RESAMPLE_U  = 1 << 2,  /* resample the empirical residuals */
    BOOT_NORMAL_U    = 1 << 3,  /* simulate normal residuals */
    BOOT_RESCALE     = 1 << 4,  /* rescale the residuals, if resampling */
    BOOT_GRAPH       = 1 << 5,  /* graph the distribution */
    BOOT_LDV         = 1 << 6,  /* model includes lagged dep var */
    BOOT_VERBOSE     = 1 << 7   /* for debugging */
};

#define resampling(b) (b->flags & BOOT_RESAMPLE_U)
#define verbose(b) (b->flags & BOOT_VERBOSE)

typedef struct boot_ boot;

struct boot_ {
    int flags;          /* option flags */
    int B;              /* number of replications */
    int k;              /* number of coefficients */
    int T;              /* number of observations used */
    int p;              /* index number of coeff to examine */
    int ldvpos;         /* col. number of lagged dep var in X matrix */
    int ldvpos0;        /* col. number of lagged dep var in full X matrix */
    gretl_matrix *y;    /* holds original, then artificial, dep. var. */
    gretl_matrix *X;    /* independent variables */
    gretl_matrix *b0;   /* coefficients used to generate dep var */
    gretl_matrix *u0;   /* original residuals for resampling */
    gretl_matrix *Xr;   /* restricted set of indep. vars (p-value case) */
    double SE;          /* original std. error of residuals */
    double point;       /* point estimate of coeff */
    double tp;          /* original t-stat for variable of interest */
    double a;           /* alpha, for confidence interval */
    char vname[VNAMELEN]; /* name of variable analysed */
};

static void boot_destroy (boot *bs)
{
    gretl_matrix_free(bs->y);
    gretl_matrix_free(bs->X);
    gretl_matrix_free(bs->b0);
    gretl_matrix_free(bs->u0);
    gretl_matrix_free(bs->Xr);

    free(bs);
}

static boot *boot_new (gretl_matrix *y,
		       gretl_matrix *X,
		       gretl_matrix *b,
		       gretl_matrix *u,
		       double a,
		       int flags)
{
    boot *bs;

    bs = malloc(sizeof *bs);
    if (bs == NULL) {
	return NULL;
    }

    bs->flags = flags;
    bs->a = a;
    bs->B = 0;
    bs->p = 0;
    bs->ldvpos = -1;
    bs->ldvpos0 = -1;
    *bs->vname = '\0';

    bs->y = y;
    bs->X = X;
    bs->b0 = b;
    bs->u0 = u;
    bs->Xr = NULL;

    bs->SE = NADBL;
    bs->point = NADBL;
    bs->tp = NADBL;

    bs->k = X->cols;
    bs->T = X->rows;

    return bs;
}

static void make_normal_y (boot *bs)
{
    gretl_matrix *X = (bs->flags & BOOT_PVAL)? bs->Xr : bs->X;
    int i, t;
    double xti;

    /* generate scaled normal errors */
    gretl_normal_dist(bs->y->val, 0, bs->T - 1);
    for (t=0; t<bs->T; t++) {
	bs->y->val[t] *= bs->SE;
    }

    /* construct y recursively */
    for (t=0; t<X->rows; t++) {
	for (i=0; i<X->cols; i++) {
	    if (t > 0 && i == bs->ldvpos) {
		gretl_matrix_set(X, t, i, bs->y->val[t-1]);
	    } 
	    xti = gretl_matrix_get(X, t, i);
	    bs->y->val[t] += bs->b0->val[i] * xti;
	}
    }  	
}

static void 
resample_vector (const gretl_matrix *u0, gretl_matrix *u,
		 double *z)
{
    int T = u->rows;
    int i, t;

    /* generate uniform random series */
    gretl_uniform_dist(z, 0, T - 1);

    /* sample from source vector based on indices */
    for (t=0; t<T; t++) {
	i = T * z[t];
	if (i > T - 1) {
	    i = T - 1;
	}
	u->val[t] = u0->val[i];
    }
}

static void 
make_resampled_y (boot *bs, double *z)
{
    gretl_matrix *X = (bs->flags & BOOT_PVAL)? bs->Xr : bs->X;
    double xti;
    int i, t;

    /* resample the residuals, into y */
    resample_vector(bs->u0, bs->y, z);

    /* construct y recursively */
    for (t=0; t<X->rows; t++) {
	for (i=0; i<X->cols; i++) {
	    if (t > 0 && i == bs->ldvpos) {
		gretl_matrix_set(X, t, i, bs->y->val[t-1]);
	    }
	    xti = gretl_matrix_get(X, t, i);
	    bs->y->val[t] += bs->b0->val[i] * xti;
	}
    }
}

/* when doing a bootstrap p-value: run a restricted regression that
   excludes the variable of interest, and save the coefficient vector
   and residuals
*/

static int do_restricted_ols (boot *bs)
{
    double s2, xti;
    int k = bs->k;
    int p = bs->p;
    int i, j, t;
    int err = 0;

    bs->Xr = gretl_matrix_alloc(bs->T, k - 1);

    if (bs->Xr == NULL) {
	return E_ALLOC;
    }

    /* make restricted X matrix, Xr */
    j = 0;
    for (i=0; i<k; i++) {
	if (i != p) {
	    for (t=0; t<bs->T; t++) {
		xti = gretl_matrix_get(bs->X, t, i);
		gretl_matrix_set(bs->Xr, t, j, xti);
	    }
	    j++;
	}
    }

    /* adjust ldv info, if needed */
    if (bs->flags & BOOT_LDV) {
	if (bs->ldvpos == p) {
	    /* excluding lagged dep var */
	    bs->flags &= ~BOOT_LDV;
	    bs->ldvpos = -1;
	} else if (bs->ldvpos > p) {
	    /* move it forward one slot */
	    bs->ldvpos -= 1;
	}
    }

    /* shrink b0 */
    gretl_matrix_reuse(bs->b0, k - 1, 1);

    /* estimate restricted model, coeffs into bs->b0 and residuals
       into bs->u0 */
    err = gretl_matrix_ols(bs->y, bs->Xr, bs->b0, NULL, bs->u0, &s2);

#if BDEBUG
    fprintf(stderr, "Restricted estimates:\n");
    for (i=0; i<k-1; i++) {
	fprintf(stderr, "b[%d] = %g\n", i, bs->b0->val[i]);
    }
    fprintf(stderr, "bs->ldvpos = %d (bs->ldvpos0 = %d)\n", bs->ldvpos,
	    bs->ldvpos0);
#endif

    if (!err) {
	bs->SE = sqrt(s2);
    }

    return err;
}

static void bs_print_result (boot *bs, double *xi, int tail, PRN *prn)
{
    pprintf(prn, _("For the coefficient on %s (point estimate %g)"), 
	    bs->vname, bs->point);
    
    pputs(prn, ":\n\n  ");

    if (bs->flags & (BOOT_CI | BOOT_GRAPH)) {
	qsort(xi, bs->B, sizeof *xi, gretl_compare_doubles);
    }

    if (bs->flags & BOOT_PVAL) {
	double pv = (double) tail / bs->B;

	pprintf(prn, "%s = %d / %d = %g", _("p-value"), tail, bs->B, pv);
    } else {
	double ql, qu;
	int i;

	/* FIXME: do something more sophisticated than simple
	   quantiles here? */

	i = bs->a * (bs->B + 1) / 2.0;
	ql = xi[i-1];

	i = bs->B - i + 1;
	qu = xi[i-1];

	pprintf(prn, "%g%% confidence interval = %g to %g", 
		100 * (1 - bs->a), ql, qu);	
    }
    
    pputs(prn, "\n\n");
    pprintf(prn, "Based on %d replications, ", bs->B);
    if (bs->flags & BOOT_RESAMPLE_U) {
	pputs(prn, "using resampled residuals");
    } else {
	pputs(prn, "with simulated normal errors");
    }

    if (bs->flags & BOOT_LDV) {
	pputc(prn, '\n');
	pputs(prn, "(recognized lagged dependent variable)");
    }

    if (bs->flags & BOOT_GRAPH) {
	int (*kdfunc) (const double *, int, const char *);
	void *handle;
	char label[48];
	int err;

	kdfunc = get_plugin_function("array_kernel_density", &handle);
	if (kdfunc == NULL) {
	    return;
	}

	if (bs->flags & BOOT_CI) {
	    strcpy(label, "bootstrap coefficient");
	} else {
	    strcpy(label, "bootstrap t-ratio");
	}

	err = (*kdfunc)(xi, bs->B, label);
	close_plugin(handle);
    }
}

/* Davidson and MacKinnon, ETM, p. 163 */

static void rescale_residuals (boot *bs)
{
    double s;
    int t, k = bs->k;

    if (bs->flags & BOOT_PVAL) {
	k--;
    }

    s = sqrt((double) bs->T / (bs->T - k));

    for (t=0; t<bs->T; t++) {
	bs->u0->val[t] *= s;
    }
}

/* do the actual bootstrap analysis: the objective is either to form a
   confidence interval or to compute a p-value; the methodology is
   either to resample the original residuals or to simulate normal
   errors with the empirically given variance.
*/

static int do_bootstrap (boot *bs, PRN *prn)
{
    gretl_matrix *XTX = NULL;   /* X'X */
    gretl_matrix *XTXI = NULL;  /* X'X^{-1} */
    gretl_matrix *b = NULL;     /* re-estimated coeffs */
    gretl_matrix *yh = NULL;    /* fitted values */
    double *z = NULL;
    double *xi = NULL;
    int k = bs->k;
    int p = bs->p;
    int tail = 0;
    int i, t, err = 0;

    if (bs->flags & BOOT_PVAL) {
	err = do_restricted_ols(bs);
	if (err) {
	    return err;
	}
    }

    b = gretl_column_vector_alloc(k);
    XTX = gretl_matrix_alloc(k, k);
    yh = gretl_column_vector_alloc(bs->T);

    if (b == NULL || XTX == NULL || yh == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (resampling(bs)) {
	/* resampling index array */
	z = malloc(bs->T * sizeof *z);
	if (z == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	if (bs->flags & BOOT_RESCALE) {
	    rescale_residuals(bs);
	}
    }

    if (bs->flags & (BOOT_CI | BOOT_GRAPH)) {
	/* array for storing results */
	xi = malloc(bs->B * sizeof *xi);
	if (xi == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }	    

    gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
			      bs->X, GRETL_MOD_NONE,
			      XTX, GRETL_MOD_NONE);

    XTXI = gretl_matrix_alloc(XTX->rows, XTX->cols);
    if (XTXI == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_matrix_cholesky_decomp(XTX);
    if (!err) {
	err = gretl_inverse_from_cholesky_decomp(XTXI, XTX);
    }

    if (verbose(bs)) {
	pprintf(prn, "%13s %13s %13s\n", "b", "se", "tval");
    }

    /* carry out B replications */

    for (i=0; i<bs->B && !err; i++) {
	double v, se, SSR, ut, tval;

#if BDEBUG > 1
	fprintf(stderr, "do_bootstrap: round %d\n", i);
#endif

	if (bs->flags & BOOT_NORMAL_U) {
	    make_normal_y(bs);
	} else {
	    make_resampled_y(bs, z); 
	} 

	if (bs->ldvpos0 >= 0) {
	    /* X matrix includes lagged dependent variable, so it has
	       to be modified */
	    for (t=1; t<bs->T; t++) {
		gretl_matrix_set(bs->X, t, bs->ldvpos0, bs->y->val[t-1]);
	    }
	    gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
				      bs->X, GRETL_MOD_NONE,
				      XTX, GRETL_MOD_NONE);
	    err = gretl_matrix_cholesky_decomp(XTX);
	    if (!err) {
		err = gretl_inverse_from_cholesky_decomp(XTXI, XTX);
	    }
	}

	if (!err) {
	    gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
				      bs->y, GRETL_MOD_NONE,
				      b, GRETL_MOD_NONE);
	}

	if (!err) {
	    err = gretl_cholesky_solve(XTX, b);
	}

	if (!err) {
	    /* form fitted values */
	    gretl_matrix_multiply(bs->X, b, yh);

	    /* coeff, standard error and t-stat */
	    SSR = 0.0;
	    for (t=0; t<bs->T; t++) {
		ut = bs->y->val[t] - yh->val[t];
		SSR += ut * ut;
	    } 
	    v = gretl_matrix_get(XTXI, p, p);
	    se = sqrt(v * SSR / (bs->T - k));
	    tval = b->val[p] / se;

	    if (verbose(bs)) {
		pprintf(prn, "%13g %13g %13g\n", b->val[p], se, tval);
	    }

	    if (bs->flags & BOOT_CI) {
		xi[i] = b->val[p];
	    } else {
		/* doing p-value */
		if (bs->flags & BOOT_GRAPH) {
		    xi[i] = tval;
		}
		if (fabs(tval) > fabs(bs->tp)) {
		    tail++;
		}
	    } 
	}
    }

    if (!err) {
	bs_print_result(bs, xi, tail, prn);
    }

 bailout:

    gretl_matrix_free(b);
    gretl_matrix_free(XTX);
    gretl_matrix_free(XTXI);
    gretl_matrix_free(yh);

    free(z);
    free(xi);
    
    return err;
}

static int 
make_model_matrices (const MODEL *pmod, const double **Z,
		     gretl_matrix **py, gretl_matrix **pX,
		     gretl_matrix **pb, gretl_matrix **pu)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *u = NULL;
    double xti;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t;

    y = gretl_column_vector_alloc(T);
    X = gretl_matrix_alloc(T, k);
    b = gretl_column_vector_alloc(k);
    u = gretl_column_vector_alloc(T);

    if (y == NULL || X == NULL || b == NULL || u == NULL) {
	gretl_matrix_free(y);
	gretl_matrix_free(X);
	gretl_matrix_free(b);
	gretl_matrix_free(u);
	return E_ALLOC;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    y->val[s] = Z[pmod->list[1]][t];
	    u->val[s] = pmod->uhat[t];
	    for (i=2; i<=pmod->list[0]; i++) {
		xti = Z[pmod->list[i]][t];
		gretl_matrix_set(X, s, i-2, xti);
	    }
	    s++;
	}
    }

    for (i=0; i<k; i++) {
	b->val[i] = pmod->coeff[i];
    }

    *py = y;
    *pX = X;
    *pb = b;
    *pu = u;

    return 0;
}

static int make_flags (gretlopt opt, int ldv)
{
    int flags = BOOT_RESCALE;

    if (opt & OPT_P) {
	flags |= BOOT_PVAL;
    } else {
	flags |= BOOT_CI;
    }

    if (opt & OPT_N) {
	flags |= BOOT_NORMAL_U;
    } else {
	flags |= BOOT_RESAMPLE_U;
    }

    if (opt & OPT_G) {
	flags |= BOOT_GRAPH;
    }

    if (ldv > 0) {
	flags |= BOOT_LDV;
    }

    return flags;
}

/* alpha * (B + 1) should be an integer, when constructing confidence
   intervals
*/

int maybe_adjust_B (int B, double a)
{
    double x = a * (B + 1);

    while (x - floor(x) > 1e-13) {
	x = a * (++B + 1);
    }

    return B;
}

/**
 * bootstrap_analysis:
 * @pmod: model to be examined.
 * @p: 0-based index number of the coefficient to analyse.
 * @B: number of replications.
 * @Z: data array.
 * @pdinfo: dataset information.
 * @opt: option flags -- may contain %OPT_P to compute p-value
 * (default is to calculate confidence interval), %OPT_N
 * to use simulated normal errors (default is to resample the
 * empirical residuals), %OPT_G to display graph.
 * @prn: printing struct.
 *
 * Calculates a bootstrap confidence interval or p-value for
 * a given coefficient in a given model, estimated via OLS.
 * If the first lag of the dependent variable is present as a
 * regressor it is handled correctly but more complex lag
 * autoregressive schemes are not (yet) handled.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int bootstrap_analysis (MODEL *pmod, int p, int B, const double **Z,
			const DATAINFO *pdinfo, gretlopt opt,
			PRN *prn)
{
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *u = NULL;
    boot *bs = NULL;
    double alpha = .05;
    int ldv, flags = 0;
    int err = 0;

    /* only OLS models for now */
    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (p < 0 || p >= pmod->ncoeff) {
	return E_DATA;
    }

    err = make_model_matrices(pmod, Z, &y, &X, &b, &u);
    if (err) {
	return err;
    }

    ldv = gretl_model_get_int(pmod, "ldepvar");

    flags = make_flags(opt, ldv);
    B = maybe_adjust_B(B, alpha);

    bs = boot_new(y, X, b, u, alpha, flags);
    if (bs == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	int v = pmod->list[p+2];

	bs->p = p;  /* coeff to examine */
	bs->B = B;  /* replications */ 
	bs->SE = pmod->sigma;
	strcpy(bs->vname, pdinfo->varname[v]);
	bs->point = pmod->coeff[p];
	bs->tp = pmod->coeff[p] / pmod->sderr[p];
	if (flags & BOOT_LDV) {
	    bs->ldvpos = bs->ldvpos0 = ldv - 2;
	}
	err = do_bootstrap(bs, prn);
    }

    if (bs != NULL) {
	boot_destroy(bs);
    } else {
	gretl_matrix_free(X);
	gretl_matrix_free(y);
	gretl_matrix_free(b);
	gretl_matrix_free(u);
    }

    return err;
}

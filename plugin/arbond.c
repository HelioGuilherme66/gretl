/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2006 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "libgretl.h"

#define ADEBUG 1

#define XPOS 5

typedef struct arbond_ arbond;

struct unit_info {
    int t1;
    int t2;
};

struct arbond_ {
    int yno;              /* ID number of dependent var */
    int p;                /* lag order for dependent variable */
    int q;                /* max lags of y used as instruments */
    int nx;               /* number of exogenous variables */
    int m;                /* number of columns in instrument matrix */
    int N;                /* number of units */
    int T;                /* total number of observations per unit */
    int maxTi;            /* maximum equations for any given unit */
    int tau;              /* maximal time-series span */
    int k;                /* number of parameters estimated */
    int nobs;             /* total observations used */
    double SSR;           /* sum of squared residuals */
    int *xlist;           /* list of independent variables */
    gretl_matrix *beta;   /* parameter estimates */
    gretl_matrix *vbeta;  /* parameter variance matrix */
    gretl_matrix *uhat;   /* residuals, differenced version */
    gretl_matrix *H;
    gretl_matrix *A;
    gretl_matrix *ZT;     /* transpose of full instrument matrix */
    gretl_matrix *Zi;     /* per-unit instrument matrix */
    gretl_matrix *dy;     /* first difference of dependent var */
    gretl_matrix *dX;     /* lagged differences of y and indep vars */
    gretl_matrix *tmp0;   /* workspace */
    gretl_matrix *tmp1;
    gretl_matrix *tmp2;
    gretl_matrix *L1;
    gretl_matrix *Lk;
    gretl_matrix *R1;
    gretl_matrix *Rk;
    struct unit_info *ui;
};

static int gretl_matrix_drop_zero_rows (gretl_matrix **a)
{
    gretl_matrix *b = NULL;
    char *drop = NULL;
    int ar = (*a)->rows;
    int ac = (*a)->cols;
    int all0, ndrop = 0;
    int i, j, err = 0;

    fprintf(stderr, "drop_zero_rows: on input a = %d x %d\n", ar, ac);

    drop = calloc(ar, 1);
    if (drop == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<ar; i++) {
	all0 = 1;
	for (j=0; j<ac; j++) {
	    if (gretl_matrix_get(*a, i, j) != 0.0) {
		all0 = 0;
		break;
	    }
	}
	if (all0) {
	    drop[i] = 1;
	    ndrop++;
	}
    }

    if (ndrop == ar) {
	/* no non-zero rows in matrix a */
	free(drop);
	return E_DATA;
    }

    if (ndrop > 0) {
	int br = ar - ndrop;

	b = gretl_matrix_alloc(br, ac);
	if (b == NULL) {
	    err = E_ALLOC;
	}
    }

    if (b != NULL) {
	double x;
	int k = 0;

	for (i=0; i<ar; i++) {
	    if (!drop[i]) {
		for (j=0; j<ac; j++) {
		    x = gretl_matrix_get(*a, i, j);
		    gretl_matrix_set(b, k, j, x);
		}
		k++;
	    }
	}
	gretl_matrix_free(*a);
	*a = b;
    }

    free(drop);

    fprintf(stderr, "drop_zero_rows: on output a = %d x %d\n\n", 
	    (*a)->rows, (*a)->cols);

    return err;
} 

static int gretl_matrix_drop_zero_cols (gretl_matrix **a)
{
    gretl_matrix *b = NULL;
    char *drop = NULL;
    int ar = (*a)->rows;
    int ac = (*a)->cols;
    int all0, ndrop = 0;
    int i, j, err = 0;

    fprintf(stderr, "drop_zero_cols: on input a = %d x %d\n", ar, ac);

    drop = calloc(ac, 1);
    if (drop == NULL) {
	return E_ALLOC;
    }

    for (j=0; j<ac; j++) {
	all0 = 1;
	for (i=0; i<ar; i++) {
	    if (gretl_matrix_get(*a, i, j) != 0.0) {
		all0 = 0;
		break;
	    }
	}
	if (all0) {
	    drop[j] = 1;
	    ndrop++;
	}
    }

    if (ndrop == ac) {
	/* no non-zero columns in matrix a */
	free(drop);
	return E_DATA;
    }

    if (ndrop > 0) {
	int bc = ac - ndrop;

	b = gretl_matrix_alloc(ar, bc);
	if (b == NULL) {
	    err = E_ALLOC;
	}
    }

    if (b != NULL) {
	double x;
	int k = 0;

	for (j=0; j<ac; j++) {
	    if (!drop[j]) {
		for (i=0; i<ar; i++) {
		    x = gretl_matrix_get(*a, i, j);
		    gretl_matrix_set(b, i, k, x);
		}
		k++;
	    }
	}
	gretl_matrix_free(*a);
	*a = b;
    }

    free(drop);

    fprintf(stderr, "drop_zero_cols: on output a = %d x %d\n\n", 
	    (*a)->rows, (*a)->cols);

    return err;
} 

static void arbond_free (arbond *ab)
{
    gretl_matrix_free(ab->beta);
    gretl_matrix_free(ab->vbeta);
    gretl_matrix_free(ab->uhat);
    gretl_matrix_free(ab->H);
    gretl_matrix_free(ab->A);
    gretl_matrix_free(ab->Zi);
    gretl_matrix_free(ab->ZT);
    gretl_matrix_free(ab->dy);
    gretl_matrix_free(ab->dX);
    gretl_matrix_free(ab->tmp0);
    gretl_matrix_free(ab->tmp1);
    gretl_matrix_free(ab->tmp2);
    gretl_matrix_free(ab->L1);
    gretl_matrix_free(ab->Lk);
    gretl_matrix_free(ab->R1);
    gretl_matrix_free(ab->Rk);

    free(ab->xlist);
    free(ab->ui);
}

static int arbond_allocate (arbond *ab)
{
    int T2 = ab->maxTi;

    ab->beta = gretl_matrix_alloc(ab->k, 1);
    ab->vbeta = gretl_matrix_alloc(ab->k, ab->k);
    ab->uhat = gretl_matrix_alloc(ab->nobs, 1);
    ab->ZT = gretl_matrix_alloc(ab->m, ab->nobs);
    ab->H = gretl_matrix_alloc(T2, T2);
    ab->A = gretl_matrix_alloc(ab->m, ab->m);
    ab->Zi = gretl_matrix_alloc(T2, ab->m);
    ab->dy = gretl_column_vector_alloc(ab->nobs);
    ab->dX = gretl_matrix_alloc(ab->nobs, ab->k);

    ab->tmp0 = gretl_matrix_alloc(T2, ab->m);
    ab->tmp1 = gretl_matrix_alloc(ab->m, ab->m);
    ab->tmp2 = gretl_matrix_alloc(ab->k, ab->m);
    ab->L1 = gretl_matrix_alloc(1, ab->m);
    ab->Lk = gretl_matrix_alloc(ab->k, ab->m);
    ab->R1 = gretl_matrix_alloc(ab->m, 1);
    ab->Rk = gretl_matrix_alloc(ab->m, ab->k);

    if (ab->ZT == NULL || ab->H == NULL || ab->A == NULL ||
	ab->Zi == NULL || ab->dy == NULL || ab->dX == NULL ||
	ab->tmp0 == NULL || ab->tmp1 == NULL || ab->tmp2 == NULL ||
	ab->L1 == NULL || ab->Lk == NULL ||
	ab->R1 == NULL || ab->Rk == NULL ||
	ab->beta == NULL || ab->vbeta == NULL ||
	ab->uhat == NULL) {
	return E_ALLOC;
    } else {
	return 0;
    }
}

static int 
arbond_init (arbond *ab, const int *list, const DATAINFO *pdinfo)
{
    if (list[0] < 4 || list[3] != LISTSEP) {
	return E_PARSE;
    }

    ab->p = list[1];
    ab->q = list[2];
    ab->yno = list[4];

    /* FIXME handling instruments list */

    if (list[0] > 4) {
	int i;

	ab->nx = list[0] - 4;
	ab->xlist = gretl_list_new(list[0] - 4);
	if (ab->xlist == NULL) {
	    return E_ALLOC;
	}
	for (i=5; i<=list[0]; i++) {
	    ab->xlist[i-4] = list[i];
	}
	printlist(ab->xlist, "ab->xlist");
    } else {
	ab->nx = 0;
	ab->xlist = NULL;
    }

    ab->N = pdinfo->paninfo->nunits;
    ab->T = pdinfo->n / ab->N;
    ab->k = ab->p + ab->nx;
    ab->maxTi = 0;
    ab->tau = 0;
    ab->m = 0;
    ab->nobs = 0;
    ab->SSR = NADBL;

#if ADEBUG
    fprintf(stderr, "yno = %d, p = %d, q = %d, nx = %d, k = %d\n",
	    ab->yno, ab->p, ab->q, ab->nx, ab->k);
#endif

    ab->beta = NULL;
    ab->vbeta = NULL;
    ab->uhat = NULL;
    ab->ZT = NULL;
    ab->H = NULL;
    ab->A = NULL;
    ab->Zi = NULL;
    ab->dy = NULL;
    ab->dX = NULL;
    ab->tmp0 = NULL;
    ab->tmp1 = NULL;
    ab->tmp2 = NULL;
    ab->L1 = NULL;
    ab->Lk = NULL;
    ab->R1 = NULL;
    ab->Rk = NULL;

    ab->ui = malloc(ab->N * sizeof *ab->ui);
    if (ab->ui == NULL) {
	free(ab->xlist);
	ab->xlist = NULL;
	return E_ALLOC;
    } else {
	return 0;
    }
}

static int anymiss (arbond *ab, const double **Z, int s)
{
    int i;

    if (na(Z[ab->yno][s])) {
	return 1;
    }

    if (ab->xlist != NULL) {
	for (i=1; i<=ab->xlist[0]; i++) {
	    if (na(Z[ab->xlist[i]][s])) {
		return 1;
	    }
	}
    }

    return 0;
}

static int 
arbond_sample_check (arbond *ab, const int *list, 
		     const double **Z, const DATAINFO *pdinfo)
{
    int i, s, t;
    int t1min = ab->T - 1;
    int t2max = 0;
    int N = ab->N;
    int err = 0;

    fprintf(stderr, "arbond_sample_check: ab->T = %d\n", ab->T);

    for (i=0; i<ab->N; i++) {
	int t1, t2 = ab->T - 1;
	int t1i = 0, Ti = 0, maxTi = 0;

	/* find the starting point and length of the longest
	   block of consecutive observations on all variables 
	*/
	for (t1=0; t1<=t2; t1++) {
	    s = i * ab->T + t1;
	    Ti = 0;
	    for (t=t1; t<=t2; t++, s++) {
		if (anymiss(ab, Z, s)) {
		    break;
		} else {
		    Ti++;
		} 
	    }
	    if (Ti > maxTi) {
		maxTi = Ti;
		t1i = t1;
	    }
	}

	t1 = t1i;
	t2 = t1 + maxTi - 1;

	/* now allow for lags of y */
	if (t1 < ab->p + 1) {
	    t1 = ab->p + 1;
	}
	s = i * ab->T + t1 - (ab->p + 1);
	for (t=t1; t<=t2; t++) {
	    if (na(Z[ab->yno][s++])) {
		t1++;
	    } else {
		break;
	    }
	}
	Ti = t2 - t1 + 1;

	if (Ti > 0) {
	    if (Ti > ab->maxTi) {
		ab->maxTi = Ti;
	    }
	    ab->nobs += Ti;
	    if (t1 < t1min) {
		t1min = t1;
	    }
	    if (t2 > t2max) {
		t2max = t2;
	    }
	    ab->ui[i].t1 = t1;
	    ab->ui[i].t2 = t2;
	    fprintf(stderr, "Unit %d: t1=%d, t2=%d, Ti = %d\n", 
		    i, t1, t2, Ti);

	} else {
	    N--;
	    ab->ui[i].t1 = -1;
	    ab->ui[i].t2 = -1;
	    fprintf(stderr, "Unit %d: not usable\n", i);
	    
	}
    }

    fprintf(stderr, "Number of units with usable observations = %d\n", N);
    fprintf(stderr, "Total usable observations = %d\n", ab->nobs);
    fprintf(stderr, "Max equations for any given unit = %d\n", ab->maxTi);
    fprintf(stderr, "Maximal overall time-series span: %d to %d = %d\n", 
	    t1min, t2max, t2max - t1min + 1);
    ab->tau = t2max - t1min + 1;

    if (ab->tau > 0) {
	ab->tau += ab->p + 1;
    }

    if (N == 0) {
	err = E_MISSDATA;
    } else {
	ab->m = (ab->tau - 2) * (ab->tau - 1) / 2 + ab->nx; /* FIXME? */
    }

    return err;
}

/* 
   Compute the variance matrix:

   N * C^{-1} * (X'*Z*A_N*\hat{V}_N*A_N*Z'*X) * C^{-1} 

   where C = X'*Z*A_N*Z'*X and
    \hat{V}_N = N^{-1} \sum Z_i'*v_i*v_i'*Z_i,
    (v_i being the step-1 residuals)
*/

static int arbond_variance (arbond *ab, gretl_matrix *den, PRN *prn)
{
    gretl_matrix *tmp = NULL;
    gretl_matrix *num = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *u = NULL;
    double x, s2;
    int T2 = ab->maxTi;
    int i, j, t, k, c;
    int err = 0;

    tmp = gretl_matrix_alloc(ab->k, ab->k);
    num = gretl_matrix_alloc(ab->k, ab->k);
    V = gretl_zero_matrix_new(ab->m, ab->m);
    u = gretl_column_vector_alloc(T2);

    if (tmp == NULL || num == NULL || V == NULL || u == NULL) {
	return E_ALLOC;
    }

    s2 = 0.0;
    c = 0;
    k = 0;

    for (i=0; i<ab->N; i++) {
	int Ti;

	if (ab->ui[i].t1 < 0) {
	    continue;
	}

	Ti = ab->ui[i].t2 - ab->ui[i].t1 + 1;

	gretl_matrix_reuse(ab->Zi, Ti, ab->m);
	gretl_matrix_reuse(u, Ti, 1);

	gretl_matrix_extract_matrix(ab->Zi, ab->ZT, 0, c,
				    GRETL_MOD_TRANSPOSE);

	for (t=0; t<Ti; t++) {
	    u->val[t] = ab->dy->val[k];
	    for (j=0; j<ab->k; j++) {
		x = gretl_matrix_get(ab->dX, k, j);
		u->val[t] -= ab->beta->val[j] * x;
	    }
	    s2 += u->val[t] * u->val[t];
	    ab->uhat->val[k] = u->val[t];
	    k++;
	}

#if ADEBUG > 1
	gretl_matrix_print(u, "ui");
#endif

	gretl_matrix_multiply_mod(u, GRETL_MOD_TRANSPOSE,
				  ab->Zi, GRETL_MOD_NONE,
				  ab->L1);

	gretl_matrix_multiply_mod(ab->L1, GRETL_MOD_TRANSPOSE,
				  ab->L1, GRETL_MOD_NONE,
				  ab->tmp1);

#if ADEBUG > 1
	gretl_matrix_print(ab->tmp1, "Z_i'*v_i*v_i'*Z_i");
#endif

	gretl_matrix_add_to(V, ab->tmp1);
	c += Ti;
    }

    gretl_matrix_divide_by_scalar(V, ab->N);

#if ADEBUG > 1
    gretl_matrix_print(V, "V");
#endif

    /* find the central term, A_N * \hat{V}_N * A_N :
       re-use V for this result */
    err += gretl_matrix_multiply(V, ab->A, ab->tmp1);
    err += gretl_matrix_multiply(ab->A, ab->tmp1, V);

    /* find the enclosing term, Z' * X */
    err += gretl_matrix_multiply(ab->ZT, ab->dX, ab->Rk);

    /* complete the large "middle bit" */
    gretl_matrix_reuse(ab->tmp2, ab->m, ab->k);
    err += gretl_matrix_multiply(V, ab->Rk, ab->tmp2);
    err += gretl_matrix_multiply_mod(ab->Rk, GRETL_MOD_TRANSPOSE,
				     ab->tmp2, GRETL_MOD_NONE,
				     num);

    /* pre- and post-multiply by C^{-1} */
    gretl_invert_symmetric_matrix(den);
    gretl_matrix_multiply(num, den, tmp);
    gretl_matrix_multiply(den, tmp, ab->vbeta);
    gretl_matrix_multiply_by_scalar(ab->vbeta, ab->N);

    gretl_matrix_print_to_prn(ab->vbeta, "Var(beta)", prn);

    /* just in case, restore original dim of tmp2 */
    gretl_matrix_reuse(ab->tmp2, ab->k, ab->m);

#if ADEBUG
    for (i=0; i<ab->k; i++) {
	x = gretl_matrix_get(ab->vbeta, i, i);
	pprintf(prn, "se(beta[%d]) = %g\n", i, sqrt(x));
    }
#endif

    ab->SSR = s2;

#if ADEBUG
    pprintf(prn, "\nRSS = %.11g\n", s2);
    s2 /= ab->nobs - ab->k; /* df correction wanted? (Ox uses it) */
    pprintf(prn, "sigma^2 = %.7g\n", s2);
    pprintf(prn, "sigma = %.7g\n", sqrt(s2));
#endif

    gretl_matrix_free(V);
    gretl_matrix_free(u);
    gretl_matrix_free(tmp);
    gretl_matrix_free(num);

    return err;
}

static void make_first_diff_matrix (gretl_matrix *m)
{
    int n = m->rows;
    double x;
    int i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    x = (i == j)? 2 : (i == j-1 || i == j+1)? -1 : 0;
	    gretl_matrix_set(m, i, j, x);
	}
    }
}

static int arbond_prepare_model (MODEL *pmod, arbond *ab,
				 const int *list,
				 const DATAINFO *pdinfo)
{
    double x;
    int i, j;
    int err = 0;

    pmod->t1 = pdinfo->t1;
    pmod->t2 = pdinfo->t2;
    pmod->dfd = ab->k;
    pmod->dfd = ab->nobs - ab->k;

    pmod->list = gretl_list_copy(list);
    if (pmod->list == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    pmod->ci = ARBOND;
    pmod->ncoeff = ab->k;
    pmod->nobs = ab->nobs;
    pmod->full_n = pdinfo->n;
    pmod->ess = ab->SSR;
    
    x = ab->SSR / pmod->dfd;
    if (x >= 0.0) {
	pmod->sigma = sqrt(x);
    }

    gretl_model_allocate_params(pmod, ab->k);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    j = 0;
    for (i=0; i<ab->p; i++) {
	/* FIXME possible varname overflow */
	sprintf(pmod->params[j++], "D%.10s(-%d)", pdinfo->varname[ab->yno], i+1);
    }

    for (i=0; i<ab->nx; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[ab->xlist[i+1]]);
    }

    err = gretl_model_allocate_storage(pmod);

    if (!err) {
	gretl_model_new_vcv(pmod, NULL);
    }

    if (!err) {
	for (i=0; i<ab->k; i++) {
	    pmod->coeff[i] = ab->beta->val[i];
	    for (j=0; j<=i; j++) {
		x = gretl_matrix_get(ab->vbeta, i, j);
		pmod->vcv[ijton(i, j, ab->k)] = x;
		if (i == j) {
		    pmod->sderr[i] = sqrt(x);
		}
	    }
	}
    }

    if (!err) {
	int t, k = 0;

	for (i=0; i<ab->N; i++) {
	    for (j=0; j<ab->tau; j++) {
		if (j >= ab->ui[i].t1 && j <= ab->ui[i].t2) {
		    t = i * ab->T + j;
		    pmod->uhat[t] = ab->uhat->val[k];
		    pmod->yhat[t] = ab->dy->val[k] - ab->uhat->val[k];
		    k++;
		}
	    }
	}
    }

    return err;
}

static int starting_col (int offset)
{
    int d = 1, c = offset;

    while (offset-- > 0) {
	c += d++;
    }

    return c;
}

/* exclude any independent variables that are zero at all
   relevant observations */

static int arbond_zero_check (arbond *ab, const double **Z)
{
    const double *x;
    int drop = 0;
    int i, j, k, t;

    for (j=0; j<ab->nx; j++) {
	int all0 = 1;

	x = Z[ab->xlist[j+1]];
	for (i=0; i<ab->N && all0; i++) {
	    for (t=ab->ui[i].t1; t<=ab->ui[i].t2 && all0; t++) {
		k = i * ab->T + t;
		if (x[k] != 0.0) {
		    all0 = 0;
		}
	    }
	}
	if (all0) {
	    ab->xlist[j+1] = -1;
	    drop++;
	}
    }

    if (drop > 0) {
	for (i=1; i<=ab->xlist[0]; i++) {
	    if (ab->xlist[i] < 0) {
		gretl_list_delete_at_pos(ab->xlist, i);
		i--;
	    }
	}
	ab->nx = ab->xlist[0];
	ab->k -= drop;
    }

    return 0;
}

static int maybe_shrink_matrices (arbond *ab)
{
    int err;

    err = gretl_matrix_drop_zero_rows(&ab->ZT);
    if (err) {
	return err;
    }

    if (ab->ZT->rows == ab->m) {
	/* nothing to be done */
	return 0;
    }

    err = gretl_matrix_drop_zero_rows(&ab->A);
    if (err) {
	return err;
    }
    
    gretl_matrix_drop_zero_cols(&ab->A);
    if (err) {
	return err;
    }

    if (ab->A->rows < ab->m) {
	ab->m = ab->A->rows;
	gretl_matrix_reuse(ab->tmp0, -1, ab->m);
	gretl_matrix_reuse(ab->tmp1, ab->m, ab->m);
	gretl_matrix_reuse(ab->tmp2, -1, ab->m);
	gretl_matrix_reuse(ab->L1, -1, ab->m);
	gretl_matrix_reuse(ab->Lk, -1, ab->m);
	gretl_matrix_reuse(ab->R1, ab->m, -1);
	gretl_matrix_reuse(ab->Rk, ab->m, -1);
    }

    return 0;
}

/* not really ready yet, just for testing.  list should be
   p q ; y xvars */

MODEL
arbond_estimate (const int *list, const double **X, 
		 const DATAINFO *pdinfo, gretlopt opt,
		 PRN *prn)
{
    MODEL mod;
    arbond ab;
    gretl_matrix *num = NULL;
    gretl_matrix *den = NULL;
    gretl_matrix *tmp = NULL;
    const double *y;
    double x;
    char zstr[8];
    int i, j, k;
    int c, s, t, v;
    int err = 0;

    gretl_model_init(&mod);

    mod.errcode = arbond_init(&ab, list, pdinfo);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_init\n", mod.errcode);
	return mod;
    }

    mod.errcode = arbond_sample_check(&ab, list, X, pdinfo);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_sample_check\n", mod.errcode);
	goto bailout;
    }

    if (ab.nx > 0) {
	arbond_zero_check(&ab, X);
    }

    mod.errcode = arbond_allocate(&ab);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_allocate\n", mod.errcode);
	goto bailout;
    }

    y = X[ab.yno];

    num = gretl_zero_matrix_new(ab.k, 1);
    den = gretl_matrix_alloc(ab.k, ab.k);
    tmp = gretl_matrix_alloc(ab.k, ab.k);

    if (num == NULL || den == NULL || tmp == NULL) {
	mod.errcode = E_ALLOC;
	goto bailout;
    }

    /* build the differenced vector \bar{y} plus the
       X matrix
    */

    s = 0;
    for (i=0; i<ab.N; i++) {
	if (ab.ui[i].t1 < 0) {
	    continue;
	}
	for (t=ab.ui[i].t1; t<=ab.ui[i].t2; t++) {
	    k = i * ab.T + t;
	    /* current difference of dependent var */
	    gretl_vector_set(ab.dy, s, y[k] - y[k-1]);
	    for (j=0; j<ab.p; j++) {
		/* lagged difference of dependent var */
		x = y[k-j-1] - y[k-j-2];
		gretl_matrix_set(ab.dX, s, j, x);
	    }
	    for (j=0; j<ab.nx; j++) {
		/* independent vars */
		v = ab.xlist[j+1];
		x = X[v][k];
		gretl_matrix_set(ab.dX, s, j + ab.p, x);
	    }
	    s++;
	}
    }

#if ADEBUG
    gretl_matrix_print(ab.dy, "dy");
    gretl_matrix_print(ab.dX, "X");
#endif

    /* build instrument matrix blocks, Z_i, and insert into
       big Z' matrix; cumulate first-stage A_N as we go */

    gretl_matrix_zero(ab.A);

    c = 0;
    for (i=0; i<ab.N; i++) {
	int Ti, xc, col = 0, offset = 0;

	if (ab.ui[i].t1 < 0) {
	    continue;
	}

	Ti = ab.ui[i].t2 - ab.ui[i].t1 + 1;
	
	gretl_matrix_reuse(ab.Zi, Ti, ab.m);
	gretl_matrix_zero(ab.Zi);

	/* column offset for initial missing obs */
	offset = ab.ui[i].t1 - (ab.p + 1);
	col = starting_col(offset);

#if ADEBUG
	fprintf(stderr, "Zi[%d]: Ti=%d, offset=%d, col=%d\n",
		i, Ti, offset, col);
#endif

	/* lagged dependent var (level) columns */
	for (t=0; t<Ti; t++) {
	    k = i * ab.T + ab.ui[i].t1 - (ab.p + 1);
	    for (j=col; j<col+t+1; j++) {
		x = y[k++]; 
		gretl_matrix_set(ab.Zi, t, j, x);
	    }
	    col = j + offset;
	}

	/* exogenous var (instr) columns */
	xc = (ab.tau - 2) * (ab.tau - 1) / 2;
	for (j=0; j<ab.nx; j++) {
	    k = i * ab.T + ab.ui[i].t1;
	    for (t=0; t<Ti; t++) {
		x = X[ab.xlist[j+1]][k];
		gretl_matrix_set(ab.Zi, t, xc, x);
		k++;
	    }
	    xc++;
	}

#if ADEBUG
	sprintf(zstr, "Z_%d", i + 1);
	gretl_matrix_print(ab.Zi, zstr);
#endif
	
	gretl_matrix_reuse(ab.H, Ti, Ti);
	make_first_diff_matrix(ab.H);
	gretl_matrix_reuse(ab.tmp0, Ti, ab.m);

	/* Add Z_i' H Z_i to A_N */
	gretl_matrix_multiply(ab.H, ab.Zi, ab.tmp0); 
	gretl_matrix_multiply_mod(ab.Zi, GRETL_MOD_TRANSPOSE,
				  ab.tmp0, GRETL_MOD_NONE,
				  ab.tmp1); 
	gretl_matrix_add_to(ab.A, ab.tmp1);

	/* Write Zi into ZT at offset 0, c */
	gretl_matrix_inscribe_matrix(ab.ZT, ab.Zi, 0, c, GRETL_MOD_TRANSPOSE);
	c += Ti;
    }

    maybe_shrink_matrices(&ab);

    gretl_matrix_divide_by_scalar(ab.A, ab.N);

#if ADEBUG
    gretl_matrix_print(ab.ZT, "ZT");
    gretl_matrix_print(ab.A, "N^{-1} * \\sum Z_i' H Z_i");
#endif

    gretl_invert_symmetric_matrix(ab.A);

#if ADEBUG
    gretl_matrix_print(ab.A, "A_N");
#endif

    /* find "Lk" = \bar{X}' Z A_N */
    gretl_matrix_multiply_mod(ab.dX, GRETL_MOD_TRANSPOSE,
			      ab.ZT, GRETL_MOD_TRANSPOSE,
			      ab.tmp2);
    gretl_matrix_multiply(ab.tmp2, ab.A, ab.Lk);

    /* calculate "numerator", \bar{X}' Z A_N Z' \bar{y} */

    gretl_matrix_multiply(ab.ZT, ab.dy, ab.R1);
    gretl_matrix_multiply(ab.Lk, ab.R1, num);

    /* calculate "denominator", \bar{X}' Z A_N Z' \bar{X} */
    
    gretl_matrix_multiply(ab.ZT, ab.dX, ab.Rk);
    gretl_matrix_multiply(ab.Lk, ab.Rk, den);

    gretl_matrix_copy_values(ab.beta, num);
    gretl_matrix_copy_values(tmp, den);

    err = gretl_LU_solve(tmp, ab.beta);

#if ADEBUG
    pputs(prn, "one-step estimates:\n\n");
    for (i=0; i<ab.k; i++) {
	pprintf(prn, "beta[%d] = %g\n", i, ab.beta->val[i]);
    }
    pputc(prn, '\n');
#endif

    err = arbond_variance(&ab, den, prn);

 bailout:

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	mod.errcode = arbond_prepare_model(&mod, &ab, list, pdinfo);
    }

    arbond_free(&ab);
    gretl_matrix_free(num);
    gretl_matrix_free(den);
    gretl_matrix_free(tmp);

    return mod;
}


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
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#include "libgretl.h" 
#include "var.h"  
#include "vartest.h"
#include "matrix_extra.h"
#include "libset.h"

int gretl_VAR_normality_test (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->E == NULL || var->S == NULL) {
	err = 1;
    } else {
	err = gretl_system_normality_test(var->E, var->S, prn);
    }

    return err;
}

int gretl_VAR_autocorrelation_test (GRETL_VAR *var, int order, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    int i, err = 0;

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
	err = autocorr_test(var->models[i], order, pZ, pdinfo,
			    OPT_Q | OPT_S, prn);
	gretl_model_test_print(var->models[i], 0, prn);
	gretl_model_destroy_tests(var->models[i]);
    }

    return err;
}

int gretl_VAR_arch_test (GRETL_VAR *var, int order, 
			 DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
	err = arch_test(var->models[i], order, pdinfo, OPT_NONE, prn);
    }

    return err;
}

int VAR_LR_lag_test (GRETL_VAR *var)
{
    double ldet;
    int err = 0;

    ldet = gretl_VAR_ldet(var, &err);

    if (!err) {
	double ll, AIC, BIC, HQC;
	int T = var->T;
	int g = var->neqns;
	int m = var->ncoeff - g;
	int k = g * m;

	var->LR = T * (ldet - var->ldet);

	ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	AIC = (-2.0 * ll + 2.0 * k) / T;
	BIC = (-2.0 * ll + k * log(T)) / T;
	HQC = (-2.0 * ll + 2.0 * k * log(log(T))) / T;
	var->Ivals[0] = AIC;
	var->Ivals[1] = BIC;
	var->Ivals[2] = HQC;
    }

    /* we're done with this set of residuals */
    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

/* make and record residuals for LR test, etc. */

int last_lag_LR_prep (GRETL_VAR *var, int ifc)
{
    int *collist = NULL;
    int g = var->ncoeff - var->neqns;
    int i, err = 0;

    if (var->F == NULL) {
	var->F = gretl_matrix_alloc(var->T, var->neqns);
	if (var->F == NULL) {
	    return E_ALLOC;
	}
    }   

    collist = gretl_list_new(var->neqns);
    if (collist == NULL) {
	return E_ALLOC;
    }

    collist[1] = ifc + var->order - 1;
    for (i=2; i<=collist[0]; i++) {
	collist[i] = collist[i-1] + var->order;
    }

    gretl_matrix_delete_columns(var->X, collist);
    gretl_matrix_reuse(var->B, g, var->neqns);
    err = gretl_matrix_multi_ols(var->Y, var->X, 
				 var->B, var->F,
				 NULL);

    free(collist);

    return err;
}

static void gretl_VAR_print_lagsel (gretl_matrix *lltab,
				    gretl_matrix *crittab,
				    int *best_row,
				    PRN *prn)
{
    int maxlag = gretl_matrix_rows(crittab);
    double x;
    int i, j;

    pprintf(prn, _("VAR system, maximum lag order %d"), maxlag);
    pputs(prn, "\n\n");

    pputs(prn, _("The asterisks below indicate the best (that is, minimized) values\n"
	  "of the respective information criteria, AIC = Akaike criterion,\n"
	  "BIC = Schwartz Bayesian criterion and HQC = Hannan-Quinn criterion."));
    pputs(prn, "\n\n");

    pputs(prn, _("lags        loglik    p(LR)       AIC          BIC          HQC"));
    pputs(prn, "\n\n");

    for (i=0; i<maxlag; i++) {
	pprintf(prn, "%4d", i + 1);
	x = gretl_matrix_get(lltab, i, 0);
	pprintf(prn, "%14.5f", x);
	if (i > 0) {
	    x = gretl_matrix_get(lltab, i, 1);
	    pprintf(prn, "%9.5f", x);
	} else {
	    pputs(prn, "         ");
	}
	for (j=0; j<N_IVALS; j++) {
	    x = gretl_matrix_get(crittab, i, j);
	    pprintf(prn, "%12.6f", x);
	    if (i == best_row[j]) {
		pputc(prn, '*');
	    } else {
		pputc(prn, ' ');
	    }
	}
	pputc(prn, '\n');
    }
}

/* apparatus for selecting the optimal lag length for a VAR */

int VAR_do_lagsel (GRETL_VAR *var, const double **Z, 
		   const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *crittab = NULL;
    gretl_matrix *lltab = NULL;

    int p = var->order;
    int r = p - 1;
    int T = var->T;
    int n = var->neqns;

    /* initialize the "best" at the longest lag */
    double best[N_IVALS] = { var->AIC, var->BIC, var->HQC };
    int best_row[N_IVALS] = { r, r, r };
    double crit[N_IVALS];
    double LRtest;
    double ldet = NADBL;
    int cols0;
    int j, m = 0;
    int err = 0;

    if (p < 2) {
	return 0;
    }

    if (var->F != NULL) {
	gretl_matrix_free(var->F);
    }

    var->F = gretl_matrix_alloc(T, n);
    if (var->F == NULL) {
	return E_ALLOC;
    }

    crittab = gretl_matrix_alloc(p, N_IVALS);
    lltab = gretl_matrix_alloc(p, 2);
    if (crittab == NULL || lltab == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* # of cols in X that are not Y lags */
    cols0 = var->ncoeff - p * n; 

    for (j=1; j<p && !err; j++) {
	int jxcols = cols0 + j * n;

	fill_VAR_X(var, j, Z, pdinfo);

	gretl_matrix_reuse(var->X, T, jxcols);
	gretl_matrix_reuse(var->B, jxcols, n);

	err = gretl_matrix_multi_ols(var->Y, var->X, var->B, 
				     var->F, NULL);

	if (!err) {
	    ldet = gretl_VAR_ldet(var, &err);
	}

	if (!err) {
	    double ll;
	    int q = var->ncoeff - (n * (p - j));
	    int c, k = n * q;

	    ll = -(n * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	    crit[0] = (-2.0 * ll + 2.0 * k) / T;               /* AIC */
	    crit[1] = (-2.0 * ll + k * log(T)) / T;            /* BIC */
	    crit[2] = (-2.0 * ll + 2.0 * k * log(log(T))) / T; /* HQC */

	    gretl_matrix_set(lltab, m, 0, ll);
	    if (j == 1) {
		gretl_matrix_set(lltab, m, 1, 0);
	    } else {
		LRtest = 2.0 * (ll - gretl_matrix_get(lltab, m-1, 0));
		gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, n * n));
	    }	
	    
	    for (c=0; c<N_IVALS; c++) {
		gretl_matrix_set(crittab, m, c, crit[c]);
		if (crit[c] < best[c]) {
		    best[c] = crit[c];
		    best_row[c] = m;
		}
	    }
	
	    m++;
	}
    }

    if (!err) {
	gretl_matrix_set(lltab, m, 0, var->ll);
	LRtest = 2.0 * (var->ll - gretl_matrix_get(lltab, m - 1, 0));
	gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, n * n));
	gretl_matrix_set(crittab, m, 0, var->AIC);
	gretl_matrix_set(crittab, m, 1, var->BIC);
	gretl_matrix_set(crittab, m, 2, var->HQC);
	gretl_VAR_print_lagsel(lltab, crittab, best_row, prn);
    }

    bailout:

    gretl_matrix_free(crittab);
    gretl_matrix_free(lltab);

    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

/* (X'X)^{-1} * X'\Omega X * (X'X)^{-1} : right now
   we're supporting only HC0 and HC1 */

static int VAR_robust_vcv (GRETL_VAR *var, gretl_matrix *V,
			   MODEL *pmod, int hcv, int k)
{
    gretl_matrix *XOX = NULL;
    double xij, xti, xtj, utk;
    int T = var->T;
    int g = var->ncoeff;
    int i, j, t;

    XOX = gretl_matrix_alloc(g, g);
    if (XOX == NULL) {
	return E_ALLOC;
    }

    /* form X' \Omega X */
    for (i=0; i<g; i++) {
	for (j=i; j<g; j++) {
	    xij = 0.0;
	    for (t=0; t<T; t++) {
		xti = gretl_matrix_get(var->X, t, i);
		xtj = gretl_matrix_get(var->X, t, j);
		utk = gretl_matrix_get(var->E, t, k);
		xij += utk * utk * xti * xtj;
	    }
	    if (hcv > 0) {
		/* cheating here, for now */
		xij *= (double) T / (T - g);
	    }
	    gretl_matrix_set(XOX, i, j, xij);
	    if (i != j) {
		gretl_matrix_set(XOX, j, i, xij);
	    }
	}
    }

    gretl_matrix_qform(var->XTX, GRETL_MOD_TRANSPOSE, XOX,
		       V, GRETL_MOD_NONE);

    gretl_model_set_int(pmod, "hc", 1);
    if (hcv > 0) {
	gretl_model_set_int(pmod, "hc_version", 1);
    }

    gretl_matrix_free(XOX);

    return 0;
}


/* Run the various per-equation omit tests (all lags of each var in
   turn, last lag of all vars) using the Wald method.  We also
   add the standard errors to the models here, since we have the
   covariance matrix to hand.
*/

int VAR_wald_omit_tests (GRETL_VAR *var, int ifc)
{
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    int hcv = get_hc_version();
    int p = var->order;
    int n = var->neqns;
    int g = var->ncoeff;
    int dim = (p > n)? p : n;
    int i, j, k, m = 0;
    int err = 0;

    if (ifc && var->robust && g - 1 > dim) {
	/* need bigger arrays for robust overall F-test */
	dim = g - 1;
    }

    V = gretl_matrix_alloc(g, g);
    C = gretl_matrix_alloc(dim, dim);
    b = gretl_column_vector_alloc(dim);

    if (V == NULL || C == NULL || b == NULL) {
	return E_ALLOC;
    }     

    for (i=0; i<n && !err; i++) {
	MODEL *pmod = var->models[i];
	int ii, jj, jpos, ipos = ifc;
	double w, vij;

	gretl_matrix_reuse(V, g, g);

	if (var->robust) {
	    err = VAR_robust_vcv(var, V, pmod, hcv, i);
	} else {
	    gretl_matrix_copy_values(V, var->XTX);
	    gretl_matrix_multiply_by_scalar(V, pmod->sigma * pmod->sigma);
	}
	
	if (!err) {
	    /* set (possibly robust) standard errors */
	    for (j=0; j<g; j++) {
		vij = gretl_matrix_get(V, j, j);
		pmod->sderr[j] = sqrt(vij);
	    }
	}

	/* exclusion of each var, all lags */

	gretl_matrix_reuse(C, p, p);
	gretl_matrix_reuse(b, p, 1);

	for (j=0; j<n && !err; j++) {
	    double w = NADBL;

	    gretl_matrix_extract_matrix(C, V, ipos, ipos, GRETL_MOD_NONE);
	    for (k=0; k<p; k++) {
		b->val[k] = pmod->coeff[k + ipos];
	    }
	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		var->Fvals[m++] = w / p;
	    }

	    ipos += p;
	}

	/* exclusion of last lag, all vars? */

	if (p > 1) {
	    gretl_matrix_reuse(C, n, n);
	    gretl_matrix_reuse(b, n, 1);

	    ipos = ifc + p - 1;
	    for (ii=0; ii<n; ii++) {
		jpos = ifc + p - 1;
		for (jj=0; jj<n; jj++) {
		    vij = gretl_matrix_get(V, ipos, jpos);
		    gretl_matrix_set(C, ii, jj, vij);
		    jpos += p;
		}
		b->val[ii] = pmod->coeff[ipos];
		ipos += p;
	    }

	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		var->Fvals[m++] = w / n;
	    }
	}

	/* exclusion of all but const? */

	if (ifc && var->robust) {
	    gretl_matrix_reuse(C, g-1, g-1);
	    gretl_matrix_reuse(b, g-1, 1);

	    gretl_matrix_extract_matrix(C, V, 1, 1, GRETL_MOD_NONE);
	    for (k=0; k<g-1; k++) {
		b->val[k] = pmod->coeff[k+1];
	    }
	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		pmod->fstt = w / (g-1);
	    }
	}
    }

    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return err;
}

#define VO_DEBUG 0

const int *gretl_VAR_get_exo_list (const GRETL_VAR *var)
{
    return var->xlist;
}

/* Based on the specification stored in the VAR struct, reconstitute
   the list that was intially passed to the gretl_VAR() function to
   set up the system.
*/

static int *rebuild_VAR_list (const GRETL_VAR *var, int *exolist, int *err)
{
    int *list = NULL;
    int lsep = (exolist[0] > 0);
    int i, j = 1;

    list = gretl_list_new(var->neqns + exolist[0] + lsep);
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<var->neqns; i++) {
	list[j++] = var->ylist[i+1];
    }

    if (lsep) {
	list[j++] = LISTSEP;
    }

    for (i=1; i<=exolist[0]; i++) {
	list[j++] = exolist[i];
    }    

    return list;
}

static int gretl_VAR_real_omit_test (const GRETL_VAR *orig,
				     const GRETL_VAR *new,
				     const DATAINFO *pdinfo,
				     PRN *prn)
{
    int *omitlist;
    double LR, pval;
    int i, df, err = 0;

#if VO_DEBUG
    fprintf(stderr, "gretl_VAR_real_omit_test: about to diff lists\n");
    printlist(orig->xlist, "orig xlist");
    printlist(new->xlist, "new xlist");
#endif

    if (new->xlist == NULL) {
	omitlist = gretl_list_copy(orig->xlist);
    } else {
	omitlist = gretl_list_diff_new(orig->xlist, new->xlist, 1);
    }

    if (omitlist == NULL) {
	return E_ALLOC;
    }

    LR = orig->T * (new->ldet - orig->ldet);
    df = orig->neqns * omitlist[0];
    pval = chisq_cdf_comp(LR, df);
    
    pputs(prn, _("\n  Null hypothesis: the regression parameters are "
		 "zero for the variables\n\n"));
    for (i=1; i<=omitlist[0]; i++) {
	pprintf(prn, "    %s\n", pdinfo->varname[omitlist[i]]);	
    }

    pprintf(prn, "\n  %s: %s(%d) = %g, ", _("Test statistic"), 
	    _("Chi-square"), df, LR);
    pprintf(prn, _("with p-value = %g\n\n"), pval);

    free(omitlist);

    return err;
}

/**
 * gretl_VAR_omit_test:
 * @omitvars: list of variables to omit from original model.
 * @orig: pointer to original VAR.
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Re-estimates a given VAR after removing the variables
 * specified in @omitvars, and reports per-equation F-tests
 * and system-wide LR tests for the null hypothesis that
 * the omitted variables have zero parameters.
 * 
 * Returns: restricted VAR on sucess, %NULL on error.
 */

GRETL_VAR *gretl_VAR_omit_test (const int *omitvars, const GRETL_VAR *orig, 
				const double **Z, DATAINFO *pdinfo, 
				PRN *prn, int *err)
{
    GRETL_VAR *var = NULL;
    gretlopt opt = OPT_NONE;
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *tmplist = NULL;
    int *varlist = NULL;
    int c1 = 0;

    *err = 0;

    if (orig == NULL || orig->xlist == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (omitvars == NULL || omitvars[0] == 0) {
	*err = E_PARSE;
	return NULL;
    }

#if VO_DEBUG
    printlist(orig->xlist, "original xlist");
#endif

    if (orig->ifc) {
	c1 = !gretl_list_const_pos(omitvars, 1, Z, pdinfo);
    } 

    /* create reduced exogenous vars list for test VAR */
    tmplist = gretl_list_omit(orig->xlist, omitvars, 1, err);
    if (tmplist == NULL) {
	goto bailout;
    }

#if VO_DEBUG
    fprintf(stderr, "c1 = %d\n", c1);
    printlist(tmplist, "exog vars list for test VAR");
#endif

    /* recreate full input VAR list for test VAR */
    varlist = rebuild_VAR_list(orig, tmplist, err);
    if (varlist == NULL) {
	goto bailout;
    }

#if VO_DEBUG
    printlist(varlist, "full list for test VAR");
#endif

    if (orig->detflags & DET_SEAS) {
	opt |= OPT_D;
    }

    if (orig->detflags & DET_TREND) {
	opt |= OPT_T;
    }

    /* If the original VAR did not include a constant, we need to
       pass OPT_N to the test VAR to suppress the constant.
       We also need to pass OPT_N in case the constant was
       present originally but is now to be omitted.
    */
    if (orig->ifc == 0 || c1 == 0) {
	opt |= OPT_N;
    }

    /* impose as sample range the estimation range of the 
       original VAR */
    pdinfo->t1 = orig->t1;
    pdinfo->t2 = orig->t2;

    var = gretl_VAR(orig->order, varlist, Z, pdinfo, opt, prn, err);

    /* now, if var is non-NULL, do the actual test(s) */
    if (var != NULL) {
	*err = gretl_VAR_real_omit_test(orig, var, pdinfo, prn);
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

 bailout:

    free(tmplist);
    free(varlist);

    return var;
}


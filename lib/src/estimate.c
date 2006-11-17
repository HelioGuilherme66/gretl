/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* estimate.c - gretl estimation procedures */

#include "libgretl.h"
#include "qr_estimate.h"
#include "gretl_panel.h"
#include "libset.h"
#include "compat.h"
#include "missing_private.h"
#include "estim_private.h"

/* There's a balancing act with 'TINY' here.  It's the minimum value
   for test that libgretl will accept before rejecting a
   data matrix as too highly collinear.  If you set it too high,
   data sets for which gretl could produce reasonable estimates will
   be rejected.  If you set it too low (and even 100 * DBL_EPSILON
   is definitely too low), gretl will produce more or less worthless
   coefficient estimates when given highly collinear data.  If you're
   tempted to change the value of TINY, check how gretl does on the
   NIST reference data sets for linear regression and ensure you're
   not getting any garbage results.  The setting of 2.1e-09 enables
   me to get decent results on the NIST nonlinear regression test
   suite, but it could be a bit too low for some contexts.
*/

#define TINY      2.1e-09 /* was 4.75e-09 (last changed 2004/07/16) */
#define SMALL     1.0e-08 /* threshold for printing a warning for collinearity */
#define YBARZERO  0.5e-14 /* threshold for treating mean of dependent
			     variable as effectively zero */
#define ESSZERO   1.0e-22 /* threshold for considering a tiny error-sum-of-
			     squares value to be effectively zero */

/* define for lots of debugging info */
#define XPX_DEBUG 0

/* private function prototypes */
static int form_xpxxpy (const int *list, int t1, int t2, 
			const double **Z, int nwt, double rho, int pwe,
			double *xpx, double *xpy, const char *mask);
static void regress (MODEL *pmod, double *xpy, double **Z, 
		     int n, double rho);
static int cholbeta (double *xpx, double *xpy, double *coeff, int nv, 
		     double *rss);
static void diaginv (double *xpx, double *xpy, double *diag, int nv);

static int hatvar (MODEL *pmod, int n, double **Z);
static void omitzero (MODEL *pmod, const double **Z, const DATAINFO *pdinfo);
static int depvar_zero (int t1, int t2, int yno, int nwt, 
			const double **Z);
static int lagdepvar (const int *list, const double **Z, const DATAINFO *pdinfo); 
static int jackknife_vcv (MODEL *pmod, const double **Z);
/* end private protos */


/* compute statistics for the dependent variable in a model */

static void model_depvar_stats (MODEL *pmod, const double **Z)
{
    double xx, sum = 0.0;
    int yno = pmod->list[1];
    int t, dwt = 0;

    if (pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy")) {
	dwt = pmod->nwt;
    }

    pmod->ybar = pmod->sdy = NADBL;

    if (pmod->nobs <= 0) {
	return;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}
	if (!model_missing(pmod, t)) {
	    sum += Z[yno][t];
	}
    }

    pmod->ybar = sum / pmod->nobs;

    sum = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}	
	if (!model_missing(pmod, t)) {
	    sum += (Z[yno][t] - pmod->ybar); 
	}
    }

    pmod->ybar = pmod->ybar + sum / pmod->nobs;

    if (fabs(pmod->ybar) < YBARZERO) {
	pmod->ybar = 0.0;
    }

    sum = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}
	if (!model_missing(pmod, t)) {
	    xx = Z[yno][t] - pmod->ybar;
	    sum += xx * xx;
	}
    }

    sum = (pmod->nobs > 1)? sum / (pmod->nobs - 1) : 0.0;

    pmod->sdy = (sum >= 0)? sqrt(sum) : NADBL;
}

/* determine the degrees of freedom for a model */

static int get_model_df (MODEL *pmod)
{
    int err = 0;

    pmod->ncoeff = pmod->list[0] - 1;
    pmod->dfd = pmod->nobs - pmod->ncoeff;

    if (pmod->dfd < 0) {
	pmod->errcode = E_DF;
        sprintf(gretl_errmsg, _("No. of obs (%d) is less than no. "
		"of parameters (%d)"), pmod->nobs, pmod->ncoeff);
	err = 1;
    } else {
	pmod->dfn = pmod->ncoeff - pmod->ifc;
    }

    return err;
}

#define LDDEBUG 0

static int
transcribe_ld_vcv (MODEL *targ, MODEL *src)
{
    int nv = targ->ncoeff;
    int nxpx = (nv * nv + nv) / 2;
    int i, j;

    if (makevcv(src, src->sigma)) {
	return 1;
    }

    if (targ->vcv == NULL) {
	targ->vcv = malloc(nxpx * sizeof *targ->vcv);
	if (targ->vcv == NULL) {
	    return 1;
	}
    }

    for (i=0; i<nv; i++) {
	for (j=i; j<nv; j++) {
	    targ->vcv[ijton(i, j, nv)] = 
		src->vcv[ijton(i, j, src->ncoeff)];
	}
    }

    return 0;
}

/* Calculate consistent standard errors (and VCV matrix) when doing
   AR(1) estimation of a model with lagged dependent variable.  See
   Ramanathan, Introductory Econometrics, 5e, p. 450.
*/

static int 
ldepvar_std_errors (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    MODEL emod;
    const double *x;

    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;

    double rho = gretl_model_get_double(pmod, "rho_in");
    int origv = pdinfo->v;
    int vnew = pmod->list[0] + 1 - pmod->ifc;

    int *list;
    int vi, vm;
    int i, t;
    int err = 0;

#if LDDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
    printlist(pmod->list, "pmod->list");
    printf("vnew = %d\n", vnew);
    printf("rho = %g\n", rho);
#endif

    list = gretl_list_new(vnew + pmod->ifc);
    if (list == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    err = dataset_add_series(vnew, pZ, pdinfo);
    if (err) {
	free(list);
	pmod->errcode = E_ALLOC;
	return 1;
    }

    vi = origv;

    /* dependent var is residual from original model */
    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[vi][t] = pmod->uhat[t];
    }    
    strcpy(pdinfo->varname[vi], "eps");
    list[1] = vi++;

    /* indep vars are rho-differenced vars from original model */
    for (i=2; i<=pmod->list[0]; i++) {
	vm = pmod->list[i];
	if (vm == 0) {
	    list[i] = 0;
	    continue;
	}
	sprintf(pdinfo->varname[vi], "%.6s_r", pdinfo->varname[vm]);
	x = (*pZ)[vm];
	for (t=0; t<pdinfo->n; t++) {
	    if (t == 0 || na(x[t]) || na(x[t-1])) {
		(*pZ)[vi][t] = NADBL;
	    } else {
		(*pZ)[vi][t] = x[t] - rho * x[t-1];
	    }
	}
	list[i] = vi++;
    }

    /* last indep var is lagged u-hat */
    for (t=0; t<pdinfo->n; t++) {
	if (t == 0) {
	    (*pZ)[vi][t] = NADBL;
	} else { 
	    (*pZ)[vi][t] = (*pZ)[pmod->list[1]][t-1];
	}
	if (na((*pZ)[vi][t])) {
	    continue;
	}
	for (i=0; i<pmod->ncoeff; i++) {
	    x = (*pZ)[pmod->list[i+2]];
	    if (na(x[t-1])) {
		(*pZ)[vi][t] = NADBL;
		break;
	    } else {
		(*pZ)[vi][t] -= pmod->coeff[i] * x[t-1];
	    }
	}
    }

    list[list[0]] = vi;
    strcpy(pdinfo->varname[vi], "uhat_1");

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    emod = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (emod.errcode) {
	err = emod.errcode;
    } else {
#if LDDEBUG
	printmodel(&emod, pdinfo, OPT_NONE, prn);
	gretl_print_destroy(prn);
#endif
	for (i=0; i<pmod->ncoeff; i++) {
	    pmod->sderr[i] = emod.sderr[i];
	}

	err = transcribe_ld_vcv(pmod, &emod);
    }
    
    clear_model(&emod);
    
    free(list);
    dataset_drop_last_variables(vnew, pZ, pdinfo);

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    if (err) {
	pmod->errcode = err;
    }

    return err;
}

/* special computation of statistics for autoregressive models */

static int compute_ar_stats (MODEL *pmod, const double **Z, double rho)
{
    int i, t, yno = pmod->list[1];
    double x, pw1 = 0.0;

    if (gretl_model_add_arinfo(pmod, 1)) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    if (pmod->ci == PWE) {
	pw1 = sqrt(1.0 - rho * rho);
    }

    pmod->arinfo->arlist[0] = pmod->arinfo->arlist[1] = 1;
    pmod->arinfo->rho[0] = rho;
    gretl_model_set_double(pmod, "rho_in", rho);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (t == pmod->t1 && pmod->ci == PWE) {
	    x = pw1 * Z[yno][t];
	    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
		x -= pmod->coeff[i] * pw1 * Z[pmod->list[i+2]][t];
	    }
	    if (pmod->ifc) {
		x -= pw1 * pmod->coeff[0];
	    }
	} else {
	    x = Z[yno][t] - rho * Z[yno][t-1];
	    for (i=0; i<pmod->ncoeff; i++) {
		x -= pmod->coeff[i] * 
		    (Z[pmod->list[i+2]][t] - 
		     rho * Z[pmod->list[i+2]][t-1]);
	    }
	}
	pmod->uhat[t] = x;
	pmod->yhat[t] = Z[yno][t] - x;
    }

    pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, Z[yno], pmod->yhat);

    pmod->adjrsq = 
	1.0 - ((1.0 - pmod->rsq) * (pmod->t2 - pmod->t1) / 
	       (double) pmod->dfd);

    return 0;
}

/* Durbin-Watson statistic for pooled model */

static void panel_dwstat (MODEL *pmod, const DATAINFO *pdinfo)
{
    double ut, us, num = 0.0;
    int s, t;

    pmod->rho = NADBL;

    for (t=pmod->t1+1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	s = t - 1;
	if (na(pmod->uhat[s])) {
	    continue;
	}
	if (pdinfo->paninfo->unit[t] != pdinfo->paninfo->unit[s]) {
	    continue;
	}
	if (pdinfo->paninfo->period[t] != pdinfo->paninfo->period[s] + 1) {
	    continue;
	}
	ut = pmod->uhat[t];
	us = pmod->uhat[s];
	num += (ut - us) * (ut - us);
    }

    pmod->dw = num / pmod->ess;
}

/* calculation of WLS stats in agreement with GNU R */

static void get_wls_stats (MODEL *pmod, const double **Z)
{
    int t, wobs = pmod->nobs, yno = pmod->list[1];
    double x, dy, wmean = 0.0, wsum = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) { 
	if (model_missing(pmod, t)) {
	    continue;
	}
	if (Z[pmod->nwt][t] == 0.0) {
	    wobs--;
	    pmod->dfd -= 1;
	} else {
	    wmean += Z[pmod->nwt][t] * Z[yno][t];
	    wsum += Z[pmod->nwt][t];
	}
    }

    wmean /= wsum;
    x = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t) || Z[pmod->nwt][t] == 0.0) {
	    continue;
	}	
	dy = Z[yno][t] - wmean;
	x += Z[pmod->nwt][t] * dy * dy;
    }

    pmod->fstt = ((x - pmod->ess) * pmod->dfd) / (pmod->dfn * pmod->ess);
    pmod->rsq = (1 - (pmod->ess / x));
    pmod->adjrsq = 1 - ((1 - pmod->rsq) * (pmod->nobs - 1)/pmod->dfd);
}

static void fix_wls_values (MODEL *pmod, double **Z)
{
    int t;

    if (gretl_model_get_int(pmod, "wt_dummy")) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    }
	}
    } else {
	double ess_orig = 0.0;
	double sw, sigma_orig;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
		pmod->nobs -= 1;
	    } else {
		sw = sqrt(Z[pmod->nwt][t]);
		pmod->yhat[t] /= sw;
		pmod->uhat[t] /= sw;
		ess_orig += pmod->uhat[t] * pmod->uhat[t];
	    }
	}

	sigma_orig = sqrt(ess_orig / pmod->dfd);
	gretl_model_set_double(pmod, "ess_orig", ess_orig);
	gretl_model_set_double(pmod, "sigma_orig", sigma_orig);
    }
}

/* drop the weight var from the list of regressors (WLS) */

static void dropwt (int *list)
{
    int i;

    list[0] -= 1;
    for (i=1; i<=list[0]; i++) {
	list[i] = list[i+1];
    }
}

static void model_stats_init (MODEL *pmod)
{
    pmod->ess = NADBL;
    pmod->sigma = NADBL;
    pmod->fstt = pmod->lnL = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
}

#define SMPL_DEBUG 0

static int 
lsq_check_for_missing_obs (MODEL *pmod, gretlopt opts,
			   DATAINFO *pdinfo, const double **Z, 
			   int *misst)
{
    int missv = 0;
    int reject_missing = 0;

    if (reference_missmask_present()) {
	int err = apply_reference_missmask(pmod);

#if SMPL_DEBUG
	fprintf(stderr, "missmask found, applied with err = %d\n",
		err);
#endif
	/* If there was a reference mask present, it was put there
	   as part of a hypothesis test on some original model, and
	   it has to be respected in estimation of this model */

	if (err) {
	    pmod->errcode = E_ALLOC;
	    return 1;
	} else {
	    return 0;
	}
    }

    /* can't do HAC VCV with missing obs in middle */
    if ((opts & OPT_R) && dataset_is_time_series(pdinfo) &&
	!get_force_hc()) {
	reject_missing = 1;
    } 

    if (opts & OPT_M) {
	reject_missing = 1;
    }

    if (reject_missing) {
	/* reject missing obs within adjusted sample */
	missv = adjust_t1t2(pmod, pmod->list, &pmod->t1, &pmod->t2,
			    pdinfo->n, Z, misst);
    } else {
	/* we'll try to compensate for missing obs */
	missv = adjust_t1t2(pmod, pmod->list, &pmod->t1, &pmod->t2,
			    pdinfo->n, Z, NULL);
    }

#if SMPL_DEBUG
    if (1) {
	char t1s[OBSLEN], t2s[OBSLEN];
	int misscount = model_missval_count(pmod);

	ntodate(t1s, pmod->t1, pdinfo);
	ntodate(t2s, pmod->t2, pdinfo);
	fprintf(stderr, "*** after adjustment, t1=%d (%s), t2=%d (%s)\n", 
		pmod->t1, t1s, pmod->t2, t2s);
	fprintf(stderr, "Valid observations in range = %d\n", 
		pmod->t2 - pmod->t1 + 1 - misscount);
    }
#endif

    return missv;
}

static void 
lagged_depvar_check (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    int ldv = lagdepvar(pmod->list, Z, pdinfo);

    if (ldv) {
	gretl_model_set_int(pmod, "ldepvar", ldv);
    } else if (gretl_model_get_int(pmod, "ldepvar")) {
	gretl_model_set_int(pmod, "ldepvar", 0);
    }
}

static void 
log_depvar_ll (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    char parent[VNAMELEN];

    if (is_log_variable(pmod->list[1], pdinfo, parent)) {
	double jll = pmod->lnL;
	int t;

	for (t=0; t<pdinfo->n; t++) {
	    if (!na(pmod->uhat[t])) {
		jll -= Z[pmod->list[1]][t];
	    }
	}
	gretl_model_set_double(pmod, "jll", jll);
	gretl_model_set_string_as_data(pmod, 
				       "log-parent", 
				       gretl_strdup(parent));
    }
}

#define COLL_DEBUG 0

int 
redundant_var (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, int **droplist)
{
    MODEL cmod;
    int targ, l0;
    int *list;
    int err = 0;
    int i, ret = 0;

    if (pmod->list[0] < 3) {
	/* shouldn't happen */
	return 0;
    }

    for (i=1; i<=pmod->list[0]; i++) {
	if (pmod->list[i] == LISTSEP) {
	    /* can't handle compound lists */
	    return 0;
	}
    }

    l0 = pmod->list[0] - 1;
    list = gretl_list_new(l0);
    if (list == NULL) {
	return 0;
    }

#if COLL_DEBUG
    fprintf(stderr, "\n*** redundant_var called ***\n");
    printlist(pmod->list, "original model list");
#endif

    list[1] = pmod->list[1];

    /* back up along the list of regressors, trying to find a single
       variable such that, when it is deleted, the exact collinearity
       problem goes away. */

    for (targ=pmod->list[0]; targ>2; targ--) {
	int j = 2;

	for (i=2; i<=pmod->list[0]; i++) {
	    if (i != targ) {
		list[j++] = pmod->list[i];
	    }
	}

#if COLL_DEBUG
	fprintf(stderr, "pass 1: target list position = %d\n", targ);
	printlist(list, "temp list for redundancy check");
#endif
	cmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_Z);

	if (cmod.errcode == 0) {
	    ret = 1;
	} else if (cmod.errcode != E_SINGULAR) {
	    /* shouldn't happen */
	    err = 1;
	}

	clear_model(&cmod);

	if (ret || err) {
	    break;
	}
    } 

    /* if that didn't work, try trimming the list of regressors more
       aggressively: look for a variable that is perfectly predicted
       by the other regressors
    */

    for (l0=pmod->list[0]; !ret && !err && l0>2; l0--) {

 	list[0] = l0 - 1;

 	for (targ=l0; targ>2; targ--) {
 	    int j = 2;
 
 	    list[1] = pmod->list[targ];
 
 	    for (i=2; i<=l0; i++) {
 		if (i != targ) {
 		    list[j++] = pmod->list[i];
 		}
 	    }
 
#if COLL_DEBUG
 	    fprintf(stderr, "pass 2: target list position = %d\n", targ);
 	    printlist(list, "temp list for redundancy check");
#endif
 	    cmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_Z);

	    if (cmod.errcode == 0) {
		if (cmod.ess == 0.0 || cmod.rsq == 1.0) {
		    ret = 1;
		}
	    } else if (cmod.errcode != E_SINGULAR) {
		err = 1;
	    }
 
 	    clear_model(&cmod);

	    if (ret || err) {
		break;
	    }
 	}
    }  

    if (ret) {
	int v = pmod->list[targ];

	/* remove var from list and reduce number of coeffs */
	gretl_list_delete_at_pos(pmod->list, targ);
	pmod->ncoeff -= 1;
	pmod->dfd = pmod->nobs - pmod->ncoeff;
	pmod->dfn = pmod->ncoeff - pmod->ifc;

	/* add redundant var to list of drops */
	gretl_list_append_term(droplist, v);
    }

    free(list);

    return ret;
}

static int check_weight_var (MODEL *pmod, const double *w, int *effobs)
{
    int t;

    if (gretl_iszero(pmod->t1, pmod->t2, w)) {
	pmod->errcode = E_WTZERO;
	return 1;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (w[t] < 0.0) {
	    pmod->errcode = E_WTNEG;
	    return 1;
	}
    }

    *effobs = gretl_isdummy(pmod->t1, pmod->t2, w);

    if (*effobs) {
	/* the weight var is a dummy, with effobs 1s */
	gretl_model_set_int(pmod, "wt_dummy", 1);
    }

    return 0;
}

void 
maybe_shift_ldepvar (MODEL *pmod, const double **Z, DATAINFO *pdinfo)
{
    if (gretl_model_get_int(pmod, "ldepvar")) {
	lagged_depvar_check(pmod, Z, pdinfo);
    }
}

static int gretl_choleski_regress (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
				   double rho, int pwe, gretlopt opt)
{
    double *xpy = NULL;
    int *droplist = NULL;
    int l0 = pmod->list[0];
    int i, k, nxpx;

 trim_var:

    if (droplist != NULL) {
	l0 = pmod->list[0];
	free(pmod->xpx);
	free(pmod->coeff);
	free(pmod->sderr);
	pmod->errcode = 0;
    }
 
    k = l0 - 1;
    nxpx = k * (k + 1) / 2;

    if (nxpx == 0) {
	fprintf(stderr, "problem: nxpx = 0 (l0 = %d)\n", l0);
	pmod->errcode = E_DATA;
	return pmod->errcode;
    }

    xpy = malloc((l0 + 1) * sizeof *xpy);
    pmod->xpx = malloc(nxpx * sizeof *pmod->xpx);
    pmod->coeff = malloc(pmod->ncoeff * sizeof *pmod->coeff);
    pmod->sderr = malloc(pmod->ncoeff * sizeof *pmod->sderr);

    if (pmod->yhat == NULL) {
	pmod->yhat = malloc(pdinfo->n * sizeof *pmod->yhat);
    } 

    if (pmod->uhat == NULL) {
	pmod->uhat = malloc(pdinfo->n * sizeof *pmod->uhat);
    }

    if (xpy == NULL || pmod->xpx == NULL || pmod->coeff == NULL ||
	pmod->sderr == NULL || pmod->yhat == NULL || pmod->uhat == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    for (i=0; i<=l0; i++) {
	xpy[i] = 0.0;
    }
    for (i=0; i<nxpx; i++) {
	pmod->xpx[i] = 0.0;
    }

    /* calculate regression results, Cholesky style */
    form_xpxxpy(pmod->list, pmod->t1, pmod->t2, (const double **) *pZ, 
		pmod->nwt, rho, pwe, pmod->xpx, xpy, pmod->missmask);

#if XPX_DEBUG
    for (i=0; i<=l0; i++) {
	fprintf(stderr, "xpy[%d] = %g\n", i, xpy[i]);
    }
    for (i=0; i<nxpx; i++) {
	fprintf(stderr, "xpx[%d] = %g\n", i, pmod->xpx[i]);
    }
    fputc('\n', stderr);
#endif

    regress(pmod, xpy, *pZ, pdinfo->n, rho);
    free(xpy);

    if (pmod->errcode == E_SINGULAR && !(opt & OPT_Z) &&
	redundant_var(pmod, pZ, pdinfo, &droplist)) {
	goto trim_var;
    }

    if (droplist != NULL) {
	/* if there's a lagged dep var, it may have moved */
	maybe_shift_ldepvar(pmod, (const double **) *pZ, pdinfo);
	gretl_model_set_list_as_data(pmod, "droplist", droplist);
    }

    return pmod->errcode;
}

/* as lsq() below, except that we allow for a non-zero value
   of the first-order quasi-differencing coefficient, rho,
   and there's no PRN.
*/

MODEL ar1_lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	       GretlCmdIndex ci, gretlopt opt, double rho)
{
    MODEL mdl;
    int effobs = 0;
    int missv = 0, misst = 0;
    int jackknife = 0;
    int use_qr = get_use_qr();
    int pwe = (ci == PWE || (opt & OPT_P));
    int yno, i;

    *gretl_errmsg = '\0';

    if (list == NULL || pZ == NULL || pdinfo == NULL) {
	fprintf(stderr, "E_DATA: lsq: list = %p, pZ = %p, pdinfo = %p\n",
		(void *) list, (void *) pZ, (void *) pdinfo);
	mdl.errcode = E_DATA;
        return mdl;
    }

    if (ci == HSK) {
	return hsk_func(list, pZ, pdinfo);
    } 

    gretl_model_init(&mdl);
    gretl_model_smpl_init(&mdl, pdinfo);
    model_stats_init(&mdl);

    if (pwe) {
	gretl_model_set_int(&mdl, "pwe", 1);
    }

    if (list[0] == 1 || pdinfo->v == 1) {
	fprintf(stderr, "E_DATA: lsq: list[0] = %d, pdinfo->v = %d\n",
		list[0], pdinfo->v);
	mdl.errcode = E_DATA;
        return mdl;
    }

    /* preserve a copy of the list supplied, for future reference */
    mdl.list = gretl_list_copy(list);
    if (mdl.list == NULL) {
        mdl.errcode = E_ALLOC;
        return mdl;
    }

    mdl.t1 = pdinfo->t1;
    mdl.t2 = pdinfo->t2;
    mdl.full_n = pdinfo->n;
    mdl.ci = ci;

    /* Doing weighted least squares? */
    if (ci == WLS) { 
	check_weight_var(&mdl, (*pZ)[mdl.list[1]], &effobs);
	if (mdl.errcode) {
	    return mdl;
	}
	mdl.nwt = mdl.list[1];
    } else {
	mdl.nwt = 0;
    }

    /* sanity check */
    if (mdl.t1 < 0 || mdl.t2 > pdinfo->n - 1) {
        mdl.errcode = E_NODATA;
        goto lsq_abort;
    }

    /* adjust sample range and check for missing obs: this
       may set the model errcode */
    missv = lsq_check_for_missing_obs(&mdl, opt, pdinfo,
				      (const double **) *pZ,
				      &misst);
    if (mdl.errcode) {
        goto lsq_abort;
    }

    /* react to presence of unhandled missing obs */
    if (missv) {
	if (dated_daily_data(pdinfo)) {
	    if (repack_missing_daily_obs(&mdl, *pZ, pdinfo)) {
		return mdl;
	    }
	} else {
	    sprintf(gretl_errmsg, _("Missing value encountered for "
		    "variable %d, obs %d"), missv, misst);
	    mdl.errcode = E_DATA;
	    return mdl;
	} 
    }

    if (ci == WLS) {
	dropwt(mdl.list);
    }

    yno = mdl.list[1];
    
    /* check for unknown vars in list */
    for (i=1; i<=mdl.list[0]; i++) {
        if (mdl.list[i] > pdinfo->v - 1) {
            mdl.errcode = E_UNKVAR;
            goto lsq_abort;
        }
    } 

    /* check for zero dependent var */
    if (depvar_zero(mdl.t1, mdl.t2, yno, mdl.nwt, (const double **) *pZ)) {  
        mdl.errcode = E_ZERO;
        goto lsq_abort; 
    } 

    /* drop any vars that are all zero and repack the list */
    omitzero(&mdl, (const double **) *pZ, pdinfo);

    /* if regressor list contains a constant, record this fact and 
       place it first among the regressors */
    mdl.ifc = reglist_check_for_const(mdl.list, (const double **) *pZ, 
				      pdinfo);

    /* Check for presence of lagged dependent variable? 
       (Don't bother if this is an auxiliary regression.) */
    if (!(opt & OPT_A)) {
	lagged_depvar_check(&mdl, (const double **) *pZ, pdinfo);
    }

    /* AR1: advance the starting observation by one? */
    if (rho != 0.0 && !pwe) {
	mdl.t1 += 1;
    }

    mdl.ncoeff = mdl.list[0] - 1; 
    if (effobs) {
	mdl.nobs = effobs; /* FIXME? */
    } else {
	mdl.nobs = mdl.t2 - mdl.t1 + 1;
	if (has_missing_obs(&mdl)) {
	    mdl.nobs -= model_missval_count(&mdl);
	}
    }

    /* check degrees of freedom */
    if (get_model_df(&mdl)) {
        goto lsq_abort; 
    }

    /* if df correction is not wanted, record this fact */
    if (opt & OPT_N) {
	gretl_model_set_int(&mdl, "no-df-corr", 1);
    }

    if (dataset_is_time_series(pdinfo)) {
	opt |= OPT_T;
    }

    if (mdl.ci == HCCM || ((opt & OPT_R) && get_hc_version() == 4)) {
	jackknife = 1;
    }

    if (!jackknife && ((opt & OPT_R) || (use_qr && !(opt & OPT_C)))) { 
	mdl.rho = rho;
	gretl_qr_regress(&mdl, pZ, pdinfo, opt);
    } else {
	gretl_choleski_regress(&mdl, pZ, pdinfo, rho, pwe, opt);
    }

    if (mdl.errcode) {
	goto lsq_abort;
    }

    /* get the mean and sd of dep. var. and make available */
    model_depvar_stats(&mdl, (const double **) *pZ);

    /* Doing an autoregressive procedure? */
    if (ci == CORC || ci == HILU || ci == PWE) {
	if (compute_ar_stats(&mdl, (const double **) *pZ, rho)) { 
	    goto lsq_abort;
	}
	if (gretl_model_get_int(&mdl, "ldepvar")) {
	    if (ldepvar_std_errors(&mdl, pZ, pdinfo)) {
		goto lsq_abort;
	    }
	}
	if (ci == HILU && (opt & OPT_B)) {
	    gretl_model_set_int(&mdl, "no-corc", 1);
	}
    }

    /* weighted least squares: fix yhat and uhat; add calculation of
       ESS and sigma based on unweighted data
    */
    if (ci == WLS) {
	get_wls_stats(&mdl, (const double **) *pZ);
	fix_wls_values(&mdl, *pZ);
    }

    if (mdl.missmask == NULL) {
	if (opt & OPT_T) {
	    mdl.rho = rhohat(1, mdl.t1, mdl.t2, mdl.uhat);
	    mdl.dw = dwstat(1, &mdl, (const double **) *pZ);
	} else if (dataset_is_panel(pdinfo)) {
	    panel_dwstat(&mdl, pdinfo);
	}
    } else {
	mdl.rho = mdl.dw = NADBL;
    }

    /* weird special case: degenerate model */
    if (mdl.ncoeff == 1 && mdl.ifc) {
	mdl.rsq = mdl.adjrsq = 0.0;
	mdl.fstt = NADBL;
    }

    /* Generate model selection statistics */
    ls_criteria(&mdl);
    if (!(opt & OPT_A) && !na(mdl.lnL)) {
	log_depvar_ll(&mdl, (const double **) *pZ, pdinfo);
    }

    /* hccm command or HC3a */
    if (jackknife) {
	mdl.errcode = jackknife_vcv(&mdl, (const double **) *pZ);
    }

 lsq_abort:

    /* If we reshuffled any missing observations, put them
       back in their right places now */
    if (gretl_model_get_int(&mdl, "daily_repack")) {
	undo_daily_repack(&mdl, *pZ, pdinfo);
    }

    if (!(opt & OPT_A)) {
	/* if it's not an auxiliary regression, set an ID number
	   on the model */
	set_model_id(&mdl);
    }

    return mdl;
}

/**
 * lsq:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: one of the command indices in #LSQ_MODEL.
 * @opt: option flags: zero or more of the following --
 *   %OPT_R compute robust standard errors;
 *   %OPT_C force use of Cholesky decomp;
 *   %OPT_A treat as auxiliary regression (don't bother checking
 *     for presence of lagged dependent var, don't augment model count);
 *   %OPT_P use Prais-Winsten for first obs;
 *   %OPT_N don't use degrees of freedom correction for standard
 *      error of regression;
 *   %OPT_M reject missing observations within sample range;
 *   %OPT_Z (internal use) suppress the automatic elimination of 
 *      perfectly collinear variables.
 *
 * Computes least squares estimates of the model specified by @list,
 * using an estimator determined by the value of @ci.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	   GretlCmdIndex ci, gretlopt opt)
{
    return ar1_lsq(list, pZ, pdinfo, ci, opt, 0.0);
}

/*
  form_xpxxpy: form the X'X matrix and X'y vector

  - if rho is non-zero, quasi-difference the data first
  - if nwt is non-zero, use that variable as weight
  - if pwe is non-zero (as well as rho) construct the
    first observation as per Prais-Winsten

    Z[v][t] = observation t on variable v
    n = number of obs in data set
    t1, t2 = starting and ending observations
    rho = quasi-differencing coefficent
    nwt = ID number of variable used as weight

    xpx = X'X matrix as a lower triangle, stacked by columns
    xpy = X'y vector
    xpy[0] = sum of y
    xpy[list[0]] = y'y
*/

static int form_xpxxpy (const int *list, int t1, int t2, 
			const double **Z, int nwt, double rho, int pwe,
			double *xpx, double *xpy, const char *mask)
{
    int i, j, t;
    int li, lj, m;
    int l0 = list[0], yno = list[1];
    double x, pw1;
    int qdiff = (rho != 0.0);

    /* Prais-Winsten term */
    if (qdiff && pwe) {
	pw1 = sqrt(1.0 - rho * rho);
    } else {
	pwe = 0;
	pw1 = 0.0;
    }

    xpy[0] = xpy[l0] = 0.0;

    for (t=t1; t<=t2; t++) {
	if (missing_masked(mask, t)) {
	    continue;
	}
	x = Z[yno][t]; 
	if (qdiff) {
	    if (pwe && t == t1) {
		x = pw1 * Z[yno][t];
	    } else {
		x -= rho * Z[yno][t-1];
	    }
	} else if (nwt) {
	    x *= sqrt(Z[nwt][t]);
	}
        xpy[0] += x;
        xpy[l0] += x * x;
    }

    if (xpy[l0] <= 0.0) {
         return yno; 
    }    

    m = 0;

    if (qdiff) {
	/* quasi-difference the data */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (pwe && t == t1) {
			x += pw1 * Z[li][t1] * pw1 * Z[lj][t];
		    } else {
			x += (Z[li][t] - rho * Z[li][t-1]) * 
			    (Z[lj][t] - rho * Z[lj][t-1]);
		    }
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		if (pwe && t == t1) {
		    x += pw1 * Z[yno][t] * pw1 * Z[li][t];
		} else {
		    x += (Z[yno][t] - rho * Z[yno][t-1]) *
			(Z[li][t] - rho * Z[li][t-1]);
		}
	    }
	    xpy[i-1] = x;
	}
    } else if (nwt) {
	/* weight the data */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!missing_masked(mask, t)) {
			x += Z[nwt][t] * Z[li][t] * Z[lj][t];
		    }
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}   
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		if (!missing_masked(mask, t)) {
		    x += Z[nwt][t] * Z[yno][t] * Z[li][t];
		}
	    }
	    xpy[i-1] = x;
	}
    } else {
	/* no quasi-differencing or weighting wanted */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!missing_masked(mask, t)) {
			x += Z[li][t] * Z[lj][t];
		    }
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		if (!missing_masked(mask, t)) {
		    x += Z[yno][t] * Z[li][t];
		}
	    }
	    xpy[i-1] = x;
	}
    }

    return 0; 
}

static int make_ess (MODEL *pmod, double **Z)
{
    int i, t, yno = pmod->list[1], l0 = pmod->list[0];
    int nwt = pmod->nwt;
    double yhat, resid;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (nwt && Z[nwt][t] == 0.0) {
	    continue;
	}
	if (model_missing(pmod, t)) {
	    continue;
	}
	yhat = 0.0;
	for (i=2; i<=l0; i++) {
	    yhat += pmod->coeff[i-2] * Z[pmod->list[i]][t];
	}
	resid = Z[yno][t] - yhat;
	if (nwt) {
	    resid *= sqrt(Z[nwt][t]);
	}
	pmod->ess += resid * resid;
    }

    return 0;
}

int check_for_effective_const (MODEL *pmod, const double *y)
{
    double x1 = 0.0, x2 = 0.0;
    double reldiff;
    int t, ret = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    x1 += pmod->yhat[t];
	    x2 += y[t];
	}
    }

    reldiff = fabs((x1 - x2) / x2);

    if (floateq(reldiff, 0.0)) {
	gretl_model_set_int(pmod, "effconst", 1);
	pmod->dfn -= 1;
	ret = 1;
    } else if (gretl_model_get_int(pmod, "effconst")) {
	gretl_model_set_int(pmod, "effconst", 0);
	pmod->dfn += 1;
    }

    return ret;
}

static void uncentered_r_squared (MODEL *pmod, const double *y)
{
    double y0 = y[pmod->t1];

    if (y0 > 0) {
	double tss = pmod->nobs * y0 * y0;

	pmod->rsq = 1 - (pmod->ess / tss);
	gretl_model_set_int(pmod, "uncentered", 1);
    }
}

static void compute_r_squared (MODEL *pmod, const double *y, int *ifc)
{
    pmod->rsq = 1.0 - (pmod->ess / pmod->tss);

    if (*ifc == 0) {
	*ifc = check_for_effective_const(pmod, y);
    }

    if (pmod->dfd > 0) {
	double den;

	if (*ifc) {
	    den = pmod->tss * pmod->dfd;
	    pmod->adjrsq = 1 - (pmod->ess * (pmod->nobs - 1) / den);
	} else {
	    int t;

	    den = 0.0;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!na(y[t])) {
		    den += y[t] * y[t];
		}
	    }
	    pmod->rsq = 1 - pmod->ess / den; /* NIST method */
	    pmod->adjrsq = 
		1.0 - ((1.0 - pmod->rsq) * (pmod->nobs - 1.0) / pmod->dfd);
	} 
    }

    if (pmod->rsq < 0.0) {
	pmod->rsq = 0.0;
    }
}

/*
  regress: takes xpx, the X'X matrix produced by form_xpxxpy(), and
  xpy (X'y), and computes ols estimates and associated statistics.

  n = total number of observations per series in data set
  ifc = 1 if constant is present in model, else = 0

  ess = error sum of squares
  sigma = standard error of regression
  fstt = F-statistic
  coeff = array of regression coefficients
  sderr = corresponding array of standard errors
*/

static void regress (MODEL *pmod, double *xpy, double **Z, 
		     int n, double rho)
{
    int v, yno = pmod->list[1];
    int ifc = pmod->ifc;
    double ysum, ypy, zz, rss = 0.0;
    double sgmasq = 0.0;
    double *diag = NULL;
    int i, err = 0;

    for (i=0; i<n; i++) {
	pmod->yhat[i] = pmod->yhat[i] = NADBL;
    }    

    ysum = xpy[0];
    ypy = xpy[pmod->ncoeff + 1];
    zz = ysum * ysum / pmod->nobs;
    pmod->tss = ypy - zz;

    /*  Cholesky-decompose X'X and find the coefficients */
    err = cholbeta(pmod->xpx, xpy, pmod->coeff, pmod->ncoeff, &rss);
    if (err) {
        pmod->errcode = err;
        return;
    }   
    
    if (rho != 0.0) {
	pmod->ess = ypy - rss;
    } else {
	make_ess(pmod, Z);
	rss = ypy - pmod->ess;
    }

    if (fabs(pmod->ess) < ESSZERO) {
	pmod->ess = 0.0;
    } else if (pmod->ess < 0.0) { 
	sprintf(gretl_errmsg, _("Error sum of squares (%g) is not > 0"),
		pmod->ess);
        return; 
    }

    if (pmod->dfd == 0) {
	pmod->sigma = 0.0;
	pmod->adjrsq = NADBL;
    } else {
	if (gretl_model_get_int(pmod, "no-df-corr")) {
	    sgmasq = pmod->ess / pmod->nobs;
	} else {
	    sgmasq = pmod->ess / pmod->dfd;
	}
	pmod->sigma = sqrt(sgmasq);
    }

    if (floatlt(pmod->tss, 0.0) || floateq(pmod->tss, 0.0)) {
	pmod->rsq = pmod->adjrsq = NADBL;
    } 

    hatvar(pmod, n, Z); 
    if (pmod->errcode) return;

    if (pmod->tss > 0.0) {
	compute_r_squared(pmod, Z[yno], &ifc);
    } else if (pmod->tss == 0.0) {
	uncentered_r_squared(pmod, Z[yno]);
    }

#if 0
    if (pmod->ifc && pmod->ncoeff == 1) {
        zz = 0.0;
        pmod->dfn = 1;
    }
#endif

    if (sgmasq <= 0.0 || pmod->dfd == 0 || pmod->dfn == 0) {
	pmod->fstt = NADBL;
    } else if (pmod->rsq == 1.0) {
	pmod->fstt = NADBL;
    } else {
	pmod->fstt = (rss - zz * ifc) / (sgmasq * pmod->dfn);
	if (pmod->fstt < 0.0) {
	    pmod->fstt = 0.0;
	}
    }

    diag = malloc(pmod->ncoeff * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    diaginv(pmod->xpx, xpy, diag, pmod->ncoeff);

    for (v=0; v<pmod->ncoeff; v++) {
	if (diag[v] >= 0.0) {
	    pmod->sderr[v] = pmod->sigma * sqrt(diag[v]);
	} else {
	    pmod->sderr[v] = 0.0;
	}
    }

    free(diag); 
}

/*
  cholbeta: does an in-place Cholesky decomposition of xpx (lower
  triangular matrix stacked in columns) and solves for the
  least-squares coefficient estimates.

  xpx = X'X on input and Cholesky decomposition on output
  xpy = the X'y vector on input and Cholesky-transformed t
        vector on output 
  coeff = array of estimated coefficients 
  nv = number of regression coefficients including the constant
  rss = location to receive regression sum of squares

  The number of floating-point operations is basically 3.5 * nv^2
  plus (nv^3) / 3.
*/

static int 
cholbeta (double *xpx, double *xpy, double *coeff, int nv, double *rss)
{
    int i, j, k, kk, l, jm1;
    double e, d, d1, d2, test, xx;

    if (xpx[0] <= 0.0) {
	fprintf(stderr, "%s %d: xpx <= 0.0\n", __FILE__, __LINE__);
	return E_NAN;
    }

    e = 1.0 / sqrt(xpx[0]);
    xpx[0] = e;
    xpy[1] *= e;
    for (i=1; i<nv; i++) {
	xpx[i] *= e;
    }

    kk = nv;

    for (j=2; j<=nv; j++) {
	/* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d1 += xx * xpy[l];
            d += xx * xx;
            k += nv-l;
        }
        d2 = xpx[kk] - d;
	test = d2 / xpx[kk];
        if (test < TINY) {
	    *rss = -1.0;
	    return E_SINGULAR;
        }
	if (test < SMALL) {
	    strcpy(gretl_msg, _("Warning: data matrix close to singularity!"));
	}
        e = 1 / sqrt(d2);
        xpx[kk] = e;
        xpy[j] = (xpy[j] - d1) * e;
        for (i=j+1; i<=nv; i++) {
	    /* off-diagonal elements */
            kk++;
            d = 0.0;
            k = j - 1;
            for (l=1; l<=jm1; l++) {
                d += xpx[k] * xpx[k-j+i];
                k += nv - l;
            }
            xpx[kk] = (xpx[kk] - d) * e;
        }
        kk++;
    }

    kk--;

    /* calculate regression sum of squares */
    d = 0.0;
    for (j=1; j<=nv; j++) {
	d += xpy[j] * xpy[j];
    }
    *rss = d;

    /* solve for the coefficients */
    for (j=0; j<nv-1; j++) {
	coeff[j] = 0.0;
    }

    coeff[nv-1] = xpy[nv] * xpx[kk];

    for (j=nv-1; j>=1; j--) {
	d = xpy[j];
	for (i=nv-1; i>=j; i--) {
	    kk--;
	    d -= coeff[i] * xpx[kk];
	}
	kk--;
	coeff[j-1] = d * xpx[kk];
    }

    for (j=0; j<nv; j++) {
	if (isnan(coeff[j])) {
	    fprintf(stderr, "%s %d: coeff %d is NaN\n", __FILE__, __LINE__, j);
	    return E_NAN;
	}
    }	

    return 0; 
}

/*
  diaginv: solves for the diagonal elements of the X'X inverse matrix.

  xpx = Cholesky-decomposed X'X matrix (input)
  xpy = X'y vector (input) used as work array
  diag = diagonal elements of X'X (output)
  nv = number of regression coefficients
*/

static void diaginv (double *xpx, double *xpy, double *diag, int nv)
{
    int kk, l, m, k, i, j;
    const int nxpx = nv * (nv + 1) / 2;
    double d, e;

    kk = 0;

    for (l=1; l<=nv-1; l++) {
        d = xpx[kk];
        xpy[l] = d;
        e = d * d;
        m = 0;
        if (l > 1) {
	    for (j=1; j<=l-1; j++) {
		m += nv - j;
	    }
	}
        for (i=l+1; i<=nv; i++) {
            d = 0.0;
            k = i + m - 1;
            for (j=l; j<=i-1; j++) {
                d += xpy[j] * xpx[k];
                k += nv - j;
            }
            d = -d * xpx[k];
            xpy[i] = d;
            e += d * d;
        }
        kk += nv + 1 - l;
        diag[l-1] = e;
    }

    diag[nv-1] = xpx[nxpx-1] * xpx[nxpx-1];
}

/**
 * makevcv:
 * @pmod: pointer to model.
 * @sigma: square root of error variance, or 1.0 to
 * produce just X'X^{-1}.
 *
 * Inverts the Cholesky-decomposed X'X matrix and computes the 
 * coefficient covariance matrix.
 * 
 * Returns: 0 on successful completion, non-zero code on error.
 */

int makevcv (MODEL *pmod, double sigma)
{
    int dec, mst, kk, i, j, kj, icnt, m, k, l = 0;
    const int nv = pmod->ncoeff;
    const int nxpx = (nv * nv + nv) / 2; 
    double d;

    if (pmod->vcv != NULL) {
	return 0;
    }

    if (pmod->xpx == NULL) {
	fprintf(stderr, "makevcv: pmod->xpx = NULL\n");
	return 1;
    }

    mst = nxpx;
    kk = nxpx - 1;

    pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nv; i++) {
	mst -= i;
	/* find diagonal element */
	d = pmod->xpx[kk];
	if (i > 0) {
	    for (j=kk+1; j<=kk+i; j++) {
		d -= pmod->xpx[j] * pmod->vcv[j];
	    }
	}
	pmod->vcv[kk] = d * pmod->xpx[kk];
	/* find off-diagonal elements indexed by kj */
	kj = kk;
	kk = kk - i - 2;
	if (i > nv - 2) {
	    continue;
	}
	for (j=i+1; j<nv; j++) {
	    icnt = i+1;
	    kj -= j;
	    d = 0.0;
	    m = mst + 1;
	    for (k=0; k<=j-1; k++) {
		if (icnt > 0) {
		    dec = 1;
		    icnt--;
		} else {
		    dec = k;
		}
		m -= dec;
		l = kj + i - k;
		d += pmod->vcv[m-1] * pmod->xpx[l];
	    }
	    pmod->vcv[kj] = -d * pmod->xpx[l-1];
	}
    }

    if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	sigma = 1.0;
    }

    if (sigma != 1.0) {
	double s2 = sigma * sigma;

	for (k=0; k<nxpx; k++) {
	    pmod->vcv[k] *= s2;
	}
    }

    return 0;
}

/**
 * dwstat:
 * @order: order of autoregression (usually 1).
 * @pmod: pointer to model.
 * @Z: data array.
 *
 * Computes the Durbin-Watson statistic for @pmod.
 * 
 * Returns: the D-W value, or #NADBL on error.
 */

double dwstat (int order, MODEL *pmod, const double **Z)
{
    double ut, u1;
    double num = 0.0;
    double den = 0.0;
    int t, t1;

    if (pmod->ess <= 0.0) {
	return NADBL;
    }

    t1 = pmod->t1 + order;

    if (pmod->nwt) {
	ut = pmod->uhat[t1 - 1];
	if (!na(ut)) {
	    den += ut * ut;
	}
    } else {
	den = pmod->ess;
    }

    for (t=t1; t<=pmod->t2; t++)  {
        ut = pmod->uhat[t];
        u1 = pmod->uhat[t-1];
        if (na(ut) || na(u1) ||
	    (pmod->nwt && (Z[pmod->nwt][t] == 0.0 || 
			   Z[pmod->nwt][t-1] == 0.0))) { 
	    continue;
	}
        num += (ut - u1) * (ut - u1);
	if (pmod->nwt) {
	    den += ut * ut;
	}
    }

    return num / den;
}

/* altrho: alternative calculation of rho */

static double altrho (int order, int t1, int t2, const double *uhat)
{
    double *ut, *u1;    
    int t, n, len = t2 - (t1 + order) + 1;
    double uht, uh1, rho;

    ut = malloc(len * sizeof *ut);
    if (ut == NULL) {
	return NADBL;
    }

    u1 = malloc(len * sizeof *u1);
    if (u1 == NULL) {
	free(ut);
	return NADBL;
    }

    n = 0;

    for (t=t1+order; t<=t2; t++) { 
        uht = uhat[t];
	uh1 = (t > 0)? uhat[t-1] : NADBL;
        if (!na(uht) && !na(uh1)) {
	    ut[n] = uht;
	    u1[n] = uh1;
	    n++;
	}
    }

    rho = gretl_corr(0, n - 1, ut, u1, NULL);

    free(ut);
    free(u1);

    return rho;
}

/**
 * rhohat:
 * @order: order of autoregression, usually 1.
 * @t1: start of sample range.
 * @t2: end of sample range.
 * @uhat: array of regression residuals.
 *
 * Computes the first order serial correlation coefficient
 * for @uhat, over the range @t1 to @t2.
 * 
 * Returns: the \hat{rho} value, or #NADBL on error.
 */

double rhohat (int order, int t1, int t2, const double *uhat)
{
    double ut, u1, uu = 0.0, xx = 0.0;
    double rho;
    int t;

    for (t=t1+order; t<=t2; t++) { 
        ut = uhat[t];
        u1 = uhat[t-1];
        if (na(ut) || na(u1)) {
	    continue;
	}
        uu += ut * u1;
        xx += u1 * u1;
    }

    if (floateq(xx, 0.0)) {
	return NADBL;
    }

    rho = uu / xx;

    if (rho > 1.0 || rho < -1.0) {
	rho = altrho(order, t1, t2, uhat);
    }

    return rho;
}

/* compute fitted values and residuals */

static int hatvar (MODEL *pmod, int n, double **Z)
{
    int xno, i, t;
    int yno = pmod->list[1];
    double x;

    for (t=0; t<n; t++) {
	pmod->yhat[t] = pmod->uhat[t] = NADBL;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}
	pmod->yhat[t] = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            xno = pmod->list[i+2];
	    x = Z[xno][t];
	    if (pmod->nwt) {
		x *= sqrt(Z[pmod->nwt][t]);
	    }
            pmod->yhat[t] += pmod->coeff[i] * x;
        }
	x = Z[yno][t];
	if (pmod->nwt) {
	    x *= sqrt(Z[pmod->nwt][t]);
	}
        pmod->uhat[t] = x - pmod->yhat[t];                
    }

    return 0;
}

static int hilu_plot (double *ssr, double *rho, int n)
{
    FILE *fp;
    int i;

    if (gnuplot_init(PLOT_REGULAR, &fp)) {
	return E_FOPEN; 
    }

    fputs("# hildreth-lu\n", fp);
    fputs("set xlabel 'rho'\n", fp);

    fprintf(fp, "set ylabel '%s'\n", I_("ESS"));

    fputs("set nokey\n", fp);
    fputs("set xrange [-1.0:1.0]\n", fp);
    fputs("plot '-' using 1:2 w impulses\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	fprintf(fp, "%g %g\n", rho[i], ssr[i]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    gnuplot_make_graph();

    return 0;
}

static double autores (MODEL *pmod, const double **Z, int ci)
{
    int t, v, t1 = pmod->t1;
    double x, num = 0.0, den = 0.0;
    double rhohat;

    if (ci == CORC || ci == HILU) {
	t1--;
    }

    for (t=t1; t<=pmod->t2; t++) {
	x = Z[pmod->list[1]][t];
	for (v=0; v<pmod->ncoeff; v++) {
	    x -= pmod->coeff[v] * Z[pmod->list[v+2]][t];
	}
	pmod->uhat[t] = x;
	if (t > t1) {
	    num += pmod->uhat[t] * pmod->uhat[t-1];
	    den += pmod->uhat[t-1] * pmod->uhat[t-1];
	}
    } 

    rhohat = num / den;

    return rhohat;
}

#define AR_DEBUG 0

/**
 * estimate_rho:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: %CORC for Cochrane-Orcutt, %HILU for Hildreth-Lu,
 *      %PWE for Prais-Winsten estimator.
 * @err: pointer for error code.
 * @opt: option flags: may include %OPT_B to suppress Cochrane-Orcutt
 *       fine-tuning of Hildreth-Lu results, %OPT_P to generate
 *       a gnuplot graph of the search in case @ci = %HILU.
 * @prn: gretl printing struct.
 *
 * Estimate the quasi-differencing coefficient for use with the
 * Cochrane-Orcutt, Hildreth-Lu or Prais-Winsten procedures for
 * handling first-order serial correlation.  Print a trace of the
 * search for rho.
 * 
 * Returns: rho estimate on successful completion, #NADBL on error.
 */

double estimate_rho (const int *list, double ***pZ, DATAINFO *pdinfo,
		     GretlCmdIndex ci, int *err, gretlopt opt, PRN *prn)
{
    double rho = 0.0, rho0 = 0.0, diff;
    double finalrho = 0.0, essmin = 1.0e8;
    double ess, ssr[199], rh[199]; 
    int iter, nn = 0;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int missv = 0, misst = 0;
    gretlopt lsqopt = OPT_A;
    int ascii = !(opt & OPT_P);
    MODEL corc_model;

    *gretl_errmsg = '\0';
    *err = 0;

    missv = adjust_t1t2(NULL, list, &pdinfo->t1, &pdinfo->t2, 
			pdinfo->n, (const double **) *pZ, &misst);
    if (missv) {
	sprintf(gretl_errmsg, _("Missing value encountered for "
				"variable %d, obs %d"), missv, misst);
	*err = E_DATA;
	goto bailout;
    }

    gretl_model_init(&corc_model);

    if (ci == PWE) {
	/* Prais-Winsten treatment of first observation */
	lsqopt |= OPT_P;
    } 

    if (ci == HILU) { 
	/* Do Hildreth-Lu first */
	for (rho = -0.990, iter = 0; rho < 1.0; rho += .01, iter++) {
	    clear_model(&corc_model);
	    corc_model = ar1_lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
	    if ((*err = corc_model.errcode)) {
		clear_model(&corc_model);
		goto bailout;
	    }
	    ess = corc_model.ess;
	    if (ascii) {
		char num[16];
		int chk;
		
		if (iter == 0) {
		    pprintf(prn, "\n RHO       %s      RHO       %s      "
			    "RHO       %s      RHO       %s     \n",
			    _("ESS"), _("ESS"), _("ESS"), _("ESS"));
		}
		sprintf(num, "%f", 100 * fabs(rho));
		chk = atoi(num);
		if (chk == 99 || chk % 10 == 0) {
		    ssr[nn] = ess;
		    rh[nn++] = rho;
		    pprintf(prn, "%5.2f %10.4g", rho, ess);
		    if (nn % 4 == 0) {
			pputc(prn, '\n');
		    } else {
			bufspace(3, prn);
		    }
		} 
	    } else {
		ssr[nn] = ess;
		rh[nn++] = rho;
	    }
	    if (iter == 0 || ess < essmin) {
		essmin = ess;
		finalrho = rho;
	    }
	} /* end of basic iteration */
	
	if (finalrho > 0.989) {
	    /* try exploring this funny region? */
	    for (rho = 0.99; rho <= 0.999; rho += .001) {
		clear_model(&corc_model);
		corc_model = ar1_lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
		if ((*err = corc_model.errcode)) {
		    clear_model(&corc_model);
		    goto bailout;
		}
		ess = corc_model.ess;
		if (ess < essmin) {
		    essmin = ess;
		    finalrho = rho;
		}
	    }
	}

	if (finalrho > 0.9989) {
	    /* this even funnier one? */
	    for (rho = 0.9991; rho <= 0.9999; rho += .0001) {
		clear_model(&corc_model);
		corc_model = ar1_lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
		if ((*err = corc_model.errcode)) {
		    clear_model(&corc_model);
		    goto bailout;
		}
		ess = corc_model.ess;
		if (ess < essmin) {
		    essmin = ess;
		    finalrho = rho;
		}
	    }
	}

	rho0 = rho = finalrho;
	pprintf(prn, _("\n\nESS is minimum for rho = %g\n\n"), rho);
	if (ascii) {
	    graphyzx(NULL, ssr, NULL, rh, nn, "ESS", "RHO", NULL, 0, prn); 
	    pputs(prn, "\n");
	} else {
	    hilu_plot(ssr, rh, nn);
	}
    } else { 
	/* Go straight to Cochrane-Orcutt (or Prais-Winsten) */
	corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A);
	if (!corc_model.errcode && corc_model.dfd == 0) {
	    corc_model.errcode = E_DF;
	}
	if ((*err = corc_model.errcode)) {
	    clear_model(&corc_model);
	    goto bailout;
	}
	rho0 = rho = corc_model.rho;
    }

    if (na(rho)) {
	*err = E_NOCONV;
	clear_model(&corc_model);
	goto bailout;
    }

    if (ci != HILU || !(opt & OPT_B)) {

	if (ci == HILU) {
	    pputs(prn, _("\nFine-tune rho using the CORC procedure...\n\n"));
	} else {
	    pputs(prn, _("\nPerforming iterative calculation of rho...\n\n"));
	}

	pputs(prn, _("                 ITER       RHO        ESS"));
	pputc(prn, '\n');

	iter = 0;
	diff = 1.0;

	while (diff > 0.001) {
	    pprintf(prn, "          %10d %12.5f", ++iter, rho);
	    clear_model(&corc_model);
	    corc_model = ar1_lsq(list, pZ, pdinfo, OLS, lsqopt, rho);
#if AR_DEBUG
	    fprintf(stderr, "corc_model: t1=%d, first two uhats: %g, %g\n",
		    corc_model.t1, 
		    corc_model.uhat[corc_model.t1],
		    corc_model.uhat[corc_model.t1+1]);
#endif
	    if ((*err = corc_model.errcode)) {
		clear_model(&corc_model);
		goto bailout;
	    }
	    pprintf(prn, "   %g\n", corc_model.ess);

	    rho = autores(&corc_model, (const double **) *pZ, ci);

#if AR_DEBUG
	    pputs(prn, "CORC model (using rho-transformed data)\n");
	    printmodel(&corc_model, pdinfo, OPT_NONE, prn);
	    pprintf(prn, "autores gives rho = %g\n", rho);
#endif

	    if (rho > .99999 || rho < -.99999) {
		*err = E_NOCONV;
		clear_model(&corc_model);
		goto bailout;
	    }

	    diff = (rho > rho0) ? rho - rho0 : rho0 - rho;
	    rho0 = rho;
	    if (iter == 30) break;
	}

	pprintf(prn, _("                final %11.5f\n\n"), rho);
    }

    clear_model(&corc_model);

 bailout:

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    if (*err) {
	rho = NADBL;
    }

    return rho;
}

/**
 * augment_regression_list:
 * @orig: list giving original regression specification.
 * @aux: either %AUX_SQ, %AUX_LOG or %AUX_WHITE.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 *
 * Augment the regression list @orig with auxiliary terms.  If @aux 
 * is %AUX_SQ add the squares of the original regressors; if @aux
 * is %AUX_WHITE add squares and cross-products, or if @aux is
 * %AUX_LOG add the natural logs of the original regressors.
 * If theye are not already present, these variables are added
 * to the data array.
 * 
 * Returns: the augmented list, or NULL on failure.
 */

int *augment_regression_list (const int *orig, int aux, 
			      double ***pZ, DATAINFO *pdinfo)
{
    int *list;
    int listlen;
    int cnum = 0;
    int i, k;

    if (aux == AUX_WHITE) {
	int cpos = gretl_list_const_pos(orig, 2, (const double **) *pZ, 
					pdinfo);
	int nt, trv = orig[0] - 1;

	if (cpos > 0) {
	    trv--;
	    cnum = orig[cpos];
	}
	nt = (trv * trv + trv) / 2;
	listlen = orig[0] + nt + 1;
    } else {
	listlen = 2 * orig[0];
    }

    list = malloc(listlen * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    /* transcribe original list */
    for (i=0; i<=orig[0]; i++) {
	list[i] = orig[i];
    }

    /* add squares, cross-products or logs of independent vars */
    k = list[0];
    for (i=2; i<=orig[0]; i++) {
	int vnew, vi = orig[i];

	if (vi == 0) {
	    continue;
	}

	if (aux == AUX_SQ || aux == AUX_WHITE) {
	    vnew = xpxgenr(vi, vi, pZ, pdinfo);
	    if (vnew > 0) {
		list[++k] = vnew;
	    }
	    if (aux == AUX_WHITE) {
		int j, vj;

		for (j=i+1; j<=orig[0]; j++) {
		    vj = orig[j];
		    if (vj == cnum) {
			continue;
		    }
		    vnew = xpxgenr(vi, vj, pZ, pdinfo);
		    if (vnew > 0) list[++k] = vnew;
		}
	    }
	} else if (aux == AUX_LOG) {
	    vnew = loggenr(vi, pZ, pdinfo);
	    if (vnew > 0) {
		list[++k] = vnew;
	    }
	}
    }

    list[0] = k;

    return list;
}

/* get_hsk_weights: take the residuals from the model pmod, square them
   and take logs; find the fitted values for this series using an
   auxiliary regression including the original independent variables
   and their squares; transform the fitted values by exponentiating
   and taking the square root; and add the resulting series to the
   data set
*/

static int get_hsk_weights (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    int oldv = pdinfo->v;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int *list = NULL;
    int err = 0, shrink = 0;
    double xx;
    MODEL aux;

    /* allocate space for an additional variable */
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    /* add transformed pmod residuals to data set */
    for (t=0; t<pdinfo->n; t++) {
	if (na(pmod->uhat[t])) {
	    (*pZ)[oldv][t] = NADBL;
	} else {
	    xx = pmod->uhat[t];
	    (*pZ)[oldv][t] = log(xx * xx);
	}
    }

    /* build regression list, adding the squares of the original
       independent vars */
    list = augment_regression_list(pmod->list, AUX_SQ, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    list[1] = oldv; /* the newly added uhat-squared */

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    aux = lsq(list, pZ, pdinfo, OLS, OPT_A);
    err = aux.errcode;
    if (err) {
	shrink = pdinfo->v - oldv;
    } else {
	/* write into the data set the required weights */
	for (t=aux.t1; t<=aux.t2; t++) {
	    if (na(aux.yhat[t])) {
		(*pZ)[oldv][t] = NADBL;
	    } else {
		xx = aux.yhat[t];
		(*pZ)[oldv][t] = 1.0 / exp(xx);
	    }
	}
	shrink = pdinfo->v - oldv - 1;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    clear_model(&aux);

    if (shrink > 0) {
	dataset_drop_last_variables(shrink, pZ, pdinfo);
    }

    free(list);

    return err;
}

/**
 * hsk_func:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using a correction for
 * heteroskedasticity.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL hsk_func (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int i, err;
    int orig_nvar = pdinfo->v;
    int *hsklist;
    MODEL hsk;

    *gretl_errmsg = '\0';

    /* run initial OLS */
    hsk = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (hsk.errcode) {
	return hsk;
    }

    /* use the residuals from the initial OLS to form weights */
    err = get_hsk_weights(&hsk, pZ, pdinfo);
    if (err) {
	hsk.errcode = err;
	return hsk;
    }

    /* allocate regression list for weighted least squares */
    hsklist = gretl_list_new(list[0] + 1);
    if (hsklist == NULL) {
	hsk.errcode = E_ALLOC;
	return hsk;
    }

    /* the last variable in the dataset will be the weight var */
    hsklist[1] = pdinfo->v - 1;

    /* put the original dependent variable in at position 2 */
    hsklist[2] = list[1];

    /* add the original independent vars into the WLS list */
    for (i=3; i<=hsklist[0]; i++) {
	hsklist[i] = list[i-1];
    }

    clear_model(&hsk);
    hsk = lsq(hsklist, pZ, pdinfo, WLS, OPT_NONE);
    hsk.ci = HSK;

    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    free(hsklist);

    return hsk;
}

static double **allocate_hccm_p (int k, int n)
{
    double **p = malloc(k * sizeof *p);
    int i;

    if (p == NULL) return NULL;

    for (i=0; i<k; i++) {
	p[i] = malloc(n * sizeof **p);
	if (p[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) {
		free(p[j]);
	    }
	    free(p);
	    p = NULL;
	    break;
	}
    }

    return p;
}

static void free_hccm_p (double **p, int m)
{
    if (p != NULL) {
	int i;

	for (i=0; i<m; i++) {
	    free(p[i]);
	}
	free(p);
    }
}

static int jackknife_vcv (MODEL *pmod, const double **Z)
{
    double *st = NULL, *ustar = NULL;
    double **p = NULL;
    int nobs, tp, nc = 0;
    int i, j, k, t;
    int t1, t2;
    double xx;
    int err = 0;

    *gretl_errmsg = '\0';

    t1 = pmod->t1;
    t2 = pmod->t2;
    nobs = pmod->nobs;
    nc = pmod->ncoeff;

    st = malloc(nc * sizeof *st);
    ustar = malloc(nobs * sizeof *ustar);
    p = allocate_hccm_p(nc, nobs);

    if (st == NULL || p == NULL || ustar == NULL) {
	err = E_ALLOC;
	goto bailout;
    }  

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
	pmod->vcv = NULL;
    }

    pmod->ci = HCCM;

    if (makevcv(pmod, 1.0)) {
	err = E_ALLOC;
	goto bailout;
    }

    /* form elements of (X'X)^{-1}X' */

    for (i=0; i<nc; i++) {
	tp = 0;
	for (t=t1; t<=t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    xx = 0.0;
	    for (j=0; j<nc; j++) {
		if (i <= j) {
		    k = ijton(i, j, nc);
		} else {
		    k = ijton(j, i, nc);
		}
		xx += pmod->vcv[k] * Z[pmod->list[j+2]][t];
	    }
	    p[i][tp++] = xx;
	}
    }

    tp = 0;
    for (t=t1; t<=t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}	
	xx = 0.0;
	for (i=0; i<nc; i++) {
	    xx += Z[pmod->list[i+2]][t] * p[i][tp];
	}
	if (floateq(xx, 1.0)) {
	    xx = 0.0;
	}
	ustar[tp++] = pmod->uhat[t] / (1.0 - xx);
    }

    for (i=0; i<nc; i++) {
	xx = 0.0;
	for (t=0; t<nobs; t++) {
	    xx += p[i][t] * ustar[t];
	}
	st[i] = xx;
    }

    for (t=0; t<nobs; t++) {
	for (i=0; i<nc; i++) {
	    p[i][t] *= ustar[t];
	}
    }

    /* MacKinnon and White, 1985, equation (13) */

    k = 0;
    for (i=0; i<nc; i++) {
	for (j=i; j<nc; j++) {
	    xx = 0.0;
	    for (t=0; t<nobs; t++) {
		xx += p[i][t] * p[j][t];
	    }
	    xx -= st[i] * st[j] / nobs;
	    /* MacKinnon and White: "It is tempting to omit the factor
	       (n - 1) / n from HC3" (1985, p. 309).  Here we leave it in
	       place, as in their simulations.
	    */
	    xx *= (nobs - 1.0) / nobs;
	    if (i == j) {
		pmod->sderr[i] = sqrt(xx);
	    }
	    pmod->vcv[k++] = xx;
	}
    }

    /* substitute robust F stat */
    if (pmod->dfd > 0 && pmod->dfn > 1) {
	pmod->fstt = robust_omit_F(NULL, pmod);
    }

    gretl_model_set_int(pmod, "robust", 1);
    gretl_model_set_int(pmod, "hc", 1);
    gretl_model_set_int(pmod, "hc_version", 4);

 bailout:

    pmod->ci = OLS;

    free(st);
    free(ustar);
    free_hccm_p(p, nc);

    return err;
}

/**
 * whites_test:
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save results to model.
 * @prn: gretl printing struct.
 *
 * Runs White's test for heteroskedasticity on the given model.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int whites_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    int lo, ncoeff, yno, t;
    int v = pdinfo->v;
    int *list = NULL;
    double zz;
    MODEL white;
    int err = 0;

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == LOGISTIC) { 
	return E_NOTIMP;
    }

    if ((err = list_members_replaced(pmod->list, pdinfo, pmod->ID))) {
	return err;
    }

    gretl_model_init(&white);

    lo = pmod->list[0];
    yno = pmod->list[1];
    ncoeff = pmod->list[0] - 1;

    /* make space in data set */
    if (dataset_add_series(1, pZ, pdinfo)) {
	err = E_ALLOC;
    }

    if (!err) {
	/* get residuals, square and add to data set */
	for (t=0; t<pdinfo->n; t++) {
	    zz = pmod->uhat[t];
	    if (na(zz)) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = zz * zz;
	    }
	}
	strcpy(pdinfo->varname[v], "uhatsq");
    }

    if (!err) {
	/* build aux regression list, adding squares and
	   cross-products of the original independent vars */
	list = augment_regression_list(pmod->list, AUX_WHITE, pZ, pdinfo);
	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    list[1] = v; /* the newly added uhat-squared */
	}
    }

    if (!err) {
	/* run auxiliary regression */
	white = lsq(list, pZ, pdinfo, OLS, OPT_A);
	err = white.errcode;
    }

    if (!err) {
	double TR2, pval;

	white.aux = AUX_WHITE;
	printmodel(&white, pdinfo, OPT_NONE, prn);

	TR2 = white.rsq * white.nobs;
	pval = chisq_cdf_comp(TR2, white.ncoeff - 1);

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_WHITES);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_TR2);
		model_test_set_dfn(test, white.ncoeff - 1);
		model_test_set_value(test, TR2);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	  
	}

	record_test_result(TR2, pval, _("White's"));
    }

    clear_model(&white);

    dataset_drop_last_variables(pdinfo->v - v, pZ, pdinfo);

    free(list);

    return err;
}

static int ar_list_max (const int *list) 
{
    int i, lmax = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] > lmax) {
	    lmax = list[i];
	}
    }

    return lmax;
}

/**
 * ar_func:
 * @list: list of lags plus dependent variable and list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list using the generalized 
 * Cochrane-Orcutt procedure for autoregressive errors.
 * 
 * Returns: #MODEL struct containing the results.
 */

MODEL ar_func (const int *list, double ***pZ, 
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    double diff, ess, tss, xx;
    int i, j, t, t1, t2, vc, yno, ryno, iter;
    int err, lag, maxlag, v = pdinfo->v;
    int *arlist = NULL, *rholist = NULL;
    int *reglist = NULL, *reglist2 = NULL;
    int pos, cpos;
    MODEL ar, rhomod;

    *gretl_errmsg = '\0';

    gretl_model_init(&ar);
    gretl_model_init(&rhomod);

    pos = gretl_list_separator_position(list);

    arlist = malloc(pos * sizeof *arlist);
    reglist = malloc((list[0] - pos + 2) * sizeof *reglist);
    reglist2 = malloc((list[0] - pos + 2) * sizeof *reglist2);
    rholist = malloc((pos + 2) * sizeof *rholist);

    if (arlist == NULL || reglist == NULL || reglist2 == NULL ||
	rholist == NULL) {
	ar.errcode = E_ALLOC;
	goto bailout;
    }
    
    arlist[0] = pos - 1;
    for (i=1; i<pos; i++) {
	arlist[i] = list[i];
    }

    reglist2[0] = reglist[0] = list[0] - pos;
    for (i=1; i<=reglist[0]; i++) {
	reglist[i] = list[i + pos];
    }

    rholist[0] = arlist[0] + 1;
    maxlag = ar_list_max(arlist);

    cpos = reglist_check_for_const(reglist, (const double **) *pZ, pdinfo);

    /* special case: ar 1 ; ... => use CORC */
    if (arlist[0] == 1 && arlist[1] == 1) {
	xx = estimate_rho(reglist, pZ, pdinfo, CORC, &err, OPT_NONE, prn);
	if (err) {
	    ar.errcode = err;
	} else {
	    ar = ar1_lsq(reglist, pZ, pdinfo, CORC, OPT_NONE, xx);
	    printmodel(&ar, pdinfo, opt, prn); 
	}
	goto bailout;
    }

    /* first pass: estimate model via OLS: use OPT_M to generate an
       error in case of missing values within sample range 
    */
    ar = lsq(reglist, pZ, pdinfo, OLS, OPT_A | OPT_M);
    if (ar.errcode) {
	goto bailout;
    }

    /* allocate space for the uhat terms and transformed data */
    if (dataset_add_series(arlist[0] + 1 + reglist[0], pZ, pdinfo)) {
	ar.errcode = E_ALLOC;
	goto bailout;
    }

    yno = reglist[1];
    t1 = ar.t1; t2 = ar.t2;
    rholist[1] = v;

    pprintf(prn, "%s\n\n", _("Generalized Cochrane-Orcutt estimation"));
    bufspace(17, prn);
    /* xgettext:no-c-format */
    pputs(prn, _("ITER             ESS           % CHANGE"));
    pputs(prn, "\n\n");

    /* now loop while ess is changing */
    diff = 1.0e6;
    ess = 0.0;
    for (iter = 1; iter <= 20 && diff > 0.005; iter++) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t > t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		/* special computation of uhat */
		xx = (*pZ)[yno][t];
		for (j=0; j<reglist[0]-1; j++) {
		    xx -= ar.coeff[j] * (*pZ)[reglist[j+2]][t];
		}
		(*pZ)[v][t] = xx;
	    }
	}		
	for (i=1; i<=arlist[0]; i++) {
	    lag = arlist[i];
	    rholist[1+i] = v + i;
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 + lag || t > t2) {
		    (*pZ)[v+i][t] = NADBL;
		} else {
		    (*pZ)[v+i][t] = (*pZ)[v][t-lag];
		}
	    }
	}

	/* now estimate the rho terms */
	if (iter > 1) {
	    clear_model(&rhomod);
	}
	rhomod = lsq(rholist, pZ, pdinfo, OLS, OPT_A);

	/* and rho-transform the data */
	ryno = vc = v + i;
	for (i=1; i<=reglist[0]; i++) {
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 + maxlag || t > t2) {
		    (*pZ)[vc][t] = NADBL;
		} else {
		    xx = (*pZ)[reglist[i]][t];
		    for (j=1; j<=arlist[0]; j++) {
			lag = arlist[j];
			xx -= rhomod.coeff[j-1] * (*pZ)[reglist[i]][t-lag];
		    }
		    (*pZ)[vc][t] = xx;
		}
	    }
	    reglist2[i] = vc++;
	}

	/* estimate the transformed model */
	clear_model(&ar);
	ar = lsq(reglist2, pZ, pdinfo, OLS, OPT_A);

        if (iter > 1) {
	    diff = 100 * (ar.ess - ess) / ess;
	}

        if (diff < 0.0) {
	    diff = -diff;
	}

	ess = ar.ess;
	pprintf(prn, "%16c%3d %20f ", ' ', iter, ess);

	if (iter > 1) {
	    pprintf(prn, "%13.3f\n", diff);
	} else {
	    pprintf(prn, "%*s\n", UTF_WIDTH(_("undefined"), 15), 
		     _("undefined")); 
	}
    } /* end "ess changing" loop */

    for (i=0; i<=reglist[0]; i++) {
	ar.list[i] = reglist[i];
    }
    if (cpos > 0) {
	ar.ifc = 1;
    }
    if (ar.ifc) {
	if (!gretl_model_get_int(&ar, "effconst")) {
	    ar.dfn -= 1;
	}
    }
    ar.ci = AR;

    /* special computation of fitted values */
    for (t=t1; t<=t2; t++) {
	xx = 0.0;
	for (j=2; j<=reglist[0]; j++) { 
	    xx += ar.coeff[j-2] * (*pZ)[reglist[j]][t];
	}
	ar.uhat[t] = (*pZ)[yno][t] - xx;
	for (j=1; j<=arlist[0]; j++) {
	    if (t - t1 >= arlist[j]) {
		xx += rhomod.coeff[j-1] * ar.uhat[t - arlist[j]];
	    }
	}
	ar.yhat[t] = xx;
    }

    for (t=t1; t<=t2; t++) { 
	ar.uhat[t] = (*pZ)[yno][t] - ar.yhat[t];
    }

    ar.rsq = gretl_corr_rsq(ar.t1, ar.t2, (*pZ)[reglist[1]], ar.yhat);
    ar.adjrsq = 1.0 - ((1.0 - ar.rsq) * (ar.nobs - 1.0) / ar.dfd);

    /* special computation of TSS */
    xx = gretl_mean(ar.t1, ar.t2, (*pZ)[ryno]);
    tss = 0.0;
    for (t=ar.t1; t<=ar.t2; t++) {
	tss += ((*pZ)[ryno][t] - xx) * ((*pZ)[ryno][t] - xx);
    }
    ar.fstt = ar.dfd * (tss - ar.ess) / (ar.dfn * ar.ess);
    ls_criteria(&ar);
    ar.dw = dwstat(maxlag, &ar, (const double **) *pZ);
    ar.rho = rhohat(maxlag, ar.t1, ar.t2, ar.uhat);

    dataset_drop_last_variables(arlist[0] + 1 + reglist[0], pZ, pdinfo);

    if (gretl_model_add_arinfo(&ar, maxlag)) {
	ar.errcode = E_ALLOC;
    } else {
	for (i=0; i<=arlist[0]; i++) { 
	    ar.arinfo->arlist[i] = arlist[i];
	    if (i >= 1) {
		ar.arinfo->rho[i-1] = rhomod.coeff[i-1];
		ar.arinfo->sderr[i-1] = rhomod.sderr[i-1];
	    }
	}
    }
    clear_model(&rhomod);

    if (!ar.errcode) {
	set_model_id(&ar);
    }  

 bailout:

    free(reglist);
    free(reglist2);
    free(rholist);
    free(arlist);

    return ar;
}

/* From 2 to end of list, omits variables with all zero observations
   and re-packs the rest of them */

static void omitzero (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    int v, lv, offset, dropmsg = 0;
    double xx = 0.0;
    char vnamebit[20];

    offset = (pmod->ci == WLS)? 3 : 2;

    for (v=offset; v<=pmod->list[0]; v++) {
        lv = pmod->list[v];
        if (gretl_iszero(pmod->t1, pmod->t2, Z[lv])) {
	    gretl_list_delete_at_pos(pmod->list, v);
	    if (pdinfo->varname[lv][0] != 0) {
		sprintf(vnamebit, "%s ", pdinfo->varname[lv]);
		strcat(gretl_msg, vnamebit);
		dropmsg = 1;
		v--;
	    }
	}
    }

    if (pmod->nwt) {
	int t, wtzero;

	for (v=offset; v<=pmod->list[0]; v++) {
	    lv = pmod->list[v];
	    wtzero = 1;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		xx = Z[lv][t] * Z[pmod->nwt][t];
		if (floatneq(xx, 0.0)) {
		    wtzero = 0;
		    break;
		}
	    }
	    if (wtzero) {
		gretl_list_delete_at_pos(pmod->list, v);
		sprintf(vnamebit, "%s ", pdinfo->varname[lv]);
		strcat(gretl_msg, vnamebit);
		dropmsg = 1;
		v--;
	    }
	}
    }

    if (dropmsg) {
	strcat(gretl_msg, _("omitted because all obs are zero."));
    }
}

static int depvar_zero (int t1, int t2, int yno, int nwt,
			const double **Z)
{
    double y;
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
	y = Z[yno][t];
	if (na(y)) {
	    continue;
	}
	if (nwt) {
	    y *= Z[nwt][t];
	}
	if (y != 0.0) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/* lagdepvar: attempt to detect presence of a lagged dependent
   variable among the regressors -- if found, return the position of
   this lagged var in the list; otherwise return 0
*/

static int 
lagdepvar (const int *list, const double **Z, const DATAINFO *pdinfo) 
{
    char depvar[VNAMELEN], othervar[VNAMELEN];
    char *p;
    int i, t, ret = 0;

    strcpy(depvar, pdinfo->varname[list[1]]);

    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    break;
	}
	strcpy(othervar, pdinfo->varname[list[i]]);
	p = strrchr(othervar, '_');
	if (p != NULL && isdigit(*(p + 1))) {
	    /* looks like a lag */
	    size_t len = strlen(othervar) - strlen(p);

	    if (!strncmp(depvar, othervar, len)) {
		int gotlag = 1;

		/* strong candidate for lagged depvar, but make sure */
		for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) {
		    if (Z[list[1]][t-1] != Z[list[i]][t]) {
			gotlag = 0;
			break;
		    }
		}
		if (gotlag) {
		    ret = i;
		    break;
		}
	    }
	}
    } 

    return ret;
}

/* if 'full' is non-zero, do the whole thing (print the ARCH test
   model, re-estimate if p-value is < .10); otherwise just do
   the ARCH test itself and print the test result
*/

static MODEL 
real_arch_test (MODEL *pmod, int order, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn, int full)
{
    MODEL archmod;
    int *wlist = NULL, *arlist = NULL;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int oldv = pdinfo->v;
    int i, t, nwt, nv, n = pdinfo->n;
    double LM, xx;
    int err = 0;

    *gretl_errmsg = '\0';

    gretl_model_init(&archmod);

    if (order == 0) {
	/* use data frequency as default lag order */
	order = pdinfo->pd;
    }

    if (order < 1 || order > T - pmod->list[0]) {
	archmod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, _("Invalid lag order for arch (%d)"), order);
	err = 1;
    }

    if (!err) {
	/* allocate workspace */
	if (dataset_add_series(order + 1, pZ, pdinfo) || 
	    (arlist = malloc((order + 3) * sizeof *arlist)) == NULL) {
	    err = archmod.errcode = E_ALLOC;
	}
    }

    if (!err) {
	/* start list for aux regression */
	arlist[0] = 2 + order;
	arlist[1] = pdinfo->v - order - 1;
	arlist[2] = 0;

	/* run OLS and get squared residuals */
	archmod = lsq(pmod->list, pZ, pdinfo, OLS, OPT_A | OPT_M);
	err = archmod.errcode;
    }

    if (!err) {
	nv = pdinfo->v - order - 1;
	strcpy(pdinfo->varname[nv], "utsq");
	for (t=0; t<n; t++) {
	    (*pZ)[nv][t] = NADBL;
	}
	for (t=archmod.t1; t<=archmod.t2; t++) {
	    xx = archmod.uhat[t];
	    (*pZ)[nv][t] = xx * xx;
	}
	/* also lags of squared resids */
	for (i=1; i<=order; i++) {
	    nv =  pdinfo->v - order + i - 1;
	    arlist[i+2] = nv;
	    sprintf(pdinfo->varname[nv], "utsq_%d", i);
	    for (t=0; t<n; t++) {
		(*pZ)[nv][t] = NADBL;
	    }
	    for (t=archmod.t1+i; t<=archmod.t2; t++) {
		(*pZ)[nv][t] = (*pZ)[arlist[1]][t-i];
	    }
	}

	/* run aux. regression */
	clear_model(&archmod);
	archmod = lsq(arlist, pZ, pdinfo, OLS, OPT_A);
	err = archmod.errcode;
    }

    if (!err) {
	archmod.aux = AUX_ARCH;
	gretl_model_set_int(&archmod, "arch_order", order);
	LM = archmod.nobs * archmod.rsq;
	xx = chisq_cdf_comp(LM, order);

	if (full) {
	    printmodel(&archmod, pdinfo, OPT_NONE, prn);
	    pprintf(prn, _("No of obs. = %d, unadjusted R^2 = %f\n"),
		    archmod.nobs, archmod.rsq);
	}

	if ((opt & OPT_S) || (opt & OPT_P)) {
	    ModelTest *test = model_test_new(GRETL_TEST_ARCH);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_TR2);
		model_test_set_order(test, order);
		model_test_set_dfn(test, order);
		model_test_set_value(test, LM);
		model_test_set_pvalue(test, xx);
		if (opt & OPT_S) {
		    maybe_add_test_to_model(pmod, test);
		} else {
		    gretl_model_test_print_direct(test, prn);
		    model_test_free(test);
		}
	    }	    
	}

	record_test_result(LM, xx, "ARCH");

	if (!full) {
	    goto arch_test_exit;
	}

	pprintf(prn, _("LM test statistic (%f) is distributed as Chi-square "
		"(%d)\nArea to the right of LM = %f  "), LM, order, xx);

	if (xx > 0.1) {
	    pprintf(prn, "\n%s.\n%s.\n",
		    _("ARCH effect is insignificant at the 10 percent level"),
		    _("Weighted estimation not done"));
	} else {
	    pprintf(prn, "\n%s.\n",
		    _("ARCH effect is significant at the 10 percent level"));
	    /* do weighted estimation */
	    wlist = gretl_list_new(pmod->list[0] + 1);
	    if (wlist == NULL) {
		archmod.errcode = E_ALLOC;
	    } else {
		nwt = wlist[1] = pdinfo->v - 1; /* weight var */
		for (i=2; i<=wlist[0]; i++) {
		    wlist[i] = pmod->list[i-1];
		}
		nv = pdinfo->v - order - 1;
		for (t=archmod.t1; t<=archmod.t2; t++) {
		    xx = archmod.yhat[t];
		    if (xx <= 0.0) {
			xx = (*pZ)[nv][t];
		    }
		    (*pZ)[nwt][t] = 1.0 / xx; /* FIXME is this right? */
		}

		strcpy(pdinfo->varname[nwt], "1/sigma");

		clear_model(&archmod);
		archmod = lsq(wlist, pZ, pdinfo, WLS, OPT_NONE);

		archmod.ci = ARCH;
		gretl_model_set_int(&archmod, "arch_order", order);
		printmodel(&archmod, pdinfo, opt, prn);
	    }
	}
    }

 arch_test_exit:

    if (arlist != NULL) free(arlist);
    if (wlist != NULL) free(wlist);

    dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo); 

    return archmod;
}

/**
 * arch_test:
 * @pmod: model to be tested.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain %OPT_O to print covariance matrix, %OPT_S
 *       to save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests @pmod for Auto-Regressive Conditional Heteroskedasticity.  
 * If this effect is significant, re-restimates the model using 
 * weighted least squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch_test (MODEL *pmod, int order, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    return real_arch_test(pmod, order, pZ, pdinfo, opt, prn, 1);
}

/**
 * arch_test_simple:
 * @pmod: model to be tested.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if %OPT_S, the test is saved to @pmod and not printed,
 * otherwise it is printed to @prn and discarded.
 * @prn: gretl printing struct.
 *
 * Tests @pmod for Auto-Regressive Conditional Heteroskedasticity.  
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int arch_test_simple (MODEL *pmod, int order, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    MODEL amod;
    int err;

    if (!(opt & OPT_S)) {
	/* if not saving, then print test */
	opt = OPT_P;
    }

    amod = real_arch_test(pmod, order, pZ, pdinfo, opt, prn, 0);
    err = amod.errcode;
    clear_model(&amod);

    return err;
}

/**
 * arch_model:
 * @list: dependent variable plus list of regressors.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_O to print covariance matrix. (?)
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list via OLS, and test for Auto-
 * Regressive Conditional Heteroskedasticity.  If the latter is
 * significant, re-restimate the model using weighted least
 * squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch_model (const int *list, int order, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    MODEL lmod, amod;

    gretl_model_init(&lmod);
    lmod.list = gretl_list_copy(list);
    if (lmod.list == NULL) {
	lmod.errcode = E_ALLOC;
	return lmod;
    } 

    /* FIXME: vcv option? */
    amod = real_arch_test(&lmod, order, pZ, pdinfo, opt, prn, 1);

    free(lmod.list);

    return amod;
}

/**
 * lad:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using the method of Least
 * Absolute Deviation (LAD).
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lad (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    MODEL lad_model;
    void *handle;
    int (*lad_driver) (MODEL *, double **, DATAINFO *);

    /* run an initial OLS to "set the model up" and check for errors.
       the lad_driver function will overwrite the coefficients etc.
    */

    lad_model = lsq(list, pZ, pdinfo, OLS, OPT_A);

    if (lad_model.errcode) {
        return lad_model;
    }

    lad_driver = get_plugin_function("lad_driver", &handle);

    if (lad_driver == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	lad_model.errcode = E_FOPEN;
	return lad_model;
    }

    (*lad_driver) (&lad_model, *pZ, pdinfo);
    close_plugin(handle);

    set_model_id(&lad_model);

    return lad_model;
}

/**
 * arma:
 * @list: dependent variable, AR and MA orders, and any exogenous
 * regressors.
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @opt: options: may include %OPT_S to suppress intercept, %OPT_V
 * for verbose results, %OPT_X to use X-12-ARIMA, %OPT_C to put
 * X-12-ARIMA into conditional maximum-likelihood mode.
 * @PRN: for printing details of iterations (or %NULL). 
 *
 * Calculate ARMA estimates, using either native gretl code or
 * by invoking X-12-ARIMA.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arma (const int *list, const double **Z, const DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn)
{
    MODEL armod;
    void *handle;
    MODEL (*arma_func) (const int *, const double **, const DATAINFO *, 
			gretlopt, PRN *);

    *gretl_errmsg = '\0';

    if (opt & OPT_X && (pdinfo->t2 - pdinfo->t1) > 719) {
	strcpy(gretl_errmsg, _("X-12-ARIMA can't handle more than 720 observations.\n"
			       "Please select a smaller sample."));
	armod.errcode = E_DATA;
	return armod;
    }	

    if (opt & OPT_X) {
	arma_func = get_plugin_function("arma_x12_model", &handle);
    } else {
	arma_func = get_plugin_function("arma_model", &handle);
    }

    if (arma_func == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	gretl_model_init(&armod);
	armod.errcode = E_FOPEN;
	return armod;
    }

    armod = (*arma_func) (list, Z, pdinfo, opt, prn);

    close_plugin(handle);
    set_model_id(&armod);

    return armod;
} 

/**
 * tobit_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @prn: printing struct for iteration info (or %NULL is this is not
 * wanted).
 *
 * Produce Tobit estimates of the model given in @list.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tobit_model (const int *list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL tmod;
    void *handle;
    MODEL (* tobit_estimate) (const int *, double ***, DATAINFO *, PRN *);

    *gretl_errmsg = '\0';

    tobit_estimate = get_plugin_function("tobit_estimate", &handle);
    if (tobit_estimate == NULL) {
	gretl_model_init(&tmod);
	tmod.errcode = E_FOPEN;
	return tmod;
    }

    tmod = (*tobit_estimate) (list, pZ, pdinfo, prn);

    close_plugin(handle);

    set_model_id(&tmod);

    return tmod;
}

static int get_offset_var (int *list)
{
    int l0 = list[0];
    int ret = 0;

    if (list[l0 - 1] == LISTSEP) {
	ret = list[l0];
	list[0] -= 2;
    }

    return ret;
}

/**
 * poisson_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: printing struct for iteration info (or %NULL is this is not
 * wanted).
 *
 * Estimate the Poisson regression model given in @list using ML.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL poisson_model (const int *list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL pmodel;
    void *handle;
    int *listcpy;
    int offvar;
    int (* poisson_estimate) (MODEL *, int, double ***, DATAINFO *, PRN *);

    *gretl_errmsg = '\0';

    gretl_model_init(&pmodel);

    listcpy = gretl_list_copy(list);
    if (listcpy == NULL) {
	pmodel.errcode = E_ALLOC;
        return pmodel;
    }

    offvar = get_offset_var(listcpy);

    /* run an initial OLS to "set the model up" and check for errors.
       the poisson_estimate_driver function will overwrite the
       coefficients etc.
    */

    pmodel = lsq(listcpy, pZ, pdinfo, OLS, OPT_A);
    free(listcpy);

    if (pmodel.errcode) {
        return pmodel;
    }

    poisson_estimate = get_plugin_function("poisson_estimate", &handle);

    if (poisson_estimate == NULL) {
	pmodel.errcode = E_FOPEN;
	return pmodel;
    }

    (*poisson_estimate) (&pmodel, offvar, pZ, pdinfo, prn);

    close_plugin(handle);

    set_model_id(&pmodel);

    return pmodel;
}

/**
 * garch:
 * @list: dependent variable plus arch and garch orders.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt: can specify robust standard errors and VCV.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Calculate GARCH estimates.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL garch (const int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn)
{
    MODEL gmod;
    void *handle;
    PRN *myprn;
    MODEL (*garch_model) (const int *, double ***, DATAINFO *, PRN *,
			  gretlopt);

    *gretl_errmsg = '\0';

    garch_model = get_plugin_function("garch_model", &handle);

    if (garch_model == NULL) {
	gretl_model_init(&gmod);
	gmod.errcode = E_FOPEN;
	return gmod;
    }

    if (opt & OPT_V) {
	myprn = prn;
    } else {
	myprn = NULL;
    }

    gmod = (*garch_model) (list, pZ, pdinfo, myprn, opt);

    close_plugin(handle);

    set_model_id(&gmod);

    return gmod;
} 

/**
 * mp_ols:
 * @list: specification of variables to use.
 * @Z: data array.
 * @pdinfo: information on the data set.
 *
 * Estimate an OLS model using multiple-precision arithmetic
 * via the GMP library.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL mp_ols (const int *list, const double **Z, DATAINFO *pdinfo)
{
    void *handle = NULL;
    int (*mplsq)(const int *, const int *, const double **, 
		 DATAINFO *, char *, MODEL *, gretlopt);
    MODEL mpmod;

    gretl_model_init(&mpmod);

    mplsq = get_plugin_function("mplsq", &handle);
    if (mplsq == NULL) {
	mpmod.errcode = 1;
	return mpmod;
    }

    if (gretl_list_has_separator(list)) {
	int *base = NULL;
	int *poly = NULL;

	gretl_list_split_on_separator(list, &base, &poly);
	if (base == NULL || poly == NULL) {
	    mpmod.errcode = E_ALLOC;
	} else {
	    mpmod.errcode = (*mplsq)(base, poly, Z, pdinfo,  
				     gretl_errmsg, &mpmod, OPT_S);
	}
	free(base);
	free(poly);
    } else {
	mpmod.errcode = (*mplsq)(list, NULL, Z, pdinfo,  
				 gretl_errmsg, &mpmod, OPT_S); 
    }

    close_plugin(handle);

    set_model_id(&mpmod);

    return mpmod;
}

static int check_panel_options (gretlopt opt)
{
    int err = 0;

    if ((opt & OPT_R) && (opt & OPT_W)) {
	/* can't specify random effects + weighted least squares */
	err = E_DATA;
    } else if ((opt & OPT_T) && !(opt & OPT_W)) {
	/* iterate option requires weighted least squares option */
	err = E_DATA;
    }

    return err;
}

/**
 * panel_model:
 * @list: regression list (dependent variable plus independent 
 * variables).
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the (panel) data set.
 * @opt: can include %OPT_Q (quiet estimation), %OPT_S
 * (silent estimation), %OPT_R (random effects model),
 * %OPT_W (weights based on the error variance for the
 * respective cross-sectional units), %OPT_T (iterate, only
 * available in conjunction with %OPT_W).
 * @prn: printing struct (or %NULL).
 *
 * Calculate estimates for a panel dataset, using fixed
 * effects (the default), random effects, or weighted
 * least squares based on the respective variances for the
 * cross-sectional units.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn)
{
    MODEL mod;

    *gretl_errmsg = '\0';

    if (check_panel_options(opt)) {
	gretl_model_init(&mod);
	mod.errcode = E_BADOPT;
    } else if (opt & OPT_W) {
	mod = panel_wls_by_unit(list, pZ, pdinfo, opt, prn);
    } else {
	mod = real_panel_model(list, pZ, pdinfo, opt, prn);
    }

    return mod;
}

/**
 * arbond_model:
 * @list: regression list.
 * @istr: may contain additional instrument specification.
 * @Z: data array.
 * @pdinfo: information on the (panel) data set.
 * @opt: to be hooked up.
 * @prn: printing struct (or %NULL).
 *
 * To be written.  This function is currently just for
 * testing.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arbond_model (const int *list, const char *istr, const double **Z, 
		    const DATAINFO *pdinfo, gretlopt opt, 
		    PRN *prn)
{
    void *handle = NULL;
    MODEL (*arbond_estimate) (const int *, const char *, const double **, 
			      const DATAINFO *, gretlopt, PRN *);
    MODEL mod;

    gretl_model_init(&mod);

    arbond_estimate = get_plugin_function("arbond_estimate", &handle);
    if (arbond_estimate == NULL) {
	mod.errcode = 1;
	return mod;
    }

    mod = (*arbond_estimate)(list, istr, Z, pdinfo, opt, prn);

    close_plugin(handle);

    if (!mod.errcode) {
	set_model_id(&mod);
    }

    return mod;    
}

/**
 * groupwise_hetero_test:
 * @pmod: pooled OLS model to be tested.
 * @pZ: pointer to data array.
 * @pdinfo: information on the (panel) data set.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Calculates iterated WLS estimates using weights based on the error
 * variance for the cross-sectional units and performs a Wald test
 * for the null hypothesis that the error variance is uniform
 * across the units.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int groupwise_hetero_test (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			   PRN *prn)
{
    MODEL wmod;
    int err;

    if (pmod->ci != OLS) {
	return E_NOTIMP;
    }

    if (!dataset_is_panel(pdinfo)) {
	strcpy(gretl_errmsg, _("This test is only available for panel data"));
	return 1;
    }

    wmod = panel_wls_by_unit(pmod->list, pZ, pdinfo, OPT_T | OPT_A, prn);
    err = wmod.errcode;

    if (!err) {
	gretl_model_set_auxiliary(&wmod, AUX_GROUPWISE);
	printmodel(&wmod, pdinfo, OPT_NONE, prn);
    }

    clear_model(&wmod);

    return err;
}


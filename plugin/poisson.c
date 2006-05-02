/* 
 * Copyright (C) 2005 Riccardo "Jack" Lucchetti
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

#undef PDEBUG 

#define POISSON_TOL 1.0e-10 
#define POISSON_MAX_ITER 100 

/* check whether a series contains nothing but non-negative
   integer values (some of which are > 1) */

static int is_count_variable (const double *x, int t1, int t2)
{
    int t, xi;
    int g1 = 0;
    int ret = 1;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	if (x[t] < 0.0) {
	    ret = 0;
	    break;
	}
	xi = x[t];
	if (x[t] != (double) xi) {
	    ret = 0;
	    break;
	}
	if (x[t] > 1.0) {
	    g1 = 1;
	}
    }

    if (g1 == 0) {
	ret = 0;
    }

    return ret;
}

static double poisson_ll (const double *y, const double *mu, 
			  int t1, int t2)
{
    double loglik = 0.0;
    double ytfact, llt;
    int t;

    for (t=t1; t<=t2; t++) {
	if (na(y[t]) || na(mu[t])) {
	    continue;
	}
	ytfact = x_factorial(y[t]);
	if (na(ytfact)) {
	    loglik = NADBL;
	    break;
	}
	llt = (-mu[t] + y[t] * log(mu[t]) - log(ytfact));
	loglik += llt;
    }  

    return loglik;
}

static int 
transcribe_poisson_results (MODEL *targ, MODEL *src, const double *y, 
			    int iter, int offvar)
{
    double sigma = src->sigma;
    double sigma2 = sigma * sigma;
    int i, j, t;
    int err = 0;

    targ->ci = POISSON;
    
    gretl_model_set_int(targ, "iters", iter);

    if (offvar > 0) {
	gretl_model_set_int(targ, "offset_var", offvar);
    }

    targ->ess = 0.0;

    for (t=targ->t1; t<=targ->t2; t++) {
	if (na(targ->yhat[t])) {
	    targ->uhat[t] = NADBL;
	} else {
	    targ->uhat[t] = y[t] - targ->yhat[t];
	    targ->ess += targ->uhat[t] * targ->uhat[t];
	}
    }

    targ->sigma = sqrt(targ->ess / targ->dfd);

    for (i=0; i<targ->ncoeff; i++) {
	targ->sderr[i] = src->sderr[i] / sigma;
    }

    targ->lnL = poisson_ll(y, targ->yhat, targ->t1, targ->t2);

#ifdef PDEBUG
    fprintf(stderr, "log-likelihood = %g\n", targ->lnL);
#endif

    mle_criteria(targ, 0); 

    /* mask invalid statistics */
    targ->rsq = NADBL;
    targ->adjrsq = NADBL;
    targ->fstt = NADBL;

    /* make and correct the variance matrix */
    if (makevcv(src)) {
	err = 1;
    } else {
	int idx;

	if (targ->vcv != NULL) {
	    free(targ->vcv);
	}
	targ->vcv = src->vcv;
	src->vcv = NULL;
	for (i=0; i<targ->ncoeff; i++) {
	    for (j=i; j<targ->ncoeff; j++) {
		idx = ijton(i, j, targ->ncoeff);
		targ->vcv[idx] /= sigma2;
	    }
	}
    }   

    return err;
}

static double *get_offset (MODEL *pmod, int offvar, double **Z,
			   double *offmean)
{
    double *offset = NULL;
    int t, err = 0;

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	} else if (na(Z[offvar][t])) {
	    err = 1;
	} else if (Z[offvar][t] < 0.0) {
	    err = 1;
	} 
    }

    if (err == 0) {
	offset = Z[offvar];
	*offmean = gretl_mean(pmod->t1, pmod->t2, offset);
    }

    return offset;
}

static int 
do_poisson (MODEL *pmod, int offvar, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int origv = pdinfo->v;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int i, t;

    int iter = 0;
    double crit = 1.0;

    double *offset = NULL;
    double offmean = NADBL;

    double *y;
    double *wgt;
    double *depvar;

    MODEL tmpmod;
    int *local_list = NULL;

    gretl_model_init(&tmpmod);

    /* set the sample to that of the initial OLS model */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    local_list = gretl_list_new(pmod->list[0] + 1);
    if (local_list == NULL) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    if (dataset_add_series(2, pZ, pdinfo)) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    if (offvar > 0) {
	offset = get_offset(pmod, offvar, *pZ, &offmean);
	if (offset == NULL) {
	    pmod->errcode = E_DATA;
	    goto bailout;
	}
    }

    /* the original dependent variable */
    y = (*pZ)[pmod->list[1]];

    /* local_list includes a weight variable */
    local_list[0] = pmod->list[0] + 1;

    /* weighting variable (first newly added var) */
    local_list[1] = origv;
    wgt = (*pZ)[origv];

    /* dependent variable for GNR (second newly added var) */
    local_list[2] = origv + 1;
    depvar = (*pZ)[origv + 1];
    
    for (i=3; i<=local_list[0]; i++) { 
	/* original independent vars */
	local_list[i] = pmod->list[i-1];
    }    

    pmod->coeff[0] = log(pmod->ybar);
    if (offvar > 0) {
	pmod->coeff[0] -= log(offmean);
    }

    for (i=1; i<pmod->ncoeff; i++) { 
	pmod->coeff[i] = 0.0;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    depvar[t] = NADBL;
	    wgt[t] = NADBL;
	} else {
	    pmod->yhat[t] = pmod->ybar;
	    if (offvar > 0) {
		pmod->yhat[t] *= offset[t] / offmean;
	    }
	    depvar[t] = y[t] / pmod->yhat[t] - 1.0;
	    wgt[t] = pmod->yhat[t];
	}
    }

    pputc(prn, '\n');

    while (iter < POISSON_MAX_ITER && crit > POISSON_TOL) {

	iter++;

	tmpmod = lsq(local_list, pZ, pdinfo, WLS, OPT_A);

	if (tmpmod.errcode) {
	    fprintf(stderr, "poisson_estimate: lsq returned %d\n", 
		    tmpmod.errcode);
	    pmod->errcode = tmpmod.errcode;
	    break;
	}

	crit = tmpmod.nobs * tmpmod.rsq;

	pprintf(prn, "%s %3d\tcrit = %g\n", _("iteration"), iter, crit);

	for (i=0; i<tmpmod.ncoeff; i++) { 
	    pmod->coeff[i] += tmpmod.coeff[i];
#ifdef PDEBUG
	    fprintf(stderr, "coeff[%d] = %g,\tgrad[%d] = %g\n", 
		    i, pmod->coeff[i], i, tmpmod.coeff[i]);
#endif
	}

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (na(pmod->uhat[t])) {
		continue;
	    }
	    pmod->yhat[t] *= exp(tmpmod.yhat[t]);
	    depvar[t] = y[t] / pmod->yhat[t] - 1;
	    wgt[t] = pmod->yhat[t];
	}

	if (crit > POISSON_TOL) {
	    clear_model(&tmpmod);
	}
    }

    pputc(prn, '\n');

    if (crit > POISSON_TOL) {
	pmod->errcode = E_NOCONV;
    } 

    if (pmod->errcode == 0) {
	transcribe_poisson_results(pmod, &tmpmod, y, iter, offvar);
    }

 bailout:

    clear_model(&tmpmod);
    free(local_list);
    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return pmod->errcode;
}

int 
poisson_estimate (MODEL *pmod, int offvar, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn) 
{
    int err = 0;

    if (!is_count_variable((*pZ)[pmod->list[1]], pmod->t1, pmod->t2)) {
	gretl_errmsg_set(_("poisson: the dependent variable must be count data"));
	err = pmod->errcode = E_DATA;
    } else {
	err = do_poisson(pmod, offvar, pZ, pdinfo, prn);
    }

    return err;
}


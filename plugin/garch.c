/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 by Allin Cottrell
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

/* GARCH plugin for gretl using the Fiorentini, Calzolari and 
   Panattoni mixed-gradient algorithm.
*/

#include "libgretl.h"
#include "libset.h"

#include "fcp.h"

static void add_garch_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				const int *list)
{
    int p = list[1];
    int q = list[2];
    int r = list[0] - 4;
    int i, j, np = 3 + p + q + r;

    free(pmod->list);
    pmod->list = copylist(list);

    pmod->params = malloc(np * sizeof pmod->params);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    pmod->nparams = np;

    for (i=0; i<np; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    pmod->errcode = E_ALLOC;
	    return;
	}
    }

    strcpy(pmod->params[0], pdinfo->varname[pmod->list[4]]);
    strcpy(pmod->params[1], pdinfo->varname[0]);

    j = 2;
    for (i=0; i<r; i++) {
	if (pmod->list[5+i] > 0) {
	    strcpy(pmod->params[j++], pdinfo->varname[pmod->list[5+i]]);
	}
    }

    strcpy(pmod->params[j++], "alpha(0)");

    for (i=0; i<q; i++) {
	sprintf(pmod->params[j++], "alpha(%d)", i + 1);
    }
    for (i=0; i<p; i++) {
	sprintf(pmod->params[j++], "beta(%d)", i + 1);
    }
}

static int make_packed_vcv (MODEL *pmod, double *vcv, int np,
			    int nc, double scale)
{
    const int nterms = np * (np + 1) / 2;
    double sfi, sfj;
    int i, j, k;

    free(pmod->vcv);
    pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return 1;  
    }

    for (i=0; i<np; i++) {
	if (i < nc) {
	    sfi = scale;
	} else if (i == nc) {
	    sfi = scale * scale;
	} else {
	    sfi = 1.0;
	}
	for (j=0; j<=i; j++) {
	    if (j < nc) {
		sfj = scale;
	    } else if (j == nc) {
		sfj = scale * scale;
	    } else {
		sfj = 1.0;
	    }
	    k = ijton(i, j, np);
	    pmod->vcv[k] = vcv[i + np * j] * sfi * sfj;
	}
    }

    return 0;
}

static int write_garch_stats (MODEL *pmod, const double **Z,
			      double scale, const DATAINFO *pdinfo,
			      const int *list, const double *theta, 
			      int nparam, int pad, const double *res,
			      const double *h)
{
    int err = 0;
    double *coeff, *sderr, *garch_h;
    int i, ynum = list[4];

    coeff = realloc(pmod->coeff, nparam * sizeof *pmod->coeff);
    sderr = realloc(pmod->sderr, nparam * sizeof *pmod->sderr);

    if (coeff == NULL || sderr == NULL) return 1;

    for (i=0; i<nparam; i++) {
	coeff[i] = theta[i+1];
	sderr[i] = theta[i+nparam+1];
    }
    
    pmod->coeff = coeff;
    pmod->sderr = sderr;
   
    pmod->ncoeff = nparam;

    pmod->ess = 0.0;
    for (i=pmod->t1; i<=pmod->t2; i++) {
	pmod->uhat[i] = res[i + pad] * scale;
	pmod->ess += pmod->uhat[i] * pmod->uhat[i];
	pmod->yhat[i] =  Z[ynum][i] - pmod->uhat[i];
    }

    pmod->sigma = NADBL;
    pmod->adjrsq = NADBL; 
    pmod->fstt = NADBL;

    pmod->criterion[C_AIC] = -2.0 * pmod->lnL + 2.0 * (pmod->ncoeff + 1);
    pmod->criterion[C_BIC] = -2.0 * pmod->lnL + (pmod->ncoeff + 1) * log(pmod->nobs);

    pmod->ci = GARCH;
    
    add_garch_varnames(pmod, pdinfo, list);

    /* add predicted error variance to model */
    garch_h = malloc(pdinfo->n * sizeof *garch_h);
    if (garch_h != NULL) {
	for (i=0; i<pdinfo->n; i++) {
	    if (i < pmod->t1 || i > pmod->t2) {
		garch_h[i] = NADBL;
	    } else {
		garch_h[i] = h[i + pad] * scale * scale;
	    }
	}
	gretl_model_set_data(pmod, "garch_h", garch_h, 
			     pdinfo->n * sizeof *garch_h);
    }

    return err;
}

static int make_garch_dataset (const int *list, double **Z,
			       int bign, int pad, int nx,
			       double **py, double ***pX)
{
    double *y = NULL, **X = NULL;
    int i, k = 0, t;
    int xnum, ynum = list[4];

    /* If pad > 0 we have to create a newly allocated, padded
       dataset.  Otherwise we can use a virtual dataset, made
       up of pointers into the original dataset, Z. 
    */

    if (pad > 0) {
	y = malloc(bign * sizeof *y);
	if (y == NULL) return 1;
    } 

    if (nx > 0) {
	X = malloc(nx * sizeof *X);
	if (X == NULL) goto bailout;

	if (pad > 0) {
	    for (i=0; i<nx; i++) {
		X[i] = malloc(bign * sizeof **X);
		if (X[i] == NULL) {
		    for (t=0; t<i; t++) {
			free(X[t]);
		    }
		    free(X);
		    goto bailout;
		}
	    } 
	}  
    }

    if (pad > 0) {
	/* build padded dataset */
	for (t=0; t<bign; t++) {
	    if (t < pad) {
		y[t] = 0.0;
		for (i=0; i<nx; i++) {
		    X[i][t] = 0.0;
		}
	    } else {
		y[t] = Z[ynum][t-pad];
		if (nx > 0) k = 5;
		for (i=0; i<nx; i++) {
		    xnum = list[k++]; 
		    if (xnum == 0) xnum = list[k++];
		    X[i][t] = Z[xnum][t-pad];
		}
	    }
	}
	*py = y;
    } else {
	/* build virtual dataset */
	*py = Z[ynum];
	if (nx > 0) k = 5;
	for (i=0; i<nx; i++) {
	    xnum = list[k++]; 
	    if (xnum == 0) xnum = list[k++];
	    X[i] = Z[xnum];
	}
    }

    *pX = X;

    return 0;

 bailout:

    free(y);
    return E_ALLOC;
}

static int get_vopt (int robust)
{
    int vopt = get_garch_vcv_version();

    /* The defaults: QML if "robust" option is in force,
       otherwise negative Hessian */
    if (vopt == VCV_UNSET) {
	if (robust) {
	    vopt = VCV_QML;
	} else {
	    vopt = VCV_HESSIAN;
	}
    }

    return vopt;
}

int do_fcp (const int *list, double **Z, double scale,
	    const DATAINFO *pdinfo, MODEL *pmod,
	    PRN *prn, gretlopt opt)
{
    int t1 = pmod->t1, t2 = pmod->t2;
    int ncoeff = pmod->ncoeff;
    int p = list[1];
    int q = list[2];
    double *y = NULL;
    double **X = NULL;
    double *h = NULL;
    double *yhat = NULL, *amax = NULL; 
    double *res = NULL, *res2 = NULL;
    double *coeff = NULL, *b = NULL;
    double *vcv = NULL;
    int err = 0, iters = 0;
    int nobs, maxlag, bign, pad = 0;
    int i, nx, nparam, vopt;

    vopt = get_vopt(opt & OPT_R);

    nx = ncoeff - 1;
    maxlag = (p > q)? p : q; 
    nparam = ncoeff + p + q + 1;

    nobs = t2 + 1; /* number of obs in full dataset */

    if (maxlag > t1) {
	/* need to pad data series at start */
	pad = maxlag - t1;
    } 

    /* length of series to pass to garch_estimate */
    bign = nobs + pad;
	
    yhat = malloc(bign * sizeof *yhat);
    res2 = malloc(bign * sizeof *res2);
    res = malloc(bign * sizeof *res);
    h = malloc(bign * sizeof *h);
    amax = malloc(bign * sizeof *amax);
    if (yhat == NULL || res2 == NULL || res == NULL || 
	amax == NULL || h == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    for (i=0; i<bign; i++) {
	yhat[i] = res2[i] = res[i] = amax[i] = 0.0;
    }   
 
    coeff = malloc(ncoeff * sizeof *coeff);
    b = malloc(ncoeff * sizeof *b);
    if (coeff == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;	
    }

    vcv = malloc((nparam * nparam) * sizeof *vcv);
    if (vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    for (i=0; i<nparam * nparam; i++) {
	vcv[i] = 0.0;
    } 

    /* create dataset for garch estimation */
    err = make_garch_dataset(list, Z, bign, pad, nx, &y, &X);
    if (err) {
	goto bailout;
    }

    /* initial coefficients from OLS */
    for (i=0; i<ncoeff; i++) {
	coeff[i] = pmod->coeff[i];
	b[i] = 0.0;
    }

    amax[0] = pmod->sigma * pmod->sigma;
    amax[1] = q;
    amax[2] = p; 
    for (i=0; i<p+q; i++) {
	/* initial alpha, beta values */
	amax[3+i] = 0.1;
    }

    err = garch_estimate(t1 + pad, t2 + pad, bign, 
			 (const double **) X, nx, yhat, coeff, ncoeff, 
			 vcv, res2, res, h, y, amax, b, scale, &iters,
			 prn, vopt);

    if (err != 0) {
	pmod->errcode = err;
    } else {
	int nparam = ncoeff + p + q + 1;

	for (i=1; i<=nparam; i++) {
	    if (i <= ncoeff) {
		amax[i] *= scale;
		amax[i + nparam] *= scale;
	    } else if (i == ncoeff + 1) {
		amax[i] *= scale * scale;
		amax[i + nparam] *= scale * scale;
	    }
	    pprintf(prn, "theta[%d]: %#14.6g (%#.6g)\n", i-1, amax[i], 
		    amax[i + nparam]);
	}
	pputc(prn, '\n');

	pmod->lnL = amax[0];
	write_garch_stats(pmod, (const double **) Z, scale, pdinfo, 
			  list, amax, nparam, pad, res, h);
	make_packed_vcv(pmod, vcv, nparam, ncoeff, scale);
	gretl_model_set_int(pmod, "iters", iters);
	gretl_model_set_int(pmod, "garch_vcv", vopt);
    }

 bailout:

    free(yhat);
    free(res2);
    free(res);
    free(h);
    free(amax);    
    free(coeff);
    free(b);
    free(vcv); 

    if (pad > 0) {
	/* don't free y if it's just a pointer into Z */
	free(y);
    }
    
    if (X != NULL) {
	if (pad > 0) {
	    /* don't free the X[i] if they're just pointers into
	       the original data matrix */
	    for (i=0; i<nx; i++) {
		free(X[i]);
	    }
	}
	free(X);
    }

    return err;
}

/* sanity/dimension check */

static int *get_garch_list (const int *list, int *err)
{
    int *ret = NULL;
    int i, p = list[1], q = list[2];
    int add0 = 1;

    *err = 0;

    /* rule out pure AR in variance */
    if (p > 0 && q == 0) {
	gretl_errmsg_set(_("Error in garch command"));
	*err = E_DATA;
	return NULL;
    }

    /* rule out excessive total GARCH terms */
    else if (p + q > 5) {
	gretl_errmsg_set(_("Error in garch command"));
	*err = E_DATA;
	return NULL;
    }

    /* insert constant if not present */
    for (i=4; i<=list[0]; i++) {
	if (list[i] == 0) {
	    add0 = 0;
	    break;
	}
    }

    ret = malloc((list[0] + 1 + add0) * sizeof *ret);
    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	ret[0] = list[0] + add0;
	for (i=1; i<=list[0]; i++) {
	    ret[i] = list[i];
	}
	if (add0) {
	    ret[i] = 0;
	}
    }

    return ret;
}

/* make regresson list for initial OLS */

static int *make_ols_list (const int *list)
{
    int *olist;
    int i;

    olist = malloc((list[0] - 2) * sizeof *olist);
    if (olist == NULL) return NULL;

    olist[0] = list[0] - 3;
    for (i=4; i<=list[0]; i++) {
	olist[i-3] = list[i];
    }

    return olist;
}

/* the driver function for the plugin */

MODEL garch_model (int *cmdlist, double ***pZ, DATAINFO *pdinfo,
		   PRN *prn, gretlopt opt) 
{
    MODEL model;
    int *list = NULL, *ols_list = NULL;
    double scale = 1.0;
    int t, err, yno = 0;

    gretl_model_init(&model);

    list = get_garch_list(cmdlist, &err);
    if (err) {
	model.errcode = err;
    }

    if (!err) {
	ols_list = make_ols_list(list);
	if (ols_list == NULL) {
	    err = model.errcode = E_ALLOC;
	}
    }

    /* run initial OLS */
    if (!err) {
	model = lsq(ols_list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (model.errcode) {
	    err = model.errcode;
	}
    }

#if 1
    if (!err) {
	yno = ols_list[1];
	scale = gretl_stddev(model.t1, model.t2, (*pZ)[yno]);
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[yno][t] /= scale;
	}
	for (t=0; t<model.ncoeff; t++) {
	    model.coeff[t] *= scale;
	}
    } 
#endif

    if (!err) {
	do_fcp(list, *pZ, scale, pdinfo, &model, prn, opt); 
    }

    if (scale != 1.0) {
	/* undo scaling of dependent variable */
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[yno][t] *= scale;
	}
    }

    free(ols_list);
    free(list);

    return model;
}

#ifdef STANDALONE

int main (void) 
{
    char *fname;
    MODEL model;
    double **Z = NULL;
    DATAINFO *datainfo;
    PRN *prn;
    int *list;
    int err;

    datainfo = datainfo_new();
    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    /* 4 1 1 999 1 */

    list = malloc(5 * sizeof *list);
    list[0] = 4;
    list[1] = list[2] = list[4] = 1;
    list[3] = 999;

    fname = malloc(128 * sizeof *fname);
    strcpy(fname, "/opt/esl/share/gretl/data/misc/b-g.gdt");

    err = get_xmldata(&Z, &datainfo, fname, 
		      NULL, DATA_NONE, prn, 0);

    if (!err) {
	model = garch_model(list, &Z, datainfo,
			    prn, OPT_NONE);
    } 

    free_Z(Z, datainfo);
    clear_model(&model, NULL);
    free_datainfo(datainfo);
    gretl_print_destroy(prn);
    free(list);
    free(fname);

    return 0;
}


#endif

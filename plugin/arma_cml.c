/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* ARMA estimation via conditional ML (BHHH method) */

#include "libgretl.h"
#include "bhhh_max.h"
#include "libset.h"
#include "arma_priv.h"

#define CML_DEBUG 0

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

/* Revising the sample: this may be necessary when we're using
   CML as initializer for exact ML. CML requires that all AR
   lags of y are available, which may dictate a later start of
   the estimation range. If we didn't adjust the t1 member of
   @ainfo we could end up calculating with missing values, or
   reading off the beginning of the dataset.
*/

static int cml_init_revise_sample (arma_info *ainfo,
				   const DATASET *dset)
{
    int pmax = ainfo->p + ainfo->P * ainfo->pd;
    int dmax = ainfo->d + ainfo->D * ainfo->pd;
    int cml_maxlag = pmax + dmax;
    int *list = ainfo->alist;
    int ypos = arma_list_y_position(ainfo);
    int t0, t1 = ainfo->t1, t2 = ainfo->t2;
    int i, vi, vlmax, k, t;
    int missing;
    int err = 0;

#if CML_DEBUG
    fprintf(stderr, "cml_init_revise_sample: initial t1=%d, T=%d\n",
	    ainfo->t1, ainfo->T);
#endif

    t0 = t1 - cml_maxlag;
    if (t0 < 0) {
	t1 -= t0;
    }

    if (arma_xdiff(ainfo)) {
	vlmax = list[0];
    } else {
	vlmax = ypos;
    }

    /* advance the starting point if need be */
    for (t=t1; t<=t2; t++) {
	missing = 0;
	for (i=ypos; i<=list[0] && !missing; i++) {
	    vi = list[i];
	    if (na(dset->Z[vi][t])) {
		missing = 1;
	    }
	    if (i <= vlmax) {
		for (k=1; k<=cml_maxlag && !missing; k++) {
		    if (na(dset->Z[vi][t-k])) {
			missing = 1;
		    }
		}
	    }
	}
	if (missing) {
	    t1++;
	} else {
	    break;
	}
    }

    /* check for missing obs within the adjusted sample range */
    for (t=t1; t<t2 && !err; t++) {
	for (i=ypos; i<=list[0]; i++) {
	    if (na(dset->Z[list[i]][t])) {
		err = E_MISSDATA;
		break;
	    }
	}
    }

    if (!err && t2 - t1 + 1 <= ainfo->nc) {
	err = E_DF;
    }

    if (!err) {
	ainfo->t1 = t1;
	ainfo->T = ainfo->t2 - ainfo->t1 + 1;
    }

#if CML_DEBUG
    fprintf(stderr, "cml_init_revise_sample: revised t1=%d, T=%d\n",
	    ainfo->t1, ainfo->T);
#endif

    return err;
}

static void do_MA_partials (double *drv,
			    arma_info *ainfo,
			    const double *theta,
			    const double *Theta,
			    int t)
{
    int i, j, k, p, s;

    k = 0;
    for (i=0; i<ainfo->q; i++) {
	if (MA_included(ainfo, i)) {
	    p = i + 1;
	    if (t - p >= 0) {
		drv[0] -= theta[k] * drv[p];
	    }
	    k++;
	}
    }

    for (j=0; j<ainfo->Q; j++) {
	p = (j + 1) * ainfo->pd;
	if (t - p >= 0) {
	    drv[0] -= Theta[j] * drv[p];
	    k = 0;
	    for (i=0; i<ainfo->q; i++) {
		if (MA_included(ainfo, i)) {
		    s = p + i + 1;
		    if (t - s >= 0) {
			drv[0] -= Theta[j] * theta[k] * drv[s];
		    }
		    k++;
		}
	    }
	}
    }
}

/* for each of the arrays of derivatives, shuffle each
   value one place up */

static void push_derivs (arma_info *ainfo, double **de, int dlen)
{
    int i, j;

    for (i=0; i<ainfo->n_aux; i++) {
	for (j=dlen-1; j>0; j--) {
	    de[i][j] = de[i][j-1];
	}
	de[i][0] = 0.0;
    }
}

static void zero_derivs (arma_info *ainfo, double **de, int dlen)
{
    int i, j;

    for (i=0; i<ainfo->n_aux; i++) {
	for (j=0; j<dlen; j++) {
	    de[i][j] = 0.0;
	}
    }
}

#if 0

/* unused at present: see comment in arma_analytical_score() */

static int adjust_score_t1 (arma_info *ainfo, const double *y)
{
    int p, pmax = ainfo->p + ainfo->pd * ainfo->P;
    int miss, t, t1 = ainfo->t1;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	miss = 0;
	for (p=1; p<=pmax; p++) {
	    if (t - p > 0 && na(y[t-p])) {
		miss = 1;
		t1++;
		break;
	    }
	}
	if (!miss) {
	    break;
	}
    }

    return t1;
}

#endif

static int arma_analytical_score (arma_info *ainfo,
				  const double *y,
				  const double **X,
				  const double *phi,
				  const double *Phi,
				  const double *theta,
				  const double *Theta,
				  double s2,
				  gretl_matrix *G)
{
    /* forecast errors */
    const double *e = ainfo->e;
    /* pointers to blocks of derivatives (workspace) */
    double **de = ainfo->aux;
    double **de_a =    de + ainfo->ifc;
    double **de_sa = de_a + ainfo->np;
    double **de_m = de_sa + ainfo->P;
    double **de_sm = de_m + ainfo->nq;
    double **de_r = de_sm + ainfo->Q;
    int dlen = 1 + ainfo->q + ainfo->pd * ainfo->Q;
    double x, Gsi;
    int t, gt, t1 = ainfo->t1;
    int i, j, k, p, s;
    int err = 0;

    zero_derivs(ainfo, de, dlen);

#if 0 /* currently this function is never called when
	 estimation is by exact ML */
    if (arma_exact_ml(ainfo)) {
	t1 = adjust_score_t1(ainfo, y);
    }
#endif

    for (t=t1, gt=0; t<=ainfo->t2 && !err; t++, gt++) {

	/* the constant term (de_0) */
	if (ainfo->ifc) {
	    de[0][0] = -1.0;
	    do_MA_partials(de[0], ainfo, theta, Theta, t);
	}

	/* non-seasonal AR terms (de_a) */
	k = 0;
	for (i=0; i<ainfo->p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    p = i + 1;
	    if (t - p >= 0) {
		de_a[k][0] = -y[t-p];
		/* cross-partial with seasonal AR */
		for (j=0; j<ainfo->P; j++) {
		    s = p + (j + 1) * ainfo->pd;
		    if (t - s >= 0) {
			de_a[k][0] += Phi[j] * y[t-s];
		    }
		}
		do_MA_partials(de_a[k], ainfo, theta, Theta, t);
	    }
	    k++;
	}

	/* seasonal AR terms (de_sa) */
	for (j=0; j<ainfo->P; j++) {
	    p = (j + 1) * ainfo->pd;
	    if (t - p >= 0) {
		de_sa[j][0] = -y[t-p];
		/* cross-partial with non-seasonal AR */
		k = 0;
		for (i=0; i<ainfo->p; i++) {
		    if (AR_included(ainfo, i)) {
			s = p + i + 1;
			if (t - s >= 0) {
			    de_sa[j][0] += phi[k] * y[t-s];
			}
			k++;
		    }
		}
		do_MA_partials(de_sa[j], ainfo, theta, Theta, t);
	    }
	}

	/* non-seasonal MA terms (de_m) */
	k = 0;
	for (i=0; i<ainfo->q; i++) {
	    if (!MA_included(ainfo, i)) {
		continue;
	    }
	    p = i + 1;
	    if (t - p >= 0) {
		de_m[k][0] = -e[t-p];
		/* cross-partial with seasonal MA */
		for (j=0; j<ainfo->Q; j++) {
		    s = p + (j + 1) * ainfo->pd;
		    if (t - s >= 0) {
			de_m[k][0] -= Theta[j] * e[t-s];
		    }
		}
		do_MA_partials(de_m[k], ainfo, theta, Theta, t);
	    }
	    k++;
	}

	/* seasonal MA terms (de_sm) */
	for (j=0; j<ainfo->Q; j++) {
	    p = (j + 1) * ainfo->pd;
	    if (t - p >= 0) {
		de_sm[j][0] = -e[t-p];
		/* cross-partial with non-seasonal MA */
		k = 0;
		for (i=0; i<ainfo->q; i++) {
		    if (MA_included(ainfo, i)) {
			s = p + i + 1;
			if (t - s >= 0) {
			    de_sm[j][0] -= theta[k] * e[t-s];
			}
			k++;
		    }
		}
		do_MA_partials(de_sm[j], ainfo, theta, Theta, t);
	    }
	}

	/* exogenous regressors (de_r) */
	for (j=0; j<ainfo->nexo; j++) {
	    de_r[j][0] = -X[j][t];
	    do_MA_partials(de_r[j], ainfo, theta, Theta, t);
	}

	/* update gradient matrix */
	x = e[t] / s2; /* sqrt(s2)? does it matter? */
	for (i=0; i<ainfo->nc; i++) {
	    Gsi = -de[i][0] * x;
	    if (na(Gsi)) {
		fprintf(stderr, "arma score, bad value at t=%d, i=%d\n", t, i);
		err = E_NAN;
		break;
	    }
	    gretl_matrix_set(G, gt, i, Gsi);
	}

	push_derivs(ainfo, de, dlen);
    }

    return err;
}

static int conditional_arma_forecast_errors (arma_info *ainfo,
					     const double *y,
					     const double **X,
					     double b0,
					     const double *phi,
					     const double *Phi,
					     const double *theta,
					     const double *Theta,
					     const double *beta,
					     double *s2)
{
    double *e = ainfo->e;
    int t1 = ainfo->t1;
    int i, j, k, s, t, p;

    *s2 = 0.0;

    for (t=t1; t<=ainfo->t2; t++) {
	e[t] = y[t];

	/* intercept */
	if (ainfo->ifc) {
	    e[t] -= b0;
	}

	/* non-seasonal AR component */
	k = 0;
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		s = t - (i + 1);
		e[t] -= phi[k++] * y[s];
	    }
	}

	/* seasonal AR component plus interactions */
	for (j=0; j<ainfo->P; j++) {
	    s = t - (j + 1) * ainfo->pd;
	    e[t] -= Phi[j] * y[s];
	    k = 0;
	    for (i=0; i<ainfo->p; i++) {
		if (AR_included(ainfo, i)) {
		    p = s - (i + 1);
		    e[t] += Phi[j] * phi[k++] * y[p];
		}
	    }
	}

	/* non-seasonal MA component */
	k = 0;
	for (i=0; i<ainfo->q; i++) {
	    if (MA_included(ainfo, i)) {
		s = t - (i + 1);
		if (s >= ainfo->t1) {
		    e[t] -= theta[k] * e[s];
		}
		k++;
	    }
	}

	/* seasonal MA component plus interactions */
	for (j=0; j<ainfo->Q; j++) {
	    s = t - (j + 1) * ainfo->pd;
	    if (s >= ainfo->t1) {
		e[t] -= Theta[j] * e[s];
		k = 0;
		for (i=0; i<ainfo->q; i++) {
		    if (MA_included(ainfo, i)) {
			p = s - (i + 1);
			if (p >= ainfo->t1) {
			    e[t] -= Theta[j] * theta[k] * e[p];
			}
			k++;
		    }
		}
	    }
	}

	/* exogenous regressors */
	for (i=0; i<ainfo->nexo; i++) {
	    e[t] -= beta[i] * X[i][t];
	}

	*s2 += e[t] * e[t];
    }

    return 0;
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a callback. */

static double bhhh_arma_callback (double *coeff,
				  gretl_matrix *G,
				  void *data,
				  int do_score,
				  int *err)
{
    arma_info *ainfo = (arma_info *) data;
    /* pointers to blocks of data */
    const double *y = ainfo->Z[0];
    const double **X = ainfo->Z + 1;
    /* pointers to blocks of coefficients */
    double *phi =   coeff + ainfo->ifc;
    double *Phi =     phi + ainfo->np;
    double *theta =   Phi + ainfo->P;
    double *Theta = theta + ainfo->nq;
    double *beta =  Theta + ainfo->Q;
    double ll, s2 = 0.0;
    int any_ma = ainfo->q > 0 || ainfo->Q > 0;

    *err = 0;

#if ARMA_DEBUG
    fprintf(stderr, "bhhh_arma_callback: do_score = %d\n", do_score);
    int i;
    for (i=0; i<ainfo->nc; i++) {
	fprintf(stderr, " coeff[%d] = %.15g\n", i, coeff[i]);
    }
#endif

    if (any_ma && maybe_correct_MA(ainfo, theta, Theta)) {
	pputs(ainfo->prn, "arma: MA estimate(s) out of bounds\n");
	fputs("bhhh_arma_callback: MA estimate(s) out of bounds\n", stderr);
	*err = E_NOCONV;
	return NADBL;
    }

    conditional_arma_forecast_errors(ainfo, y, X, coeff[0],
				     phi, Phi, theta, Theta,
				     beta, &s2);

    /* error variance and log-likelihood */
    s2 /= (double) ainfo->T;
    ll = -ainfo->T * (0.5 * log(s2) + LN_SQRT_2_PI_P5);

    if (isnan(ll)) {
	*err = E_NAN;
    }

    if (do_score) {
	ainfo->ll = ll;
	arma_analytical_score(ainfo, y, X,
			      phi, Phi, theta, Theta,
			      s2, G);
    }

    return ll;
}

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

static const double **make_arma_Z (arma_info *ainfo,
				   const DATASET *dset)
{
    const double **aZ;
    int *list = ainfo->alist;
    int ypos = arma_list_y_position(ainfo);
    int nx = list[0] - ypos;
    int v, i;

#if ARMA_DEBUG
    fprintf(stderr, "make_arma_Z: allocating %d series pointers\n",
	    nx + 1);
#endif

    aZ = malloc((nx + 1) * sizeof *aZ);
    if (aZ == NULL) {
	return NULL;
    }

    /* the dependent variable */
    if (ainfo->y != NULL) {
	aZ[0] = ainfo->y;
    } else {
	aZ[0] = dset->Z[list[ypos]];
    }

    /* the independent variables */
    for (i=1; i<=nx; i++) {
	v = list[i + ypos];
	aZ[i] = dset->Z[v];
    }

    return aZ;
}

/* add extra OPG-related stuff to the arma info struct */

static int set_up_arma_OPG_info (arma_info *ainfo,
				 const DATASET *dset)
{
    /* array length needed for derivatives */
    int nd = 1 + ainfo->q + ainfo->pd * ainfo->Q;
    /* number of derivatives */
    int k = ainfo->nc;
    int err = 0;

    /* construct virtual dataset for dep var, real regressors */
    ainfo->Z = make_arma_Z(ainfo, dset);
    if (ainfo->Z == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	/* allocate gradient matrix */
	ainfo->G = gretl_zero_matrix_new(ainfo->T, k);
	if (ainfo->G == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && !arma_exact_ml(ainfo)) {
	/* allocate covariance matrix */
	ainfo->V = gretl_matrix_alloc(k, k);
	if (ainfo->V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* forecast errors array */
	ainfo->e = malloc((ainfo->t2 + 1) * sizeof *ainfo->e);
	if (ainfo->e == NULL) {
	    err = E_ALLOC;
	} else {
	    int t;

	    for (t=0; t<=ainfo->t2; t++) {
		ainfo->e[t] = 0.0;
	    }
	}
    }

    if (!err) {
	/* derivatives arrays */
	ainfo->aux = doubles_array_new0(k, nd);
	if (ainfo->aux == NULL) {
	    err = E_ALLOC;
	} else {
	    ainfo->n_aux = k;
	}
    }

    return err;
}

/* retrieve results specific to bhhh procedure */

static int
conditional_arma_model_prep (MODEL *pmod, arma_info *ainfo,
			     double *theta)
{
    int i, t, err;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->nobs = pmod->t2 - pmod->t1 + 1;
    pmod->ncoeff = ainfo->nc;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	return err;
    }

    pmod->lnL = ainfo->ll;
    pmod->sigma = NADBL; /* will be replaced */

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = theta[i];
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = ainfo->e[t];
    }

    err = gretl_model_write_vcv(pmod, ainfo->V);

    return err;
}

int bhhh_arma (double *theta, const DATASET *dset,
	       arma_info *ainfo, MODEL *pmod,
	       gretlopt opt)
{
    gretlopt bhhh_opt = OPT_NONE;
    double tol = libset_get_double(BHHH_TOLER);
    int err = 0;

    err = set_up_arma_OPG_info(ainfo, dset);
    if (err) {
	pmod->errcode = err;
	return err;
    }

    if (opt & OPT_V) {
	bhhh_opt |= OPT_V;
    }

    err = bhhh_max(theta, ainfo->nc, ainfo->G,
		   bhhh_arma_callback, tol,
		   &ainfo->fncount, &ainfo->grcount,
		   ainfo, ainfo->V, bhhh_opt, ainfo->prn);

    if (err) {
	fprintf(stderr, "arma: bhhh_max returned %d\n", err);
    } else {
	pmod->full_n = dset->n;
	err = conditional_arma_model_prep(pmod, ainfo, theta);
    }

    if (!err) {
	gretl_model_set_int(pmod, "fncount", ainfo->fncount);
	gretl_model_set_int(pmod, "grcount", ainfo->grcount);
	write_arma_model_stats(pmod, ainfo, dset);
	arma_model_add_roots(pmod, ainfo, theta);
    }

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }

    return pmod->errcode;
}

static int bhhh_arma_simple (double *theta, const DATASET *dset,
			     arma_info *ainfo, gretlopt opt)
{
    int err = set_up_arma_OPG_info(ainfo, dset);

    if (err) {
	return err;
    } else {
	gretlopt bhhh_opt = OPT_I; /* initializing */
	double tol = 1.0e-4;
	PRN *prn = NULL;

	if (opt & OPT_V) {
	    bhhh_opt |= OPT_V;
	    prn = ainfo->prn;
	    pprintf(prn, "%s\n\n", _("BHHH iteration"));
	}

	err = bhhh_max(theta, ainfo->nc, ainfo->G,
		       bhhh_arma_callback, tol,
		       &ainfo->fncount, &ainfo->grcount,
		       ainfo, NULL, bhhh_opt, prn);
    }

    return err;
}

/* In case of exact ML estimation, run some CML iterations
   to jump-start the process */

int cml_arma_init (double *theta, const DATASET *dset,
		   arma_info *ainfo, gretlopt opt)
{
    double *tmp = copyvec(theta, ainfo->nc);
    int save_t1 = ainfo->t1;
    int save_T = ainfo->T;
    int i, err = 0;

    /* temporarily remove the "exact" flag */
    ainfo->flags &= ~ARMA_EXACT;

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	if (ainfo->p > 0 || ainfo->P > 0) {
	    err = cml_init_revise_sample(ainfo, dset);
	}
	if (!err) {
	    err = bhhh_arma_simple(tmp, dset, ainfo, opt);
	}
	if (!err) {
	    for (i=0; i<ainfo->nc; i++) {
		theta[i] = tmp[i];
#if CML_DEBUG
		fprintf(stderr, "cml_init pre %g\n", theta[i]);
#endif
	    }
	}
	free(tmp);
    }

    /* Whether this succeeded or not, we should convert the
       constant (if present) back to the form wanted by
       exact ML estimation. In addition, if we're scaling
       the y data for exact ML we should rescale the constant
       correspondingly.
    */
    if (ainfo->ifc) {
	transform_arma_const(theta, ainfo);
	if (ainfo->yscale != 1.0) {
	    theta[0] *= ainfo->yscale;
	}
#if CML_DEBUG
	fprintf(stderr, "adjusted const %g\n", theta[0]);
#endif
    }

    ainfo->t1 = save_t1;
    ainfo->T = save_T;
    /* reinstate the "exact" flag */
    ainfo->flags |= ARMA_EXACT;

    /* for now we'll disregard any errors in here */

    return 0;
}

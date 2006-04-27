/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
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

/* Native gretl code for ARMA estimation.  Much of the code here
   was contributed by Riccardo "Jack" Lucchetti, the rest is due
   to Allin Cottrell; thanks also to Stephen Moshier for cephes.
*/

#include "../cephes/polrt.c"

#include "libgretl.h"
#include "bhhh_max.h"
#include "libset.h"
#include "kalman.h"

#define ARMA_DEBUG 0

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

#include "arma_common.c"

static PRN *errprn;

/* check whether the MA estimates have gone out of bounds in the
   course of iteration */

static int 
ma_out_of_bounds (struct arma_info *ainfo, const double *theta,
		  const double *Theta)
{
    static double *temp;
    static double *tmp2;
    static cmplx *roots;
    static int qmax;

    double re, im, rt;
    int i, j, k;
    int err = 0, allzero = 1;

    if (ainfo == NULL) {
	/* signal for cleanup */
	free(temp);
	temp = NULL;
	free(tmp2);
	tmp2 = NULL;
	free(roots);
	roots = NULL;
	qmax = 0;
	return 0;
    }

    for (i=0; i<ainfo->q && allzero; i++) {
	if (theta[i] != 0.0) {
	    allzero = 0;
	}    
    }  

    for (i=0; i<ainfo->Q && allzero; i++) {
	if (Theta[i] != 0.0) {
	    allzero = 0;
	}    
    }  
    
    if (allzero) {
	return 0;
    }

    if (temp == NULL) {
	/* not allocated yet */
	qmax = ainfo->q + ainfo->Q * ainfo->pd;

	temp  = malloc((qmax + 1) * sizeof *temp);
	tmp2  = malloc((qmax + 1) * sizeof *tmp2);
	roots = malloc(qmax * sizeof *roots);

	if (temp == NULL || tmp2 == NULL || roots == NULL) {
	    free(temp);
	    temp = NULL;
	    free(tmp2);
	    tmp2 = NULL;
	    free(roots);
	    roots = NULL;
	    return 1;
	}
    }

    temp[0] = 1.0;

    /* initialize to non-seasonal MA or zero */
    for (i=0; i<qmax; i++) {
        if (i < ainfo->q) {
            temp[i+1] = theta[i];
        } else {
            temp[i+1] = 0.0;
        }
    }

    /* add seasonal MA and interaction */
    for (i=0; i<ainfo->Q; i++) {
	k = (i + 1) * ainfo->pd;
	temp[k] += Theta[i];
	for (j=0; j<ainfo->q; j++) {
	    int m = k + j + 1;

	    temp[m] += Theta[i] * theta[j];
	}
    }

    polrt(temp, tmp2, qmax, roots);

    for (i=0; i<qmax; i++) {
	re = roots[i].r;
	im = roots[i].i;
	rt = re * re + im * im;
	if (rt > DBL_EPSILON && rt <= 1.0) {
	    pprintf(errprn, "MA root %d = %g\n", i, rt);
	    err = 1;
	    break;
	}
    }

    return err;
}

static void do_MA_partials (double *drv,
			    struct arma_info *ainfo,
			    const double *theta,
			    const double *Theta,
			    int t)
{
    int i, j, s, p;

    for (i=0; i<ainfo->q; i++) {
	s = t - (i + 1);
	if (s >= 0) {
	    drv[t] -= theta[i] * drv[s];
	}
    }

    for (i=0; i<ainfo->Q; i++) {
	s = t - ainfo->pd * (i + 1);
	if (s >= 0) {
	    drv[t] -= Theta[i] * drv[s];
	    for (j=0; j<ainfo->q; j++) {
		p = s - (j + 1);
		if (p >= 0) {
		    drv[t] -= Theta[i] * theta[j] * drv[p];
		}
	    }
	}
    }
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a "callback". */

static int arma_ll (double *coeff, 
		    const double **bhX, double **Z, 
		    model_info *arma,
		    int do_score)
{
    int i, j, s, t;
    int t1 = model_info_get_t1(arma);
    int t2 = model_info_get_t2(arma);
    int n = t2 - t1 + 1;

    const double *y = bhX[0];
    const double **X = bhX + 1;
    double **series = model_info_get_series(arma);
    double *e = series[0];
    double **de = series + 1;
    double **de_a, **de_sa, **de_m, **de_sm, **de_r;
    const double *phi, *Phi;
    const double *theta, *Theta;
    const double *beta;

    struct arma_info *ainfo;

    double ll, s2 = 0.0;
    int err = 0;

    /* retrieve ARMA-specific information */
    ainfo = model_info_get_extra_info(arma);

    /* pointers to blocks of coefficients */
    phi = coeff + ainfo->ifc;
    Phi = phi + ainfo->p;
    theta = Phi + ainfo->P;
    Theta = theta + ainfo->q;
    beta = Theta + ainfo->Q;

    /* pointers to blocks of derivatives */
    de_a = de + ainfo->ifc;
    de_sa = de_a + ainfo->p;
    de_m = de_sa + ainfo->P;
    de_sm = de_m + ainfo->q;
    de_r = de_sm + ainfo->Q;

#if ARMA_DEBUG
    fprintf(stderr, "arma_ll: p=%d, q=%d, P=%d, Q=%d, pd=%d\n",
	    ainfo->p, ainfo->q, ainfo->P, ainfo->Q, ainfo->pd);
#endif

    if (ma_out_of_bounds(ainfo, theta, Theta)) {
	pputs(errprn, "arma: MA estimate(s) out of bounds\n");
	fputs("arma: MA estimate(s) out of bounds\n", stderr);
	return 1;
    }

    /* update forecast errors */

    for (t=t1; t<=t2; t++) {
	int p;

	e[t] = y[t];

	/* intercept */
	if (ainfo->ifc) {
	    e[t] -= coeff[0];
	} 

	/* non-seasonal AR component */
	for (i=0; i<ainfo->p; i++) {
	    s = t - (i + 1);
	    e[t] -= phi[i] * y[s];
	}

	/* seasonal AR component plus interactions */
	for (i=0; i<ainfo->P; i++) {
	    s = t - ainfo->pd * (i + 1);
	    e[t] -= Phi[i] * y[s];
	    for (j=0; j<ainfo->p; j++) {
		p = s - (j + 1);
		e[t] += Phi[i] * phi[j] * y[p];
	    }
	}

	/* non-seasonal MA component */
	for (i=0; i<ainfo->q; i++) {
	    s = t - (i + 1);
	    if (s >= t1) {
		e[t] -= theta[i] * e[s];
	    }
	}

	/* seasonal MA component plus interactions */
	for (i=0; i<ainfo->Q; i++) {
	    s = t - ainfo->pd * (i + 1);
	    if (s >= t1) {
		e[t] -= Theta[i] * e[s];
		for (j=0; j<ainfo->q; j++) {
		    p = s - (j + 1);
		    if (p >= t1) {
			e[t] -= Theta[i] * theta[j] * e[p];
		    }
		}
	    }
	}

	/* exogenous regressors */
	for (i=0; i<ainfo->nexo; i++) {
	    e[t] -= beta[i] * X[i][t];
	}

	s2 += e[t] * e[t];
    }

    /* get error variance and log-likelihood */

    s2 /= (double) n;

    ll = -n * (0.5 * log(s2) + LN_SQRT_2_PI_P5);
    model_info_set_ll(arma, ll, do_score);

    if (do_score) {
	int lag, xlag;
	double x;

	for (t=t1; t<=t2; t++) {

	    /* the constant term (de_0) */
	    if (ainfo->ifc) {
		de[0][t] = -1.0;
		do_MA_partials(de[0], ainfo, theta, Theta, t);
	    }

	    /* non-seasonal AR terms (de_a) */
	    for (j=0; j<ainfo->p; j++) {
		lag = j + 1;
		if (t >= lag) {
		    de_a[j][t] = -y[t-lag];
		    /* cross-partial with seasonal AR */
		    for (i=0; i<ainfo->P; i++) {
			xlag = lag + ainfo->pd * (i + 1);
			if (t >= xlag) {
			    de_a[j][t] += Phi[i] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_a[j], ainfo, theta, Theta, t);
		}
	    }

	    /* seasonal AR terms (de_sa) */
	    for (j=0; j<ainfo->P; j++) {
		lag = ainfo->pd * (j + 1);
		if (t >= lag) {
		    de_sa[j][t] = -y[t-lag];
		    /* cross-partial with non-seasonal AR */
		    for (i=0; i<ainfo->p; i++) {
			xlag = lag + (i + 1);
			if (t >= xlag) {
			    de_sa[j][t] += phi[i] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_sa[j], ainfo, theta, Theta, t);
		}
	    }

	    /* non-seasonal MA terms (de_m) */
	    for (j=0; j<ainfo->q; j++) {
		lag = j + 1;
		if (t >= lag) {
		    de_m[j][t] = -e[t-lag];
		    /* cross-partial with seasonal MA */
		    for (i=0; i<ainfo->Q; i++) {
			xlag = lag + ainfo->pd * (i + 1);
			if (t >= xlag) {
			    de_m[j][t] -= Theta[i] * e[t-xlag];
			}
		    }
		    do_MA_partials(de_m[j], ainfo, theta, Theta, t);
		}
	    }

	    /* seasonal MA terms (de_sm) */
	    for (j=0; j<ainfo->Q; j++) {
		lag = ainfo->pd * (j + 1);
		if (t >= lag) {
		    de_sm[j][t] = -e[t-lag];
		    /* cross-partial with non-seasonal MA */
		    for (i=0; i<ainfo->q; i++) {
			xlag = lag + (i + 1);
			if (t >= xlag) {
			    de_sm[j][t] -= theta[i] * e[t-xlag];
			}
		    }
		    do_MA_partials(de_sm[j], ainfo, theta, Theta, t);
		}
	    }

	    /* exogenous regressors (de_r) */
	    for (j=0; j<ainfo->nexo; j++) {
		de_r[j][t] = -X[j][t]; 
		do_MA_partials(de_r[j], ainfo, theta, Theta, t);
	    }

	    /* update OPG data set */
	    x = e[t] / s2; /* sqrt(s2)? does it matter? */
	    for (i=0; i<ainfo->nc; i++) {
		Z[i+1][t] = -de[i][t] * x;
	    }
	}
    }

    if (isnan(ll)) {
	err = 1;
    }

    return err;
}

/*
  Given an ARMA process $A(L)B(L) y_t = C(L)D(L) \epsilon_t$, finds the 
  roots of the four polynomials -- or just two polynomials if seasonal
  AR and MA effects, B(L) and D(L) are not present -- and attaches
  this information to the ARMA model.

  pmod: MODEL pointer to which the roots info should be attached.

  ainfo: gives various pieces of information on the ARMA model,
  including seasonal and non-seasonal AR and MA orders.

  coeff: ifc + p + q + P + Q vector of coefficients (if an intercept
  is present it is element 0 and is ignored)

  returns: zero on success, non-zero on failure
*/

static int arma_model_add_roots (MODEL *pmod, struct arma_info *ainfo,
				 const double *coeff)
{
    const double *phi = coeff + ainfo->ifc;
    const double *Phi = phi + ainfo->p;
    const double *theta = Phi + ainfo->P;
    const double *Theta = theta + ainfo->q;

    int nr = ainfo->p + ainfo->P + ainfo->q + ainfo->Q;
    int pmax, qmax, lmax;
    double *temp = NULL, *temp2 = NULL;
    cmplx *rptr, *roots = NULL;
    int i;

    pmax = (ainfo->p > ainfo->P)? ainfo->p : ainfo->P;
    qmax = (ainfo->q > ainfo->Q)? ainfo->q : ainfo->Q;
    lmax = (pmax > qmax)? pmax : qmax;

    if (pmax == 0 && qmax == 0) {
	return 0;
    }

    temp  = malloc((lmax + 1) * sizeof *temp);
    temp2 = malloc((lmax + 1) * sizeof *temp2);
    roots = malloc(nr * sizeof *roots);

    if (temp == NULL || temp2 == NULL || roots == NULL) {
	free(temp);
	free(temp2);
	free(roots);
	return E_ALLOC;
    }

    temp[0] = 1.0;
    rptr = roots;

    if (ainfo->p > 0) {
	/* A(L), non-seasonal */
	for (i=0; i<ainfo->p; i++) {
	    temp[i+1] = -phi[i];
	}
	polrt(temp, temp2, ainfo->p, rptr);
	rptr += ainfo->p;
    }

    if (ainfo->P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo->P; i++) {
	    temp[i+1] = -Phi[i];
	}    
	polrt(temp, temp2, ainfo->P, rptr);
	rptr += ainfo->P;
    }

    if (ainfo->q > 0) {
	/* C(L), non-seasonal */
	for (i=0; i<ainfo->q; i++) {
	    temp[i+1] = theta[i];
	}  
	polrt(temp, temp2, ainfo->q, rptr);
	rptr += ainfo->q;
    }

    if (ainfo->Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo->Q; i++) {
	    temp[i+1] = Theta[i];
	}  
	polrt(temp, temp2, ainfo->Q, rptr);
    }
    
    free(temp);
    free(temp2);

    gretl_model_set_data(pmod, "roots", roots, MODEL_DATA_CMPLX_ARRAY,
			 nr * sizeof *roots);

    return 0;
}

/* below: exact ML using Kalman filter apparatus */

static gretl_matrix *S = NULL;
static gretl_matrix *P = NULL;
static gretl_matrix *F = NULL;
static gretl_matrix *A = NULL;
static gretl_matrix *H = NULL;
static gretl_matrix *Q = NULL;

static gretl_matrix *Tmp;
static gretl_matrix *vecP;
static gretl_matrix *vecQ;

static double *ac;
static double *mc;

static struct arma_info *kainfo;    

static int ainfo_get_r (struct arma_info *ainfo)
{
    int pmax = ainfo->p + ainfo->pd * ainfo->P;
    int qmax = ainfo->q + ainfo->pd * ainfo->Q;

    return (pmax > qmax + 1)? pmax : qmax + 1;
}

static int allocate_ac_mc (struct arma_info *ainfo)
{
    if (ainfo->P > 0) {
	int pmax = ainfo->p + ainfo->pd * ainfo->P;

	ac = malloc((pmax + 1) * sizeof *ac);
	if (ac == NULL) {
	    return E_ALLOC;
	}
    }

    if (ainfo->Q > 0) {
	int qmax = ainfo->q + ainfo->pd * ainfo->Q;

	mc = malloc((qmax + 1) * sizeof *mc);
	if (mc == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static void free_ac_mc (void)
{
    if (ac != NULL) free(ac);
    if (mc != NULL) free(mc);
}

static void write_big_phi (const double *phi, 
			   const double *Phi,
			   struct arma_info *ainfo,
			   gretl_matrix *F)
{
    int pmax = ainfo->p + ainfo->pd * ainfo->P;
    double x, y;
    int i, j, k;

    for (i=0; i<=pmax; i++) {
	ac[i] = 0.0;
    }

    for (i=0; i<=ainfo->P; i++) {
	x = (i == 0)? -1 : Phi[i-1];
	for (j=0; j<=ainfo->p; j++) {
	    y = (j == 0)? -1 : phi[j-1];
	    k = j + ainfo->pd * i;
	    ac[k] -= x * y;
	}
    }

    for (i=0; i<pmax; i++) {
	gretl_matrix_set(F, 0, i, ac[i+1]);
    }
}

static void write_big_theta (const double *theta, 
			     const double *Theta,
			     struct arma_info *ainfo,
			     gretl_matrix *H)
{
    int qmax = ainfo->q + ainfo->pd * ainfo->Q;
    double x, y;
    int i, j, k;

    for (i=0; i<=qmax; i++) {
	mc[i] = 0.0;
    }

    for (i=0; i<=ainfo->Q; i++) {
	x = (i == 0)? -1 : Theta[i-1];
	for (j=0; j<=ainfo->q; j++) {
	    y = (j == 0)? -1 : theta[j-1];
	    k = j + ainfo->pd * i;
	    mc[k] -= x * y;
	}
    }

    for (i=1; i<=qmax; i++) {
	gretl_vector_set(H, i, mc[i]);
    }    
}

static int write_kalman_matrices (const double *b)
{
    const double *phi = b + kainfo->ifc;
    const double *Phi = phi + kainfo->p;
    const double *theta = Phi + kainfo->P;
    const double *Theta = theta + kainfo->q;
    const double *beta = Theta + kainfo->Q;
    double s2 = *(beta + kainfo->nexo);
    double mu = (kainfo->ifc)? b[0] : 0.0;
    int i, r, err = 0;

    gretl_matrix_zero(S);
    gretl_matrix_zero(P);
    gretl_matrix_zero(F);
    gretl_matrix_zero(H);
    gretl_matrix_zero(Q);

    /* See Hamilton, Time Series Analysis, ch 13, p. 375 */

    r = gretl_matrix_rows(F);

    /* form the F matrix using phi and/or Phi */
    if (kainfo->P > 0) {
	write_big_phi(phi, Phi, kainfo, F);
    } else {
	for (i=0; i<kainfo->p; i++) {
	    gretl_matrix_set(F, 0, i, phi[i]);
	}
    } 
    gretl_matrix_inscribe_I(F, 1, 0, r - 1);

#if ARMA_DEBUG
    gretl_matrix_print(F, "F");
#endif

    /* form the H vector using theta and/or Theta */
    gretl_vector_set(H, 0, 1.0);
    if (kainfo->Q > 0) {
	write_big_theta(theta, Theta, kainfo, H);
    } else {
	for (i=0; i<kainfo->q; i++) {
	    gretl_vector_set(H, i + 1, theta[i]);
	}
    }

#if ARMA_DEBUG
    gretl_matrix_print(H, "H");
#endif

    gretl_matrix_set(Q, 0, 0, s2);
    gretl_matrix_vectorize(vecQ, Q);

#if ARMA_DEBUG
    gretl_matrix_print(Q, "Q");
#endif

    gretl_vector_set(A, 0, mu);
    for (i=0; i<kainfo->nexo; i++) {
	gretl_vector_set(A, i + 1, beta[i]);
    }

#if ARMA_DEBUG
    gretl_matrix_print(A, "A");
#endif

    /* form $P_{1|0}$ (MSE) matrix, as per Hamilton, ch 13, p. 378.
       Is there a cheaper way of doing this?
    */
    gretl_matrix_kronecker_product(F, F, Tmp);
    gretl_matrix_I_minus(Tmp);
    err = gretl_invert_general_matrix(Tmp);
    if (!err) {
	gretl_matrix_multiply(Tmp, vecQ, vecP);
	gretl_matrix_unvectorize(P, vecP);
    }  

#if ARMA_DEBUG
    gretl_matrix_print(P, "P");
#endif

    return err;
}

static int rewrite_kalman_matrices (kalman *K, const double *b)
{
    int err = write_kalman_matrices(b);

    if (!err) {
	kalman_set_initial_state_vector(K, S);
	kalman_set_initial_MSE_matrix(K, P);
    }

    return err;
}

/* add innovations based on Kalman forecast, and also covariance
   matrix and standard errors based on Outer Product of the
   estimated innovations */

static int arma_OPG_stderrs (MODEL *pmod, kalman *K, double *b, int m, int T)
{
    gretl_matrix *E = NULL;
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    const double eps = 1.0e-8;
    double g0, g1;
    double x = 0.0;
    int k = m - 1;
    int i, j, s, t;
    int err = 0;

    E = gretl_column_vector_alloc(T);
    G = gretl_matrix_alloc(k, T);
    V = gretl_matrix_alloc(k, k);
    if (E == NULL || G == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* get estimate of innovations */
    kalman_forecast(K, E);
    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = gretl_vector_get(E, s++);
    }

    /* construct approximation to G matrix based
       on finite differences */
    for (i=0; i<k && !err; i++) {
	b[i] -= eps;
	rewrite_kalman_matrices(K, b);
	err = kalman_forecast(K, E);
	for (t=0; t<T; t++) {
	    g0 = gretl_vector_get(E, t);
	    gretl_matrix_set(G, i, t, g0);
	}
	b[i] += 2.0 * eps;
	rewrite_kalman_matrices(K, b);
	err = kalman_forecast(K, E);
	for (t=0; t<T; t++) {
	    g1 = gretl_vector_get(E, t);
	    g0 = gretl_matrix_get(G, i, t);
	    gretl_matrix_set(G, i, t, (g1 - g0) / (2.0 * eps));
	}
	b[i] -= eps;
    }

    gretl_matrix_multiply_mod(G, GRETL_MOD_NONE,
			      G, GRETL_MOD_TRANSPOSE,
			      V);

    err = gretl_invert_symmetric_matrix(V);

    if (!err) {
	int idx;

	for (i=0; i<k; i++) {
	    for (j=0; j<=i; j++) {
		idx = ijton(i, j, k);
		x = b[k] * gretl_matrix_get(V, i, j);
		pmod->vcv[idx] = x;
		if (i == j) {
		    pmod->sderr[i] = sqrt(x);
		}
	    }
	}
    }

 bailout:

    gretl_matrix_free(E);
    gretl_matrix_free(G);
    gretl_matrix_free(V);
    
    return err;
}

/* in Kalman case the basic model struct is empty, so we have
   to allocate for coefficients, residuals and so on */

static int kalman_arma_model_allocate (MODEL *pmod, int k, int T)
{
    int t;

    pmod->ncoeff = k;

    pmod->coeff = malloc(k * sizeof *pmod->coeff);
    if (pmod->coeff == NULL) {
	return E_ALLOC;
    }

    pmod->sderr = malloc(k * sizeof *pmod->sderr);
    if (pmod->sderr == NULL) {
	return E_ALLOC;
    }

    if (gretl_model_new_vcv(pmod, NULL)) {
	return E_ALLOC;
    }

    pmod->uhat = malloc(T * sizeof *pmod->uhat);
    if (pmod->uhat == NULL) {
	return E_ALLOC;
    }    

    pmod->yhat = malloc(T * sizeof *pmod->yhat);
    if (pmod->yhat == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<T; t++) {
	pmod->uhat[t] = pmod->yhat[t] = NADBL;
    }

    return 0;
}

static int kalman_arma_finish (MODEL *pmod, const int *alist,
			       struct arma_info *ainfo,
			       const double **Z, const DATAINFO *pdinfo, 
			       kalman *K, double *b, int m, int T)
{
    int i, k = m - 1;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->nobs = T;

    pmod->errcode = kalman_arma_model_allocate(pmod, k, ainfo->T);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    for (i=0; i<k; i++) {
	pmod->coeff[i] = b[i];
    }

    pmod->sigma = sqrt(b[k]);
    pmod->errcode = arma_OPG_stderrs(pmod, K, b, m, T);

    pmod->lnL = kalman_get_loglik(K);

    write_arma_model_stats(pmod, alist, ainfo, Z, pdinfo);
    arma_model_add_roots(pmod, ainfo, b);

    gretl_model_set_int(pmod, "arma_flags", ARMA_EXACT);

    return pmod->errcode;
}

static double kalman_arma_ll (const double *b, void *p)
{
    int offset = kainfo->ifc + kainfo->p + kainfo->P;
    const double *theta = b + offset;
    const double *Theta = theta + kainfo->q;
    double ll = NADBL;
    kalman *K;
    int err = 0;

#if ARMA_DEBUG
    int i;
    fprintf(stderr, "kalman_arma_ll():\n");
    for (i=0; i<kainfo->q; i++) {
	fprintf(stderr, "theta[%d] = %g\n", i, theta[i]);
    }
    for (i=0; i<kainfo->Q; i++) {
	fprintf(stderr, "Theta[%d] = %g\n", i, Theta[i]);
    }    
#endif

    if (ma_out_of_bounds(kainfo, theta, Theta)) {
	pputs(errprn, "arma: MA estimate(s) out of bounds\n");
	fputs("arma: MA estimate(s) out of bounds\n", stderr);
	return NADBL;
    }

    K = (kalman *) p;
    err = rewrite_kalman_matrices(K, b);
    if (!err) {
	err = kalman_forecast(K, NULL);
	ll = kalman_get_loglik(K);
    }

#if ARMA_DEBUG
    fprintf(stderr, "loglik = %g\n", ll);
#endif

    return ll;
}

static int kalman_arma_gradient (double *b, double *g, void *p)
{
    kalman *K = (kalman *) p;
    double eps = 1.0e-8;
    double bi0, ll1, ll2;
    int i, k, err = 0;

#if ARMA_DEBUG
    fprintf(stderr, "kalman_gradient_ll():\n\n");
#endif

    k = kalman_get_ncoeff(K);

    for (i=0; i<k && !err; i++) {
	ll1 = ll2 = 0.0;
	bi0 = b[i];
	b[i] -= eps;
#if ARMA_DEBUG
	fprintf(stderr, "trying b[%d] = %.15g\n", i, b[i]);
#endif
	err = rewrite_kalman_matrices(K, b);
	if (!err) {
	    err = kalman_forecast(K, NULL);
	    ll1 = kalman_get_loglik(K);
	}
	if (err) {
	    b[i] = bi0;
	    break;
	}
	b[i] = bi0 + eps;
#if ARMA_DEBUG
	fprintf(stderr, "trying b[%d] = %.15g\n", i, b[i]);
#endif
	err = rewrite_kalman_matrices(K, b);
	if (!err) {
	    err = kalman_forecast(K, NULL);
	    ll2 = kalman_get_loglik(K);
	}
	if (err) {
	    b[i] = bi0;
	    break;
	}
	b[i] = bi0; /* reset to original value */
	g[i] = (ll2 - ll1) / (2.0 * eps); /* sign? */
#if ARMA_DEBUG
	fprintf(stderr, "ll2 = %g, ll1 = %g, ll2 - ll1 = %g, g[%d] = %.8g\n", 
		ll2, ll1, ll2 - ll1, i, g[i]);
#endif
    }

    return err;
}

static gretl_matrix *form_arma_y_vector (const int *alist, 
					 const double **Z,
					 struct arma_info *ainfo,
					 int *err)
{
    gretl_matrix *yvec;
    const double *y;
    int s, t, T;

#if ARMA_DEBUG
    fprintf(stderr, "ainfo->t1 = %d, ainfo->t2 = %d\n",
	    ainfo->t1, ainfo->t2);
#endif

    if (ainfo->dy != NULL) {
	y = ainfo->dy;
	/* should be handled earlier? */
	for (t=ainfo->t1; t<=ainfo->t2; t++) {
	    if (na(y[t])) {
		ainfo->t1 += 1;
	    } else {
		break;
	    }
	}
    } else {
	y = Z[ainfo->yno];
    }

    T = ainfo->t2 - ainfo->t1 + 1;
    if (T == 0) {
	*err = E_DATA;
	return NULL;
    }

    yvec = gretl_column_vector_alloc(T);
    if (yvec == NULL) {
	*err = E_ALLOC;
	return NULL;
    }    

    s = 0;
    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	if (na(y[t])) {
	    *err = E_DATA;
	}
	gretl_vector_set(yvec, s++, y[t]);
    }

#if ARMA_DEBUG
    gretl_matrix_print(yvec, "y");
    fprintf(stderr, "y has %d rows\n", gretl_matrix_rows(yvec));
#endif

    if (*err) {
	gretl_matrix_free(yvec);
	yvec = NULL;
    }
    
    return yvec;
}

static gretl_matrix *form_arma_x_matrix (const int *alist, 
					 const double **Z,
					 struct arma_info *ainfo)
{
    gretl_matrix *x;
    int i, xstart;
    int *xlist;

    xlist = gretl_list_new(ainfo->nexo);
    if (xlist == NULL) {
	return NULL;
    }

    xstart = arma_list_y_position(ainfo) + 1;
    for (i=xstart; i<=alist[0]; i++) {
	xlist[i - xstart + 1] = alist[i];
    }

#if ARMA_DEBUG
    printlist(alist, "alist (arma list)");
    printlist(xlist, "xlist (exog vars)");
#endif

    x = gretl_matrix_data_subset(xlist, Z, ainfo->t1, ainfo->t2, NULL);
    if (x == NULL) {
	free(xlist);
	return NULL;
    }

#if ARMA_DEBUG
    gretl_matrix_print(x, "x");
    fprintf(stderr, "x has %d rows\n", gretl_matrix_rows(x));
#endif

    free(xlist);
    
    return x;
}

/* Given an estimate of the ARMA constant via OLS, convert to the form
   wanted for initializing the Kalman filter
*/

static void transform_arma_const (double *b, struct arma_info *ainfo)
{
    const double *phi = b + 1;
    const double *Phi = phi + ainfo->p;
    double narfac = 1.0;
    double sarfac = 1.0;
    int i;

    for (i=0; i<ainfo->p; i++) {
	narfac -= phi[i];
    }

    for (i=0; i<ainfo->P; i++) {
	sarfac -= Phi[i];
    }

    b[0] /= (narfac * sarfac);
}

static int kalman_arma (const int *alist, double *coeff, double s2,
			const double **Z, const DATAINFO *pdinfo,
			struct arma_info *ainfo, MODEL *pmod,
			PRN *prn)
{
    kalman *K = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *x = NULL;

    int k = 1 + ainfo->nexo; /* number of exog vars plus space for const */
    int r;

    /* BFGS apparatus */
    int ncoeff = ainfo->nc + 1;
    int maxit = 1000;
    double reltol = 1.0e-12;
    int fncount = 0;
    int grcount = 0;

    double *b;
    int i, T, err = 0;

    b = malloc(ncoeff * sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<ainfo->nc; i++) {
	b[i] = coeff[i];
    }
    b[ainfo->nc] = s2;

#if ARMA_DEBUG
    for (i=0; i<ncoeff; i++) {
	fprintf(stderr, "initial b[%d] = %g\n", i, b[i]);
    }
#endif

    y = form_arma_y_vector(alist, Z, ainfo, &err);
    if (y == NULL) {
	free(b);
	return err;
    }

    if (ainfo->nexo > 0) {
	x = form_arma_x_matrix(alist, Z, ainfo);
	if (x == NULL) {
	    free(b);
	    gretl_matrix_free(y);
	    return E_ALLOC;
	}
    }

    if (allocate_ac_mc(ainfo)) {
	free(b);
	gretl_matrix_free(y);
	gretl_matrix_free(x);
	return E_ALLOC;
    }

    r = ainfo_get_r(ainfo);
    T = gretl_matrix_rows(y);

    S = gretl_column_vector_alloc(r);
    P = gretl_matrix_alloc(r, r);
    F = gretl_matrix_alloc(r, r);
    A = gretl_column_vector_alloc(k);
    H = gretl_column_vector_alloc(r);
    Q = gretl_matrix_alloc(r, r);

    Tmp = gretl_matrix_alloc(r * r, r * r);
    vecQ = gretl_column_vector_alloc(r * r);
    vecP = gretl_column_vector_alloc(r * r);

    if (S == NULL || P == NULL || F == NULL || A == NULL ||
	H == NULL || Q == NULL || Tmp == NULL || 
	vecP == NULL || vecQ == NULL) {
	free(b);
	gretl_matrix_free(y);
	gretl_matrix_free(x);
	return E_ALLOC;
    }

#if ARMA_DEBUG
    fprintf(stderr, "ready to estimate: ainfo specs:\n"
	    "p=%d, P=%d, q=%d, Q=%d, ifc=%d, nexo=%d, t1=%d, t2=%d\n", 
	    ainfo->p, ainfo->P, ainfo->q, ainfo->Q, ainfo->ifc, 
	    ainfo->nexo, ainfo->t1, ainfo->t2);
    fprintf(stderr, "Kalman dims: r = %d, k = %d, T = %d, ncoeff=%d\n", 
	    r, k, T, ncoeff);
#endif

    /* publish ainfo */
    kainfo = ainfo;

    K = kalman_new(S, P, F, A, H, Q, NULL, y, x, ncoeff, ainfo->ifc, &err);

    if (err) {
	fprintf(stderr, "kalman_new(): err = %d\n", err);
    } else {
	err = BFGS_max(ncoeff, b, maxit, reltol, &fncount, &grcount,
		       kalman_arma_ll, kalman_arma_gradient, K,
		       (prn != NULL)? OPT_V : OPT_NONE, prn);
	fprintf(stderr, "BFGS_max returned %d\n", err);
    }

    if (err) {
	pmod->errcode = err;
    } else {
	kalman_arma_finish(pmod, alist, ainfo, Z, pdinfo, 
			   K, b, ncoeff, T);
    } 

    kalman_free(K);

    gretl_matrix_free(S);
    gretl_matrix_free(P);
    gretl_matrix_free(F);
    gretl_matrix_free(A);
    gretl_matrix_free(H);
    gretl_matrix_free(Q);

    gretl_matrix_free(y);
    gretl_matrix_free(x);

    gretl_matrix_free(Tmp);
    gretl_matrix_free(vecP);
    gretl_matrix_free(vecQ);

    free(b);
    free_ac_mc();

    /* unpublish ainfo */
    kainfo = NULL;

    return err;
}

/* end of Kalman-specific material */

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

static const double **
make_armax_X (const int *list, struct arma_info *ainfo, const double **Z)
{
    const double **X;
    int ypos, nx;
    int v, i;

    ypos = arma_list_y_position(ainfo);
    nx = list[0] - ypos;

#if ARMA_DEBUG
    fprintf(stderr, "make_armax_X: allocating %d series pointers\n",
	    nx + 1);
#endif    

    X = malloc((nx + 1) * sizeof *X);
    if (X == NULL) {
	return NULL;
    }

    /* the dependent variable */
    if (ainfo->dy != NULL) {
	X[0] = ainfo->dy;
    } else {
	X[0] = Z[list[ypos]];
    }

    /* the independent variables */
    for (i=1; i<=nx; i++) {
	v = list[i + ypos];
	X[i] = Z[v];
    }

    return X;
}

static void y_Xb_at_lag (char *spec, struct arma_info *ainfo, 
			 int narmax, int lag)
{
    char term[32];
    int i, nt;

    if (narmax == 0) {
	sprintf(term, "y_%d", lag);
	strcat(spec, term);
	return;
    }

    nt = ainfo->ifc + narmax;

    sprintf(term, "(y_%d-", lag);
    strcat(spec, term);

    if (nt > 1) {
	strcat(spec, "(");
    }

    if (ainfo->ifc) {
	strcat(spec, "b0");
    }

    for (i=1; i<=narmax; i++) {
	if (ainfo->ifc || i > 1) {
	    strcat(spec, "+");
	} 
	sprintf(term, "b%d*x%d_%d", i, i, lag);
	strcat(spec, term);
    }

    if (nt > 1) {
	strcat(spec, "))");
    } else {
	strcat(spec, ")");
    }
}

static int arma_get_nls_model (MODEL *amod, struct arma_info *ainfo,
			       int narmax, double ***pZ, DATAINFO *pdinfo) 
{
#if ARMA_DEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
    gretlopt opt = OPT_A | OPT_V;
#else
    PRN *prn = NULL;
    gretlopt opt = OPT_A;
#endif
    char fnstr[MAXLINE];
    char term[32];
    nls_spec *spec;
    int *plist = NULL;
    int v, oldv = pdinfo->v;
    int nparam;
    int i, j, k, err = 0;

    spec = nls_spec_new(NLS, pdinfo);
    if (spec == NULL) {
	return E_ALLOC;
    }

    nparam = ainfo->ifc + ainfo->p + ainfo->P + ainfo->nexo;

    plist = gretl_list_new(nparam);
    if (plist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = dataset_add_scalars(nparam, pZ, pdinfo); 
    if (err) {
	goto bailout;
    }

    /* FIXME initialization of NLS parameters? */

    /* construct names for the parameters, and param list */
    v = oldv;
    k = 1;
    if (ainfo->ifc) {
	(*pZ)[v][0] = gretl_mean(0, pdinfo->n - 1, (*pZ)[1]);
	strcpy(pdinfo->varname[v], "b0");
	plist[k++] = v++;
    }
    for (i=1; i<=ainfo->p; i++) {
	if (i == 1) {
	    (*pZ)[v][0] = 0.1; /* ?? */
	}
	sprintf(pdinfo->varname[v], "phi_%d", i);
	plist[k++] = v++;
    }
    for (i=1; i<=ainfo->P; i++) {
	sprintf(pdinfo->varname[v], "Phi_%d", i);
	plist[k++] = v++;
    }
    for (i=1; i<=ainfo->nexo; i++) {
	sprintf(pdinfo->varname[v], "b%d", i);
	plist[k++] = v++;
    }

    /* construct NLS specification */

    strcpy(fnstr, "y=");

    if (ainfo->ifc) {
	strcat(fnstr, "b0");
    } else {
	strcat(fnstr, "0");
    } 

    for (i=0; i<=ainfo->p; i++) {
	if (i > 0) {
	    sprintf(term, "+phi_%d*", i);
	    strcat(fnstr, term);
	    y_Xb_at_lag(fnstr, ainfo, narmax, i);
	}
	for (j=0; j<=ainfo->P; j++) {
	    if (i == 0 && j > 0) {
		sprintf(term, "+Phi_%d*", j);
		strcat(fnstr, term);
		y_Xb_at_lag(fnstr, ainfo, narmax, j * ainfo->pd);
	    }
	    if (i > 0 && j > 0) {
		sprintf(term, "-phi_%d*Phi_%d*", i, j);
		strcat(fnstr, term);
		y_Xb_at_lag(fnstr, ainfo, narmax, j * ainfo->pd + i);
	    }
	}
    } 

    for (i=1; i<=ainfo->nexo; i++) {
	sprintf(term, "+b%d*x%d", i, i);
	strcat(fnstr, term);
    }	

    err = nls_spec_set_regression_function(spec, fnstr, pdinfo);

    if (!err) {
	err = nls_spec_add_param_list(spec, plist, (const double **) *pZ,
				      pdinfo);
    }

    if (!err) {
	*amod = model_from_nls_spec(spec, pZ, pdinfo, opt, prn);
	err = amod->errcode;
#if ARMA_DEBUG
	if (!err) {
	    printmodel(amod, pdinfo, OPT_NONE, prn);
	}
	gretl_print_destroy(prn);
#endif
    }

    bailout:

    nls_spec_destroy(spec);
    free(plist);

    return err;
}

static int *make_ar_ols_list (struct arma_info *ainfo, int av, int ptotal)
{
    int * alist = gretl_list_new(av);
    int i, offset;

    if (alist == NULL) {
	return NULL;
    }

    alist[1] = 1;

    if (ainfo->ifc) {
	alist[2] = 0;
	offset = 2;
    } else {
	alist[0] -= 1;
	offset = 1;
    }

    for (i=1; i<=ainfo->p; i++) {
	alist[offset + i] = 1 + i;
    }

    for (i=1; i<=ainfo->P; i++) {
	alist[offset + ainfo->p + i] = ainfo->p + 1 + i;
    }

    offset += ptotal;  

    for (i=1; i<=ainfo->nexo; i++) {
	alist[offset + i] = ptotal + 1 + i;
    }

    return alist;
}

/* Run a least squares model to get initial values for the AR
   coefficients */

static int ar_init_by_ls (const int *list, double *coeff, double *s2,
			  const double **Z, const DATAINFO *pdinfo,
			  struct arma_info *ainfo)
{
    int an = pdinfo->t2 - ainfo->t1 + 1;
    int np = ainfo->p, nq = ainfo->q;
    int nmixed = ainfo->p * ainfo->P;
    int ptotal = ainfo->p + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    const double *y;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int nonlin = 0;
    int narmax = 0;
    int xstart, lag;
    int axi = 0, ayi = 0;
    int i, j, k, t;
    int err = 0;

    /* dependent variable */
    if (ainfo->dy != NULL) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    if (ptotal == 0 && ainfo->nexo == 0 && !ainfo->ifc) {
	/* special case of pure MA model */
	for (i=0; i<nq; i++) {
	    coeff[i] = 0.0; 
	} 
	if (s2 != NULL) {
	    *s2 = gretl_variance(ainfo->t1, ainfo->t2, y);
	}
	return 0;
    }

    if (arma_exact_ml(ainfo)) {
	narmax = ainfo->nexo;
	if (narmax > 0) {
	    /* ARMAX-induced lags of exog vars */
	    av += ainfo->nexo * ptotal;
	}
    }

    gretl_model_init(&armod); 

    adinfo = create_new_dataset(&aZ, av, an, 0);
    if (adinfo == NULL) {
	return 1;
    }

    if (narmax > 0 || nmixed > 0) {
	/* have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	alist = make_ar_ols_list(ainfo, av, ptotal);
    }

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    /* construct the variable names */

    strcpy(adinfo->varname[1], "y");

    axi = ptotal + ainfo->nexo + 2;

    for (i=1; i<=ainfo->p; i++) {
	sprintf(adinfo->varname[i+1], "y_%d", i);
	for (j=1; j<=narmax; j++) {
	    sprintf(adinfo->varname[axi++], "x%d_%d", j, i);
	}
    }

    ayi = ainfo->p + ainfo->P + 2;

    for (i=1; i<=ainfo->P; i++) {
	k = ainfo->p + 1 + i;
	sprintf(adinfo->varname[k], "y_%d", ainfo->pd * i);
	for (j=1; j<=narmax; j++) {
	    sprintf(adinfo->varname[axi++], "x%d_%d", j, ainfo->pd * i);
	}
	for (j=1; j<=ainfo->p; j++) {
	    sprintf(adinfo->varname[ayi++], "y_%d", ainfo->pd * i + j);
	    for (k=1; k<=narmax; k++) {
		sprintf(adinfo->varname[axi++], "x%d_%d", k, ainfo->pd * i + j);
	    }
	}
    }

    axi = ptotal + 2;

    for (i=1; i<=ainfo->nexo; i++) {
	sprintf(adinfo->varname[axi++], "x%d", i);
    }

    /* build temporary dataset including lagged vars */

    for (t=0; t<an; t++) {
	int s, m;

	aZ[1][t] = y[t + ainfo->t1];

	axi = ptotal + ainfo->nexo + 2;

	for (i=1; i<=ainfo->p; i++) {
	    s = t + ainfo->t1 - i;
	    aZ[i+1][t] = (s >= 0)? y[s] : NADBL;
	    for (j=1; j<=narmax; j++) {
		m = list[xstart + j - 1];
		aZ[axi++][t] = (s >= 0)? Z[m][s] : NADBL;
	    }
	}

	ayi = ainfo->p + ainfo->P + 2;

	for (i=1; i<=ainfo->P; i++) {
	    lag = ainfo->pd * i;
	    s = t + ainfo->t1 - lag;
	    k = ainfo->p + 1 + i;
	    aZ[k][t] = (s >= 0)? y[s] : NADBL;
	    for (j=1; j<=ainfo->p; j++) {
		lag = ainfo->pd * i + j;
		s = t + ainfo->t1 - lag;
		aZ[ayi++][t] = (s >= 0)? y[s] : NADBL;
		for (k=1; k<=narmax; k++) {
		    m = list[xstart + k - 1];
		    aZ[axi++][t] = (s >= 0)? Z[m][s] : NADBL;
		}
	    }
	}

	s = t + ainfo->t1;
	axi = ptotal + 2;

	for (i=1; i<=ainfo->nexo; i++) {
	    m = list[xstart + i - 1];
	    aZ[axi++][t] = Z[m][s];
	}
    }

    if (arma_has_seasonal(ainfo)) {
	np += ainfo->P;
	nq += ainfo->Q;
    }

#if ARMA_DEBUG
    printlist(alist, "'alist' in ar_init_by_ls");
#endif

    if (nonlin) {
	err = arma_get_nls_model(&armod, ainfo, narmax, &aZ, adinfo);
    } else {
	/* just use OLS */
	armod = lsq(alist, &aZ, adinfo, OLS, OPT_A | OPT_Z);
	err = armod.errcode;
    }

    if (!err) {
	j = 0;
	for (i=0; i<armod.ncoeff; i++) {
	    if (i == np + ainfo->ifc) {
		j += nq; /* reserve space for MA coeffs */
	    }
	    coeff[j++] = armod.coeff[i];
	}
	for (i=0; i<nq; i++) {
	    /* insert zeros for MA coeffs */
	    coeff[i + np + ainfo->ifc] = 0.0;
	} 
	if (s2 != NULL) {
	    *s2 = armod.sigma * armod.sigma;
	}
    }

    if (!err && arma_exact_ml(ainfo) && ainfo->ifc) {
	if (!nonlin || ainfo->nexo == 0) {
	    transform_arma_const(coeff, ainfo);
	}
    }

#if ARMA_DEBUG
    if (!err) {
	fprintf(stderr, "LS init: ncoeff = %d, nobs = %d\n", 
		armod.ncoeff, armod.nobs);
	for (i=0; i<armod.ncoeff; i++) {
	    fprintf(stderr, " coeff[%d] = %g\n", i, armod.coeff[i]);
	}
    } else {
	fprintf(stderr, "LS init: armod.errcode = %d\n", err);
    }
#endif

    /* clean up */
    free(alist);
    destroy_dataset(aZ, adinfo);
    clear_model(&armod);

    return err;
}

/* set up a model_info struct for passing to bhhh_max */

static model_info *
set_up_arma_model_info (struct arma_info *ainfo)
{
    double tol = get_bhhh_toler();
    model_info *arma;

    if (na(tol)) {
	tol = 1.0e-6;
    }

    arma = model_info_new(ainfo->nc, ainfo->t1, ainfo->t2, ainfo->T, tol);

    if (arma == NULL) return NULL;

    model_info_set_opts(arma, PRESERVE_OPG_MODEL);
    model_info_set_n_series(arma, ainfo->nc + 1);

    /* add pointer to ARMA-specific details */
    model_info_set_extra_info(arma, ainfo);

    return arma;
}

/* retrieve results specific to bhhh procedure */

static void 
conditional_arma_model_prep (MODEL *pmod, model_info *minfo,
			     double *theta)
{
    double **series;
    int i, t;

    pmod->lnL = model_info_get_ll(minfo);

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = theta[i];
    }

    series = model_info_get_series(minfo);
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = series[0][t];
    }

    pmod->sigma = NADBL; /* will be replaced */
}

static int bhhh_arma (const int *alist, double *coeff, 
		      const double **Z, const DATAINFO *pdinfo,
		      struct arma_info *ainfo, MODEL *pmod,
		      PRN *prn)
{
    model_info *minfo = NULL;
    const double **X = NULL;
    int err = 0;

    /* construct virtual dataset for dep var, real regressors */
    X = make_armax_X(alist, ainfo, Z);
    if (X == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    /* create model_info struct to feed to bhhh_max() */
    minfo = set_up_arma_model_info(ainfo);
    if (minfo == NULL) {
	pmod->errcode = E_ALLOC;
	free(X);
	return pmod->errcode;
    }

    /* call BHHH conditional ML function (OPG regression) */
    err = bhhh_max(arma_ll, X, coeff, minfo, prn);
    
    if (err) {
	fprintf(stderr, "arma: bhhh_max returned %d\n", err);
	pmod->errcode = E_NOCONV;
    } else {
	MODEL *omod = model_info_capture_OPG_model(minfo);
	double *theta = model_info_get_theta(minfo);

	conditional_arma_model_prep(omod, minfo, theta);
	write_arma_model_stats(omod, alist, ainfo, Z, pdinfo);
	arma_model_add_roots(omod, ainfo, theta);
	*pmod = *omod;
	free(omod);
    }

    free(X);
    model_info_free(minfo);

    return pmod->errcode;
}

MODEL arma_model (const int *list, const double **Z, const DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    double *coeff = NULL;
    int *alist = NULL;
    PRN *aprn = NULL;
    MODEL armod;
    double s2 = 0.0;
    struct arma_info ainfo;
    char flags = 0;
    int err = 0;

    if (!(opt & OPT_C)) {
	flags = ARMA_EXACT;
    }

    if (opt & OPT_V) {
	aprn = prn;
	errprn = prn;
    } else {
	errprn = NULL;
    }

    arma_info_init(&ainfo, flags, pdinfo);
    gretl_model_init(&armod); 
    gretl_model_smpl_init(&armod, pdinfo);

    alist = gretl_list_copy(list);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    err = arma_check_list(alist, opt, Z, pdinfo, &ainfo);
    if (err) {
	armod.errcode = err;
	goto bailout;
    } 

    /* calculate maximum lag */
    calc_max_lag(&ainfo);

    /* adjust sample range if need be */
    if (arma_adjust_sample(pdinfo, Z, alist, &ainfo)) {
        armod.errcode = E_DATA;
	goto bailout;
    }

    /* allocate initial coefficient vector */
    coeff = malloc(ainfo.nc * sizeof *coeff);
    if (coeff == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    /* create differenced series if needed */
    if (ainfo.d > 0 || ainfo.D > 0) {
	err = arima_difference(Z[ainfo.yno], &ainfo);
    }

    /* initialize the coefficients: AR and regression part by least
       squares, MA at 0 */
    err = ar_init_by_ls(alist, coeff, &s2, Z, pdinfo, &ainfo);
    if (err) {
	armod.errcode = err;
	goto bailout;
    }

    if (flags & ARMA_EXACT) {
	kalman_arma(alist, coeff, s2, Z, pdinfo, &ainfo, &armod, aprn);
    } else {
	bhhh_arma(alist, coeff, Z, pdinfo, &ainfo, &armod, aprn);
    }

 bailout:

    free(alist);
    free(coeff);
    free(ainfo.dy);

    /* cleanup in MA roots checker */
    ma_out_of_bounds(NULL, NULL, NULL);

    errprn = NULL;

    return armod;
}

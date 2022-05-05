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

#define FULL_XML_HEADERS 1

#include "libgretl.h"
#include "uservar.h"
#include "gretl_func.h"
#include "gretl_xml.h"
#include "matrix_extra.h"
#include "libset.h"
#include "kalman.h"

/**
 * SECTION:kalman
 * @short_description: The Kalman filter
 * @title: Kalman
 * @include: gretl/libgretl.h, gretl/kalman.h
 *
 */

#define KDEBUG 0
#define EXACT_DEBUG 0
#define EXACT_SMDIST 1

/* try using incomplete observations? */
#define USE_INCOMPLETE_OBS 1

#define K_TINY 1.0e-15

enum {
    KALMAN_USER    = 1 << 0, /* user-defined filter? */
    KALMAN_DIFFUSE = 1 << 1, /* using diffuse P_{1|0} */
    KALMAN_FORWARD = 1 << 2, /* running forward filtering pass */
    KALMAN_SMOOTH  = 1 << 3, /* preparing for smoothing pass */
    KALMAN_SIM     = 1 << 4, /* running simulation */
    KALMAN_CHECK   = 1 << 5, /* checking user-defined matrices */
    KALMAN_BUNDLE  = 1 << 6, /* kalman is inside a bundle */
    KALMAN_SSFSIM  = 1 << 7, /* on simulation, emulate SsfPack */
    KALMAN_ARMA_LL = 1 << 8, /* filtering for ARMA estimation */
    KALMAN_SM_AM   = 1 << 9, /* Anderson-Moore smoothing */
    KALMAN_UNI     = 1 << 10 /* using univariate representation */
};

enum {
    STD_VAR, /* standard, basic variance representation */
    DJ_VAR,  /* as per DeJong, supporting cross correlation */
    DK_VAR   /* as per Durbin and Koopman */
};

typedef struct stepinfo_ stepinfo;

struct stepinfo_ {
    gretl_matrix *T;   /* N x (r * r) */
    gretl_matrix *ZT;  /* N x (r * n) */
};

struct kalman_ {
    int flags;   /* for recording any options */
    int exact;   /* exact initial iterations? */
    int vartype; /* STD_VAR, DJ_VAR or DK_VAR */

    int r;   /* rows of a = number of elements in state */
    int n;   /* columns of y = number of observables */
    int k;   /* columns of B = number of exogenous vars in obs eqn */
    int p;   /* length of combined disturbance vector */
    int q;   /* length of state disturbance vector */
    int N;   /* rows of y = number of observations */
    int okN; /* N - number of missing observations */
    int t;   /* current time step, when filtering */
    int d;   /* time-step at which standard iterations start */
    int j;   /* observable at which standard iters start */

    int ifc; /* boolean: obs equation includes an implicit constant? */

    double SSRw;    /* \sum_{t=1}^N v_t^{\prime} F_t^{-1} v_t */
    double loglik;  /* log-likelihood */
    double s2;      /* = SSRw / k */

    /* continuously updated matrices */
    gretl_matrix *a0; /* r x 1: state vector, before updating */
    gretl_matrix *a1; /* r x 1: state vector, after updating */
    gretl_matrix *P0; /* r x r: MSE matrix, before updating */
    gretl_matrix *P1; /* r x r: MSE matrix, after updating */
    gretl_matrix *vt; /* n x 1: one-step forecast error(s), time t */

    /* for the exact diffuse case */
    gretl_matrix *Pk0; /* P∞,t */
    gretl_matrix *Pk1; /* P∞,t update */
    gretl_matrix *Fk;  /* F∞,t */
    gretl_matrix *Ck;  /* C∞,t */
    gretl_matrix *PK;  /* P∞, all t */

    /* input data matrices: note that the order matters for various
       functions, including matrix_is_varying()
    */
    gretl_matrix *T;  /* r x r: state transition matrix */
    gretl_matrix *BT; /* k x n: B', coeffs on exogenous vars, obs eqn */
    gretl_matrix *ZT; /* r x n: Z', coeffs on state variables, obs eqn */
    gretl_matrix *VS; /* r x r: contemp covariance matrix, state eqn */
    gretl_matrix *VY; /* n x n: contemp covariance matrix, obs eqn */
    gretl_matrix *mu; /* r x 1: constant term in state transition */
    gretl_matrix *y;  /* N x n: dependent variable vector (or matrix) */
    gretl_matrix *x;  /* N x k: independent variables matrix */
    gretl_matrix *aini; /* r x 1: a_0 */
    gretl_matrix *Pini; /* r x r: P_0 */

    /* user inputs for cross-correlated disturbances */
    gretl_matrix *H;  /* r x p: state var factor */
    gretl_matrix *G;  /* n x p: obs var factor */
    /* and the cross-matrix itself */
    gretl_matrix *HG; /* r x n */

    /* OR: Durbin-Koopman R and Q (as in statevar = RQR') */
    gretl_matrix *Q;   /* q x q, symmetric */
    gretl_matrix *R;   /* r x q, selection matrix */
    gretl_matrix *QRT; /* Q*R' */

    /* apparatus for registering time-variation of matrices */
    char *matcall;
    char *varying;

    /* optional matrices for recording extra info */
    gretl_matrix *LL;  /* N x 1: loglikelihood, all time-steps */

    /* optional run-time export matrices */
    gretl_matrix *V;   /* N x n: forecast errors, all time-steps */
    gretl_matrix *F;   /* N x nn: MSE for observables, all time-steps */
    gretl_matrix *A;   /* N x r: state vector, all time-steps */
    gretl_matrix *P;   /* N x rr: MSE for state, all time-steps */
    gretl_matrix *K;   /* N x rn: gain matrix, all time-steps */
    gretl_matrix *U;   /* N x ??: smoothed disturbances */
    gretl_matrix *Vsd; /* Variance of smoothed disturbance */

    /* structure needed only when smoothing in the time-varying case */
    stepinfo *step;

    /* workspace matrices */
    gretl_matrix_block *Blk; /* holder for the following */
    gretl_matrix *PZ;
    gretl_matrix *Ft;
    gretl_matrix *iFt;
    gretl_matrix *Kt;
    gretl_matrix *Mt;
    gretl_matrix *Ct;

    /* backups for matrices that are shrunk in case of
       incomplete observations, plus recorder array for
       use with the smoother (if needed)
    */
    gretl_matrix *saveZT;
    gretl_matrix *saveVY;
    gretl_matrix *saveG;
    guint16 *nt;

    gretl_bundle *b; /* the bundle of which this struct is a member */
    void *data;      /* handle for attaching additional info */
    PRN *prn;        /* verbose printer */
};

enum {
    K_V,
    K_F,
    K_A,
    K_BIG_P,
    K_K,
    K_LL,
};

/* max number of time-varying matrices: T, B, Z, VS, VY, mu */
#define K_N_MATCALLS 6

#define set_kalman_running(K) (K->flags |= KALMAN_FORWARD)
#define set_kalman_stopped(K) (K->flags &= ~KALMAN_FORWARD)
#define kalman_is_running(K)  (K->flags & KALMAN_FORWARD)
#define kalman_simulating(K)  (K->flags & KALMAN_SIM)
#define kalman_checking(K)    (K->flags & KALMAN_CHECK)
#define kalman_ssfsim(K)      (K->flags & KALMAN_SSFSIM)
#define kalman_diffuse(K)     (K->flags & KALMAN_DIFFUSE)
#define kalman_smoothing(K)   (K->flags & KALMAN_SMOOTH)
#define kalman_arma_ll(K)     (K->flags & KALMAN_ARMA_LL)
#define basic_smoothing(K)    (K->flags & KALMAN_SM_AM)
#define kalman_univariate(K)  (K->flags & KALMAN_UNI)

#define kalman_xcorr(K)       (K->vartype == DJ_VAR)
#define kalman_dkvar(K)       (K->vartype == DK_VAR)

#define filter_is_varying(K) (K->matcall != NULL)

static const char *kalman_matrix_name (int sym);
static int kalman_revise_variance (kalman *K);
static int check_for_matrix_updates (kalman *K, ufunc *uf);
static int construct_Pini (kalman *K);
static int kalman_set_diffuse (kalman *K, int d);
static int kfilter_univariate (kalman *K, PRN *prn);
static int ksmooth_univariate (kalman *K, int dist);

/* symbolic identifiers for input matrices: note that potentially
   time-varying matrices must appear first in the enumeration, and
   the order must match the order of the "input data matrices" in
   the Kalman struct (above).
*/

enum {
    K_T = 0,
    K_BT,
    K_ZT,
    K_VS,
    K_VY,
    K_m,
    K_y,
    K_x,
    K_a,
    K_P,
    K_R,
    K_MMAX /* sentinel */
};

enum {
    SM_TYPE_NONE, /* not doing smoothing */
    SM_STATE_STD, /* regular state smoother */
    SM_STATE_INI, /* exact initial state smoother */
    SM_DIST_BKWD, /* disturbance smoother, backward pass */
    SM_DIST_FRWD, /* disturbance smoother, forward pass */
};

void free_stepinfo (kalman *K)
{
    if (K->step != NULL) {
        gretl_matrix_free(K->step->T);
        gretl_matrix_free(K->step->ZT);
        free(K->step);
        K->step = NULL;
    }
}

void kalman_free (kalman *K)
{
    if (K == NULL) {
        return;
    }

    /* internally allocated matrices */
    gretl_matrix_free(K->a0);
    gretl_matrix_free(K->a1);
    gretl_matrix_free(K->P0);
    gretl_matrix_free(K->P1);
    gretl_matrix_free(K->vt);
    gretl_matrix_free(K->LL);
    gretl_matrix_free(K->Pk0);
    gretl_matrix_free(K->Pk1);
    gretl_matrix_free(K->Fk);
    gretl_matrix_free(K->Ck);

    /* internally allocated workspace */
    gretl_matrix_block_destroy(K->Blk);

    if (K->flags & KALMAN_BUNDLE) {
        gretl_matrix **mptr[] = {
            &K->T, &K->BT, &K->ZT, &K->VS, &K->VY,
            &K->mu, &K->y, &K->x, &K->aini, &K->Pini,
	    &K->R
        };
        int i;

        for (i=0; i<K_MMAX; i++) {
            gretl_matrix_free(*mptr[i]);
        }

        /* @K also owns these "export" matrices */
        gretl_matrix_free(K->V);
        gretl_matrix_free(K->F);
        gretl_matrix_free(K->A);
        gretl_matrix_free(K->P);
        gretl_matrix_free(K->K);
        gretl_matrix_free(K->U);
        gretl_matrix_free(K->Vsd);
    }

    /* time-variation info */
    free(K->matcall);
    free(K->varying);

    if (kalman_xcorr(K)) {
        /* correlated errors info */
        gretl_matrix_free(K->H);
        gretl_matrix_free(K->G);
        gretl_matrix_free(K->HG);
    } else if (kalman_dkvar(K)) {
	/* Durbin-Koopman errors info */
	gretl_matrix_free(K->Q);
	gretl_matrix_free(K->QRT);
    }

    if (K->step != NULL) {
        free_stepinfo(K);
    }

    free(K);
}

static kalman *kalman_new_empty (int flags)
{
    kalman *K = malloc(sizeof *K);

    if (K != NULL) {
        K->exact = 0;
	K->vartype = STD_VAR;
        K->aini = K->Pini = NULL;
        K->a0 = K->a1 = NULL;
        K->P0 = K->P1 = NULL;
        K->Pk0 = K->Pk1 = NULL;
        K->PK = NULL;
        K->Fk = K->Ck = NULL;
        K->LL = NULL;
        K->vt = NULL;
        K->Blk = NULL;
        K->T = K->BT = K->ZT = NULL;
        K->VS = K->VY = NULL;
        K->H = K->G = K->HG = NULL;  /* DeJong variance */
	K->Q = K->R = K->QRT = NULL; /* Durbin-Koopman variance */
        K->V = K->F = K->A = K->K = K->P = NULL;
        K->y = K->x = NULL;
        K->mu = NULL;
        K->U = NULL;
        K->Vsd = NULL;
        K->matcall = NULL;
        K->varying = NULL;
        K->step = NULL;
        K->flags = flags;
        K->t = 0;
        K->prn = NULL;
        K->data = NULL;
        K->b = NULL;
        K->d = 0;
    }

    return K;
}

static void kalman_set_dimensions (kalman *K)
{
    K->r = gretl_matrix_rows(K->T);  /* T->rows defines r */
    K->k = gretl_matrix_rows(K->BT); /* BT->rows defines k */
    K->n = gretl_matrix_cols(K->y);  /* y->cols defines n */

    if (!kalman_simulating(K)) {
        /* y->rows defines N, except when simulating */
        K->N = gretl_matrix_rows(K->y);
    }

    K->okN = K->N;

    /* K->p is non-zero only under cross-correlation; in that case the
       matrix given as 'Q' in Kalman set-up in fact represents H (as
       in v_t = H \varepsilon_t) and it must be r x p, where p is the
       number of elements in the "combined" disturbance vector
       \varepsilon_t.
    */
    K->p = (K->H != NULL)? gretl_matrix_cols(K->H): 0;

    /* K->q is non-zero only when using the Durbin-Koopman
       representation of the variance of the state disturbance as
       R*Q*R'
    */
    K->q = (K->Q != NULL)? gretl_matrix_rows(K->Q): 0;
}

static int missing_matrix_error (const char *name)
{
    if (name == NULL) {
        gretl_errmsg_set(_("kalman: a required matrix is missing"));
    } else {
        gretl_errmsg_sprintf(_("kalman: required matrix %s is missing"),
                             name);
    }
    return E_DATA;
}

static int check_matrix_dims (kalman *K, const gretl_matrix *m, int i)
{
    int r = 0, c = 0, symm = (i == K_VS || i == K_VY);
    int err = 0;

    if (i == K_T || i == K_VS || i == K_P) {
        r = c = K->r;
    } else if (i == K_BT) {
        r = K->k;
        c = K->n;
    } else if (i == K_ZT) {
        r = K->r;
        c = K->n;
    } else if (i == K_VY)  {
        r = c = K->n;
    } else if (i == K_a || i == K_m) {
        r = K->r;
        c = 1;
    }

    if (m->rows != r || m->cols != c) {
        gretl_errmsg_sprintf("kalman: %s is %d x %d, should be %d x %d\n",
                             kalman_matrix_name(i), m->rows, m->cols, r, c);
        err = E_NONCONF;
    } else if (symm && !gretl_matrix_is_symmetric(m)) {
        gretl_errmsg_sprintf("kalman: %s is not symmetric\n",
                kalman_matrix_name(i));
        err = E_NONCONF;
    }

    return err;
}

static int maybe_resize_export_matrix (kalman *K, gretl_matrix *m, int i)
{
    int rows = K->N, cols = 0;
    int err = 0;

    if (i == K_V) {
        cols = K->n;
    } else if (i == K_F) {
	if (K->flags & KALMAN_UNI) {
	    cols = K->n; /* diagonal only */
	} else {
	    cols = (K->n * K->n + K->n) / 2;
	}
    } else if (i == K_A) {
        cols = K->r;
    } else if (i == K_BIG_P) {
        cols = (K->r * K->r + K->r) / 2;
    } else if (i == K_LL) {
        cols = 1;
    } else if (i == K_K) {
        cols = K->r * K->n;
    } else {
        err = E_DATA;
    }

    if (!err && (m->rows != rows || m->cols != cols)) {
        err = gretl_matrix_realloc(m, rows, cols);
        if (!err) {
            gretl_matrix_zero(m);
        }
    }

    return err;
}

static int kalman_check_dimensions (kalman *K)
{
    int err = 0;

    if (K->r < 1 || K->n < 1 || K->N < 2) {
        /* the state and observation vectors must have at least one
           element, and there must be at least two observations
        */
        err = E_DATA;
    }

    /* T is mandatory, should be r x r */
    if (!err) {
        err = check_matrix_dims(K, K->T, K_T);
    }

    /* Z' is mandatory, should be r x n */
    if (!err) {
        err = check_matrix_dims(K, K->ZT, K_ZT);
    }

    if (K->H != NULL) {
        /* cross-correlated disturbances */
        if (K->G == NULL) {
            err = missing_matrix_error("obsymat");
        } else if (K->H->rows != K->r || K->H->cols != K->p ||
                   K->G->rows != K->n || K->G->cols != K->p) {
            fprintf(stderr, "H is %d x %d, G is %d x %d, p=%d, r=%d, n=%d\n",
                    K->H->rows, K->H->cols, K->G->rows, K->G->cols,
                    K->p, K->r, K->n);
            err = E_NONCONF;
        }
    } else if (K->R != NULL) {
	/* state disturbances as per Durbin-Koopman */
	if (K->Q == NULL) {
	    err = missing_matrix_error("statevar");
	} else if (K->R->rows != K->r || K->R->cols != K->Q->rows) {
            fprintf(stderr, "R is %d x %d, Q is %d x %d, r=%d\n",
                    K->R->rows, K->R->cols, K->Q->rows, K->Q->cols,
                    K->r);
            err = E_NONCONF;
        }
    } else {
        /* VS is mandatory, should be r x r and symmetric */
        if (!err) {
            err = check_matrix_dims(K, K->VS, K_VS);
        }
        /* VY should be n x n and symmetric, if present */
        if (!err && K->VY != NULL) {
            err = check_matrix_dims(K, K->VY, K_VY);
        }
    }

    /* initial a should be r x 1, if present */
    if (!err && K->aini != NULL) {
        err = check_matrix_dims(K, K->aini, K_a);
    }

    /* initial P should be r x r, if present */
    if (!err && K->Pini != NULL) {
        err = check_matrix_dims(K, K->Pini, K_P);
    }

    /* B' should be k x n, if present (BT->rows defines k) */
    if (!err) {
        if (K->BT != NULL) {
            err = check_matrix_dims(K, K->BT, K_BT);
        } else if (K->x != NULL) {
            /* B is NULL => can't have a non-NULL x */
            err = E_NONCONF;
        }
    }

    /* mu should be r x 1, if present */
    if (!err && K->mu != NULL) {
        err = check_matrix_dims(K, K->mu, K_m);
    }

    if (err) {
        goto bailout;
    }

    K->ifc = 0;

    /* x should have N rows to match y; and it should have either k or k - 1
       columns (the latter case indicating an implicit const) */
    if (K->x != NULL) {
        if (K->x->rows < K->N) {
            fprintf(stderr, "kalman: %s has %d rows, should have %d\n",
                    kalman_matrix_name(K_x), K->x->rows, K->N);
            return E_NONCONF;
        } else if (K->x->cols != K->k && K->x->cols != K->k - 1) {
            fprintf(stderr, "kalman: %s has %d columns, should have %d or %d\n",
                    kalman_matrix_name(K_x), K->x->cols, K->k, K->k - 1);
            return E_NONCONF;
        } else if (K->x->cols == K->k - 1) {
            /* register the implicit const */
            K->ifc = 1;
        }
    } else if (K->k == 1) {
        /* B has one row but there's no x => implicit const */
        K->ifc = 1;
    } else if (K->k > 1) {
        /* B has more than one row but there's no x => error */
        return missing_matrix_error("obsxmat");
    }

    /* Below we have the optional "export" matrices for shipping out
       results. If these are present but not sized correctly we'll
       try to fix them up -- but note that they are not used in a
       simulation run.
    */

    if (kalman_simulating(K)) {
        return err;
    }

    /* big V should be N x n */
    if (K->V != NULL) {
        err = maybe_resize_export_matrix(K, K->V, K_V);
    }

    /* big F should be N x nn */
    if (!err && K->F != NULL) {
        err = maybe_resize_export_matrix(K, K->F, K_F);
    }

    /* big A should be N x r */
    if (!err && K->A != NULL) {
        err = maybe_resize_export_matrix(K, K->A, K_A);
    }

    /* big P should be N x nr */
    if (!err && K->P != NULL) {
        err = maybe_resize_export_matrix(K, K->P, K_BIG_P);
    }

    /* LL should be N x 1 */
    if (!err && K->LL != NULL) {
        err = maybe_resize_export_matrix(K, K->LL, K_LL);
    }

    /* K (gain) should be N x (r * n) */
    if (!err && K->K != NULL) {
        err = maybe_resize_export_matrix(K, K->K, K_K);
    }

 bailout:

    if (err) {
        fprintf(stderr, "kalman_check_dimensions: err = %d\n", err);
    }

    return err;
}

static int kalman_init (kalman *K)
{
    int err = 0;

    K->SSRw = NADBL;
    K->loglik = NADBL;
    K->s2 = NADBL;

    clear_gretl_matrix_err();

    if (K->aini != NULL) {
        K->a0 = gretl_matrix_copy(K->aini);
        K->a1 = gretl_matrix_copy(K->aini);
    } else {
        K->a0 = gretl_zero_matrix_new(K->r, 1);
        K->a1 = gretl_zero_matrix_new(K->r, 1);
    }

    if (K->Pini != NULL) {
        K->P0 = gretl_matrix_copy(K->Pini);
        K->P1 = gretl_matrix_copy(K->Pini);
    } else {
        K->P0 = gretl_zero_matrix_new(K->r, K->r);
        K->P1 = gretl_zero_matrix_new(K->r, K->r);
    }

    /* forecast error vector, per observation */
    K->vt = gretl_matrix_alloc(K->n, 1);

    err = get_gretl_matrix_err();
    if (err) {
        return err;
    }

    K->Blk = gretl_matrix_block_new(&K->PZ,  K->r, K->n, /* P*Z */
                                    &K->Ft,  K->n, K->n, /* (Z'*P*Z + R)^{-1} */
                                    &K->iFt, K->n, K->n, /* Ft-inverse */
                                    &K->Kt,  K->r, K->n, /* gain at t */
                                    &K->Mt,  K->r, K->n, /* intermediate term */
                                    &K->Ct,  K->r, K->r, /* intermediate term */
                                    NULL);

    if (K->Blk == NULL) {
        err = E_ALLOC;
    }

    if (!err && K->Pini == NULL && !(K->flags & KALMAN_USER)) {
        /* in the "user" case we do this later */
        err = construct_Pini(K);
    }

    return err;
}

/* The following section includes functions that support the "plain C"
   Kalman API, as opposed to the bundle-based interface that subserves
   userland state-space functionality. We use this API in gretl's arma
   plugin in the special case of ARIMA with non-zero order of
   integration along with missing values. It's also possible that
   third-party users of libgretl might wish to use it.
*/

/**
 * kalman_new:
 * @a: r x 1 initial state vector.
 * @P: r x r initial precision matrix.
 * @T: r x r state transition matrix.
 * @BT: n x k matrix of coefficients on exogenous variables in the
 * observation equation, transposed.
 * @ZT: n x r matrix of coefficients on the state variables in the
 * observation equation, transposed.
 * @VS: r x r contemporaneous covariance matrix for the errors in the
 * state equation.
 * @VY: n x n contemporaneous covariance matrix for the errors in the
 * observation equation (or NULL if this is not applicable).
 * @y: T x n matrix of observable variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be NULL if there
 * are no exogenous variables, or if there's only a constant.
 * @mu: r x 1 vector of constants in the state transition, or NULL.
 * @V: T x n matrix in which to record forecast errors (or NULL if
 * this is not required).
 * @err: location to receive error code.
 *
 * Allocates and initializes a Kalman struct, which can subsequently
 * be used for forecasting with kalman_forecast().
 *
 * Returns: pointer to allocated struct, or NULL on failure, in
 * which case @err will receive a non-zero code.
 */

kalman *kalman_new (gretl_matrix *a, gretl_matrix *P,
                    gretl_matrix *T, gretl_matrix *BT,
                    gretl_matrix *ZT, gretl_matrix *VS,
                    gretl_matrix *VY, gretl_matrix *y,
                    gretl_matrix *x, gretl_matrix *mu,
                    gretl_matrix *V, int *err)
{
    kalman *K;

    *err = 0;

    if (y == NULL || T == NULL || ZT == NULL || VS == NULL) {
        fprintf(stderr, "kalman_new: y=%p, T=%p, ZT=%p, VS=%p\n",
                (void *) y, (void *) T, (void *) ZT, (void *) VS);
        *err = missing_matrix_error(NULL);
        return NULL;
    }

    K = kalman_new_empty(0);
    if (K == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    /* use pointers for input matrices, don't copy */
    K->T = T;
    K->BT = BT;
    K->ZT = ZT;
    K->VS = VS;
    K->VY = VY;
    K->y = y;
    K->x = x;
    K->aini = a;
    K->Pini = P;
    K->mu = mu;

    /* output, but again use external pointer */
    K->V = V;

    kalman_set_dimensions(K);

    *err = kalman_check_dimensions(K);
    if (*err) {
        free(K);
        return NULL;
    }

    *err = kalman_init(K);

    if (*err) {
        kalman_free(K);
        K = NULL;
    } else {
        gretl_matrix_zero(K->vt);
    }

    return K;
}

/**
 * kalman_get_loglik:
 * @K: pointer to Kalman struct.
 *
 * Retrieves the log-likelhood calculated via a run of
 * kalman_forecast().
 *
 * Returns: ll value, or #NADBL on failure.
 */

double kalman_get_loglik (const kalman *K)
{
    return K->loglik;
}

/**
 * kalman_get_arma_variance:
 * @K: pointer to Kalman struct.
 *
 * Retrieves the estimated variance for an ARMA model
 * estimated using the Kalman filter.
 *
 * Returns: sigma-squared value, or #NADBL on failure.
 */

double kalman_get_arma_variance (const kalman *K)
{
    if (na(K->SSRw)) {
        return NADBL;
    } else {
        return K->SSRw / K->okN;
    }
}

/**
 * kalman_set_initial_state_vector:
 * @K: pointer to Kalman struct.
 * @a: vector of values to set.
 *
 * Resets the initial value of the state vector in a Kalman
 * struct, using the values from @a.  See also kalman_new().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_state_vector (kalman *K, const gretl_vector *a)
{
    return gretl_matrix_copy_values(K->a0, a);
}

/**
 * kalman_set_initial_MSE_matrix:
 * @K: pointer to Kalman struct.
 * @P: matrix of values to set.
 *
 * Resets the initial value of the MSE matrix in a Kalman
 * struct, using the values from @P.  See also kalman_new().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P)
{
    return gretl_matrix_copy_values(K->P0, P);
}

void kalman_attach_data (kalman *K, void *data)
{
    if (K != NULL) {
        K->data = data;
    }
}

void *kalman_get_data (const kalman *K)
{
    return (K != NULL)? K->data : NULL;
}

void kalman_attach_printer (kalman *K, PRN *prn)
{
    if (K != NULL) {
        K->prn = prn;
    }
}

PRN *kalman_get_printer (const kalman *K)
{
    return (K != NULL)? K->prn : NULL;
}

void kalman_set_arma_ll (kalman *K)
{
    K->flags |= KALMAN_ARMA_LL;
}

/* end of functions dedicated to non-bundle C API */

#define kappa 1.0e7 /* 1.0e6? */

static int diffuse_Pini (kalman *K)
{
    gretl_matrix_zero(K->P0);

    if (K->exact) {
        if (K->Pk0 == NULL) {
            K->Pk0 = gretl_identity_matrix_new(K->r);
            K->Pk1 = gretl_matrix_alloc(K->r, K->r);
            K->Fk  = gretl_matrix_alloc(K->n, K->n);
            K->Ck  = gretl_matrix_alloc(K->r, K->r);
            if (K->Pk0 == NULL || K->Pk1 == NULL ||
                K->Fk == NULL || K->Ck == NULL) {
                return E_ALLOC;
            }
        } else {
            gretl_matrix_inscribe_I(K->Pk0, 0, 0, K->r);
        }
    } else {
        /* old-style */
        int i;

        for (i=0; i<K->r; i++) {
            gretl_matrix_set(K->P0, i, i, kappa);
        }
    }

    return 0;
}

static int statemat_out_of_bounds (kalman *K)
{
    gretl_matrix *evals;
    double r, c, x;
    int i, err = 0;

    evals = gretl_general_matrix_eigenvals(K->T, &err);

    for (i=0; i<evals->rows && !err; i++) {
        r = gretl_matrix_get(evals, i, 0);
        c = gretl_matrix_get(evals, i, 1);
        x = sqrt(r*r + c*c);
        if (x >= 1.0) {
            fprintf(stderr, "T: modulus of eigenvalue %d = %g\n", i+1, x);
            err = E_SINGULAR;
        }
    }

    gretl_matrix_free(evals);

    return err;
}

/* If the user has not given an initial value for P_{1|0}, compute
   this automatically as per Hamilton, ch 13, p. 378.  This works only
   if the eigenvalues of K->T lie inside the unit circle.  Failing
   that, or if the --diffuse option is given for the user Kalman
   filter, we apply a diffuse initialization.
*/

static int construct_Pini (kalman *K)
{
    gretl_matrix *Svar;
    gretl_matrix *vQ;
    int r2, err = 0;

    if (K->flags & KALMAN_DIFFUSE) {
        return diffuse_Pini(K);
    }

    r2 = K->r * K->r;

    Svar = gretl_matrix_alloc(r2, r2);
    vQ = gretl_column_vector_alloc(r2);

    if (Svar == NULL || vQ == NULL) {
        gretl_matrix_free(Svar);
        gretl_matrix_free(vQ);
        return E_ALLOC;
    }

    gretl_matrix_kronecker_product(K->T, K->T, Svar);
    gretl_matrix_I_minus(Svar);
    gretl_matrix_vectorize(vQ, K->VS);

    err = gretl_LU_solve(Svar, vQ);
    if (err) {
        /* failed: are some of the eigenvalues out of bounds? */
        err = statemat_out_of_bounds(K);
        if (err == E_SINGULAR) {
            err = diffuse_Pini(K);
            K->flags |= KALMAN_DIFFUSE;
        }
    } else {
        gretl_matrix_unvectorize(K->P0, vQ);
    }

    gretl_matrix_free(Svar);
    gretl_matrix_free(vQ);

    return err;
}

/* variant of gretl_matrix_copy_values for us when we already
   know that the matrices are non-NULL and conformable
*/

static inline void fast_copy_values (gretl_matrix *B, const gretl_matrix *A)
{
    memcpy(B->val, A->val, B->rows * B->cols * sizeof(double));
}

/* Write the vech of @src into row @t of @targ */

static void load_to_vech (gretl_matrix *targ,
                          const gretl_matrix *src,
                          int n, int t)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            x = gretl_matrix_get(src, i, j);
            gretl_matrix_set(targ, t, m++, x);
        }
    }
}

/* Write the vech of @src into column @t of @targ */

static void load_to_vechT (gretl_matrix *targ,
                           const gretl_matrix *src,
                           int n, int t)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            x = gretl_matrix_get(src, i, j);
            gretl_matrix_set(targ, m++, t, x);
        }
    }
}

/* Write the vec of @src into row @t of @targ */

static void load_to_vec (gretl_matrix *targ,
                         const gretl_matrix *src,
                         int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* Write the square root of diagonal of square matrix @src
   into row @t of @targ, starting at column offset @j
*/

static void load_to_diag (gretl_matrix *targ,
                          const gretl_matrix *src,
                          int t, int j)
{
    int i, n = gretl_vector_get_length(src);
    double x;

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, i, i);
        if (x <= 0.0) {
            gretl_matrix_set(targ, t, i+j, 0.0);
        } else {
            gretl_matrix_set(targ, t, i+j, sqrt(x));
        }
    }
}

/* copy from vector @src into row @t of @targ */

static void load_to_row (gretl_matrix *targ,
                         const gretl_vector *src,
                         int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* copy from vector @src into row @t of @targ,
   allowing for insertion of missing values
*/

static void load_to_row_special (gretl_matrix *targ,
                                 const gretl_vector *src,
                                 int t, char *missmap)
{
    int i, k;

    for (i=0, k=0; i<targ->cols; i++) {
        if (missmap[i]) {
            gretl_matrix_set(targ, t, i, NADBL);
        } else {
            gretl_matrix_set(targ, t, i, src->val[k++]);
        }
    }
}

/* copy from vector @src into row @t of @targ,
   starting at column offset @j in @targ */

static void load_to_row_offset (gretl_matrix *targ,
                                const gretl_vector *src,
                                int t, int j)
{
    int i, n = gretl_vector_get_length(src);

    for (i=0; i<n; i++) {
        gretl_matrix_set(targ, t, i+j, src->val[i]);
    }
}

static void set_row_to_value (gretl_matrix *targ, int t, double x)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, x);
    }
}

/* supports hansl function for creating a named Kalman bundle */

kalman *kalman_new_minimal (gretl_matrix *M[], int copy[],
                            int nmat, int dkvar, int *err)
{
    gretl_matrix **targ[5];
    kalman *K;
    int i;

    *err = 0;

    if (M[0] == NULL || M[1] == NULL || M[2] == NULL || M[3] == NULL) {
        fprintf(stderr, "kalman_new_minimal: nmat=%d, y=%p, Z=%p, T=%p, Q=%p\n",
                nmat, (void *) M[0], (void *) M[1], (void *) M[2], (void *) M[3]);
        *err = missing_matrix_error(NULL);
        return NULL;
    }

    K = kalman_new_empty(KALMAN_USER | KALMAN_BUNDLE);
    if (K == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    targ[0] = &K->y;
    targ[1] = &K->ZT;
    targ[2] = &K->T;

    if (nmat == 5) {
	if (dkvar) {
	    K->vartype = DK_VAR;
	    targ[3] = &K->Q;
	    targ[4] = &K->R;
	} else {
	    K->vartype = DJ_VAR;
	    targ[3] = &K->H;
	    targ[4] = &K->G;
	}
    } else {
        targ[3] = &K->VS; /* statevar */
    }

    for (i=0; i<nmat; i++) {
        if (copy[i]) {
            *targ[i] = gretl_matrix_copy(M[i]);
        } else {
            *targ[i] = M[i];
        }
    }

    kalman_set_dimensions(K);

    if (K->vartype != STD_VAR) {
        *err = kalman_revise_variance(K);
    }
    if (!*err) {
        *err = kalman_check_dimensions(K);
    }
    if (!*err) {
        *err = kalman_init(K);
    }

    if (*err) {
        kalman_free(K);
        K = NULL;
    } else {
        gretl_matrix_zero(K->vt);
    }

    return K;
}

static int matrix_is_varying (kalman *K, int i)
{
    if (K->matcall != NULL) {
        if (K->varying == NULL) {
            check_for_matrix_updates(K, NULL);
        }
        if (K->varying != NULL) {
            return K->varying[i];
        }
    }

    return 0;
}

enum {
    UPDATE_INIT, /* initialization of matrices */
    UPDATE_STEP  /* refreshing matrices per time-step */
};

/* Ensure we have suitable matrices into which to write
   VS, VY and HG'
*/

static int ensure_cross_covariance_matrices (kalman *K)
{
    int r[] = {K->H->rows, K->G->rows, K->H->rows};
    int c[] = {K->H->rows, K->G->rows, K->G->rows};
    gretl_matrix **V[] = {&K->VS, &K->VY, &K->HG};
    gretl_matrix **targ;
    int i, err = 0;

    for (i=0; i<3 && !err; i++) {
        targ = V[i];
        if (*targ == NULL) {
            *targ = gretl_matrix_alloc(r[i], c[i]);
        } else {
            gretl_matrix *Vi = *targ;

            if (Vi->rows != r[i] || Vi->cols != c[i]) {
                gretl_matrix_realloc(Vi, r[i], c[i]);
            }
        }
        if (*targ == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

static int ensure_DK_covariance_matrices (kalman *K)
{
    int r[] = {K->r, K->q};
    int c[] = {K->r, K->r};
    gretl_matrix **V[] = {&K->VS, &K->QRT};
    gretl_matrix **targ;
    int i, err = 0;

    for (i=0; i<2 && !err; i++) {
        targ = V[i];
        if (*targ == NULL) {
            *targ = gretl_matrix_alloc(r[i], c[i]);
        } else {
            gretl_matrix *Vi = *targ;

            if (Vi->rows != r[i] || Vi->cols != c[i]) {
                gretl_matrix_realloc(Vi, r[i], c[i]);
            }
        }
        if (*targ == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

/* After reading H and G from the user, either at (re-)initialization
   or at a given time-step in the case where either of these matrices
   is time-varying, record the user input in K->H and K->G and form
   the 'real' VS and VY.  But note that in the time-step case it may be
   that only one of VS, VY needs to be treated in this way (if only one
   is time-varying, only one will have been redefined via a function
   call).
*/

static int kalman_update_crossinfo (kalman *K, int mode)
{
    int err = 0;

    /* Note that H and G may be needed as such for simulation */

    if (mode == UPDATE_INIT) {
        err = ensure_cross_covariance_matrices(K);
        if (err) {
            return err;
        }
    }

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_VS)) {
        /* (re)create VS using modified H */
        err = gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                        K->H, GRETL_MOD_TRANSPOSE,
                                        K->VS, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_VY))) {
        /* (re)create VY using modified G */
        err = gretl_matrix_multiply_mod(K->G, GRETL_MOD_NONE,
                                        K->G, GRETL_MOD_TRANSPOSE,
                                        K->VY, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_VS) ||
                 matrix_is_varying(K, K_VY))) {
        /* (re)create HG' using modified H and/or G */
        err = gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                        K->G, GRETL_MOD_TRANSPOSE,
                                        K->HG, GRETL_MOD_NONE);
    }

    return err;
}

static int kalman_update_dkvar (kalman *K, int mode)
{
    int err = 0;

    if (mode == UPDATE_INIT) {
        err = ensure_DK_covariance_matrices(K);
        if (err) {
            return err;
        }
    }

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_VS)) {
        /* (re)create VS and QR' using modified Q */
	//fprintf(stderr, "*** kalman_update_dkvar, mode %d\n", mode);
	gretl_matrix_qform(K->R, GRETL_MOD_NONE, K->Q,
			   K->VS, GRETL_MOD_NONE);
        err = gretl_matrix_multiply_mod(K->Q, GRETL_MOD_NONE,
                                        K->R, GRETL_MOD_TRANSPOSE,
                                        K->QRT, GRETL_MOD_NONE);
     }

    return err;
}

/* kalman_revise_DJ_variance: the user has actually given H in place
   of VS and G in place of VY; so we have to form VS, VY, and HG'
   (cross-correlated disturbances).

   This function is called in the course of initial set-up of a
   filter, and also when the Kalman matrices are being re-checked at
   the start of filtering, smoothing or simulation.
*/

static int kalman_revise_DJ_variance (kalman *K)
{
    int err = 0;

    if (K->H == NULL || K->G == NULL) {
	fprintf(stderr, "K->H %p, K->G %p\n", (void *) K->H, (void *) K->G);
        return missing_matrix_error("'statevar' or 'obsvar'");
    }

    err = kalman_update_crossinfo(K, UPDATE_INIT);

    if (err) {
        fprintf(stderr, "kalman_revise_cross_variance: err = %d\n", err);
    }

    return err;
}

static int kalman_revise_DK_variance (kalman *K)
{
    int err = 0;

    if (K->Q == NULL || K->R == NULL) {
        return missing_matrix_error("'statevar' or 'statsel'");
    }

    err = kalman_update_dkvar(K, UPDATE_INIT);

    if (err) {
        fprintf(stderr, "kalman_revise_DK_variance: err = %d\n", err);
    }

    return err;
}

static int kalman_revise_variance (kalman *K)
{
    if (K->vartype == DJ_VAR) {
	return kalman_revise_DJ_variance(K);
    } else {
	return kalman_revise_DK_variance(K);
    }
}

#if KDEBUG > 1
static void kalman_print_state (kalman *K)
{
   int j;

    /* if (t > 5) return; */

    fprintf(stderr, "t = %d:\n", K->t);

    for (j=0; j<K->n; j++) {
        fprintf(stderr, "y[%d] = %.10g, err[%d] = %.10g\n", j,
                gretl_matrix_get(K->y, K->t, j),
                j, gretl_vector_get(K->vt, j));
    }

    gretl_matrix_print(K->a0, "K->a0");
    gretl_matrix_print(K->P0, "K->P0");
}
#endif

/* On filtering: record the state and/or its variance, as
   wanted. Plus in the "exact initial" case record P∞,t
   in K->PK.
*/

static int kalman_record_state (kalman *K)
{
    int err = 0;

    if (K->A != NULL) {
        load_to_row(K->A, K->a0, K->t);
    }
    if (K->P != NULL) {
        load_to_vech(K->P, K->P0, K->r, K->t);
    }

    if (K->exact && K->PK != NULL) {
        if (K->t >= K->PK->cols) {
            err = gretl_matrix_realloc(K->PK, K->PK->rows, K->t + 1);
        }
        if (!err) {
            load_to_vechT(K->PK, K->Pk0, K->r, K->t);
        }
    }

    return err;
}

/* Read from the appropriate row of x (N x k) and multiply by B' to
   form B'x_t.  Note: the flag K->ifc is used to indicate that the
   observation equation has an implicit constant, with an entry in
   the B matrix (the first) but no explicit entry in the x matrix.

   The case where x is NULL and B an n-vector (implicit constant)
   is also handled.

   There's no need to store B'x_t: we either want to add this to,
   or subtract it from, the n-vector @targ (the operator being
   indicated by the @mod argument).

   Returns: the number of missing values encountered.
*/

static int kalman_do_Bx (kalman *K, gretl_matrix *targ,
                         GretlMatrixMod mod)
{
    double xjt, bji, bxi;
    int i, j, missvals = 0;

    for (i=0; i<K->n; i++) {
        if (K->x == NULL) {
            /* the implicit constant case */
            bxi = K->BT->val[i];
        } else {
            bxi = 0;
            for (j=0; j<K->k; j++) {
                if (K->ifc) {
                    xjt = (j == 0)? 1.0 : gretl_matrix_get(K->x, K->t, j-1);
                } else {
                    xjt = gretl_matrix_get(K->x, K->t, j);
                }
                bji = gretl_matrix_get(K->BT, j, i);
                if (bji != 0) {
                    /* here we'll implicitly take 0 * NA as 0 */
                    if (na(xjt)) {
                        missvals++;
                    } else {
                        bxi += xjt * bji;
                    }
                }
            }
        }
        if (mod == GRETL_MOD_DECREMENT) {
            targ->val[i] -= bxi;
        } else {
            targ->val[i] += bxi;
        }
    }

    return missvals;
}

/* Apparatus for handling the case where the number of observables,
   K->n, is > 1 and we hit a partially missing obs at time t. Instead
   of discarding the entire observation it's optimal to make use of
   the components of y_t for which we have valid values. This requires
   shrinking the observation-related matrices ZT and VY, plus the
   vector v_t and all other matrices whose dimensions involve K->n.
   (See the brief discussion in Koopman's 1997 JASA paper.)

   In addition, these matrices must be restored to their original sizes
   when we've finished dealing with period t.
*/

/* For use in the smoother: retrieve the effective number of
   observables at time t (recorded in the forward pass).
*/

static int get_effective_n (kalman *K, int t)
{
    if (K->nt == NULL || K->nt[t] == 0) {
        /* no missing values */
        return K->n;
    } else {
        return K->nt[t];
    }
}

/* We have two cases to handle here, the more complex one being when
   the disturbances are correlated across the observation and
   state-transition equations.
*/

static void shrink_obsvar (kalman *K, int nt)
{
    double gij;
    int i, j, k;

    if (K->G == NULL) {
        /* the relatively simple case */
        K->saveVY = K->VY;
        K->VY = gretl_matrix_alloc(nt, nt);

        for (j=0, k=0; j<K->n; j++) {
            if (!na(K->vt->val[j])) {
                for (i=0; i<K->n; i++) {
                    if (!na(K->vt->val[i])) {
                        gij = gretl_matrix_get(K->saveVY, i, j);
                        K->VY->val[k++] = gij;
                    }
                }
            }
        }
    } else {
        /* the more complex case */
        K->saveG = K->G;
        K->G = gretl_matrix_alloc(nt, K->p);

        for (i=0, k=0; i<K->n; i++) {
            if (!na(K->vt->val[i])) {
                for (j=0; j<K->p; j++) {
                    gij = gretl_matrix_get(K->saveG, i, j);
                    gretl_matrix_set(K->G, k, j, gij);
                }
                k++;
            }
        }
        gretl_matrix_multiply_mod(K->G, GRETL_MOD_NONE,
                                  K->G, GRETL_MOD_TRANSPOSE,
                                  K->VY, GRETL_MOD_NONE);
        if (K->HG != NULL) {
            gretl_matrix_reuse(K->HG, K->r, nt);
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                      K->G, GRETL_MOD_TRANSPOSE,
                                      K->HG, GRETL_MOD_NONE);
        }
    }
}

static int obsvar_resize_needed (kalman *K, int smtype)
{
    if (K->VY != NULL) {
        return smtype == SM_TYPE_NONE ||
            smtype == SM_DIST_BKWD ||
            smtype == SM_DIST_FRWD ||
            smtype == SM_STATE_INI;
    } else {
        return 0;
    }
}

/* The inverse operation of shrink_obsvar() above */

static void unshrink_obsvar (kalman *K)
{
    if (K->G == NULL) {
        gretl_matrix_free(K->VY);
        K->VY = K->saveVY;
        K->saveVY = NULL;
    } else {
        gretl_matrix_free(K->G);
        K->G = K->saveG;
        K->saveG = NULL;
        gretl_matrix_reuse(K->VY, K->n, K->n);
        gretl_matrix_multiply_mod(K->G, GRETL_MOD_NONE,
                                  K->G, GRETL_MOD_TRANSPOSE,
                                  K->VY, GRETL_MOD_NONE);
        if (K->HG != NULL) {
            gretl_matrix_reuse(K->HG, K->r, K->n);
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                      K->G, GRETL_MOD_TRANSPOSE,
                                      K->HG, GRETL_MOD_NONE);
        }
    }
}

static void shrink_vt (kalman *K, int nt, int smtype)
{
    double *tmp = malloc(nt * sizeof *tmp);
    int i, k;

    if (smtype < SM_DIST_BKWD) {
        /* We need to shrink ZT too */
        int j;

        K->saveZT = K->ZT; /* save the original matrix */
        K->ZT = gretl_matrix_alloc(K->r, nt);

        for (j=0, k=0; j<K->n; j++) {
            if (!na(K->vt->val[j])) {
                for (i=0; i<K->r; i++) {
                    K->ZT->val[k++] = gretl_matrix_get(K->saveZT, i, j);
                }
            }
        }
    }

    if (obsvar_resize_needed(K, smtype)) {
        shrink_obsvar(K, nt);
    }

    /* K-member workspace matrices */
    gretl_matrix_reuse(K->Ft, nt, nt);
    gretl_matrix_reuse(K->iFt, nt, nt);
    gretl_matrix_reuse(K->Mt, K->r, nt);
    gretl_matrix_reuse(K->Kt, K->r, nt);
    gretl_matrix_reuse(K->PZ, K->r, nt);

    /* remove missing rows of vt */
    for (i=0, k=0; i<K->n; i++) {
        if (!na(K->vt->val[i])) {
            tmp[k++] = K->vt->val[i];
        }
    }
    memcpy(K->vt->val, tmp, nt * sizeof(double));
    gretl_matrix_reuse(K->vt, nt, 1);
    free(tmp);
}

static void unshrink_vt (kalman *K, int smtype)
{
    if (smtype < SM_DIST_BKWD) {
        /* We need to restore ZT too */
        gretl_matrix_free(K->ZT);
        K->ZT = K->saveZT;
        K->saveZT = NULL;
    }
    gretl_matrix_reuse(K->vt, K->n, 1);

    /* K-member workspace matrices */
    gretl_matrix_reuse(K->Ft, K->n, K->n);
    gretl_matrix_reuse(K->iFt, K->n, K->n);
    gretl_matrix_reuse(K->Mt, K->r, K->n);
    gretl_matrix_reuse(K->Kt, K->r, K->n);
    gretl_matrix_reuse(K->PZ, K->r, K->n);

    if (obsvar_resize_needed(K, smtype)) {
        unshrink_obsvar(K);
    }
}

/* Shrink the observation vector and associated matrices in the case
   where we encounter a partially missing observation. Used only on a
   filtering pass.
*/

static void shrink_obs_to_ok (kalman *K, int nt)
{
    /* includes back up and replacement of ZT */
    shrink_vt(K, nt, SM_TYPE_NONE);

    if (basic_smoothing(K)) {
        /* make a record for the smoother */
        if (K->nt == NULL) {
            K->nt = calloc(K->N, sizeof *K->nt);
        }
        if (K->nt != NULL) {
            K->nt[K->t] = nt;
        }
    }
}

/* Compute the one-step ahead forecast error:

   v_t = y_t - B_t*x_t - Z_t*a_t

   if @missmap is non-NULL, return the number of observables for which
   there's a valid observation at time t; otherwise return 0 if any
   elements of y_t are missing.
*/

static int compute_forecast_error (kalman *K, char *missmap)
{
    int i, nt = K->n;

    if (missmap != NULL) {
        memset(missmap, 0, K->n);
    }

    /* initialize v_t to y_t */
    for (i=0; i<K->n; i++) {
        K->vt->val[i] = gretl_matrix_get(K->y, K->t, i);
    }

    if (K->BT != NULL) {
        /* subtract effect of exogenous terms, if any */
        kalman_do_Bx(K, K->vt, GRETL_MOD_DECREMENT);
    }

    /* check for any missing values in v_t */
    for (i=0; i<K->n; i++) {
        if (na(K->vt->val[i])) {
            if (missmap != NULL) {
                missmap[i] = 1;
            }
            nt--;
        }
    }

    if (nt > 0 && nt < K->n) {
        if (missmap != NULL) {
            shrink_obs_to_ok(K, nt);
        } else {
            nt = 0; /* skip this observation */
        }
    }

    if (nt > 0) {
        /* subtract contribution from state */
        gretl_matrix_multiply_mod(K->ZT,  GRETL_MOD_TRANSPOSE,
                                  K->a0, GRETL_MOD_NONE,
                                  K->vt, GRETL_MOD_DECREMENT);
    }

    return nt;
}

/* Given a unified function to update one or more of the potentially
   time-varying matrices, try to figure out which matrix or matrices
   are actually modified by this function.
*/

static int check_for_matrix_updates (kalman *K, ufunc *uf)
{
    char **lines;
    int nlines = 0;

    if (K->varying != NULL) {
        free(K->varying);
        K->varying = NULL;
    }

    if (uf == NULL) {
        uf = get_user_function_by_name(K->matcall);
        if (uf == NULL) {
            gretl_errmsg_sprintf("Couldn't find function '%s'", K->matcall);
            return E_DATA;
        }
    }

    K->varying = calloc(K_N_MATCALLS, 1);

    lines = gretl_function_retrieve_code(uf, &nlines);

    if (lines != NULL) {
        const char *bname = fn_param_name(uf, 0);
        char test[VNAMELEN+1];
        const char *s;
        int n = strlen(bname) + 1;
        int i, j;

        sprintf(test, "%s.", bname);
        for (i=0; i<nlines; i++) {
            if (!strncmp(lines[i], test, n)) {
                for (j=K_T; j<=K_m; j++) {
                    s = kalman_matrix_name(j);
                    if (!strncmp(lines[i] + n, s, strlen(s))) {
                        fprintf(stderr, "matrix %s is varying\n", s);
                        K->varying[j] = 1;
                        break;
                    }
                }
            }
        }
        free(lines);
    }

    return 0;
}

/* Function to update any time-varying matrices, for use with a kalman
   bundle. Bypasses the regular "genr" apparatus, passing the attached
   bundle directly to the given user function after it has been found
   by name.
*/

static int kalman_update_matrices (kalman *K, PRN *prn)
{
    ufunc *uf;
    fncall *fc;
    int err = 0;

    uf = get_user_function_by_name(K->matcall);
    if (uf == NULL) {
        gretl_errmsg_sprintf("Couldn't find function '%s'", K->matcall);
        return E_DATA;
    }

    if (K->varying == NULL) {
        check_for_matrix_updates(K, uf);
    }

    fc = fncall_new(uf, 0);
    err = push_anon_function_arg(fc, GRETL_TYPE_BUNDLE_REF, K->b);

    if (!err) {
        err = gretl_function_exec(fc, GRETL_TYPE_NONE, NULL, NULL,
                                  NULL, prn);
    }
    if (err) {
        fprintf(stderr, "kalman_update_matrices: call='%s', err=%d\n",
                K->matcall, err);
    }

    return err;
}

/* If we have any time-varying coefficient matrices, refresh these for
   the current time step. This is called on a forward filtering pass.
*/

static int kalman_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
        &K->T, &K->BT, &K->ZT, &K->VS, &K->VY, &K->mu
    };
    int cross_update = 0;
    int dkvar_update = 0;
    int i, err = 0;

    if (kalman_xcorr(K)) {
        mptr[3] = &K->H;
        mptr[4] = &K->G;
    } else if (kalman_dkvar(K)) {
	mptr[3] = &K->Q;
	mptr[4] = &K->R;
    }

    if (K->matcall != NULL) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<K_N_MATCALLS && !err; i++) {
        if (matrix_is_varying(K, i)) {
            if (kalman_xcorr(K) && (i == K_VS || i == K_VY)) {
                /* handle revised H and/or G */
                cross_update = 1;
	    } else if (kalman_dkvar(K) && i == K_VS) {
		/* handle revised Q and/or R */
		dkvar_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], i);
            }
            if (err) {
                fprintf(stderr, "kalman_refresh_matrices: err = %d at t = %d\n",
                        err, K->t);
            }
        }
    }

    if (!err && K->step != NULL) {
        /* keep a record of T and/or Z' at the given time step */
        if (K->step->T != NULL) {
            load_to_vec(K->step->T, K->T, K->t);
        }
        if (K->step->ZT != NULL) {
            load_to_vec(K->step->ZT, K->ZT, K->t);
        }
    }

    if (!err) {
	if (cross_update) {
	    /* cross-correlated case */
	    err = kalman_update_crossinfo(K, UPDATE_STEP);
	} else if (dkvar_update) {
	    /* Durbin-Koopman case */
	    err = kalman_update_dkvar(K, UPDATE_STEP);
	}
    }

    return err;
}

/* Variant of the above for use when Koopman-smoothing */

static int ksmooth_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
        &K->VS, &K->VY
    };
    int idx[] = {
        K_VS, K_VY
    };
    int var_update = 0;
    int i, ii, err = 0;

    if (kalman_xcorr(K)) {
        mptr[0] = &K->H;
        mptr[1] = &K->G;
    } else if (kalman_dkvar(K)) {
        mptr[0] = &K->Q;
        mptr[1] = &K->R;
    }

    if (K->matcall != NULL) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<2 && !err; i++) {
        ii = idx[i];
        if (matrix_is_varying(K, ii)) {
            if (kalman_xcorr(K) && (ii == K_VS || ii == K_VY)) {
                /* handle revised H and/or G */
                var_update = 1;
	    } else if (kalman_dkvar(K) && ii == K_VS) {
		var_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], ii);
            }
            if (err) {
                fprintf(stderr, "ksmooth_refresh_matrices: err = %d at t = %d\n",
                        err, K->t);
            }
        }
    }

    if (!err && var_update) {
	if (K->vartype == DJ_VAR) {
	    err = kalman_update_crossinfo(K, UPDATE_STEP);
	} else {
	    err = kalman_update_dkvar(K, UPDATE_STEP);
	}
    }

    return err;
}

/* exact initial iteration for multivariate y_t */

static int koopman_exact_general (kalman *K,
                                  double *ldet,
                                  double *qt)
{
    gretl_matrix_block *B;
    gretl_matrix *Mk;
    gretl_matrix *Fmk;
    gretl_matrix *Fmt;
    gretl_matrix *MFk;
    gretl_matrix *PkZ;
    gretl_matrix *V, *J;
    gretl_matrix *l = NULL;
    double xij, rlj;
    int i, j, n, ns;
    size_t sz;
    int err = 0;

    n = K->vt->rows; /* may be < K->n */

    B = gretl_matrix_block_new(&Mk,  K->r, n,
                               &Fmk, n, n,
                               &Fmt, n, n,
                               &MFk, K->r, n,
                               &PkZ, K->r, n,
                               &V,   n, n,
                               &J,   n, n,
                               NULL);

    if (B == NULL) {
        err = E_ALLOC;
    } else {
        /* Mk = T * P∞ * Z' */
        gretl_matrix_multiply(K->Pk0, K->ZT, PkZ);
        gretl_matrix_multiply(K->T, PkZ, Mk);
        /* eigenanalysis */
        l = gretl_gensymm_eigenvals(K->Fk, K->Ft, V, &err);
    }

    if (err) {
        goto bailout;
    }

    /* @ns counts the "zero" eigenvalues */
    ns = 0;
    for (i=0; i<n; i++) {
        if (l->val[i] < 1.0e-12) {
            ns++;
        } else {
            break;
        }
    }

    /* J = V[,1:ns] */
    gretl_matrix_reuse(J, n, ns);
    sz = n * ns * sizeof(double);
    memcpy(J->val, V->val, sz);
    /* Fmt = F_{*,t}^{-} = J*J' */
    gretl_matrix_multiply_mod(J, GRETL_MOD_NONE,
                              J, GRETL_MOD_TRANSPOSE,
                              Fmt, GRETL_MOD_NONE);

    /* J = V[,ns+1:] ./ sqrt(l[ns+1:]') */
    gretl_matrix_reuse(J, n, n - ns);
    sz = K->n * J->cols * sizeof(double);
    memcpy(J->val, V->val + n*ns, sz);
    for (j=0; j<J->cols; j++) {
        rlj = sqrt(l->val[ns+j]);
        for (i=0; i<n; i++) {
            xij = gretl_matrix_get(J, i, j);
            gretl_matrix_set(J, i, j, xij / rlj);
        }
    }
    /* Fmk = F_{∞,t}^{-} = J*J' */
    gretl_matrix_multiply_mod(J, GRETL_MOD_NONE,
                              J, GRETL_MOD_TRANSPOSE,
                              Fmk, GRETL_MOD_NONE);
    /* copy to iFt for recording (?) */
    fast_copy_values(K->iFt, Fmk);

    /* MFk = M∞ * Fm∞ */
    gretl_matrix_multiply(Mk, Fmk, MFk);

    /* K★ = M★ * Fmt + M∞ * F∞ */
    fast_copy_values(K->Kt, MFk);
    gretl_matrix_multiply_mod(K->Mt, GRETL_MOD_NONE,
                              Fmt, GRETL_MOD_NONE,
                              K->Kt, GRETL_MOD_CUMULATE);

    /* C★ = M★ Kt' (first component) */
    gretl_matrix_multiply_mod(K->Mt, GRETL_MOD_NONE,
                              K->Kt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_NONE);

    /* Mt <- M★ - M∞ * F∞ * Fmt (Wrong? Should be F★?) */
#if 0
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              Fmt, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_DECREMENT);
#else
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              K->Ft, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_DECREMENT);
#endif

    /* Ct += M∞ * F∞ * (M★ - M∞ * F∞ * Fmt)' (second component) */
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_CUMULATE);

    /* C∞ = M∞ * F∞ * M∞' */
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              Mk,  GRETL_MOD_TRANSPOSE,
                              K->Ck, GRETL_MOD_NONE);

    gretl_matrix_add_to(Fmk, Fmt);
    *ldet = gretl_matrix_log_determinant(Fmk, &err);
    *qt = gretl_scalar_qform(K->vt, Fmt, &err);

 bailout:

    gretl_matrix_block_destroy(B);

    return err;
}

/* exact initial iteration for univariate y_t,
   when F_\infty ("Fk") > 0
*/

static int koopman_exact_nonsingular (kalman *K)
{
    gretl_matrix *Mk;
    gretl_matrix *PkZ;
    int n, err = 0;

    n = K->vt->rows; /* may be < K->n */
    Mk  = gretl_matrix_alloc(K->r, n);
    PkZ = gretl_matrix_alloc(K->r, n);

    /* M∞ = T * P∞ * Z' */
    gretl_matrix_multiply(K->Pk0, K->ZT, PkZ);
    gretl_matrix_multiply(K->T, PkZ, Mk);

    /* Kt = M∞ ./ F∞[1] */
    fast_copy_values(K->Kt, Mk);
    gretl_matrix_divide_by_scalar(K->Kt, K->Fk->val[0]);

    /* Ct = Mt * Kt' */
    gretl_matrix_multiply_mod(K->Mt, GRETL_MOD_NONE,
                              K->Kt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_NONE);

    /* Mt <- Mt - Kt * Ft */
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                              K->Ft, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_DECREMENT);

    /* Ct += Kt * (Mt - Kt * Ft)' */
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_CUMULATE);

    /* Ck = Kt * M∞' */
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                              Mk, GRETL_MOD_TRANSPOSE,
                              K->Ck, GRETL_MOD_NONE);

    gretl_matrix_free(Mk);
    gretl_matrix_free(PkZ);

    return err;
}

static double max_val (const gretl_matrix *m)
{
    int i, n = m->rows * m->cols;
    double ret = 0;

    for (i=0; i<n; i++) {
        if (m->val[i] > ret) {
            ret = m->val[i];
        }
    }

    return ret;
}

static void P_infty_update (kalman *K, int C_zero)
{
    gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->Pk0,
                       K->Pk1, GRETL_MOD_NONE);
    if (!C_zero) {
        gretl_matrix_subtract_from(K->Pk1, K->Ck);
    }
    fast_copy_values(K->Pk0, K->Pk1);

    if (K->d == 0 && max_val(K->Pk0) < K_TINY) {
        K->d = K->t + 1;
    }
}

/* Handling of missing y_t or x_t, for the case of
   univariate y_t (or all y_t values missing).
*/

static void handle_missing_obs (kalman *K)
{
    /* state update: a1 = T*a0 */
    gretl_matrix_multiply(K->T, K->a0, K->a1);
    /* handle stconst if present */
    if (K->mu != NULL) {
        gretl_matrix_add_to(K->a1, K->mu);
    }
    fast_copy_values(K->a0, K->a1);

    /* var(state) update: P1 = T*P0*T' + VS (C = 0) */
    fast_copy_values(K->P1, K->VS);
    gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->P0,
                       K->P1, GRETL_MOD_CUMULATE);
    fast_copy_values(K->P0, K->P1);

    if (K->exact && max_val(K->Pk0) > K_TINY) {
        /* update P∞ (what about K->d?) */
        gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->Pk0,
                           K->Pk1, GRETL_MOD_NONE);
        fast_copy_values(K->Pk0, K->Pk1);
    }

    /* record stuff if wanted */
    if (K->F != NULL) {
        if (kalman_smoothing(K)) {
            set_row_to_value(K->F, K->t, 0.0);
        } else {
            set_row_to_value(K->F, K->t, NADBL); /* 0 ? */
        }
    }
    if (K->LL != NULL) {
        gretl_vector_set(K->LL, K->t, 0.0); /* ? */
    }
    if (K->K != NULL) {
        set_row_to_value(K->K, K->t, 0.0);
    }
    if (K->V != NULL) {
        if (kalman_smoothing(K)) {
            set_row_to_value(K->V, K->t, 0.0);
        } else {
            set_row_to_value(K->V, K->t, NADBL); /* 0 ? */
        }
    }
}

/**
 * kalman_forecast:
 * @K: pointer to Kalman struct.
 * @prn: printing apparatus (or NULL).
 *
 * Generates a series of one-step ahead forecasts for y, based on
 * information in the kalman struct @K.
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_forecast (kalman *K, PRN *prn)
{
    double ll0 = K->n * LN_2_PI;
    double sumldet = 0;
    char *missmap = NULL;
    int err = 0;

#if KDEBUG
    fprintf(stderr, "kalman_forecast: N = %d\n", K->N);
#endif

#if USE_INCOMPLETE_OBS
    if (!kalman_smoothing(K) || basic_smoothing(K)) {
        /* we won't do this for disturbance smoothing, yet */
        missmap = calloc(K->n, 1);
    }
#endif

    if (kalman_diffuse(K) && !K->exact) {
        K->d = K->r;
    } else {
        K->d = 0;
    }
    K->SSRw = K->loglik = 0.0;
    K->s2 = NADBL;
    K->okN = K->N;
    set_kalman_running(K);

    for (K->t = 0; K->t < K->N && !err; K->t += 1) {
        int Kt_done = 0;
        int nt = K->n;
        double llt = NADBL;
        double ldet = 0, qt = 0;

#if KDEBUG > 1
        kalman_print_state(K);
#endif

        if (K->A != NULL || K->P != NULL) {
            kalman_record_state(K);
        }

        if (filter_is_varying(K)) {
            /* we have time-varying coefficients */
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                K->loglik = NADBL;
                break;
            }
        }

        /* calculate v_t, checking for missing values */
        nt = compute_forecast_error(K, missmap);
        if (nt == 0) {
            /* skip this observation */
            K->okN -= 1;
            handle_missing_obs(K);
            continue;
        }

        /* calculate F_t = ZPZ' [+ VY] */
        if (K->VY != NULL) {
            fast_copy_values(K->Ft, K->VY);
            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
                               K->Ft, GRETL_MOD_CUMULATE);
        } else {
            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
                               K->Ft, GRETL_MOD_NONE);
        }

        /* calculate M_t = TPZ' [+ HG'] */
        gretl_matrix_multiply(K->P0, K->ZT, K->PZ);
        gretl_matrix_multiply(K->T, K->PZ, K->Mt);
        if (K->HG != NULL) {
            gretl_matrix_add_to(K->Mt, K->HG);
        }

        if (K->exact && max_val(K->Pk0) > K_TINY) {
            /* handle initial exact iterations */
            int C_zero = 0;

            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE,
                               K->Pk0, K->Fk, GRETL_MOD_NONE);
            if (K->n == 1 && K->Fk->val[0] < K_TINY) {
                /* univariate F∞ = 0 */
                ldet = log(K->Ft->val[0]);
                K->iFt->val[0] = 1.0 / K->Ft->val[0];
                qt = K->vt->val[0] * K->vt->val[0] * K->iFt->val[0];
                C_zero = 1;
            } else if (K->n == 1) {
                /* univariate F∞ nonsingular */
                ldet = log(K->Fk->val[0]);
                K->iFt->val[0] = 0;
                qt = 0;
                err = koopman_exact_nonsingular(K);
                Kt_done = 1;
            } else {
                /* multivariate case */
                err = koopman_exact_general(K, &ldet, &qt);
                Kt_done = 1;
            }
            P_infty_update(K, C_zero);
        } else {
            /* standard Kalman procedure */
            fast_copy_values(K->iFt, K->Ft);
            err = gretl_invert_symmetric_matrix2(K->iFt, &ldet);
            if (err) {
                fprintf(stderr, "kalman_forecast: failed to invert Ft\n");
            } else {
                qt = gretl_scalar_qform(K->vt, K->iFt, &err);
            }
        }

        if (K->F != NULL) {
            /* we're recording F_t for all t */
            if (kalman_smoothing(K)) {
                /* record inverse */
                load_to_vech(K->F, K->iFt, nt, K->t);
            } else {
                /* record F_t itself */
                load_to_vech(K->F, K->Ft, nt, K->t);
            }
        }

        /* determine and record loglikelihood */
        if (err) {
            K->loglik = NADBL;
            break;
        } else {
            llt = -0.5 * (ll0 + ldet + qt);
            if (na(llt)) {
                K->loglik = NADBL;
                break;
            }
            K->SSRw += qt;
            K->loglik += llt;
            sumldet += ldet;
        }
        if (K->LL != NULL) {
            gretl_vector_set(K->LL, K->t, llt);
        }

        if (!Kt_done) {
            /* Calculate gain K_t = M_t F_t^{-1}, and C matrix.
               Note that this will have been done correctly already
               in two of the "exact initial" cases above.
            */
            gretl_matrix_multiply(K->Mt, K->iFt, K->Kt);
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                      K->Mt, GRETL_MOD_TRANSPOSE,
                                      K->Ct, GRETL_MOD_NONE);
        }

        if (!err && K->K != NULL) {
            /* record the gain */
            load_to_vec(K->K, K->Kt, K->t);
        }

        /* update state: a1 = T a0 + K_t v_t [+ mu] */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->vt, GRETL_MOD_NONE,
                                  K->a1, GRETL_MOD_CUMULATE);
        fast_copy_values(K->a0, K->a1);

        /* update var(state): P1 = TPT' + VS - C */
        fast_copy_values(K->P1, K->VS);
        gretl_matrix_qform(K->T, GRETL_MOD_NONE,
                           K->P0, K->P1, GRETL_MOD_CUMULATE);
        gretl_matrix_subtract_from(K->P1, K->Ct);
        fast_copy_values(K->P0, K->P1);

        /* record forecast errors if wanted */
        if (!err && K->V != NULL) {
            if (missmap != NULL) {
                load_to_row_special(K->V, K->vt, K->t, missmap);
            } else {
                load_to_row(K->V, K->vt, K->t);
            }
        }

        if (nt > 0 && nt < K->n) {
            /* restore certain matrices to full size */
            unshrink_vt(K, SM_TYPE_NONE);
        }
    }

    set_kalman_stopped(K);

    if (missmap != NULL) {
	free(missmap);
    }

    if (na(K->loglik)) {
        err = E_NAN;
    } else if (kalman_arma_ll(K)) {
        double ll1 = 1.0 + LN_2_PI + log(K->SSRw / K->okN);

        K->loglik = -0.5 * (K->okN * ll1 + sumldet);
    } else {
        K->s2 = K->SSRw / (K->n * K->okN - K->d);
    }

#if KDEBUG
    fprintf(stderr, "kalman_forecast: err=%d, ll=%#.12g, d=%d\n",
            err, K->loglik, K->d);
#endif

    return err;
}

struct K_input_mat {
    int sym;
    const char *name;
};

/* mapping to names used in setting elements of kalman bundle */

struct K_input_mat K_input_mats[] = {
    { K_y,  "obsy" },
    { K_ZT, "obsymat" },
    { K_x,  "obsx" },
    { K_BT, "obsxmat" },
    { K_VY, "obsvar" },
    { K_T,  "statemat" },
    { K_VS, "statevar" },
    { K_m,  "stconst" },
    { K_a,  "inistate" },
    { K_P,  "inivar" },
    { K_R,  "statesel" }
};

int extra_mats[] = {
    K_BT,
    K_VY,
    K_m,
    K_x,
    K_a,
    K_P,
    K_R
};

static int n_extra_mats = G_N_ELEMENTS(extra_mats);

static int obsy_check (kalman *K)
{
    if (K->y == NULL) {
        return missing_matrix_error("obsy");
    } else if (K->y->rows != K->N || K->y->cols != K->n) {
        fprintf(stderr, "obsy_check: K->y should be %d x %d, is %d x %d\n",
                K->N, K->n, gretl_matrix_rows(K->y),
                gretl_matrix_cols(K->y));
        return E_NONCONF;
    } else {
        return 0;
    }
}

static int kalman_bundle_recheck_matrices (kalman *K, PRN *prn)
{
    int err = 0;

    K->flags |= KALMAN_CHECK;

    if (filter_is_varying(K)) {
        err = kalman_update_matrices(K, prn);
    }

    K->flags ^= KALMAN_CHECK;

    if (!err && (K->ZT == NULL || K->T == NULL || K->VS == NULL)) {
        fprintf(stderr, "kalman_bundle_kalman_recheck_matrices: Z=%p, T=%p, Q=%p\n",
                K->ZT, K->T, K->VS);
        err = missing_matrix_error(NULL);
    }

    if (err) {
        return err;
    }

    /* redundant? */
    kalman_set_dimensions(K);

    if (gretl_matrix_rows(K->T) != K->r ||
        gretl_matrix_rows(K->BT) != K->k) {
        err = E_NONCONF;
    } else if (!kalman_simulating(K)) {
        err = obsy_check(K);
    }

    if (!err && K->vartype != STD_VAR) {
        err = kalman_revise_variance(K);
    }

    if (!err) {
        err = kalman_check_dimensions(K);
    }

    if (!err) {
        if (K->aini != NULL) {
            gretl_matrix_copy_values(K->a0, K->aini);
        } else {
            gretl_matrix_zero(K->a0);
        }
        if (K->Pini != NULL) {
            gretl_matrix_copy_values(K->P0, K->Pini);
        } else {
            err = construct_Pini(K);
        }
    }

    return err;
}

static const char *kalman_matrix_name (int sym)
{
    int i;

    for (i=0; i<K_MMAX; i++) {
        if (K_input_mats[i].sym == sym) {
            return K_input_mats[i].name;
        }
    }

    /* failed */
    return "matrix";
}

static int kalman_ensure_output_matrices (kalman *K)
{
    int err = 0;

    if (K->V == NULL) {
        K->V = gretl_null_matrix_new();
    }
    if (K->F == NULL) {
        K->F = gretl_null_matrix_new();
    }
    if (K->A == NULL) {
        K->A = gretl_null_matrix_new();
    }
    if (K->P == NULL) {
        K->P = gretl_null_matrix_new();
    }
    if (K->K == NULL) {
        K->K = gretl_null_matrix_new();
    }

    if (K->V == NULL || K->F == NULL || K->A == NULL ||
        K->P == NULL || K->K == NULL) {
        err = E_ALLOC;
    }

    return err;
}

int kalman_run (kalman *K, PRN *prn, int *errp)
{
    int err = kalman_ensure_output_matrices(K);

    if (!err) {
        gretl_matrix_zero(K->vt);
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err && K->LL == NULL) {
        K->LL = gretl_matrix_alloc(K->N, 1);
        if (K->LL == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
	if (K->flags & KALMAN_UNI) {
	    err = kfilter_univariate(K, prn);
	} else {
	    err = kalman_forecast(K, prn);
	}
    }

    if (err != E_NAN) {
        *errp = err;
    } else {
        /* we'll flag E_NAN with a return value of 1 but
           won't count it as a 'true' error */
        *errp = 0;
    }

    return err;
}

/* Implements the userland kfilter() function */

int kalman_bundle_run (gretl_bundle *b, PRN *prn, int *errp)
{
    kalman *K = gretl_bundle_get_private_data(b);

    if (getenv("KALMAN_UNI") != NULL) {
	K->flags |= KALMAN_UNI;
    }

    K->b = b; /* attach bundle pointer */

    return kalman_run(K, prn, errp);
}

/* Copy row @t from @src into @targ; or add row @t of @src to @targ;
   or subtract row @t of @src from @targ.  We allow the possibility
   that the length of vector @targ is less than the number of columns
   in @src, but not the converse.
*/

static int load_from_row (gretl_vector *targ,
                          const gretl_matrix *src,
                          int t, GretlMatrixMod mod)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols) {
        fprintf(stderr, "load_from_row: targ length = %d, but src "
                "has %d columns\n", n, src->cols);
        return E_NONCONF;
    }

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, t, i);
        if (mod == GRETL_MOD_CUMULATE) {
            targ->val[i] += x;
        } else if (mod == GRETL_MOD_DECREMENT) {
            targ->val[i] -= x;
        } else {
            targ->val[i] = x;
        }
    }

    return 0;
}

/* As load_from_row(), except that a column offset, @j, is supported
   for the reading of a row from @src, and we don't support the @mod
   option.
*/

static int load_from_row_offset (gretl_vector *targ,
                                 const gretl_matrix *src,
                                 int t, int j)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols - j) {
        return E_NONCONF;
    }

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, t, i + j);
        targ->val[i] = x;
    }

    return 0;
}

/* Row @t of @src represents the vech of an n x n matrix: extract the
   row and apply the inverse operation of vech to reconstitute the
   matrix in @targ -- or subtract the newly reconstituted matrix
   from @targ.
*/

static void load_from_vech (gretl_matrix *targ, const gretl_matrix *src,
                            int n, int t, int mod)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            x = gretl_matrix_get(src, t, m++);
            if (mod == GRETL_MOD_DECREMENT) {
                x = gretl_matrix_get(targ, i, j) - x;
            }
            gretl_matrix_set(targ, i, j, x);
            if (i != j) {
                gretl_matrix_set(targ, j, i, x);
            }
        }
    }
}

/* Column @t of @src represents the vech of an n x n matrix: extract the
   the column and apply the inverse operation of vech to reconstitute the
   matrix in @targ.
*/

static void load_from_vechT (gretl_matrix *targ, const gretl_matrix *src,
                             int n, int t)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            x = gretl_matrix_get(src, m++, t);
            gretl_matrix_set(targ, i, j, x);
            if (i != j) {
                gretl_matrix_set(targ, j, i, x);
            }
        }
    }
}

/* Row @t of @src represents the vec of a certain matrix: extract the
   row and reconstitute the matrix in @targ.
*/

static int load_from_vec (gretl_matrix *targ,
                          const gretl_matrix *src,
                          int t)
{
    int i, k = targ->rows * targ->cols;

    for (i=0; i<k; i++) {
        targ->val[i] = gretl_matrix_get(src, t, i);
    }

    return 0;
}

/* For disturbance smoothing: ensure we have on hand matrices that are
   correctly sized to hold estimates of the variance of the
   disturbance(s) in the state and (if applicable) observation
   equations.
*/

static int maybe_resize_dist_mse (kalman *K,
                                  gretl_matrix **Vwt,
                                  gretl_matrix **Vut)
{
    int n = K->VY == NULL ? 0 : K->n;
    int k, err = 0;

    /* combined results: how many columns do we need? */
    k = K->r + n;

    if (K->Vsd == NULL) {
        K->Vsd = gretl_matrix_alloc(K->N, k);
        if (K->Vsd == NULL) {
            err = E_ALLOC;
        }
    } else if (K->Vsd->rows != K->N || K->Vsd->cols != k) {
        err = gretl_matrix_realloc(K->Vsd, K->N, k);
    }

    if (!err) {
        /* step-t square state matrix */
        *Vwt = gretl_matrix_alloc(K->r, K->r);
        if (*Vwt == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && n > 0) {
        /* step-t square obs matrix */
        *Vut = gretl_matrix_alloc(K->n, K->n);
        if (*Vut == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

static int retrieve_Tt (kalman *K, int t)
{
    if (K->step == NULL || K->step->T == NULL) {
        return E_DATA;
    } else {
        return load_from_vec((gretl_matrix *) K->T, K->step->T, t);
    }
}

static int retrieve_Zt (kalman *K, int t)
{
    if (K->step == NULL || K->step->ZT == NULL) {
        return E_DATA;
    } else {
        return load_from_vec((gretl_matrix *) K->ZT, K->step->ZT, t);
    }
}

/* For use with smoothing: load what we need for step @t, from the
   record that was kept on the prior forecasting pass. What we need
   from the forward pass depends on what exactly we're computing,
   which is conveyed by @smtype. If @pnt is non-NULL, use it to pass
   back the effective size of the observables vector.

   We need to handle some complications here. In particular we must
   update anything that's time-varying, and (unless we're doing
   disturbance smoothing) deal with any incomplete observations, which
   have implications for the dimensions of various matrices.
*/

static int load_filter_data (kalman *K, int t, int *pnt, int smtype)
{
    int nt, err = 0;

    /* load the forecast error */
    load_from_row(K->vt, K->V, t, GRETL_MOD_NONE);

    if (smtype == SM_DIST_BKWD) {
        /* disturbances */
        if (filter_is_varying(K)) {
            K->t = t;
            ksmooth_refresh_matrices(K, NULL);
        }
    } else {
        /* state: get T_t and/or Z_t if need be */
        if (matrix_is_varying(K, K_T)) {
            err = retrieve_Tt(K, t);
        }
        if (!err && matrix_is_varying(K, K_ZT)) {
            err = retrieve_Zt(K, t);
        }
        if (err) {
            return err;
        }
    }

    /* check for an incomplete observation */
    nt = get_effective_n(K, t);
    if (nt < K->n) {
        shrink_vt(K, nt, smtype);
    }
    if (pnt != NULL) {
	*pnt = nt;
    }

    if (smtype < SM_DIST_BKWD) {
        /* load the state and its MSE */
        load_from_row(K->a0, K->A, t, GRETL_MOD_NONE);
        load_from_vech(K->P0, K->P, K->r, t, GRETL_MOD_NONE);
    }

    if (smtype == SM_STATE_INI) {
        /* load P∞ */
        load_from_vechT(K->Pk0, K->PK, K->r, t);
    } else if (smtype == SM_STATE_STD || smtype == SM_DIST_BKWD) {
        /* load the gain and F^{-1} */
        load_from_vec(K->Kt, K->K, t);
        load_from_vech(K->iFt, K->F, nt, t, GRETL_MOD_NONE);
    }

    return err;
}

/* wrapper for extra info that's needed for exact intial
   disturbance smoothing */

typedef struct dsinfo_ {
    gretl_matrix *R;
    gretl_matrix *D;
    gretl_matrix *u;
    gretl_matrix *n1;
    gretl_matrix *r1;
    gretl_matrix *Vwt;
    gretl_matrix *Vut;
    int DKstyle;
} dsinfo;

/* Calculate the variance of the smoothed disturbances for the
   cross-correlated case. See Koopman, Shephard and Doornik (1998),
   page 19, var(\varepsilon_t|Y_n).
*/

static int combined_dist_variance (kalman *K,
                                   gretl_matrix *D,
                                   gretl_matrix *Nt,
                                   gretl_matrix *Vwt,
                                   gretl_matrix *Vut,
                                   gretl_matrix_block *BX,
                                   int DKstyle)
{
    gretl_matrix *DG, *KN, *Veps, *NH, *NK;

    DG   = gretl_matrix_block_get_matrix(BX, 0);
    KN   = gretl_matrix_block_get_matrix(BX, 1);
    Veps = gretl_matrix_block_get_matrix(BX, 2);
    NH   = gretl_matrix_block_get_matrix(BX, 3);

    /* First chunk of Veps in Koopman's notation:
       G_t'(D_t*G_t - K_t'*N_t*H_t)
    */
    KN = gretl_matrix_reuse(KN, K->n, K->r);
    gretl_matrix_multiply(D, K->G, DG);
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
                              Nt, GRETL_MOD_NONE,
                              KN, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(KN, GRETL_MOD_NONE,
                              K->H, GRETL_MOD_NONE,
                              DG, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
                              DG, GRETL_MOD_NONE,
                              Veps, GRETL_MOD_NONE);

    /* Second chunk of Veps, to be added to the above
       H_t'(N_t*H_t - N_t*K_t*G_t)
    */
    NK = gretl_matrix_reuse(KN, K->r, K->n);
    gretl_matrix_multiply(Nt, K->H, NH);
    gretl_matrix_multiply_mod(Nt, GRETL_MOD_TRANSPOSE,
                              K->Kt, GRETL_MOD_NONE,
                              NK, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(NK, GRETL_MOD_NONE,
                              K->G, GRETL_MOD_NONE,
                              NH, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
                              NH, GRETL_MOD_NONE,
                              Veps, GRETL_MOD_CUMULATE);

    if (DKstyle) {
        /* Veps = I_p - Veps */
        double vii;
        int i;

        gretl_matrix_multiply_by_scalar(Veps, -1.0);
        for (i=0; i<K->p; i++) {
            vii = gretl_matrix_get(Veps, i, i);
            gretl_matrix_set(Veps, i, i, 1.0 + vii);
        }
    }

    /* Veps (p x p) holds the variance of \epsilon_t
       conditional on Y_n: now form the per-equation
       disturbance variance matrices, @Vwt and @Vut,
       for this time-step.
    */
    gretl_matrix_qform(K->H, GRETL_MOD_NONE, Veps,
                       Vwt, GRETL_MOD_NONE);
    gretl_matrix_qform(K->G, GRETL_MOD_NONE, Veps,
                       Vut, GRETL_MOD_NONE);

    return 0;
}

static int dist_variance (kalman *K,
			  gretl_matrix *D,
			  gretl_matrix *Vwt,
			  gretl_matrix *Vut,
			  gretl_matrix *Nt,
			  gretl_matrix_block *BX,
			  int t, int DKstyle)
{
    int err = 0;

    if (D != NULL) {
	/* needed only in presence of obs disturbance */
	if (K->exact && t < K->d) {
	    /* D_t = K_t' N_t K_t */
	    gretl_matrix_qform(K->Kt, GRETL_MOD_TRANSPOSE,
			       Nt, D, GRETL_MOD_NONE);
	} else {
	    /* D_t = F_t^{-1} + K_t' N_t K_t */
	    fast_copy_values(D, K->iFt);
	    if (t < K->N - 1) {
		gretl_matrix_qform(K->Kt, GRETL_MOD_TRANSPOSE,
				   Nt, D, GRETL_MOD_CUMULATE);
	    }
	}
    }

    if (K->p == 0) {
	/* variance of state disturbance */
	if (DKstyle) {
	    /* VS - VS N_t VS */
	    fast_copy_values(Vwt, K->VS);
	    gretl_matrix_qform(K->VS, GRETL_MOD_TRANSPOSE,
			       Nt, Vwt, GRETL_MOD_DECREMENT);
	} else {
	    /* VS N_t VS */
	    gretl_matrix_qform(K->VS, GRETL_MOD_TRANSPOSE,
			       Nt, Vwt, GRETL_MOD_NONE);
	}
	load_to_diag(K->Vsd, Vwt, t, 0);

        /* variance of obs disturbance */
        if (DKstyle) {
            /* VY - VY D_t VY */
            fast_copy_values(Vut, K->VY);
            gretl_matrix_qform(K->VY, GRETL_MOD_TRANSPOSE,
                               D, Vut, GRETL_MOD_DECREMENT);
        } else {
            /* VY D_t VY */
            gretl_matrix_qform(K->VY, GRETL_MOD_TRANSPOSE,
                               D, Vut, GRETL_MOD_NONE);
        }
        load_to_diag(K->Vsd, Vut, t, K->r);
    } else {
        /* cross-correlated disturbance variance */
        err = combined_dist_variance(K, D, Nt, Vwt, Vut, BX,
                                     DKstyle);
        if (!err) {
            load_to_diag(K->Vsd, Vwt, t, 0);
            load_to_diag(K->Vsd, Vut, t, K->r);
        }
    }

    return err;
}

/* Initial smoothed state: a + P*r0. Note that this does
   not apply in the case of exact initial smoothing.
*/

static void koopman_calc_a0 (kalman *K, gretl_matrix *r0)
{
    if (K->Pini != NULL) {
        gretl_matrix_multiply(K->Pini, r0, K->a0);
    } else {
        construct_Pini(K);
        gretl_matrix_multiply(K->P0, r0, K->a0);
    }
    if (K->aini != NULL) {
        gretl_matrix_add_to(K->a0, K->aini);
    }
    load_to_row(K->A, K->a0, 0);
}

static void transcribe_r0_N0 (kalman *K,
			      gretl_matrix *r0,
			      gretl_matrix *rdag,
			      gretl_matrix *N0,
			      gretl_matrix *Ndag)
{
    double nij;
    int i, j;

    for (j=0; j<K->r; j++) {
	r0->val[j] = rdag->val[j];
	for (i=0; i<K->r; i++) {
	    nij = gretl_matrix_get(Ndag, i, j);
	    gretl_matrix_set(N0, i, j, nij);
	}
    }
}

/* Implement exact initial state smoothing for t <= d:
   if @dsi is non-zero we're doing disturbance smoothing,
   otherwise it's state smoothing. Note: this approach
   does not work for a multivariate observable.
*/

static int exact_initial_smooth (kalman *K,
                                 gretl_matrix *r0,
                                 gretl_matrix *L0,
                                 gretl_matrix *N0,
                                 gretl_matrix *N1,
				 dsinfo *dsi)
{
    gretl_matrix_block *B;
    gretl_matrix *Mk;
    gretl_matrix *F1, *F2;
    gretl_matrix *K0, *K1;
    gretl_matrix *L1, *L2;
    gretl_matrix *Pdag, *Ldag;
    gretl_matrix *rdag, *rdag_;
    gretl_matrix *Ndag, *Ndag_;
    gretl_matrix *tmp = K->PZ;
    gretl_matrix *rbit;
    int dist = (dsi != NULL);
    int offset;
    int rr = 2 * K->r;
    int nt = K->n;
    int i, t, err = 0;

    B = gretl_matrix_block_new(&Mk, K->r, K->n,
                               &F1, K->n, K->n,
                               &F2, K->n, K->n,
                               &K0, K->r, K->n,
                               &K1, K->r, K->n,
                               &L1, K->r, K->r,
                               &L2, K->r, K->r,
                               &Pdag, K->r, rr,
                               &Ldag, rr, rr,
                               &rdag, rr, 1,
                               &rdag_, rr, 1,
                               &rbit, K->r, 1,
                               &Ndag, rr, rr,
                               &Ndag_, rr, rr,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    /* initial r† = [rd 0]' */
    for (i=0; i<rr; i++) {
        rdag->val[i] = (i < K->r)? r0->val[i] : 0;
    }

    gretl_matrix_zero(Ndag);
    gretl_matrix_inscribe_matrix(Ndag, N0, 0, 0, GRETL_MOD_NONE);
    gretl_matrix_zero(Ndag_);

#if EXACT_DEBUG
    fprintf(stderr, "*** exact_initial_smooth ***\n");
    gretl_matrix_print(r0, "rd");
    gretl_matrix_print(N0, "Nd");
#endif

    /* Two cases must be handled below:

       Case 1: F∞ is nonsingular (the only case discussed in
               Durbin and Koopman, 2012)
       Case 2: F∞ is singular (see Koopman and Durbin, 2003)

       D and K describe the second case as one in which F∞ = 0,
       but it may be that F∞ is singular without being a zero
       matrix, for example if the observable is of higher
       dimension than the state.
    */

    for (t=K->d-1; t>=0; t--) {
	int singular = 0;

        err = load_filter_data(K, t, &nt, SM_STATE_INI);
        if (err) {
            break;
        } else if (nt < K->n) {
            gretl_matrix_reuse(Mk, K->r, nt);
            gretl_matrix_reuse(F1, nt, nt);
            gretl_matrix_reuse(F2, nt, nt);
            gretl_matrix_reuse(K0, K->r, nt);
            gretl_matrix_reuse(K1, K->r, nt);
            gretl_matrix_reuse(tmp, K->r, nt);
        }

	if (dist) {
	    /* correct? */
	    load_to_row(dsi->R, r0, t);
	}

        /* F∞ = Z * P∞ * Z' */
        gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->Pk0,
                           F1, GRETL_MOD_NONE);
	/* F1 = inv(F∞), if possible */
	if (gretl_invert_symmetric_matrix(F1) != 0) {
	    singular = 1;
	}

        /* F★ = Z * P★ * Z' [+ VY] */
        if (K->VY != NULL) {
            fast_copy_values(K->Ft, K->VY);
            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
                               K->Ft, GRETL_MOD_CUMULATE);
        } else {
            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
                               K->Ft, GRETL_MOD_NONE);
        }

	if (singular) {
	    /* we need inv(F★) under the name F1 */
	    fast_copy_values(F1, K->Ft);
	    if (gretl_invert_symmetric_matrix(F1)) {
		fprintf(stderr, "F_star is singular: can't continue\n");
		err = E_NOTPD;
		break;
	    }
	}

	/* M★ = P★ * Z' */
	gretl_matrix_multiply(K->P0, K->ZT, K->Mt);

	if (singular) {
	    /* K0 = T * M★ * F1 */
	    gretl_matrix_multiply(K->T, K->Mt, tmp);
	    gretl_matrix_multiply(tmp, F1, K0);
	} else {
	    /* M∞ = P∞ * Z' */
	    gretl_matrix_multiply(K->Pk0, K->ZT, Mk);

	    /* F2 = -iF∞ * F★ * iF∞ */
	    gretl_matrix_zero(F2);
	    gretl_matrix_qform(F1, GRETL_MOD_TRANSPOSE, K->Ft,
			       F2, GRETL_MOD_DECREMENT);

	    /* K0 = T * M∞ * F1 */
	    gretl_matrix_multiply(K->T, Mk, tmp);
	    gretl_matrix_multiply(tmp, F1, K0);
	}

	if (dist) {
	    /* state disturbance: VS r0_t */
	    gretl_matrix_multiply(K->VS, r0, dsi->r1);
	    load_to_row_offset(K->U, dsi->r1, t, 0);
	    if (K->VY != NULL) {
		/* obs disturbance: -VY K0_t' r0_t */
		gretl_matrix_multiply_mod(K0, GRETL_MOD_TRANSPOSE,
					  r0, GRETL_MOD_NONE,
					  dsi->n1, GRETL_MOD_NONE);
		gretl_matrix_multiply(K->VY, dsi->n1, dsi->u);
		gretl_matrix_multiply_by_scalar(dsi->u, -1.0);
		load_to_row_offset(K->U, dsi->u, t, K->r);
	    }
	    load_to_row(K->Kt, K0, t);
	    dist_variance(K, dsi->D, dsi->Vwt, dsi->Vut, N0,
			  NULL, t, dsi->DKstyle);
	}

	/* L0 = T - K0 * Z */
	fast_copy_values(L0, K->T);
	gretl_matrix_multiply_mod(K0, GRETL_MOD_NONE,
				  K->ZT, GRETL_MOD_TRANSPOSE,
				  L0, GRETL_MOD_DECREMENT);

	if (singular) {
	    /* the relevant big-L matrix is fairly simple */
	    gretl_matrix_zero(Ldag);
	    gretl_matrix_inscribe_matrix(Ldag, L0, 0, 0, GRETL_MOD_NONE);
	    gretl_matrix_inscribe_matrix(Ldag, K->T, K->r, K->r, GRETL_MOD_NONE);
	} else {
	    /* K1 = T * M★ * F1 + T * M∞ * F2 */
	    gretl_matrix_multiply(K->T, K->Mt, tmp);
	    gretl_matrix_multiply(tmp, F1, K1);
	    gretl_matrix_multiply(K->T, Mk, tmp);
	    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
				      F2, GRETL_MOD_NONE,
				      K1, GRETL_MOD_CUMULATE);
	    /* L1 = -K1 * Z */
	    gretl_matrix_zero(L1);
	    gretl_matrix_multiply_mod(K1, GRETL_MOD_NONE,
				      K->ZT, GRETL_MOD_TRANSPOSE,
				      L1, GRETL_MOD_DECREMENT);
	    /* L† = (L0 ~ L1) | (0 ~ L0) */
	    gretl_matrix_zero(Ldag);
	    gretl_matrix_inscribe_matrix(Ldag, L0, 0, 0, GRETL_MOD_NONE);
	    gretl_matrix_inscribe_matrix(Ldag, L1, 0, K->r, GRETL_MOD_NONE);
	    gretl_matrix_inscribe_matrix(Ldag, L0, K->r, K->r, GRETL_MOD_NONE);
	}
#if EXACT_DEBUG
	gretl_matrix_print(Ldag, "Ldag");
#endif

	/* P† = P★ ~ P∞ */
	gretl_matrix_inscribe_matrix(Pdag, K->P0, 0, 0,
				     GRETL_MOD_NONE);
	gretl_matrix_inscribe_matrix(Pdag, K->Pk0, 0, K->r,
				     GRETL_MOD_NONE);

        /* r†_{t-1} = rPart1 + L†_t * r†_t */
        gretl_matrix_zero(rdag_);
        gretl_matrix_multiply(K->ZT, F1, tmp);
        gretl_matrix_multiply(tmp, K->vt, rbit);
	offset = singular ? 0 : K->r;
	gretl_matrix_inscribe_matrix(rdag_, rbit, offset, 0,
				     GRETL_MOD_NONE);
        gretl_matrix_multiply_mod(Ldag, GRETL_MOD_TRANSPOSE,
                                  rdag, GRETL_MOD_NONE,
                                  rdag_, GRETL_MOD_CUMULATE);
#if EXACT_DEBUG
	gretl_matrix_print(Pdag, "Pdag");
        gretl_matrix_print(rdag, "rdag");
        gretl_matrix_print(rdag_, "rdag_minus");
#endif

	/* \hat{\alpha}_t = a_t + P†_t * r†_{t-1} */
	fast_copy_values(K->a1, K->a0);
	gretl_matrix_multiply_mod(Pdag, GRETL_MOD_NONE,
				  rdag_, GRETL_MOD_NONE,
				  K->a1, GRETL_MOD_CUMULATE);
	load_to_row(K->A, K->a1, t);
#if EXACT_DEBUG
        gretl_matrix_print(K->a1, "ahat");
#endif

        fast_copy_values(rdag, rdag_);

	/* using P1 as workspace */
        gretl_matrix_qform(K->ZT, GRETL_MOD_NONE, F1,
                           K->P1, GRETL_MOD_NONE);

	if (singular) {
	    /* N†_{t-1} = (Z'*F1*Z ~ 0) | 0 + L†'*N†*L† */
	    gretl_matrix_inscribe_matrix(Ndag_, K->P1, 0, 0,
					 GRETL_MOD_NONE);
	} else {
	    /* N†_{t-1} = (0 ~ Z'*F1*Z) | (Z'*F1*Z ~ Z'*F2*Z) + L†'*N†*L† */
	    gretl_matrix_inscribe_matrix(Ndag_, K->P1, 0, K->r,
					 GRETL_MOD_NONE);
	    gretl_matrix_inscribe_matrix(Ndag_, K->P1, K->r, 0,
					 GRETL_MOD_NONE);
	    gretl_matrix_qform(K->ZT, GRETL_MOD_NONE, F2,
			       K->P1, GRETL_MOD_NONE);
	    gretl_matrix_inscribe_matrix(Ndag_, K->P1, K->r, K->r,
					 GRETL_MOD_NONE);
	}
        gretl_matrix_qform(Ldag, GRETL_MOD_TRANSPOSE, Ndag,
                           Ndag_, GRETL_MOD_CUMULATE);

	if (!dist) {
	    /* Vt = P★ - P† * N†_{t-1} * P†' */
	    fast_copy_values(K->P1, K->P0);
	    gretl_matrix_qform(Pdag, GRETL_MOD_NONE, Ndag_,
			       K->P1, GRETL_MOD_DECREMENT);
	    load_to_vech(K->P, K->P1, K->r, t);
#if EXACT_DEBUG
	    gretl_matrix_print(Ndag, "Ndag");
	    gretl_matrix_print(Ndag_, "Ndag_minus");
	    gretl_matrix_print(K->P1, "Vt");
#endif
	}

        fast_copy_values(Ndag, Ndag_);

	if (dist && t > 0) {
	    /* transcribe for next step */
	    transcribe_r0_N0(K, r0, rdag, N0, Ndag);
	}

        if (nt < K->n) {
            unshrink_vt(K, SM_STATE_INI);
            gretl_matrix_reuse(Mk, K->r, K->n);
            gretl_matrix_reuse(F1, K->n, K->n);
            gretl_matrix_reuse(F2, K->n, K->n);
            gretl_matrix_reuse(K0, K->r, K->n);
            gretl_matrix_reuse(K1, K->r, K->n);
            gretl_matrix_reuse(tmp, K->r, K->n);
        }
    }

    gretl_matrix_block_destroy(B);

    return err;
}

/* This iteration is in common between the state smoother
   (Anderson-Moore) and the disturbance smoother (Koopman).
*/

static void LrN_iteration (kalman *K,
			   gretl_matrix *L,
			   gretl_matrix *n1,
			   gretl_matrix *r0,
			   gretl_matrix *r1,
			   gretl_matrix *N0,
			   gretl_matrix *N1,
			   int t)
{
    if (t < K->N - 1) {
	/* L_t = T_t - K_t Z_t */
	fast_copy_values(L, K->T);
	gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				  K->ZT, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_DECREMENT);
    }

    /* r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t */
    gretl_matrix_multiply(K->iFt, K->vt, n1);
    if (t == K->N - 1) {
	gretl_matrix_multiply(K->ZT, n1, r0);
    } else {
	gretl_matrix_multiply(K->ZT, n1, r1);
	gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE,
				  r0, GRETL_MOD_NONE,
				  r1, GRETL_MOD_CUMULATE);
	fast_copy_values(r0, r1);
    }

    /* N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t */
    if (t == K->N - 1) {
	gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
			   K->iFt, N0, GRETL_MOD_NONE);
    } else {
	gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
			   K->iFt, N1, GRETL_MOD_NONE);
	gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
			   N0, N1, GRETL_MOD_CUMULATE);
	fast_copy_values(N0, N1);
    }
}

/* Disturbance smoothing -- see Koopman, Shephard and Doornik (SsfPack
   doc), section 4.4; also Durbin and Koopman, 2012.

   As of 2022-02-13 we're experimenting with the exact diffuse initial
   variant of this. But proper testing is needed before it can go
   public.
*/

static int koopman_smooth (kalman *K, int DKstyle)
{
    gretl_matrix_block *B, *BX = NULL;
    gretl_matrix *u, *L, *R;
    gretl_matrix *r0, *r1, *N0, *N1, *n1, *tr;
    gretl_matrix *D = NULL;
    gretl_matrix *Vwt = NULL;
    gretl_matrix *Vut = NULL;
    gretl_matrix *DG = NULL;
    gretl_matrix *KN = NULL;
    gretl_matrix *RZS = NULL;
    gretl_matrix *NH = NULL;
    gretl_matrix *Ut = NULL;
    dsinfo dsi = {0};
    int ft_min = 0;
    int t, err = 0;

    B = gretl_matrix_block_new(&u,  K->n, 1,
                               &L,  K->r, K->r,
                               &R,  K->N, K->r,
                               &r0, K->r, 1,
                               &r1, K->r, 1,
                               &N0, K->r, K->r,
                               &N1, K->r, K->r,
                               &n1, K->n, 1,
                               &tr, K->r, 1,
                               NULL);

    if (B == NULL) {
        return E_ALLOC;
    }

    /* for variance of smoothed disturbances */
    err = maybe_resize_dist_mse(K, &Vwt, &Vut);

    if (K->VY != NULL) {
	/* for variance of observable */
	D = gretl_matrix_alloc(K->n, K->n);
    }

    if (K->p > 0) {
	/* cross-correlated disturbances */
        BX = gretl_matrix_block_new(&DG,  K->n, K->p,
                                    &KN,  K->n, K->r,
                                    &RZS, K->p, K->p,
                                    &NH,  K->r, K->p,
                                    &Ut,  K->p, 1,
                                    NULL);
        if (BX == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        gretl_matrix_block_destroy(B);
        gretl_matrix_block_destroy(BX);
        gretl_matrix_free(Vwt);
        gretl_matrix_free(Vut);
	gretl_matrix_free(D);
        return err;
    }

    if (K->exact) {
	/* wrap up matrices, etc that we'll need */
	dsi.R = R;
	dsi.D = D;
	dsi.u = u;
	dsi.n1 = n1;
	dsi.r1 = r1;
	dsi.Vwt = Vwt;
	dsi.Vut = Vut;
	dsi.DKstyle = DKstyle;
	ft_min = K->d;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    /* The backward recursion */

    for (t=K->N-1; t>=0 && !err; t--) {
        if (K->exact && t == K->d - 1) {
            err = exact_initial_smooth(K, r0, L, N0, N1, &dsi);
            break;
        }

        err = load_filter_data(K, t, NULL, SM_DIST_BKWD);
        if (err) {
            break;
        }

        /* u_t = F_t^{-1} v_t - K_t' r_t */
        gretl_matrix_multiply(K->iFt, K->vt, u);
        if (t < K->N - 1) {
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
                                      r0, GRETL_MOD_NONE,
                                      u, GRETL_MOD_DECREMENT);
        }
        /* Store u_t values in K->V: these are needed in
           the forward pass to compute the smoothed
           disturbances. Also save r_t in R.
        */
        load_to_row(K->V, u, t);
	load_to_row(R, r0, t);

	/* compute variance of disturbances */
	err = dist_variance(K, D, Vwt, Vut, N0, BX, t, DKstyle);
	if (err) {
	    break;
	}

	/* compute r_{t-1}, N_{t-1} */
	LrN_iteration(K, L, n1, r0, r1, N0, N1, t);

	if (t == 0) {
	    /* compute initial smoothed state */
	    koopman_calc_a0(K, r0);
        }
    }

    /* Forward iteration for smoothed disturbances, all time steps,
       plus smoothed state from t = 1 (or t = d in the exact initial
       case) onward.
    */
    for (t=ft_min; t<K->N; t++) {
        err = load_filter_data(K, t, NULL, SM_DIST_FRWD);
        if (err) {
            break;
        }

	/* state disturbance */
        load_from_row(r0, R, t, GRETL_MOD_NONE);
        if (K->p > 0) {
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
                                      r0, GRETL_MOD_NONE,
                                      Ut, GRETL_MOD_NONE);
            load_from_row(K->vt, K->V, t, GRETL_MOD_NONE);
            gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
                                      K->vt, GRETL_MOD_NONE,
                                      Ut, GRETL_MOD_CUMULATE);
            gretl_matrix_multiply(K->H, Ut, r1);
        } else {
            gretl_matrix_multiply(K->VS, r0, r1);
        }
        load_to_row(R, r1, t);
	load_to_row_offset(K->U, r1, t, 0);

	/* observation disturbance */
	if (K->VY != NULL) {
	    if (K->p > 0) {
		gretl_matrix_multiply(K->G, Ut, n1);
	    } else {
		gretl_matrix_multiply(K->VY, K->vt, n1);
	    }
	    load_to_row_offset(K->U, n1, t, K->r);
	}

	if (t >= K->d) {
	    /* state: a_{t+1} = T a_t + w_t (or + H*eps_t) */
	    load_from_row(K->a0, K->A, t-1, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->T, K->a0, K->a1);
	    if (K->exact && t == K->d) {
		/* pick up prior smoothed state disturbance */
		load_from_row(K->a1, K->U, t-1, GRETL_MOD_CUMULATE);
	    } else {
		load_from_row(K->a1, R, t-1, GRETL_MOD_CUMULATE);
	    }
	    if (K->mu != NULL) {
		gretl_matrix_add_to(K->a1, K->mu);
	    }
	    load_to_row(K->A, K->a1, t);
	}
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(BX);
    gretl_matrix_free(Vwt);
    gretl_matrix_free(Vut);
    gretl_matrix_free(D);

    return err;
}

/* Anderson-Moore Kalman smoothing.  This method uses a_{t|t-1} and
   P_{t|t-1} for all t, but we overwrite these with the smoothed
   values as we go. We also need stored values for the prediction
   error, its MSE, and the gain at each time step.  Note that r_t and
   N_t are set to zero for t = N - 1.
*/

static int anderson_moore_smooth (kalman *K)
{
    gretl_matrix_block *B;
    gretl_matrix *r0, *r1, *N0, *N1, *n1, *L;
    int nt = K->n;
    int t, err = 0;

    B = gretl_matrix_block_new(&r0,  K->r, 1,
                               &r1,  K->r, 1,
                               &N0,  K->r, K->r,
                               &N1,  K->r, K->r,
                               &n1,  K->n, 1,
                               &L,   K->r, K->r,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    for (t=K->N-1; t>=0 && !err; t--) {
        if (K->exact && t == K->d - 1) {
            err = exact_initial_smooth(K, r0, L, N0, N1, NULL);
            break;
        }

        err = load_filter_data(K, t, &nt, SM_STATE_STD);
        if (err) {
            break;
        } else if (nt < K->n) {
            gretl_matrix_reuse(n1, nt, 1);
        }

	/* compute r_{t-1}, N_{t-1} */
	LrN_iteration(K, L, n1, r0, r1, N0, N1, t);

        /* a_{t|T} = a_{t|t-1} + P_{t|t-1} r_{t-1} */
        fast_copy_values(K->a1, K->a0);
        gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
                                  r0, GRETL_MOD_NONE,
                                  K->a1, GRETL_MOD_CUMULATE);
        load_to_row(K->A, K->a1, t);

        /* P_{t|T} = P_{t|t-1} - P_{t|t-1} N_{t-1} P_{t|t-1} */
        fast_copy_values(K->P1, K->P0);
        gretl_matrix_qform(K->P0, GRETL_MOD_NONE, N0,
                           K->P1, GRETL_MOD_DECREMENT);
        load_to_vech(K->P, K->P1, K->r, t);

        if (nt < K->n) {
            unshrink_vt(K, SM_STATE_STD);
            gretl_matrix_reuse(n1, K->n, 1);
        }
    }

    gretl_matrix_block_destroy(B);

    return err;
}

/* If we're doing smoothing for a system that has time-varying
   coefficients in K->T or K->ZT we'll record the vec of the
   coefficient matrices for each time-step on the forward pass.
   Here we allocate the required storage.
*/

static int kalman_add_stepinfo (kalman *K)
{
    int err = 0;

    K->step = malloc(sizeof *K->step);

    if (K->step == NULL) {
        return E_ALLOC;
    }

    K->step->T = K->step->ZT = NULL;

    if (matrix_is_varying(K, K_T)) {
        K->step->T = gretl_matrix_alloc(K->N, K->r * K->r);
        if (K->step->T == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && matrix_is_varying(K, K_ZT)) {
        K->step->ZT = gretl_matrix_alloc(K->N, K->r * K->n);
        if (K->step->ZT == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        free_stepinfo(K);
    }

    return err;
}

/* optional matrix to hold smoothed disturbances a la Koopman */

static int ensure_U_matrix (kalman *K)
{
    int Ucols = K->r;
    int Urows = K->N;
    int err = 0;

    if (K->VY != NULL) {
        Ucols += K->n;
    }

    if (K->U == NULL) {
        K->U = gretl_matrix_alloc(Urows, Ucols);
        if (K->U == NULL) {
            err = E_ALLOC;
        }
    } else if (K->U->rows != Urows || K->U->cols != Ucols) {
        err = gretl_matrix_realloc(K->U, Urows, Ucols);
    }

    return err;
}

static int real_kalman_smooth (kalman *K, int dist, PRN *prn)
{
    int err = kalman_ensure_output_matrices(K);
    int save_exact = K->exact;
    gretl_matrix *save_P = NULL;

    if (!err && dist && !kalman_univariate(K)) {
        err = ensure_U_matrix(K);
    }

    if (err) {
        return err;
    }

    if (matrix_is_varying(K, K_T) || matrix_is_varying(K, K_ZT)) {
        /* add recorder for T_t and/or Z_t */
        err = kalman_add_stepinfo(K);
        if (err) {
            goto bailout;
        }
    }

#if EXACT_SMDIST
    if (dist > 0 && K->exact && K->p > 0) {
	/* we don't know how to do exact initial disturbance
	   smoothing in the cross-correlated case
	*/
	pputs(prn, "Warning: exact initial disturbance smoothing is not "
	      "supported for cross-correlated errors\n\n");
	kalman_set_diffuse(K, 1);
    }
#else
    if (dist > 0 && K->exact && !kalman_univariate(K)) {
	/* exact initial disturbance smoother not ready? */
	pputs(prn, "Warning: exact initial disturbance smoothing is not "
	      "supported\n\n");
	kalman_set_diffuse(K, 1);
    }
#endif

    if (K->exact) {
        int rr = (K->r * K->r + K->r) / 2;

        /* Note: PK will need more than K->r cols if d > K->r */
        K->PK = gretl_zero_matrix_new(rr, K->r);
    }

    if (!err) {
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    /* recorder for partially missing observations */
    K->nt = NULL;

    if (!err) {
        /* prior forward pass */
        K->flags |= KALMAN_SMOOTH;
        if (dist == 0) {
            /* Anderson-Moore */
            K->flags |= KALMAN_SM_AM;
        } else if (0 /* !K->exact */) {
	    /* don't overwrite P (not needed for the disturbance
	       smoother, other than in the exact initial case)
	    */
            save_P = K->P;
            K->P = NULL;
        }
	if (kalman_univariate(K)) {
	    err = kfilter_univariate(K, NULL);
	} else {
	    err = kalman_forecast(K, NULL);
	}
        K->flags &= ~KALMAN_SMOOTH;
        K->flags &= ~KALMAN_SM_AM;
    }

    K->t = 0;

    if (!err) {
	if (kalman_univariate(K)) {
	    err = ksmooth_univariate(K, dist);
	} else if (dist > 1) {
            err = koopman_smooth(K, 1);
        } else if (dist == 1) {
            err = koopman_smooth(K, 0);
        } else {
            err = anderson_moore_smooth(K);
        }
    }

    if (save_exact) {
        /* re-establish 'exact' status */
        kalman_set_diffuse(K, 2);
    }
    if (save_P != NULL) {
        K->P = save_P;
    }

    /* free "special case" storage */
    if (K->nt != NULL) {
        free(K->nt);
        K->nt = NULL;
    }
    if (K->PK != NULL) {
        gretl_matrix_free(K->PK);
        K->PK = NULL;
    }

 bailout:

    /* trash the "stepinfo" storage */
    free_stepinfo(K);

    return err;
}

/* For use with userland bundle-based API */

int kalman_bundle_smooth (gretl_bundle *b, int dist, PRN *prn)
{
    kalman *K = gretl_bundle_get_private_data(b);

    if (K == NULL) {
        fprintf(stderr, "kalman_bundle_smooth: K is NULL\n");
        return E_DATA;
    }

    K->b = b; /* attach bundle pointer */

    return real_kalman_smooth(K, dist, prn);
}

/* For use with "plain C" API. TODO: support disturbance
   smoothing via @opt -- right now only state smoothing
   is offered.
*/

gretl_matrix *kalman_smooth (kalman *K, gretlopt opt,
                             PRN *prn, int *err)
{
    gretl_matrix *S = NULL;

    *err = real_kalman_smooth(K, 0, prn);

    if (K->A != NULL && K->P != NULL) {
        int r = K->A->rows;
        int c = K->A->cols;
        size_t Asize = r * c * sizeof(double);
        size_t Psize = 0;

        if (opt & OPT_M) {
            /* include MSE of state */
            c += K->P->cols;
            Psize = r * K->P->cols * sizeof(double);
        }
        if (r > 0 && c > 0) {
            S = gretl_matrix_alloc(r, c);
            if (S == NULL) {
                *err = E_ALLOC;
            }
        } else {
            *err = E_DATA;
        }
        if (S != NULL) {
            memcpy(S->val, K->A->val, Asize);
            if (opt & OPT_M) {
                memcpy(S->val + r * K->A->cols, K->P->val, Psize);
            }
        }
    }

    return S;
}

static gretl_matrix *extract_Q (kalman *K,
                                const gretl_matrix *Sim0)
{
    gretl_matrix *Q;
    double x;
    int i, j;

    Q = gretl_matrix_alloc(K->r, K->r);

    if (Q != NULL) {
        for (i=0; i<K->r; i++) {
            for (j=0; j<K->r; j++) {
                x = gretl_matrix_get(Sim0, i, j);
                gretl_matrix_set(Q, i, j, x);
            }
        }
    }

    return Q;
}

/* See the account in Koopman, Shephard and Doornik, Econometrics
   Journal, 1999 (volume 2, pp. 113-166), section 4.2, regarding the
   initialization of the state under simulation.
*/

static int sim_state_0 (kalman *K, const gretl_matrix *U,
                        const gretl_matrix *Sim0)
{
    gretl_matrix *Q, *v0 = NULL, *bv = NULL;
    int getroot = 1;
    int i, err = 0;

    if (!kalman_ssfsim(K)) {
        if (Sim0 != NULL) {
            /* Sim0 contains the state for t = 1 */
            err = gretl_matrix_copy_values(K->a0, Sim0);
        }
        /* error or not, we're done */
        return err;
    }

    /* now we're in the "ssfsim" case, emulating ssfpack */

    if (Sim0 != NULL) {
        /* Sim0 contains state variance factor
           plus the state for t = 0
        */
        Q = extract_Q(K, Sim0);
        getroot = 0;
    } else {
        Q = gretl_matrix_copy(K->P0);
    }

    if (Q == NULL) {
        err = E_ALLOC;
    } else if (getroot) {
        err = gretl_matrix_psd_root(Q, 0);
    }

    if (!err) {
        int vlen = K->p > 0 ? K->p : K->r;

        v0 = gretl_matrix_alloc(vlen, 1);
        if (v0 == NULL) {
            err = E_ALLOC;
        }
    }

    if (K->p > 0) {
        bv = gretl_matrix_alloc(K->r, 1);
        if (bv == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && Sim0 != NULL) {
        /* set a0 from last row of Sim0 */
        for (i=0; i<K->r; i++) {
            K->a0->val[i] = gretl_matrix_get(Sim0, K->r, i);
        }
    }

    if (!err) {
        /* handle the t = 0 disturbance */
        load_from_row(v0, U, 0, GRETL_MOD_NONE);
        if (K->p > 0) {
            /* cross-correlated */
            gretl_matrix_multiply(K->H, v0, bv);
            gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                                      bv, GRETL_MOD_NONE,
                                      K->a0, GRETL_MOD_CUMULATE);
        } else {
            gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                                      v0, GRETL_MOD_NONE,
                                      K->a0, GRETL_MOD_CUMULATE);
        }
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(v0);
    gretl_matrix_free(bv);

    return err;
}

/* note: it's OK for @S to be NULL (if the simulated state is not
   wanted), so watch out for that!
*/

static int kalman_simulate (kalman *K,
                            const gretl_matrix *U,
                            const gretl_matrix *Sim0,
                            gretl_matrix *Y,
                            gretl_matrix *S,
                            PRN *prn)
{
    gretl_matrix *yt, *et = NULL;
    int obs_offset = 0;
    int obsdist = 0;
    int tmin = 0;
    int err = 0;

    yt = gretl_zero_matrix_new(K->n, 1);
    if (yt == NULL) {
        return E_ALLOC;
    }

    if (K->p > 0) {
        et = gretl_matrix_alloc(K->p, 1);
        if (et == NULL) {
            gretl_matrix_free(yt);
            return E_ALLOC;
        }
    }

    if (Y->cols == K->r + K->n) {
        /* combined (state, obs) in @Y */
        S = Y;
        obs_offset = K->r;
    }

    err = sim_state_0(K, U, Sim0);

    if (!err && kalman_ssfsim(K)) {
        if (S != NULL) {
            load_to_row_offset(S, K->a0, 0, 0);
        }
        load_to_row_offset(Y, yt, 0, obs_offset);
        /* the first row of output is handled */
        tmin = 1;
    }

    if (K->p == 0 && K->VY != NULL) {
        /* we want to read observation disturbances */
        obsdist = 1;
    }

    for (K->t = tmin; K->t < K->N && !err; K->t += 1) {
        if (filter_is_varying(K)) {
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                break;
            }
        }

        /* y_t = B'*x_t + Z'*a_t + w_t */
        gretl_matrix_multiply_mod(K->ZT, GRETL_MOD_TRANSPOSE,
                                  K->a0, GRETL_MOD_NONE,
                                  yt, GRETL_MOD_NONE);
        if (K->BT != NULL) {
            /* handle missing values? */
            kalman_do_Bx(K, yt, GRETL_MOD_CUMULATE);
        }
        if (K->p > 0) {
            /* G \varepsilon_t */
            load_from_row(et, U, K->t, GRETL_MOD_NONE);
            gretl_matrix_multiply(K->G, et, K->vt);
        } else if (obsdist) {
            load_from_row_offset(K->vt, U, K->t, K->r);
        }
        gretl_matrix_add_to(yt, K->vt);

        /* record the t-dated observables */
        load_to_row_offset(Y, yt, K->t, obs_offset);

        /* record the t-dated state? */
        if (S != NULL && tmin == 0) {
            load_to_row_offset(S, K->a0, K->t, 0);
        }

        /* a_{t+1} = T*a_t + v_t */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->p > 0) {
            /* H \varepsilon_t */
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                      et, GRETL_MOD_NONE,
                                      K->a1, GRETL_MOD_CUMULATE);
        } else {
            load_from_row(K->a1, U, K->t, GRETL_MOD_CUMULATE);
        }

        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }

        /* record the (t+1)-dated state? */
        if (S != NULL && tmin == 1) {
            load_to_row_offset(S, K->a1, K->t, 0);
        }

        fast_copy_values(K->a0, K->a1);
    }

    gretl_matrix_free(yt);
    gretl_matrix_free(et);

    return err;
}

static int check_simul_inputs (kalman *K,
                               const gretl_matrix *U,
                               const gretl_matrix *Sim0,
                               const gretl_matrix *SimX,
                               int ssfsim,
                               PRN *prn)
{
    int err = 0;

    if (U == NULL) {
        err = missing_matrix_error("U");
    } else {
        int ncols;

        if (K->p > 0) {
            /* cross-correlated */
            ncols = K->p;
        } else {
            ncols = K->VY == NULL ? K->r : K->r + K->n;
        }

        if (U->cols != ncols) {
            pprintf(prn, "U should have %d columns but has %d\n",
                    ncols, U->cols);
            err = E_NONCONF;
        }
    }

    if (!err && Sim0 != NULL) {
        int r = ssfsim ? K->r + 1 : K->r;
        int c = ssfsim ? K->r : 1;

        if (Sim0->rows != r || Sim0->cols != c) {
            pprintf(prn, "simstart should be %d x %d, is %d x %d\n",
                    r, c, Sim0->rows, Sim0->cols);
        }
    }

    if (!err && K->x != NULL) {
        /* do we have enough "obsx" data? */
        const gretl_matrix *X = SimX != NULL ? SimX : K->x;

        if (X->rows < U->rows) {
            pprintf(prn, "obsx should have %d rows but has %d\n",
                    U->rows, X->rows);
            err = E_NONCONF;
        }
    }

    return err;
}

gretl_matrix *kalman_bundle_simulate (gretl_bundle *b,
                                      const gretl_matrix *U,
                                      int get_state,
                                      PRN *prn, int *err)
{
    kalman *K = gretl_bundle_get_private_data(b);
    const gretl_matrix *Sim0 = NULL;
    const gretl_matrix *SimX = NULL;
    gretl_matrix *Ret = NULL;
    gretl_matrix *savex = NULL;
    double ssfx;
    int ssfsim = 0;
    int saveN;

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    /* try accessing auxiliary info from the bundle */
    Sim0 = gretl_bundle_get_matrix(b, "simstart", NULL);
    ssfx = gretl_bundle_get_scalar(b, "ssfsim", NULL);
    if (K->x != NULL) {
        SimX = gretl_bundle_get_matrix(b, "simx", NULL);
    }

    ssfsim = !na(ssfx) && ssfx != 0;

    *err = check_simul_inputs(K, U, Sim0, SimX, ssfsim, prn);
    if (*err) {
        return NULL;
    }

    K->b = b; /* attach bundle pointer */

    saveN = K->N;
    savex = K->x;

    /* we let U temporarily define the sample length */
    K->N = U->rows;

    /* and we allow temporary replacement of K->x */
    if (SimX != NULL) {
        K->x = (gretl_matrix *) SimX;
    }

    /* set state */
    if (ssfsim) {
        K->flags |= (KALMAN_SIM | KALMAN_SSFSIM);
    } else {
        K->flags |= KALMAN_SIM;
    }

    /* now, are the other needed matrices in place? */
    *err = kalman_bundle_recheck_matrices(K, prn);

    /* matrix to hold simulated observables, and state
       if wanted */
    if (!*err) {
        int ncols = get_state ? (K->r + K->n) : K->n;

        Ret = gretl_matrix_alloc(K->N, ncols);
        if (Ret == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        *err = kalman_simulate(K, U, Sim0, Ret, NULL, prn);
    }

    if (*err) {
        gretl_matrix_free(Ret);
        Ret = NULL;
    }

    /* restore state */
    K->flags &= ~KALMAN_SIM;
    K->flags &= ~KALMAN_SSFSIM;
    K->N = saveN;
    K->x = savex;

    return Ret;
}

static int matrix_is_diagonal (const gretl_matrix *m)
{
    int i, j;

    for (j=0; j<m->cols; j++) {
        for (i=0; i<m->rows; i++) {
            if (i != j && gretl_matrix_get(m, i, j) != 0.0) {
		return 0;
	    }
        }
    }

    return 1;
}

static int simdata_refresh_QR (kalman *K, PRN *prn)
{
    int err = 0;

    if (matrix_is_varying(K, K_VS) || matrix_is_varying(K, K_VY)) {
        err = kalman_update_matrices(K, prn);
    }

    return err;
}

/* Return a matrix in which the standard normal variates
   in @U are scaled according to K->VS, and K->VY if present.
*/

gretl_matrix *kalman_bundle_simdata (gretl_bundle *b,
                                     const gretl_matrix *U,
                                     PRN *prn, int *err)
{
    kalman *K = gretl_bundle_get_private_data(b);
    gretl_matrix *E = NULL;

    if (K == NULL || U == NULL) {
        *err = E_DATA;
        return NULL;
    }

    if (K->p > 0) {
        if (U->cols != K->p) {
            *err = E_DATA;
            return NULL;
        } else {
            /* nothing to be done, really */
            E = gretl_matrix_copy(U);
        }
    } else {
        int n = K->VY == NULL ? 0 : K->n;
        int t, j, rn = K->r + n;
        int T = U->rows;
        int varying = 0;
        double vjj, utj;

        if (U->cols != rn) {
            *err = E_DATA;
            return NULL;
        }

        E = gretl_matrix_alloc(U->rows, rn);
        if (E == NULL) {
            *err = E_ALLOC;
            return NULL;
        }

        if (matrix_is_varying(K, K_VS) || matrix_is_varying(K, K_VY)) {
            varying = 1;
        }

        K->b = b;
        set_kalman_running(K);

        if (matrix_is_diagonal(K->VS) &&
            (K->VY == NULL || matrix_is_diagonal(K->VY))) {
            for (t=0; t<T && !*err; t++) {
                if (varying) {
                    K->t = t;
                    *err = simdata_refresh_QR(K, prn);
                }
                for (j=0; j<rn && !*err; j++) {
                    if (j < K->r) {
                        vjj = gretl_matrix_get(K->VS, j, j);
                    } else {
                        vjj = gretl_matrix_get(K->VY, j-K->r, j-K->r);
                    }
                    utj = gretl_matrix_get(U, t, j);
                    gretl_matrix_set(E, t, j, sqrt(vjj) * utj);
                }
            }
        } else {
            gretl_matrix *V = gretl_zero_matrix_new(rn, rn);
            gretl_matrix *Ut = NULL;
            gretl_matrix *Et = NULL;

            if (V == NULL) {
                *err = E_ALLOC;
                goto bailout;
            }

            if (varying) {
                Ut = gretl_matrix_alloc(1, rn);
                Et = gretl_matrix_alloc(1, rn);
                if (Ut == NULL || Et == NULL) {
                    gretl_matrix_free(V);
                    *err = E_ALLOC;
                    goto bailout;
                }
            }

            if (varying) {
                for (t=0; t<T && !*err; t++) {
                    K->t = t;
                    *err = simdata_refresh_QR(K, prn);
                    if (!*err) {
                        gretl_matrix_inscribe_matrix(V, K->VS, 0, 0,
                                                     GRETL_MOD_NONE);
                        if (n > 0) {
                            gretl_matrix_inscribe_matrix(V, K->VY, K->r, K->r,
                                                         GRETL_MOD_NONE);
                        }
                        *err = gretl_matrix_psd_root(V, 0);
                        if (*err) {
                            gretl_errmsg_set("Failed to compute factor of Omega_t");
                        } else {
                            load_from_row(Ut, U, t, GRETL_MOD_NONE);
                            gretl_matrix_multiply_mod(Ut, GRETL_MOD_NONE,
                                                      V, GRETL_MOD_TRANSPOSE,
                                                      Et, GRETL_MOD_NONE);
                            load_to_row(E, Et, t);
                        }
                    }
                }

                gretl_matrix_free(Ut);
                gretl_matrix_free(Et);
            } else {
                gretl_matrix_inscribe_matrix(V, K->VS, 0, 0,
                                             GRETL_MOD_NONE);
                if (n > 0) {
                    gretl_matrix_inscribe_matrix(V, K->VY, K->r, K->r,
                                                 GRETL_MOD_NONE);
                }
                *err = gretl_matrix_psd_root(V, 0);
                if (*err) {
                    gretl_errmsg_set("Failed to compute factor of Omega");
                } else {
                    gretl_matrix_multiply_mod(U, GRETL_MOD_NONE,
                                              V, GRETL_MOD_TRANSPOSE,
                                              E, GRETL_MOD_NONE);
                }
            }

            gretl_matrix_free(V);
        }
    }

 bailout:

    K->t = 0;
    set_kalman_stopped(K);

    if (E == NULL && !*err) {
        *err = E_ALLOC;
    } else if (E != NULL && *err) {
        gretl_matrix_free(E);
        E = NULL;
    }

    return E;
}

/*
   below: functions to support the "wrapping" of a Kalman
   struct in a gretl bundle
*/

static int check_replacement_dims (const gretl_matrix *orig,
                                   const gretl_matrix *repl,
                                   int id)
{
    int err = 0;

    if (id == K_ZT || id == K_VY || id == K_T || id == K_VS ||
        id == K_BT || id == K_m || id == K_a || id == K_P || id == K_R) {
        if (repl->rows != orig->rows || repl->cols != orig->cols) {
            err = E_DATA;
        }
    } else if (id == K_y || id == K_x) {
        if (repl->cols != orig->cols) {
            err = E_DATA;
        }
    }

    if (err) {
        gretl_errmsg_set("You cannot resize a state-space "
                         "system matrix");
    }

    return err;
}

/* On input @targ is a pointer to the matrix to be replaced, if
   there's already a matrix corresponding to @id in place; @src is
   the new matrix, to be copied in if @copy is non-zero or
   otherwise just attached. On output @targ points to the new
   matrix.
*/

static int add_or_replace_k_matrix (kalman *K,
                                    gretl_matrix **targ,
                                    gretl_matrix *src,
                                    int id, int copy)
{
    int err = 0;

    if (*targ != src) {
        if (*targ != NULL) {
            err = check_replacement_dims(*targ, src, id);
            if (err) {
                return err;
            }
            /* destroy old Kalman-owned matrix */
            gretl_matrix_free(*targ);
        }
        if (copy) {
            *targ = gretl_matrix_copy(src);
            if (*targ == NULL) {
                err = E_ALLOC;
            }
        } else {
            *targ = src;
        }
    }

    return err;
}

static gretl_matrix **get_input_matrix_target_by_id (kalman *K, int i)
{
    gretl_matrix **targ = NULL;

    if (i == K_T) {
        targ = &K->T;
    } else if (i == K_BT) {
        targ = &K->BT;
    } else if (i == K_ZT) {
        targ = &K->ZT;
    } else if (i == K_VS) {
	/* variance of state */
        if (kalman_xcorr(K)) {
            targ = &K->H;
	} else if (kalman_dkvar(K)) {
	    targ = &K->Q;
        } else {
            targ = &K->VS;
        }
    } else if (i == K_VY) {
	/* variance of observable */
        if (kalman_xcorr(K)) {
            targ = &K->G;
        } else {
            targ = &K->VY;
        }
    } else if (i == K_m) {
        targ = &K->mu;
    } else if (i == K_y) {
        targ = &K->y;
    } else if (i == K_x) {
        targ = &K->x;
    } else if (i == K_a) {
        targ = &K->aini;
    } else if (i == K_P) {
        targ = &K->Pini;
    } else if (i == K_R) {
	targ = &K->R;
    }

    return targ;
}

/* Try attaching a matrix to a Kalman bundle: similar to
   kalman_bundle_set_matrix() below, but allowing for the possibility
   that the @data input of type @vtype has to be converted first.
*/

static int
kalman_bundle_try_set_matrix (kalman *K, void *data,
                              GretlType vtype, int id,
                              int copy)
{
    gretl_matrix **targ;
    int err = 0;

    /* determine the location for this matrix */
    targ = get_input_matrix_target_by_id(K, id);

    if (targ == NULL) {
        err = E_DATA;
    } else {
        gretl_matrix *m;

        if (vtype == GRETL_TYPE_MATRIX) {
            m = data;
            err = add_or_replace_k_matrix(K, targ, m, id, copy);
        } else if (vtype == GRETL_TYPE_DOUBLE) {
            m = gretl_matrix_from_scalar(*(double *) data);
            if (m == NULL) {
                err = E_ALLOC;
            } else {
                err = add_or_replace_k_matrix(K, targ, m, id, 0);
            }
        } else {
            err = E_TYPES;
        }
    }

    return err;
}

/* Called by kalman_deserialize() when reconstructing a kalman bundle
   from XML, and also by kalman_bundle_copy() when duplicating such a
   bundle. In these contexts we know we have a gretl_matrix on input.
*/

static int
kalman_bundle_set_matrix (kalman *K, gretl_matrix *m,
                          int i, int copy)
{
    gretl_matrix **targ;
    int err = 0;

    targ = get_input_matrix_target_by_id(K, i);

    if (targ == NULL) {
        err = E_DATA;
    } else {
        err = add_or_replace_k_matrix(K, targ, m, i, copy);
    }

    return err;
}

static gretl_matrix **kalman_output_matrix (kalman *K,
                                            const char *key)
{
    gretl_matrix **pm = NULL;

    if (!strcmp(key, "prederr")) {
        pm = &K->V;
    } else if (!strcmp(key, "pevar")) {
        pm = &K->F;
    } else if (!strcmp(key, "state")) {
        pm = &K->A;
    } else if (!strcmp(key, "stvar")) {
        pm = &K->P;
    } else if (!strcmp(key, "gain")) {
        pm = &K->K;
    } else if (!strcmp(key, "llt")) {
        pm = &K->LL;
    } else if (!strcmp(key, "smdist")) {
        pm = &K->U;
    } else if (!strcmp(key, "smdisterr")) {
        pm = &K->Vsd;
    } else if (!strcmp(key, "uhat")) {
        pm = &K->vt;
    }

    return pm;
}

#define K_N_OUTPUTS 9

static const char *kalman_output_matrix_names[K_N_OUTPUTS] = {
    "prederr",
    "pevar",
    "state",
    "stvar",
    "gain",
    "llt",
    "smdist",
    "smdisterr",
    "uhat"
};

static int output_matrix_slot (const char *s)
{
    int i;

    for (i=0; i<K_N_OUTPUTS; i++) {
        if (!strcmp(s, kalman_output_matrix_names[i])) {
            return i;
        }
    }

    return -1;
}

#define K_N_SCALARS 14

enum {
    Ks_t = 0,
    Ks_DIFFUSE,
    Ks_EXACT,
    Ks_CROSS,
    Ks_DKVAR,
    Ks_UNI,
    Ks_S2,
    Ks_LNL,
    Ks_r,
    Ks_n,
    Ks_N,
    Ks_p,
    Ks_d,
    Ks_j
};

static const char *kalman_output_scalar_names[K_N_SCALARS] = {
    "t",
    "diffuse",
    "exact",
    "cross",
    "dkvar",
    "univariate",
    "s2",
    "lnl",
    "r",
    "n",
    "N",
    "p",
    "d",
    "j"
};

static double *kalman_output_scalar (kalman *K,
                                     const char *key)
{
    /* static storage for on-the-fly scalars */
    static double retval[K_N_SCALARS];
    int i, idx = -1;

    for (i=0; i<K_N_SCALARS; i++) {
        if (!strcmp(key, kalman_output_scalar_names[i])) {
            idx = i;
            break;
        }
    }

    if (idx < 0 && !strcmp(key, "T")) {
        /* backward compatibility */
        idx = Ks_N;
    }

    if (idx < 0) {
        return NULL;
    }

    switch (idx) {
    case Ks_t:
        if (kalman_is_running(K)) {
            retval[idx] = K->t + 1;
        } else {
            retval[idx] = kalman_checking(K) ? 1 : 0;
        }
        break;
    case Ks_DIFFUSE:
        retval[idx] = (K->flags & KALMAN_DIFFUSE)? 1 : 0;
        break;
    case Ks_EXACT:
        retval[idx] = K->exact;
        break;
    case Ks_CROSS:
        retval[idx] = (K->vartype == DJ_VAR);
        break;
    case Ks_DKVAR:
        retval[idx] = (K->vartype == DK_VAR);
        break;
    case Ks_UNI:
        retval[idx] = (K->flags & KALMAN_UNI)? 1 : 0;
        break;
    case Ks_S2:
        retval[idx] = K->s2;
        break;
    case Ks_LNL:
        retval[idx] = K->loglik;
        break;
    case Ks_r:
        retval[idx] = K->r;
        break;
    case Ks_n:
        retval[idx] = K->n;
        break;
    case Ks_N:
        retval[idx] = K->N;
        break;
    case Ks_p:
        retval[idx] = K->p;
        break;
    case Ks_d:
        retval[idx] = K->d;
        break;
    case Ks_j:
        retval[idx] = K->j;
        break;
    default:
        break;
    }

    return &retval[idx];
}

static const gretl_matrix *k_input_matrix_by_id (kalman *K, int i)
{
    const gretl_matrix *m = NULL;

    if (i == K_T) {
        m = K->T;
    } else if (i == K_BT) {
        m = K->BT;
    } else if (i == K_ZT) {
        m = K->ZT;
    } else if (i == K_VS) {
        if (kalman_xcorr(K)) {
            m = K->H;
	} else if (kalman_dkvar(K)) {
	    m = K->Q;
        } else {
            m = K->VS;
        }
    } else if (i == K_VY) {
        if (kalman_xcorr(K)) {
            m = K->G;
        } else {
            m = K->VY;
        }
    } else if (i == K_m) {
        m = K->mu;
    } else if (i == K_y) {
        m = K->y;
    } else if (i == K_x) {
        m = K->x;
    } else if (i == K_a) {
        m = K->aini;
    } else if (i == K_P) {
        m = K->Pini;
    } else if (i == K_R) {
	m = K->R;
    }

    return m;
}

static int input_matrix_slot (const char *s)
{
    int i;

    for (i=0; i<K_MMAX; i++) {
        if (!strcmp(s, K_input_mats[i].name)) {
            return K_input_mats[i].sym;
        }
    }

    return -1;
};

static GretlType kalman_extra_type (const char *key)
{
    if (!strcmp(key, "ssfsim")) {
        return GRETL_TYPE_DOUBLE;
    } else if (!strcmp(key, "simstart") ||
               !strcmp(key, "simx")) {
        return GRETL_TYPE_MATRIX;
    } else {
        return GRETL_TYPE_NONE;
    }
}

/* respond to "diffuse" setting */

static int kalman_set_diffuse (kalman *K, int d)
{
    if (d != 0 && d != 1 && d != 2) {
	return E_INVARG;
    } else if (d) {
        K->exact = (d == 2);
        K->flags |= KALMAN_DIFFUSE;
        return diffuse_Pini(K);
    } else {
        K->exact = 0;
        K->flags &= ~KALMAN_DIFFUSE;
        if (K->Pk0 != NULL) {
            gretl_matrix_free(K->Pk0);
            gretl_matrix_free(K->Pk1);
            gretl_matrix_free(K->Fk);
            K->Pk0 = K->Pk1 = K->Fk = NULL;
        }
        return 0;
    }
}

/* respond to "univariate" setting */

static int kalman_set_univariate (kalman *K, int u)
{
    if (u && K->vartype == DJ_VAR) {
	gretl_errmsg_set("kalman: the 'univariate' setting is not compatible with\n"
			 "cross-correlated disturbances");
	return E_INVARG;
    }

    if (u) {
	K->flags |= KALMAN_UNI;
    } else {
	K->flags &= ~KALMAN_UNI;
    }

    return 0;
}

/* Called by real_bundle_set_data() in gretl_bundle.c.  The return
   value indicates whether the putative setting was handled (1) or not
   (0). Not being handled here is not necessarily an error.

   The @copy flag here is inherited from the specific caller of
   real_bundle_set_data(): @copy = 1 if the caller was
   gretl_bundle_set_data(), 0 if it was gretl_bundle_donate_data().
   Either way the kalman struct takes ownership.
*/

int maybe_set_kalman_element (void *kptr,
                              const char *key,
                              void *vptr,
                              GretlType vtype,
                              int copy,
                              int *err)
{
    GretlType targtype;
    kalman *K = kptr;
    int fncall = 0;
    int i, id = -1;
    int done = 0;

    if (K == NULL) {
        *err = E_DATA;
        return 0;
    }

    /* Check for optional "extra" kalman items that
       live outside of the kalman struct itself.
    */
    targtype = kalman_extra_type(key);
    if (targtype != GRETL_TYPE_NONE) {
        if (vtype != targtype) {
            *err = E_TYPES;
        }
        return 0;
    }

    if (!strcmp(key, "diffuse") || !strcmp(key, "univariate")) {
	/* scalar config settings */
        if (vtype == GRETL_TYPE_DOUBLE) {
            double v = *(double *) vptr;

	    if (!strcmp(key, "diffuse")) {
		*err = kalman_set_diffuse(K, (int) v);
	    } else {
		*err = kalman_set_univariate(K, (int) v);
	    }
        } else {
            *err = E_TYPES;
        }
        return 1; /* done, error or no */
    }

    if (!strcmp(key, "timevar_call")) {
        /* try for a function call specifier (string) */
        if (vtype == GRETL_TYPE_STRING) {
            fncall = 1;
        } else {
            *err = E_TYPES;
        }
    } else {
        /* try for a matrix specifier */
        i = input_matrix_slot(key);
        if (i >= 0) {
            if (vtype == GRETL_TYPE_MATRIX ||
                vtype == GRETL_TYPE_DOUBLE) {
                id = i;
            } else {
                *err = E_TYPES;
            }
        }
    }

    if (*err) {
        return 0;
    } else if (fncall) {
        if (copy) {
            K->matcall = gretl_strdup((char *) vptr);
        } else {
            K->matcall = (char *) vptr;
        }
        /* re-evaluate what's actually varying */
        *err = check_for_matrix_updates(K, NULL);
        if (!*err) {
            done = 1;
        }
    } else if (id < 0) {
        if (kalman_output_matrix(K, key) != NULL ||
            kalman_output_scalar(K, key) != NULL) {
            *err = E_DATA;
            gretl_errmsg_sprintf("The member %s is read-only", key);
        }
    } else {
        *err = kalman_bundle_try_set_matrix(K, vptr, vtype, id, copy);
        done = (*err == 0);
    }

    return done;
}

int maybe_delete_kalman_element (void *kptr,
                                 const char *key,
                                 int *err)
{
    kalman *K = kptr;
    gretl_matrix **pm;
    int done = 0;

    if (K == NULL) {
        return 0;
    }

    if (kalman_output_scalar(K, key) != NULL ||
        input_matrix_slot(key) >= 0 || !strcmp(key, "uhat")) {
        /* note: the matrix under the key "uhat" is part of
           the internal kalman apparatus */
        gretl_errmsg_sprintf("%s: cannot be deleted", key);
        *err = E_DATA;
    } else if ((pm = kalman_output_matrix(K, key)) != NULL) {
        /* OK to delete a user-output matrix */
        gretl_matrix_free(*pm);
        *pm = NULL;
    } else if (!strcmp(key, "timevar_call")) {
        /* OK to delete time-variation call */
        if (K->matcall != NULL) {
            free(K->matcall);
            K->matcall = NULL;
            free(K->varying);
            K->varying = NULL;
            done = 1;
        } else {
            *err = E_DATA;
        }
    }

    return done;
}

void *maybe_retrieve_kalman_element (void *kptr,
                                     const char *key,
                                     GretlType *type,
                                     int *reserved,
                                     int *err)
{
    kalman *K = kptr;
    void *ret = NULL;
    int i, id = -1;

    *type = GRETL_TYPE_NONE;

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    if (!strcmp(key, "timevar_call")) {
        /* function call specifier? */
        *reserved = 1;
        if (K->matcall != NULL) {
            ret = K->matcall;
            *type = GRETL_TYPE_STRING;
        }
    } else {
        /* try for an input matrix specifier */
        for (i=0; i<K_MMAX; i++) {
            if (!strcmp(key, K_input_mats[i].name)) {
                id = K_input_mats[i].sym;
                ret = (gretl_matrix *) k_input_matrix_by_id(K, id);
                if (ret != NULL) {
                    *type = GRETL_TYPE_MATRIX;
                }
                break;
            }
        }
        if (id < 0) {
            /* try for an output matrix */
            gretl_matrix **pm = kalman_output_matrix(K, key);

            if (pm != NULL) {
                *reserved = 1;
                if (*pm != NULL) {
                    ret = *pm;
                    *type = GRETL_TYPE_MATRIX;
                }
            }
        }
        if (id < 0 && *reserved == 0) {
            /* nothing matched yet: try scalar member */
            ret = kalman_output_scalar(K, key);
            if (ret != NULL) {
                *type = GRETL_TYPE_DOUBLE;
            }
        }
    }

    if (id >= 0 && *reserved == 0) {
        /* flag the fact that @key was a kalman-reserved
           identifier */
        *reserved = 1;
    }

    if (*reserved && ret == NULL) {
        gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
        *err = E_DATA;
    }

    return ret;
}

static int output_matrix_count (kalman *K)
{
    gretl_matrix **pm;
    int i, n = 0;

    for (i=0; i<K_N_OUTPUTS; i++) {
        pm = kalman_output_matrix(K, kalman_output_matrix_names[i]);
        n += (pm != NULL && *pm != NULL);
    }

    return n;
}

int print_kalman_bundle_info (void *kptr, PRN *prn)
{
    kalman *K = kptr;
    int err = 0;

    if (K == NULL) {
        pputs(prn, "Kalman struct: empty\n");
        err = E_DATA;
    } else {
        const gretl_matrix *m;
        gretl_matrix **pm;
        double *px;
        const char *name;
        int i, id;

        pputs(prn, "\nKalman input matrices\n");

        for (i=0; i<K_MMAX; i++) {
            id = K_input_mats[i].sym;
            m = k_input_matrix_by_id(K, id);
            if (m != NULL) {
                pprintf(prn, " %s: ", K_input_mats[i].name);
                pprintf(prn, "%d x %d\n", m->rows, m->cols);
            }
        }

        if (output_matrix_count(K) > 0) {
            pputs(prn, "\nKalman output matrices\n");
            for (i=0; i<K_N_OUTPUTS; i++) {
                name = kalman_output_matrix_names[i];
                pm = kalman_output_matrix(K, name);
                if (pm != NULL && *pm != NULL) {
                    m = *pm;
                    pprintf(prn, " %s: ", name);
                    pprintf(prn, "%d x %d\n", m->rows, m->cols);
                }
            }
        }

        pputs(prn, "\nKalman scalars\n");

        for (i=0; i<K_N_SCALARS; i++) {
            name = kalman_output_scalar_names[i];
            pprintf(prn, " %s: ", name);
            px = kalman_output_scalar(K, name);
            if (px == NULL || na(*px)) {
                pputs(prn, "NA\n");
            } else {
                pprintf(prn, "%g\n", *px);
            }
        }

        if (K->matcall != NULL) {
            pputs(prn, "\nKalman strings\n");
            pprintf(prn, " timevar_call: %s\n", K->matcall);
        }
    }

    return err;
}

/* For use in context of a kalman bundle: serialize the information in
   the kalman struct to XML
*/

int kalman_serialize (void *kptr, PRN *prn)
{
    kalman *K = kptr;
    const gretl_matrix *m;
    gretl_matrix **pm;
    double *px;
    const char *name;
    int i, err = 0;

    if (K == NULL) {
        fputs("kalman_serialize: got NULL\n", stderr);
        return E_DATA;
    }

    pputs(prn, "<gretl-kalman>\n");

    for (i=0; i<K_MMAX; i++) {
        m = k_input_matrix_by_id(K, K_input_mats[i].sym);
        if (m != NULL) {
            gretl_matrix_serialize(m, K_input_mats[i].name, prn);
        }
    }

    for (i=0; i<K_N_OUTPUTS; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != NULL && *pm != NULL) {
            gretl_matrix_serialize(*pm, name, prn);
        }
    }

    for (i=0; i<=Ks_LNL; i++) {
        name = kalman_output_scalar_names[i];
        px = kalman_output_scalar(K, name);
        if (px != NULL && !na(*px)) {
            gretl_finite_scalar_serialize(*px, name, prn);
        }
    }

    if (K->matcall != NULL) {
        gretl_string_serialize(K->matcall, "timevar_call", prn);
    }

    pputs(prn, "</gretl-kalman>\n");

    return err;
}

static int required_matrix_slot (const char *s)
{
    if (!strcmp(s, "obsy"))     return 0;
    if (!strcmp(s, "obsymat"))  return 1;
    if (!strcmp(s, "statemat")) return 2;
    if (!strcmp(s, "statevar")) return 3;
    if (!strcmp(s, "obsvar"))   return 4;
    return -1;
};

/* For use in context of a kalman bundle: deserialize the kalman
   struct from XML
*/

gretl_bundle *kalman_deserialize (void *p1, void *p2, int *err)
{
    xmlNodePtr cur, node = p1;
    xmlDocPtr doc = p2;
    gretl_matrix *Mreq[5] = {NULL};
    gretl_matrix *Mopt[K_MMAX] = {NULL};
    gretl_matrix *Mout[K_N_OUTPUTS] = {NULL};
    char *tvcall = NULL;
    double s2 = NADBL;
    double lnl = NADBL;
    int copy[5] = {0};
    int i, nmats = 0;
    int Kflags = 0;
    int vtype = 0;
    gretl_matrix *m;
    double x;
    char *key, *strv;
    gretl_bundle *b = NULL;

    while (node != NULL && !*err) {
        if (!xmlStrcmp(node->name, (XUC) "gretl-kalman")) {
            cur = node->xmlChildrenNode;
            while (cur != NULL && !*err) {
                key = (char *) xmlGetProp(cur, (XUC) "name");
                if (!xmlStrcmp(cur->name, (XUC) "gretl-matrix")) {
                    /* pick up kalman matrices */
                    m = gretl_xml_get_matrix(cur, doc, err);
                    if ((i = required_matrix_slot(key)) >= 0) {
                        nmats++;
                        Mreq[i] = m;
                    } else if ((i = input_matrix_slot(key)) >= 0) {
                        Mopt[i] = m;
                    } else if ((i = output_matrix_slot(key)) >= 0) {
                        Mout[i] = m;
                    }
                } else if (!xmlStrcmp(cur->name, (XUC) "scalar")) {
                    /* pick up kalman scalars */
                    if (gretl_xml_get_prop_as_double(cur, "value", &x)) {
                        if (!strcmp(key, "diffuse") && x > 0) {
                            Kflags |= KALMAN_DIFFUSE;
                        } else if (!strcmp(key, "cross") && x > 0) {
			    vtype = DJ_VAR;
			} else if (!strcmp(key, "dkvar") && x > 0) {
			    vtype = DK_VAR;
                        } else if (!strcmp(key, "s2")) {
                            s2 = x;
                        } else if (!strcmp(key, "lnl")) {
                            lnl = x;
                        }
                    }
                } else if (!xmlStrcmp(cur->name, (XUC) "string")) {
                    /* pick up kalman strings */
                    if (!strcmp(key, "timevar_call") &&
                        gretl_xml_get_prop_as_string(cur, "value", &strv)) {
                            tvcall = strv;
                    }
                }
                free(key);
                cur = cur->next;
            }
            break;
        }
        node = node->next;
    }

    if (vtype == DK_VAR) {
	if (Mreq[K_VS] == NULL || Mopt[K_R] == NULL) {
	    *err = E_DATA;
	    goto bailout;
	} else {
	    Mopt[K_VY] = Mreq[4];
	    Mreq[4] = Mopt[K_R];
	    Mopt[K_R] = NULL;
	}
    } else if (vtype == STD_VAR && nmats == 5) {
        /* drop obsvar from initialization */
        Mopt[K_VY] = Mreq[4];
        Mreq[4] = NULL;
        nmats--;
    }

    if ((vtype > 0 && nmats != 5) || (vtype == 0 && nmats != 4)) {
        *err = E_DATA;
    } else {
        b = kalman_bundle_new(Mreq, copy, nmats, vtype == DK_VAR, err);
        if (b != NULL) {
            kalman *K = gretl_bundle_get_private_data(b);
            gretl_matrix **pm;
            const char *name;

            K->flags = Kflags;
	    K->vartype = vtype;
            K->s2 = s2;
            K->loglik = lnl;

            for (i=0; i<K_MMAX; i++) {
                if (Mopt[i] != NULL) {
                    kalman_bundle_set_matrix(K, Mopt[i], i, 0);
                }
            }
            for (i=0; i<K_N_OUTPUTS; i++) {
                if (Mout[i] != NULL) {
                    name = kalman_output_matrix_names[i];
                    pm = kalman_output_matrix(K, name);
                    *pm = Mout[i];
                }
            }
            K->matcall = tvcall;
        }
    }

 bailout:

    if (*err) {
        /* clean up */
        for (i=0; i<5; i++) {
            gretl_matrix_free(Mreq[i]);
        }
        for (i=0; i<K_MMAX; i++) {
            gretl_matrix_free(Mopt[i]);
        }
        for (i=0; i<K_N_OUTPUTS; i++) {
            gretl_matrix_free(Mout[i]);
        }
        free(tvcall);
    }

    return b;
}

/* Called from gretl_bundle.c to meet the case where the user calls
   for a kalman bundle to be copied: here we create a new kalman
   struct and copy across the required elements (since they are
   not regular bundle members).
*/

gretl_bundle *kalman_bundle_copy (const gretl_bundle *src, int *err)
{
    kalman *K, *Knew;
    gretl_bundle *b = NULL;
    gretl_matrix *M[5] = {NULL};
    int copy[5] = {1, 1, 1, 1, 1};
    gretl_matrix *m, **pm, **pm1;
    const char *name;
    int dkopt = 0;
    int i, id, k = 4;

    K = gretl_bundle_get_private_data((gretl_bundle *) src);

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    /* set pointers to required Kalman matrices */
    M[0] = K->y;
    M[1] = K->ZT;
    M[2] = K->T;

    /* variants dependent on presence/absence of cross-correlation */
    if (kalman_xcorr(K)) {
        M[3] = K->H;
        M[4] = K->G;
        k = 5;
    } else if (kalman_dkvar(K)) {
        M[3] = K->Q;
        M[4] = K->R;
        k = 5;
	dkopt = 1;
    } else {
        M[3] = K->VS;
    }

    b = kalman_bundle_new(M, copy, k, dkopt, err);

    if (*err) {
        return b;
    }

    Knew = gretl_bundle_get_private_data(b);
    Knew->flags = K->flags;

    /* add any "extra" matrices, beyond the required ones */
    for (i=0; i<n_extra_mats && !*err; i++) {
        id = extra_mats[i];
        m = (gretl_matrix *) k_input_matrix_by_id(K, id);
        if (m != NULL) {
            *err = kalman_bundle_set_matrix(Knew, m, id, 1);
        }
    }

    for (i=0; i<K_N_OUTPUTS && !*err; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != NULL && *pm != NULL) {
            pm1 = kalman_output_matrix(Knew, name);
            *pm1 = gretl_matrix_copy(*pm);
        }
    }

    Knew->ifc = K->ifc;
    Knew->s2 = K->s2;
    Knew->loglik = K->loglik;
    Knew->vartype = K->vartype;

    if (K->flags & KALMAN_DIFFUSE) {
        Knew->flags |= KALMAN_DIFFUSE;
    }

    if (K->matcall != NULL) {
        Knew->matcall = gretl_strdup(K->matcall);
    }

    return b;
}

/* for use in constructing GUI bundle save menu */

char **kalman_bundle_get_matrix_names (kalman *K, int *ns)
{
    char **S = NULL;
    gretl_matrix **pm;
    const char *name;
    int i, id, err = 0;

    *ns = 0;

    for (i=0; i<K_MMAX && !err; i++) {
        id = K_input_mats[i].sym;
        if (k_input_matrix_by_id(K, id) != NULL) {
            err = strings_array_add(&S, ns, K_input_mats[i].name);
        }
    }

    for (i=0; i<K_N_OUTPUTS && !err; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != NULL && *pm != NULL) {
            err = strings_array_add(&S, ns, name);
        }
    }

    return S;
}

/* also for use in constructing GUI bundle save menu */

char **kalman_bundle_get_scalar_names (kalman *K, int *ns)
{
    char **S;

    *ns = K_N_SCALARS - na(K->s2) - na(K->loglik);
    S = strings_array_new(*ns);

    if (S != NULL) {
	int i, j = 0;

	for (i=0; i<K_N_SCALARS; i++) {
	    if (i == Ks_S2 && na(K->s2)) {
		continue;
	    } else if (i == Ks_LNL && na(K->loglik)) {
		continue;
	    } else {
		S[j++] = gretl_strdup(kalman_output_scalar_names[i]);
	    }
	}
    }

    return S;
}

/* to support the nelem() function for kalman bundles */

int kalman_bundle_n_members (gretl_bundle *b)
{
    kalman *K = gretl_bundle_get_private_data(b);
    int n = 0;

    if (K != NULL) {
        int i, id;

        n = K_N_SCALARS;

        for (i=0; i<K_MMAX; i++) {
            id = K_input_mats[i].sym;
            n += (k_input_matrix_by_id(K, id) != NULL);
        }

        n += output_matrix_count(K);
    }

    return n;
}

/* not great placement, but for now... here comes some univariate code */

static void dj_from_Finf (const gretl_matrix *Fk, int *d, int *j)
{
    int rmax = Fk->rows - 1;
    int cmax = Fk->cols - 1;
    int i, k;

    *d = *j = 0;

    for (i=rmax; i >= 0 && *d == 0; i--) {
	for (k=cmax; k>=0; k--) {
	    if (gretl_matrix_get(Fk, i, k) > 0) {
		*d = i + 1;
		*j = k + 1;
		break;
	    }
	}
    }
}

static gretl_matrix *first_n_rows (gretl_matrix *m, int n)
{
    gretl_matrix *ret = gretl_matrix_alloc(n, m->cols);
    size_t colsize = n * sizeof(double);
    double *dest = ret->val;
    double *src = m->val;
    int j;

    for (j=0; j<m->cols; j++) {
	memcpy(dest, src, colsize);
	dest += n;
	src += m->rows;
    }

    gretl_matrix_free(m);

    return ret;
}

static int bundle_add_matrix (gretl_bundle *b,
			      const char *key,
			      gretl_matrix *m)
{
    return gretl_bundle_donate_data(b, key, m, GRETL_TYPE_MATRIX, 0);
}

/* write row @i of matrix @src into matrix @targ */

static void mat_from_row (gretl_matrix *targ, const gretl_matrix *src, int i)
{
    int j;

    for (j=0; j<src->cols; j++) {
	targ->val[j] = gretl_matrix_get(src, i, j);
    }
}

/* write vec of matrix @src into row @t of @targ */

void row_from_mat (gretl_matrix *targ, const gretl_matrix *src, int t)
{
    int i, n = src->rows * src->cols;

    for (i=0; i<n; i++) {
	gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* write content of vector @src into column @j of @targ */

static void col_from_vec (gretl_matrix *targ, const gretl_vector *src, int j)
{
    int i;

    for (i=0; i<targ->rows; i++) {
	gretl_matrix_set(targ, i, j, src->val[i]);
    }
}

/* write column @j of src into vector @targ */

static void vec_from_col (gretl_vector *targ, const gretl_matrix *src, int j)
{
    int i;

    for (i=0; i<targ->rows; i++) {
	targ->val[i] = gretl_matrix_get(src, i, j);
    }
}

/* simple wrapper ignoring error check */

static double dotprod (const gretl_vector *a, const gretl_vector *b)
{
    int err = 0;
    return gretl_vector_dot_product(a, b, &err);
}

static void increment_state (gretl_matrix *a,
			     const gretl_matrix *b,
			     gretl_matrix *tmp,
			     double x)
{
    gretl_matrix_copy_values(tmp, b);
    gretl_matrix_multiply_by_scalar(tmp, x);
    gretl_matrix_add_to(a, tmp);
}

static void decrement_state_var (gretl_matrix *P,
				 const gretl_matrix *K,
				 gretl_matrix *tmp,
				 double x)
{
    gretl_matrix_multiply_mod(K, GRETL_MOD_NONE,
			      K, GRETL_MOD_TRANSPOSE,
			      tmp, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(tmp, x);
    gretl_matrix_subtract_from(P, tmp);
}

static void state_cross_update (gretl_matrix *Pti,
				const gretl_matrix *Kti,
				const gretl_matrix *Kki,
				gretl_matrix *tmp,
				double Fti, double Fkinv)
{
    gretl_matrix_multiply_mod(Kki, GRETL_MOD_NONE,
			      Kki, GRETL_MOD_TRANSPOSE,
			      tmp, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(tmp, Fti * Fkinv * Fkinv);
    gretl_matrix_add_to(Pti, tmp);
    gretl_matrix_multiply_mod(Kti, GRETL_MOD_NONE,
			      Kki, GRETL_MOD_TRANSPOSE,
			      tmp, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(Kki, GRETL_MOD_NONE,
			      Kti, GRETL_MOD_TRANSPOSE,
			      tmp, GRETL_MOD_CUMULATE);
    gretl_matrix_multiply_by_scalar(tmp, Fkinv);
    gretl_matrix_subtract_from(Pti, tmp);
}

/* smoother-related matrices */

struct smo_recorder {
    gretl_matrix *PK;
    gretl_matrix *Fk;
    gretl_matrix *Kinf;
};

static int allocate_smo_recorders (struct smo_recorder *smr,
				   int m, int p,
				   int N, int Nd)
{
    smr->PK   = gretl_zero_matrix_new(Nd, m*m);
    smr->Fk   = gretl_zero_matrix_new(Nd, p);
    smr->Kinf = gretl_zero_matrix_new(Nd, p*m);

    return 0;
}

static gretl_matrix *kalman_ldl (const gretl_matrix *V,
				 gretl_matrix **pLinv,
				 int *err)
{
    gretl_matrix *d, *L;
    double dj, lij;
    int n = V->rows;
    int i, j;

    L = gretl_matrix_copy(V);
    *err = gretl_matrix_cholesky_decomp(L);
    if (*err) {
	return NULL;
    }

    d = gretl_matrix_alloc(n, 1);
    for (j=0; j<n; j++) {
	dj = gretl_matrix_get(L, j, j);
	d->val[j] = dj * dj;
	for (i=j; i<n; i++) {
	    lij = gretl_matrix_get(L, i, j);
	    gretl_matrix_set(L, i, j, lij / dj);
	}
    }

    *err = gretl_invert_triangular_matrix(L, 'L');
    if (*err) {
	gretl_matrix_free(d);
	gretl_matrix_free(L);
	d = *pLinv = NULL;
    } else {
	*pLinv = L;
    }

    return d;
}

static gretl_vector *kalman_diagonalize (kalman *K,
					 gretl_matrix *Z,
					 int *err)
{
    gretl_matrix *d, *Linv = NULL;
    int p = Z->rows;
    int m = Z->cols;

    d = kalman_ldl(K->VY, &Linv, err);

    if (!*err) {
	gretl_matrix *tmp = gretl_matrix_alloc(p, m);
	gretl_matrix *yt = gretl_matrix_alloc(p, 1);
	int t;

	gretl_matrix_zero(K->VY);
	gretl_matrix_set_diagonal(K->VY, d, 0);
	gretl_matrix_multiply(Linv, Z, tmp);
	gretl_matrix_copy_values(Z, tmp);
	gretl_matrix_reuse(tmp, p, 1);

	for (t=0; t<K->N; t++) {
	    mat_from_row(tmp, K->y, t);
	    gretl_matrix_multiply(Linv, tmp, yt);
	    row_from_mat(K->y, yt, t);
	}

	bundle_add_matrix(K->b, "Linv", Linv);

	gretl_matrix_free(tmp);
	gretl_matrix_free(yt);
    }

    return d;
}

static int handle_matrix_transforms (kalman *K,
				     gretl_matrix **pg,
				     gretl_matrix **pZ,
				     int *gtrans,
				     int *Ztrans)
{
    int err = 0;

    if (K->n == 1 && K->r == 1) {
	/* nothing to be done */
	*pg = K->VY;
	*pZ = K->ZT;
	return 0;
    }

    /* first see if any required transformations have
       already been done */
    if (gretl_bundle_has_key(K->b, "g") &&
	gretl_bundle_has_key(K->b, "Z")) {
	*pg = gretl_bundle_get_matrix(K->b, "g", &err);
	*pZ = gretl_bundle_get_matrix(K->b, "Z", &err);
	return err;
    }

    /* ZT (at least) needs modification */
    *pZ = gretl_matrix_alloc(K->n, K->r);
    gretl_matrix_transpose(*pZ, K->ZT);
    *Ztrans = 1;

    if (K->n > 1 && K->VY != NULL) {
	/* obsvar needs transformation */
	if (matrix_is_diagonal(K->VY)) {
	    *pg = gretl_matrix_get_diagonal(K->VY, &err);
	} else {
	    printf("kfilter: diagonalizing\n");
	    *pg = kalman_diagonalize(K, *pZ, &err);
	}
	*gtrans = 1;
    } else if (K->VY != NULL) {
	*pg = K->VY;
    }

    return err;
}

static int kfilter_univariate (kalman *K, PRN *prn)
{
    struct smo_recorder sm = {0};
    gretl_matrix *Kkt = K->Mt;
    gretl_matrix *at = K->a0;
    gretl_matrix *Pt = K->P0;
    gretl_matrix *ati = K->a1;
    gretl_matrix *Pti = K->P1;
    gretl_matrix *Pk = K->Pk0;
    gretl_matrix *g = NULL;
    gretl_matrix *Z = NULL;

    /* workspace matrices */
    gretl_matrix_block *B;
    gretl_matrix *Zi, *Kti;
    gretl_matrix *Pki = NULL;
    gretl_matrix *Kki = NULL;
    gretl_matrix *m1, *mm;

    int p = K->n; /* # of observables */
    int m = K->r; /* length of state vector */
    int t, i, err = 0;

    int smo = kalman_smoothing(K);
    int all_kappa = 0; /* FIXME conditionality */
    int TI = gretl_is_identity_matrix(K->T);

    double yti, vti, Fti, llct;
    double Ftinv, Fkinv, Fki = 0;
    double qt, ldt;
    double l2pi = log(M_2PI);
    int rankPk = kalman_diffuse(K)? m : 0;
    int d = K->exact ? 0 : -1;
    int j = K->exact ? 0 : -1;
    int gtrans = 0;
    int Ztrans = 0;
    int Nd = 0;

    /* extra matrices to be returned via k->b */
    gretl_matrix *P = gretl_zero_matrix_new(K->N, m*m);
    gretl_matrix *att = gretl_matrix_alloc(K->N, m);
    gretl_matrix *Ptt = gretl_matrix_alloc(K->N, m*m);

    if (P == NULL) {
	err = E_ALLOC;
    } else if (K->exact) {
	B = gretl_matrix_block_new(&Pki, m, m,
				   &Zi, 1, m,
				   &Kti, m, 1,
				   &Kki, m, 1,
				   &m1, m, 1, &mm, m, m,
				   NULL);
	if (B == NULL) err = E_ALLOC;
    } else {
	B = gretl_matrix_block_new(&Zi, 1, m,
				   &Kti, m, 1,
				   &m1, m, 1, &mm, m, m,
				   NULL);
	if (B == NULL) err = E_ALLOC;
    }

    if (!err) {
	err = handle_matrix_transforms(K, &g, &Z, &gtrans, &Ztrans);
    }

    if (err) {
	return err;
    }

    if (smo && K->exact) {
	Nd = all_kappa ? K->N : 4*m;
	allocate_smo_recorders(&sm, m, p, K->N, Nd);
    }

    K->SSRw = 0;
    K->loglik = 0;

    printf("\nCCC kfilt, m=%d, p=%d CCC\n", m, p);

    for (t=0; t<K->N; t++) {
	load_to_vec(K->A, at, t);
	row_from_mat(P, Pt, t);
	fast_copy_values(ati, at);
	fast_copy_values(Pti, Pt);
	if (Pki != NULL) {
	    fast_copy_values(Pki, Pk);
	}
	qt = 0;
	ldt = 0;
	llct = 0;
	if (smo && d == 0 && t < Nd) {
	    row_from_mat(sm.PK, Pk, t);
	}
	for (i=0; i<p; i++) {
	    if (t < 2 && K->exact) {
		printf("t,i = %d,%d, rankPk = %d\n", t, i, rankPk);
	    }
	    yti = gretl_matrix_get(K->y, t, i);
	    mat_from_row(Zi, Z, i);
	    gretl_matrix_multiply_mod(Pti, GRETL_MOD_NONE,
	    			      Zi, GRETL_MOD_TRANSPOSE,
	    			      Kti, GRETL_MOD_NONE);
	    Fti = dotprod(Zi, Kti);
	    if (g != NULL) {
		Fti += g->val[i];
	    }
	    if (smo) {
		gretl_matrix_set(K->F, t, i, Fti);
		col_from_vec(K->Kt, Kti, i);
	    }
	    if (d == 0) {
		/* still initial */
		gretl_matrix_multiply_mod(Pki, GRETL_MOD_NONE,
					  Zi, GRETL_MOD_TRANSPOSE,
					  Kki, GRETL_MOD_NONE);
		Fki = dotprod(Zi, Kki);
		if (smo && t < Nd) {
		    gretl_matrix_set(sm.Fk, t, i, Fki);
		    col_from_vec(Kkt, Kki, i);
		}
	    }

	    if (na(yti)) {
		gretl_matrix_set(K->V, t, i, NADBL);
		continue;
	    } else {
		vti = yti - dotprod(Zi, ati);
		gretl_matrix_set(K->V, t, i, vti);
	    }

	    if (rankPk > 0) {
		if (Fki > 0) {
		    Fkinv = 1.0 / Fki;
		    ldt += log(Fki);
		    increment_state(ati, Kki, m1, Fkinv * vti);
		    state_cross_update(Pti, Kti, Kki, mm, Fti, Fkinv);
		    decrement_state_var(Pki, Kki, mm, Fkinv);
		    --rankPk;
		} else if (Fti > 0) {
		    Ftinv = 1.0 / Fti;
		    qt += vti * vti * Ftinv;
		    ldt += log(Fti);
		    llct += l2pi;
		    increment_state(ati, Kti, m1, Ftinv * vti);
		    decrement_state_var(Pti, Kti, mm, Ftinv);
		}
	    } else {
		/* rankPk = 0 */
		if (j == 0) {
		    j = (p == 1)? 1 : i;
		    printf("*** filter: got j = %d (i=%d)\n", j, i);
		}
		Ftinv = 1/Fti;
		qt += vti * vti * Ftinv;
		ldt += log(Fti);
		llct += l2pi;
		increment_state(ati, Kti, m1, Ftinv * vti);
		decrement_state_var(Pti, Kti, mm, Ftinv);
	    }
	} /* observables at time t */

	if (smo) {
	    row_from_mat(K->K, K->Kt, t);
	    if (d == 0 && t < Nd) {
		row_from_mat(sm.Kinf, Kkt, t);
	    }
	}
	if (j > 0 && d == 0) {
	    d = t + 1; /* note: 1-based */
	    printf("*** filter: got d = %d\n", d);
	}

	/* update for t+1 */
	if (TI) {
	    /* shortcut in case T = I_m */
	    gretl_matrix_copy_values(at, ati);
	    gretl_matrix_copy_values(Pt, Pti);
	    gretl_matrix_add_to(Pt, K->VS);
	    gretl_matrix_copy_values(Pk, Pki);
	} else {
	    gretl_matrix_multiply(K->T, ati, at);
	    fast_copy_values(Pt, K->VS);
	    gretl_matrix_qform(K->T, GRETL_MOD_NONE, Pti,
			       Pt, GRETL_MOD_CUMULATE);
	    gretl_matrix_qform(K->T, GRETL_MOD_NONE, Pki,
			       Pk, GRETL_MOD_NONE);
	}

	if (K->LL != NULL) {
	    K->LL->val[t] = -0.5 * (llct + ldt + qt);
	    K->loglik += K->LL->val[t];
	} else {
	    K->loglik -= 0.5 * (llct + ldt + qt);
	}
	K->SSRw += qt;
	if (att != NULL) {
	    load_to_vec(att, ati, t);
	}
	if (Ptt != NULL) {
	    row_from_mat(Ptt, Pti, t);
	}
    } /* end of time loop */

    if (0) { /* Hmm, conditionality? */
	K->d = K->N;
	K->j = p;
    } else if (smo && K->exact) {
	dj_from_Finf(sm.Fk, &d, &j);
	printf("after filtering: d=%d, j=%d\n", d, j);
	if (d > 0 && d < sm.PK->rows) {
	    sm.PK = first_n_rows(sm.PK, d);
	    sm.Fk = first_n_rows(sm.Fk, d);
	    sm.Kinf = first_n_rows(sm.Kinf, d);
	}
	K->d = d > 0 ? d : 0;
	K->j = j > 0 ? j : 0;
    }

    printf("univariate: loglik = %#.8g\n", K->loglik);

    bundle_add_matrix(K->b, "P", P);
    bundle_add_matrix(K->b, "att", att);
    bundle_add_matrix(K->b, "Ptt", Ptt);
    if (gtrans) {
	bundle_add_matrix(K->b, "g", g);
    }
    if (Ztrans) {
	bundle_add_matrix(K->b, "Z", Z);
    }

    if (smo && K->exact) {
	/* donate the extra matrices needed for smoothing */
	bundle_add_matrix(K->b, "PK", sm.PK);
	bundle_add_matrix(K->b, "Fk", sm.Fk);
	bundle_add_matrix(K->b, "Kk", sm.Kinf);
    }

    /* free stuff! */
    gretl_matrix_block_destroy(B);

    return 0;
}

/* univariate smoothing */

/* cumulant matrices */

struct cumulants {
    gretl_matrix *r0;
    gretl_matrix *N0;
    gretl_matrix *L0;
    gretl_matrix *r1;
    gretl_matrix *N1;
    gretl_matrix *N2;
    gretl_matrix *mm;
};

static int allocate_cumulants (struct cumulants *c,
			       kalman *K)
{
    int m = K->r;

    c->r0 = gretl_zero_matrix_new(m, 1);
    c->N0 = gretl_zero_matrix_new(m, m);
    c->L0 = gretl_matrix_alloc(m, m);
    c->mm = gretl_matrix_alloc(m, m);
    if (K->exact) {
	c->r1 = gretl_zero_matrix_new(m, 1);
	c->N1 = gretl_zero_matrix_new(m, m);
	c->N2 = gretl_zero_matrix_new(m, m);
    }

    return 0;
}

static void clear_cumulants (struct cumulants *c,
			     kalman *K)
{
    gretl_matrix_free(c->r0);
    gretl_matrix_free(c->N0);
    gretl_matrix_free(c->L0);
    gretl_matrix_free(c->mm);
    if (K->exact) {
	gretl_matrix_free(c->r1);
	gretl_matrix_free(c->N1);
	gretl_matrix_free(c->N2);
    }
}

/* complex smoothing iteration when F_\infty > 0 */

static void fkpos (double fkinv,
		   double Fti,
		   double vti,
		   const gretl_matrix *Kki,
		   const gretl_matrix *Kti,
		   const gretl_matrix *Zti,
		   struct cumulants *c)
{
    gretl_matrix *Linf, *Lmid;
    gretl_matrix *m1, *mm2;
    int m = Zti->cols;

    mm2 = gretl_matrix_alloc(m, m);

    /* Linf = I(m) - Kki * Zti * fkinv */
    Linf = gretl_identity_matrix_new(m);
    gretl_matrix_multiply(Kki, Zti, c->mm);
    gretl_matrix_multiply_by_scalar(c->mm, fkinv);
    gretl_matrix_subtract_from(Linf, c->mm);

    /* Lmid = Kki * Fti * fkinv - Kti */
    Lmid = gretl_matrix_copy(Kki);
    gretl_matrix_multiply_by_scalar(Lmid, Fti * fkinv);
    gretl_matrix_subtract_from(Lmid, Kti);

    /* L0 = Lmid * Zti * fkinv */
    gretl_matrix_multiply(Lmid, Zti, c->L0);
    gretl_matrix_multiply_by_scalar(c->L0, fkinv);

    /* r1 = Linf' * r1 + L0' * r0 */
    m1 = gretl_matrix_copy(c->r1);
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
			      m1, GRETL_MOD_NONE,
			      c->r1, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
			      c->r0, GRETL_MOD_NONE,
			      c->r1, GRETL_MOD_CUMULATE);

    /* r1 += vti * fkinv * Zti' */
    gretl_matrix_copy_values_shaped(m1, Zti);
    gretl_matrix_multiply_by_scalar(m1, vti * fkinv);
    gretl_matrix_add_to(c->r1, m1);

    /* r0 = Linf' * r0 */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
			      c->r0, GRETL_MOD_NONE,
			      m1, GRETL_MOD_NONE);
    gretl_matrix_copy_values(c->r0, m1);

    /* N2 = Linf' * N2 * Linf */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
			      c->N2, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, Linf, c->N2);

    /* N2 -= Fti * fkinv * fkinv * Zti' * Zti */
    gretl_matrix_multiply_mod(Zti, GRETL_MOD_TRANSPOSE,
			      Zti, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(c->mm, Fti * fkinv * fkinv);
    gretl_matrix_subtract_from(c->N2, c->mm);

    /* N2 += L0' * N0 * L0 */
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
			      c->N0, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, c->L0, mm2);
    gretl_matrix_add_to(c->N2, mm2);

    /* mm = Linf' * N1 * L0 */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
			      c->N1, GRETL_MOD_NONE,
			      mm2, GRETL_MOD_NONE);
    gretl_matrix_multiply(mm2, c->L0, c->mm);

    /* N2 += mm + mm' */
    gretl_matrix_add_to(c->N2, c->mm);
    gretl_matrix_add_transpose_to(c->N2, c->mm);

    /* N1 = Linf' * N1 * Linf */
    gretl_matrix_multiply(mm2, Linf, c->N1);

    /* N1 += fkinv * Zti' * Zti */
    gretl_matrix_multiply_mod(Zti, GRETL_MOD_TRANSPOSE,
			      Zti, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(c->mm, fkinv);
    gretl_matrix_add_to(c->N1, c->mm);

    /* N1 += L0' * N0 * Linf */
    gretl_matrix_multiply(c->N0, Linf, c->mm);
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
			      c->mm, GRETL_MOD_NONE,
			      c->N1, GRETL_MOD_CUMULATE);

    /* N0 = Linf' * N0 * Linf */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
			      c->N0, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, Linf, c->N0);

    gretl_matrix_free(m1);
    gretl_matrix_free(mm2);
    gretl_matrix_free(Linf);
    gretl_matrix_free(Lmid);
}

/* relatively simple smoothing iteration when F_\infty = 0 */

static void fkzero (double ftinv,
		    double vti,
		    const gretl_matrix *Kti,
		    const gretl_matrix *Zti,
		    struct cumulants *c)
{
    gretl_matrix *m1;
    int m = Zti->cols;

    m1 = gretl_matrix_alloc(m, 1);

    /* L0 = I(m) - ftinv * Kti * Zti */
    gretl_matrix_inscribe_I(c->L0, 0, 0, m);
    gretl_matrix_multiply(Kti, Zti, c->mm);
    gretl_matrix_multiply_by_scalar(c->mm, ftinv);
    gretl_matrix_subtract_from(c->L0, c->mm);

    /* r0 = L0' * r0 + vti * ftinv * Zti' */
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
			      c->r0, GRETL_MOD_NONE,
			      m1, GRETL_MOD_NONE);
    gretl_matrix_copy_values(c->r0, m1);
    gretl_matrix_copy_values_shaped(m1, Zti);
    gretl_matrix_multiply_by_scalar(m1, vti * ftinv);
    gretl_matrix_add_to(c->r0, m1);

    if (c->r1 != NULL) {
	/* r1 = L0' * r1 */
	gretl_matrix_copy_values(m1, c->r1);
	gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
				  m1, GRETL_MOD_NONE,
				  c->r1, GRETL_MOD_NONE);
    }

    /* N0 = L0' * N0 * L0 + ftinv * Z[i,]' * Z[i,] */
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
			      c->N0, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, c->L0, c->N0);
    gretl_matrix_multiply_mod(Zti, GRETL_MOD_TRANSPOSE,
			      Zti, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(c->mm, ftinv);
    gretl_matrix_add_to(c->N0, c->mm);

    if (c->N1 != NULL && c->N2 != NULL) {
	/* N1 = N1 * L0 */
	gretl_matrix_multiply(c->N1, c->L0, c->mm);
	gretl_matrix_copy_values(c->N1, c->mm);
	/* N2 = N2 * L0 */
	gretl_matrix_multiply(c->N2, c->L0, c->mm);
	gretl_matrix_copy_values(c->N2, c->mm);
    }

    gretl_matrix_free(m1);
}

struct filter_mats {
    gretl_matrix *Z;
    gretl_matrix *A;
    gretl_matrix *P;
    gretl_matrix *g;
    gretl_matrix *FK;
    gretl_matrix *KK;
    gretl_matrix *PK;
};

/* smoothing update of Ahat and Vhat */

static void dagger_calc (struct filter_mats *f,
			 struct cumulants *c,
			 gretl_matrix *Pt,
			 gretl_matrix *Pk,
			 gretl_matrix *Ahat,
			 gretl_matrix *Vhat,
			 int t)
{
    gretl_matrix *rdag, *Ndag, *Pdag;
    gretl_matrix *at, *tmp;
    int m = c->N0->rows;
    int m2 = m * 2;
    int mm = m * m;
    size_t rsize = m * sizeof(double);
    size_t Nsize = mm * sizeof(double);

    /* workspace */
    at = gretl_matrix_alloc(m, 1);
    tmp = gretl_matrix_alloc(m, m2);
    rdag = gretl_matrix_alloc(m2, 1);
    Ndag = gretl_matrix_alloc(m2, m2);
    Pdag = gretl_matrix_alloc(m, 2*m);

    /* rdag = r0 | r1 */
    memcpy(rdag->val, c->r0->val, rsize);
    memcpy(rdag->val + m, c->r1->val, rsize);

    /* Ndag = (N0 ~ N1') | (N1 ~ N2) */
    gretl_matrix_inscribe_matrix(Ndag, c->N0, 0, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(Ndag, c->N1, 0, m, GRETL_MOD_TRANSPOSE);
    gretl_matrix_inscribe_matrix(Ndag, c->N1, m, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(Ndag, c->N2, m, m, GRETL_MOD_NONE);

    mat_from_row(at, f->A, t);
    mat_from_row(Pt, f->P, t);
    mat_from_row(Pk, f->PK, t);

    /* Pdag = Pt ~ Pk */
    memcpy(Pdag->val, Pt->val, Nsize);
    memcpy(Pdag->val + mm, Pk->val, Nsize);

    /* Ahat[t,] = (at + Pdag * rdag)' */
    gretl_matrix_multiply_mod(Pdag, GRETL_MOD_NONE,
			      rdag, GRETL_MOD_NONE,
			      at, GRETL_MOD_CUMULATE);
    row_from_mat(Ahat, at, t);

    /* Vhat[t,] = vec(Pt - Pdag * Ndag * Pdag')' */
    gretl_matrix_multiply(Pdag, Ndag, tmp);
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
			      Pdag, GRETL_MOD_TRANSPOSE,
			      Pt, GRETL_MOD_DECREMENT);
    row_from_mat(Vhat, Pt, t);

    gretl_matrix_free(at);
    gretl_matrix_free(tmp);
    gretl_matrix_free(rdag);
    gretl_matrix_free(Ndag);
    gretl_matrix_free(Pdag);
}

static void eps_smooth_Ft (double ftinv, double gi,
			   double vti,
			   const gretl_matrix *Kti,
			   struct cumulants *c,
			   gretl_matrix *epshat,
			   gretl_matrix *veps,
			   int t, int i)
{
    gretl_matrix *tmp = gretl_matrix_alloc(Kti->rows, 1);
    double x;

    x = dotprod(Kti, c->r0);
    gretl_matrix_set(epshat, t, i, gi * ftinv * (vti - x));

    gretl_matrix_multiply(c->N0, Kti, tmp);
    gretl_matrix_multiply_by_scalar(tmp, ftinv * ftinv);
    x = gi - gi*gi * (ftinv + dotprod(Kti, tmp));
    gretl_matrix_set(veps, t, i, x);

    gretl_matrix_free(tmp);
}

static void eps_smooth_Fk (double fkinv, double gi,
			   const gretl_matrix *Kki,
			   struct cumulants *c,
			   gretl_matrix *epshat,
			   gretl_matrix *veps,
			   int t, int i)
{
    gretl_matrix *tmp = gretl_matrix_alloc(Kki->rows, 1);
    double x;

    x = dotprod(Kki, c->r0);
    gretl_matrix_set(epshat, t, i, -gi * x * fkinv);

    gretl_matrix_multiply(c->N0, Kki, tmp);
    x = gi - gi*gi * fkinv * fkinv * dotprod(Kki, tmp);
    gretl_matrix_set(veps, t, i, x);

    gretl_matrix_free(tmp);
}

static void eta_smooth (const gretl_matrix *Q,
			const gretl_matrix *QRT,
			const gretl_matrix *r,
			const gretl_matrix *N,
			gretl_matrix *etahat,
			gretl_matrix *veta,
			int t)
{
    gretl_matrix *tmp = gretl_matrix_copy(Q);

    /* FIXME is qform safe here? */

    /* Q - QRT * N * QRT' */
    gretl_matrix_qform(QRT, GRETL_MOD_NONE,
		       N, tmp, GRETL_MOD_DECREMENT);
    row_from_mat(veta, tmp, t);

    if (t > 0) {
	/* QRT * r */
	gretl_matrix_reuse(tmp, -1, 1);
	gretl_matrix_multiply(QRT, r, tmp);
	row_from_mat(etahat, tmp, t-1);
    }

    gretl_matrix_free(tmp);
}

static void state_smooth (const gretl_matrix *at,
			  const gretl_matrix *Pt,
			  const gretl_matrix *r0,
			  const gretl_matrix *N0,
			  gretl_matrix *Ahat,
			  gretl_matrix *Vhat,
			  int t)
{
    gretl_matrix *tmp = gretl_matrix_copy(Pt);

    /* Pt - Pt * N0 * Pt' */
    gretl_matrix_qform(Pt, GRETL_MOD_NONE,
		       N0, tmp, GRETL_MOD_DECREMENT);
    row_from_mat(Vhat, tmp, t);

    /* at + Pt * r0 */
    gretl_matrix_reuse(tmp, -1, 1);
    gretl_matrix_multiply(Pt, r0, tmp);
    gretl_matrix_add_to(tmp, at);
    row_from_mat(Ahat, tmp, t);

    gretl_matrix_free(tmp);
}

static void regular_backdate (struct cumulants *c,
			      gretl_matrix *T)
{
    int m = c->mm->rows;

    /* N0 = T' * N0 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->N0, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, T, c->N0);

    /* r0 = T' * r0 */
    gretl_matrix_reuse(c->mm, m, 1);
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->r0, GRETL_MOD_NONE,
			      c->mm, GRETL_MOD_NONE);
    gretl_matrix_copy_values(c->r0, c->mm);
    gretl_matrix_reuse(c->mm, m, m);
}

static void diffuse_backdate (struct cumulants *c,
			      const gretl_matrix *T,
			      int m)
{
    gretl_matrix *mm = gretl_matrix_alloc(m, m);

    /* N0 = T' * N0 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->N0, GRETL_MOD_NONE,
			      mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(mm, T, c->N0);

    /* N1 = T' * N1 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->N1, GRETL_MOD_NONE,
			      mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(mm, T, c->N1);

    /* N2 = T' * N2 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->N2, GRETL_MOD_NONE,
			      mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(mm, T, c->N2);

    gretl_matrix_reuse(mm, m, 1);
    /* r0 = T' * r0 */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->r0, GRETL_MOD_NONE,
			      mm, GRETL_MOD_NONE);
    gretl_matrix_copy_values(c->r0, mm);
    /* r1 = T' * r1 */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
			      c->r1, GRETL_MOD_NONE,
			      mm, GRETL_MOD_NONE);
    gretl_matrix_copy_values(c->r1, mm);

    gretl_matrix_free(mm);
}

static void load_filter_matrices (struct filter_mats *fm,
				  kalman *K)
{
    fm->Z = gretl_bundle_get_matrix(K->b, "Z", NULL);
    fm->P = gretl_bundle_get_matrix(K->b, "P", NULL);
    fm->g = gretl_bundle_get_matrix(K->b, "g", NULL);
    if (K->exact) {
	fm->FK = gretl_bundle_get_matrix(K->b, "Fk", NULL);
	fm->KK = gretl_bundle_get_matrix(K->b, "Kk", NULL);
	fm->PK = gretl_bundle_get_matrix(K->b, "PK", NULL);
    }
    /* these two may be untransformed */
    if (fm->Z == NULL) fm->Z = K->ZT;
    if (fm->g == NULL) fm->g = K->VY;

    fm->A = K->A;
}

static int state_dist_setup (kalman *K,
			     gretl_matrix **Q,
			     gretl_matrix **QRT,
			     gretl_matrix **etahat,
			     gretl_matrix **veta)
{
    int q = K->r;

    if (K->Q != NULL) {
	q = K->Q->rows;
	*Q = K->Q;
	*QRT = K->QRT;
    } else {
	*Q = K->VS;
	*QRT = K->VS; /* FIXME? */
    }

    *etahat = gretl_zero_matrix_new(K->N, q);
    *veta = gretl_zero_matrix_new(K->N, q*q);

    return 1;
}

static int ksmooth_univariate (kalman *K, int dist)
{
    gretl_matrix_block *B = NULL;
    struct filter_mats f = {0};
    struct cumulants c = {0};
    gretl_matrix *Q = NULL;
    gretl_matrix *QRT = NULL;
    gretl_matrix *g = NULL;

    load_filter_matrices(&f, K);
    g = f.g;    /* convenience pointer */

    int i, t, N = K->N;
    int d = K->d;
    int j = K->j;
    int TI = gretl_is_identity_matrix(K->T);
    int p = K->n; /* # observables */
    int m = K->r; /* length of state vector */
    int eps_smo = 0;
    int eta_smo = 0;

    /* matrices to be returned in bundle */
    gretl_matrix *Ahat = gretl_zero_matrix_new(N, m);
    gretl_matrix *Vhat = gretl_zero_matrix_new(N, m*m);
    gretl_matrix *epshat = NULL;
    gretl_matrix *veps = NULL;
    gretl_matrix *etahat = NULL;
    gretl_matrix *veta = NULL;

    /* block-workspace matrices */
    gretl_matrix *Nt, *Pt, *Pk, *vt;
    gretl_matrix *Ft, *Kt, *at;
    gretl_matrix *Fk, *Kk;
    gretl_matrix *Zti, *Kti, *Kki;

    if (K->exact) {
	B = gretl_matrix_block_new(&Nt, m, m, &Kt, m, p,
				   &Pt, m, m, &Ft, p, 1,
				   &vt, p, 1, &at, m, 1,
				   &Fk, p, 1, &Kk, m, m,
				   &Pk, m, m, &Zti, 1, m,
				   &Kti, m, 1, &Kki, m, 1,
				   NULL);
    } else {
	B = gretl_matrix_block_new(&Nt, m, m, &Kt, m, p,
				   &Pt, m, m, &Ft, p, 1,
				   &vt, p, 1, &at, m, 1,
				   &Zti, 1, m, &Kti, m, 1,
				   NULL);
    }
    gretl_matrix_zero(Nt);

    allocate_cumulants(&c, K);

    /* workspace scalars */
    double vti, Fti, Fki;
    double ftinv, fkinv;

    if (dist) {
	if (g != NULL) {
	    /* smoothing of obs disturbances wanted */
	    epshat = gretl_zero_matrix_new(N, p);
	    veps = gretl_zero_matrix_new(N, p);
	    eps_smo = 1;
	}
	/* smoothing of state disturbances wanted */
	eta_smo = state_dist_setup(K, &Q, &QRT, &etahat, &veta);
    }

    printf("\nCCC ksmo: d = %d, j = %d CCC\n\n", d, j);

    for (t=K->N-1; t>=d; t--) {
	mat_from_row(vt, K->V, t);
	mat_from_row(Ft, K->F, t);
	mat_from_row(Kt, K->K, t);
	mat_from_row(Pt, f.P, t);
	mat_from_row(at, K->A, t);

	/* loop across observables */
	for (i=p-1; i>=0; i--) {
	    vti = vt->val[i];
	    Fti = Ft->val[i];
	    if (na(vti) || Fti == 0) {
		if (eps_smo) {
		    gretl_matrix_set(veps, t, i, g->val[i]);
		}
		continue;
	    }
	    ftinv = 1.0 / Fti;
	    mat_from_row(Zti, f.Z, i);
	    vec_from_col(Kti, Kt, i);
	    if (eps_smo) {
		eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
			      epshat, veps, t, i);
	    }
	    fkzero(ftinv, vti, Kti, Zti, &c);
	}

	if (eta_smo) {
	    eta_smooth(Q, QRT, c.r0, Nt, etahat, veta, t);
	}

	/* smoothed state and its variance */
	state_smooth(at, Pt, c.r0, c.N0, Ahat, Vhat, t);

	/* regular backdate for t-1 */
	if (t > 0) {
	    gretl_matrix_copy_values(Nt, c.N0);
	    if (!TI) {
		regular_backdate(&c, K->T);
	    }
	}
    } /* end first time loop */

    if (d > 0) {
	t = d - 1; /* first time-step in diffuse phase */
	mat_from_row(Kt, K->K, t);
	mat_from_row(Ft, K->F, t);
	mat_from_row(vt, K->V, t);

	/* countdown to observable @j */
	for (i=p-1; i>=j; i--) {
	    vti = vt->val[i];
	    Fti = Ft->val[i];
	    if (na(vti) || Fti == 0) {
		if (eps_smo) {
		    gretl_matrix_set(veps, t, i, g->val[i]);
		}
		continue;
	    }
	    ftinv = 1.0 / Fti;
	    mat_from_row(Zti, f.Z, i);
	    vec_from_col(Kti, Kt, i);
	    if (eps_smo) {
		eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
			      epshat, veps, t, i);
	    }
	    fkzero(ftinv, vti, Kti, Zti, &c);
	}

	/* retrieve the diffuse values */
	mat_from_row(Fk, f.FK, t);
	mat_from_row(Kk, f.KK, t);
	mat_from_row(Pk, f.PK, t);

	/* further countdown to first observable */
	for (i=j-1; i>=0; i--) {
	    Fki = Fk->val[i];
	    Fti = Ft->val[i];
	    vti = vt->val[i];
	    printf("t,i = %d,%d, i-countdown (2) j to 0, Fki = %g\n", t, i, Fki);
	    if (na(vti)) {
		if (eps_smo) {
		    gretl_matrix_set(veps, t, i, g->val[i]);
		}
		continue;
	    }
	    if (Fki > 0) {
		fkinv = 1.0 / Fki;
		mat_from_row(Zti, f.Z, i);
		vec_from_col(Kki, Kk, i);
		vec_from_col(Kti, Kt, i);
		if (eps_smo) {
		    eps_smooth_Fk(fkinv, g->val[i], Kki, &c,
				  epshat, veps, t, i);
		}
		fkpos(fkinv, Fti, vti, Kki, Kti, Zti, &c);
	    } else if (Fti > 0) {
		ftinv = 1.0 / Fti;
		mat_from_row(Zti, f.Z, i);
		vec_from_col(Kti, Kt, i);
		if (eps_smo) {
		    eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
				  epshat, veps, t, i);
		}
		fkzero(ftinv, vti, Kti, Zti, &c);
	    }
	}

	if (eta_smo) {
	    eta_smooth(Q, QRT, c.r0, Nt, etahat, veta, t);
	}

	printf("dagger calc, call 1 (t=%d)\n", t);
	dagger_calc(&f, &c, Pt, Pk, Ahat, Vhat, t);

	if (t > 0) {
	    printf("call diffuse_backdate\n");
	    gretl_matrix_copy_values(Nt, c.N0);
	    if (!TI) {
		diffuse_backdate(&c, K->T, m);
	    }
	}

	for (t=d-2; t>=0; t--) {
	    // printf("t,i = %d,%d, t-countdown d-2 to 0\n", t, i);
	    mat_from_row(vt, K->V, t);
	    mat_from_row(Ft, K->F, t);
	    mat_from_row(Fk, f.FK, t);
	    mat_from_row(Kt, K->K, t);
	    mat_from_row(Kk, f.KK, t);

	    for (i=p-1; i>=0; i--) {
		printf("t,i = %d,%d, i-countdown p to 0\n", t, i);
		vti = vt->val[i];
		Fti = Ft->val[i];
		Fki = Fk->val[i];
		if (na(vti)) {
		    if (eps_smo) {
			gretl_matrix_set(veps, t, i, g->val[i]);
		    }
		    continue;
		}
		if (Fki > 0) {
		    fkinv = 1.0 / Fki;
		    mat_from_row(Zti, f.Z, i);
		    vec_from_col(Kki, Kk, i);
		    vec_from_col(Kti, Kt, i);
		    if (eps_smo) {
			eps_smooth_Fk(fkinv, g->val[i], Kki, &c,
				      epshat, veps, t, i);
		    }
		    fkpos(fkinv, Fti, vti, Kki, Kti, Zti, &c);
		} else if (Fti > 0) {
		    ftinv = 1.0 / Fti;
		    ftinv = 1.0 / Fti;
		    mat_from_row(Zti, f.Z, i);
		    vec_from_col(Kti, Kt, i);
		    if (eps_smo) {
			eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
				      epshat, veps, t, i);
		    }
		    fkzero(ftinv, vti, Kti, Zti, &c);
		}
	    }

	    printf("dagger calc, call 2 (t=%d)\n", t);
	    dagger_calc(&f, &c, Pt, Pk, Ahat, Vhat, t);

	    if (eta_smo) {
		eta_smooth(Q, QRT, c.r0, Nt, etahat, veta, t);
	    }

	    if (t > 0) {
		gretl_matrix_copy_values(Nt, c.N0);
		if (!TI) {
		    diffuse_backdate(&c, K->T, m);
		}
	    }
	} /* end t=d-2..0 */
    }

    printf("\n");

    /* add smoothed quantities to @b */
    bundle_add_matrix(K->b, "Ahat", Ahat);
    bundle_add_matrix(K->b, "Vhat", Vhat);
    if (eps_smo) {
	bundle_add_matrix(K->b, "epshat", epshat);
	bundle_add_matrix(K->b, "veps",  veps);
    }
    if (eta_smo) {
	bundle_add_matrix(K->b, "etahat", etahat);
	bundle_add_matrix(K->b, "veta", veta);
    }

    gretl_matrix_block_destroy(B);
    clear_cumulants(&c, K);

    return 0;
}

/* end of univariate smoother functions */

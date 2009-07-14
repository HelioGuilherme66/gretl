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

/* libset.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "usermat.h"
#include "matrix_extra.h"
#include "gretl_string_table.h"
#include "gretl_func.h"
#include "gretl_scalar.h"

#include <unistd.h>
#include <errno.h>

#define PDEBUG 0

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE,
    AUTO_LAG_NEWEYWEST
};

enum {
    STATE_USE_CWD         = 1 << 0,  /* store: use current dir as default */
    STATE_ECHO_ON         = 1 << 1,  /* echoing commands or not */
    STATE_MSGS_ON         = 1 << 2,  /* emitting non-error messages or not */
    STATE_FORCE_DECPOINT  = 1 << 3,  /* override locale decimal separator */
    STATE_USE_PCSE        = 1 << 4,  /* Beck-Katz panel-corrected std errs */
    STATE_USE_SVD         = 1 << 5,  /* SVD decomposition is matrix OLS default */
    STATE_PREWHITEN       = 1 << 6,  /* HAC pre-whitening? */
    STATE_FORCE_HC        = 1 << 7,  /* don't use HAC for time series */
    STATE_HALT_ON_ERR     = 1 << 8,  /* errors fatal in batch mode */
    STATE_USE_LBFGS       = 1 << 9,  /* prefer LBFGS to BFGS? */
    STATE_SHELL_OK        = 1 << 10, /* "shell" facility is approved? */
    STATE_MAX_VERBOSE     = 1 << 11, /* verbose output from maximizer? */
    STATE_USE_FCP         = 1 << 12, /* use FCP garch code */
    STATE_WARN_ON         = 1 << 13, /* print numerical warning messages */
    STATE_VERBOSE_INCLUDE = 1 << 14, /* verbose include */
    STATE_SKIP_MISSING    = 1 << 15, /* skip NAs when building matrix from series */
    STATE_LOOPING         = 1 << 16, /* loop is in progress at this level */
    STATE_LOOP_QUIET      = 1 << 17  /* loop commands should be quiet */
};    

/* for values that really want a non-negative integer */
#define UNSET_INT -1
#define is_unset(i) (i == UNSET_INT)

typedef struct set_vars_ set_vars;

struct robust_opts {
    int auto_lag;
    int user_lag;
    int hc_version;
    int hkern;
    double qsband;
};

struct set_vars_ {
    int flags;
    unsigned int seed;          /* for PRNG */
    double hp_lambda;           /* for Hodrick-Prescott filter */
    int horizon;                /* for VAR impulse responses */ 
    int bootrep;                /* bootstrap replications */
    double nls_toler;           /* NLS convergence criterion */
    int loop_maxiter;           /* max no. of iterations in non-for loops */
    char delim;                 /* delimiter for CSV data export */
    int longdigits;             /* digits for printing data in long form */
    int vecm_norm;              /* VECM beta normalization */
    int bfgs_maxiter;           /* max iterations, BFGS */         
    double bfgs_toler;          /* convergence tolerance, BFGS */
    int bhhh_maxiter;           /* max iterations, BHHH */          
    double bhhh_toler;          /* convergence tolerance, BHHH */
    int garch_vcv;              /* GARCH vcv variant */
    int garch_robust_vcv;       /* GARCH vcv variant, robust estimation */
    int arma_vcv;               /* ARMA vcv variant */
    int bkbp_k;                 /* Baxter-King k */
    int bkbp_p0;                /* Baxter-King periods[0] */
    int bkbp_p1;                /* Baxter-King periods[1] */
    int rq_maxiter;             /* max iterations for quantreg, simplex */
    gretl_matrix *initvals;     /* parameter initializer */
    struct robust_opts ropts;   /* robust standard error options */
    char shelldir[MAXLEN];      /* working dir for shell commands */
};

#define ECHO "echo"
#define MESSAGES "messages"
#define WARNINGS "warnings"
#define GRETL_DEBUG "debug"
#define BLAS_NMK_MIN "blas_nmk_min"
#define PROTECT_LISTS "protect_lists"

#define libset_boolvar(s) (!strcmp(s, ECHO) || \
                           !strcmp(s, MESSAGES) || \
                           !strcmp(s, WARNINGS) || \
                           !strcmp(s, FORCE_DECP) || \
			   !strcmp(s, FORCE_HC) || \
			   !strcmp(s, HALT_ON_ERR) || \
                           !strcmp(s, MAX_VERBOSE) || \
			   !strcmp(s, USE_LBFGS) || \
			   !strcmp(s, PCSE) || \
			   !strcmp(s, PREWHITEN) || \
			   !strcmp(s, USE_SVD) || \
			   !strcmp(s, SHELL_OK) || \
			   !strcmp(s, USE_CWD) || \
			   !strcmp(s, USE_FCP) || \
                           !strcmp(s, PROTECT_LISTS) || \
                           !strcmp(s, VERBOSE_INCLUDE) || \
                           !strcmp(s, SKIP_MISSING) || \
			   !strcmp(s, R_FUNCTIONS) || \
			   !strcmp(s, R_LIB))

#define libset_double(s) (!strcmp(s, BFGS_TOLER) || \
			  !strcmp(s, BHHH_TOLER) || \
			  !strcmp(s, HP_LAMBDA) || \
			  !strcmp(s, NLS_TOLER) || \
			  !strcmp(s, QS_BANDWIDTH))

#define libset_int(s) (!strcmp(s, BFGS_MAXITER) || \
		       !strcmp(s, BHHH_MAXITER) || \
                       !strcmp(s, BKBP_K) || \
		       !strcmp(s, BOOTREP) || \
		       !strcmp(s, HAC_KERNEL) || \
                       !strcmp(s, HC_VERSION) || \
		       !strcmp(s, HORIZON) || \
		       !strcmp(s, LONGDIGITS) || \
		       !strcmp(s, LOOP_MAXITER) || \
                       !strcmp(s, RQ_MAXITER) || \
		       !strcmp(s, VECM_NORM) || \
		       !strcmp(s, GRETL_DEBUG) || \
		       !strcmp(s, BLAS_NMK_MIN))

/* global state */
set_vars *state;
static int gretl_debug;
static int protect_lists = 1;
static int user_mp_bits;
static int R_functions;
static int R_lib;

static int boolvar_get_flag (const char *s);
static const char *hac_lag_string (void);

static void robust_opts_init (struct robust_opts *r)
{
    r->auto_lag = AUTO_LAG_STOCK_WATSON;
    r->user_lag = UNSET_INT;
    r->hc_version = 0;
    r->hkern = KERNEL_BARTLETT;
    r->qsband = NADBL;
}

static void robust_opts_copy (struct robust_opts *r)
{
    r->auto_lag = state->ropts.auto_lag;
    r->user_lag = state->ropts.user_lag;
    r->hc_version = state->ropts.hc_version;
    r->hkern = state->ropts.hkern; 
    r->qsband = state->ropts.qsband;
}

static const char *csv_delim_args[] = {
    "comma",
    "space",
    "tab",
    "semicolon",
    NULL
};

static const char *garch_vcv_strs[] = {
    "unset",
    "hessian",
    "im",
    "op",
    "qml",
    "bw",
    NULL
};

static const char *arma_vcv_strs[] = {
    "hessian",
    "op",
    NULL
};

static const char *hac_kernel_strs[] = {
    "bartlett", 
    "parzen", 
    "qs",
    NULL
};

static const char *hc_version_strs[] = {
    "0", "1", "2", "3", "3a", NULL
};

static const char *vecm_norm_strs[] = {
    "phillips",
    "diag",
    "first",
    "none",
    NULL
};

static const char **libset_option_strings (const char *s)
{
    if (!strcmp(s, GARCH_VCV)) {
	return garch_vcv_strs;
    } else if (!strcmp(s, ARMA_VCV)) {
	return arma_vcv_strs;
    } else if (!strcmp(s, HAC_KERNEL)) {
	return hac_kernel_strs;
    } else if (!strcmp(s, HC_VERSION)) {
	return hc_version_strs;
    } else if (!strcmp(s, VECM_NORM)) {
	return vecm_norm_strs;
    } else if (!strcmp(s, "csv_delim")) {
	return csv_delim_args;
    } else {
	return NULL;
    }
}

static void coded_var_show_opts (const char *s, PRN *prn)
{
    const char **S = libset_option_strings(s);

    if (S != NULL) {
	pputs(prn, "valid settings:");
	while (*S != NULL) {
	    pprintf(prn, " %s", *S);
	    S++;
	}
	pputc(prn, '\n');
    }
}

static const char *get_arma_vcv_str (int v)
{
    if (v == VCV_HESSIAN) {
	return arma_vcv_strs[0];
    } else if (v == VCV_OP) {
	return arma_vcv_strs[1];
    } else {
	return "unknown";
    }
}

static const char *libset_option_string (const char *s)
{
    if (!strcmp(s, HAC_LAG)) {
	return hac_lag_string(); /* special */
    } else if (!strcmp(s, GARCH_VCV)) {
	return garch_vcv_strs[state->garch_vcv];
    } else if (!strcmp(s, ARMA_VCV)) {
	return get_arma_vcv_str(state->arma_vcv);
    } else if (!strcmp(s, HAC_KERNEL)) {
	return hac_kernel_strs[state->ropts.hkern];
    } else if (!strcmp(s, HC_VERSION)) {
	return hc_version_strs[state->ropts.hc_version];
    } else if (!strcmp(s, VECM_NORM)) {
	return vecm_norm_strs[state->vecm_norm];
    } else {
	return "?";
    }
}

static void print_initvals (const gretl_matrix *ivals, PRN *prn)
{
    if (ivals == NULL) {
	pputs(prn, " initvals = auto\n");
    } else {
	gretl_matrix_print_to_prn(ivals, " initvals =", prn);
    }
}

/* check_for_state() returns non-zero if the program options
   state is not readable */

static int check_for_state (void) 
{
    if (state == NULL) {
	return libset_init();
    } else {
#if PDEBUG > 1
	fprintf(stderr, "check_for_state: state = %p\n", (void *) state);
#endif
	return 0;
    }
}

static int flag_to_bool (set_vars *sv, int flag)
{
    if (!sv) {
	return 0;
    } else {
	return (sv->flags & flag)? 1 : 0;
    }
}

static void state_vars_copy (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_copy() called\n");
#endif
    sv->flags = state->flags;
    /* We're not (yet) looping at the current level of execution (but
       note that the STATE_LOOP_QUIET flag should be inherited).
    */
    sv->flags &= ~STATE_LOOPING;

    sv->seed = state->seed;
    sv->hp_lambda = state->hp_lambda;
    sv->horizon = state->horizon;
    sv->bootrep = state->bootrep;
    sv->loop_maxiter = state->loop_maxiter;
    sv->rq_maxiter = state->rq_maxiter;
    sv->nls_toler = state->nls_toler;
    sv->delim = state->delim; 
    sv->longdigits = state->longdigits; 
    sv->vecm_norm = state->vecm_norm;
    sv->bfgs_maxiter = state->bfgs_maxiter;
    sv->bfgs_toler = state->bfgs_toler;
    sv->bhhh_maxiter = state->bhhh_maxiter;
    sv->bhhh_toler = state->bhhh_toler;
    sv->garch_vcv = state->garch_vcv;
    sv->garch_robust_vcv = state->garch_robust_vcv;
    sv->bkbp_k = state->bkbp_k;
    sv->bkbp_p0 = state->bkbp_p0;
    sv->bkbp_p1 = state->bkbp_p1;

    sv->initvals = gretl_matrix_copy(state->initvals);
    strcpy(sv->shelldir, state->shelldir);

    robust_opts_copy(&sv->ropts);
}

static void state_vars_init (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_init called\n");
#endif
    sv->flags = STATE_ECHO_ON | STATE_MSGS_ON | STATE_WARN_ON | 
	STATE_HALT_ON_ERR | STATE_SKIP_MISSING;
    sv->seed = 0;
    sv->hp_lambda = NADBL;
    sv->horizon = UNSET_INT;
    sv->bootrep = 1000;
    sv->nls_toler = NADBL;
    sv->loop_maxiter = 250;
    sv->rq_maxiter = 1000;
    sv->delim = UNSET_INT;
    sv->longdigits = 10;
    sv->vecm_norm = NORM_PHILLIPS;
    sv->initvals = NULL;

    sv->bfgs_maxiter = -1;
    sv->bfgs_toler = NADBL;
    sv->bhhh_maxiter = 500;
    sv->bhhh_toler = NADBL;
    sv->garch_vcv = VCV_UNSET;
    sv->arma_vcv = VCV_HESSIAN;
    sv->garch_robust_vcv = VCV_UNSET;

    sv->bkbp_k = UNSET_INT;
    sv->bkbp_p0 = UNSET_INT;
    sv->bkbp_p1 = UNSET_INT;

    *sv->shelldir = '\0';

    robust_opts_init(&sv->ropts);
}

int get_bkbp_k (const DATAINFO *pdinfo)
{
    if (check_for_state()) {
	return 0;
    }

    if (is_unset(state->bkbp_k)) {
	if (pdinfo->pd == 1) {
	    return 3;
	} else if (pdinfo->pd == 4) {
	    return 12;
	} else if (pdinfo->pd == 12) {
	    return 36;
	} else {
	    return 3;
	}
    } else {
	return state->bkbp_k;
    }
}

int set_bkbp_k (int k)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (k > 0) {
	state->bkbp_k = k;
	return 0;
    } else {
	return 1;
    }
}

void unset_bkbp_k (void)
{
    if (check_for_state()) {
	return;
    }

    state->bkbp_k = UNSET_INT;
}

void get_bkbp_periods (const DATAINFO *pdinfo, int *l, int *u)
{
    if (check_for_state()) {
	return;
    }

    if (is_unset(state->bkbp_p0)) {
	*l = (pdinfo->pd == 4)? 6 :
	    (pdinfo->pd == 12)? 18 : 2;
    } else {
	*l = state->bkbp_p0;
    }

    if (is_unset(state->bkbp_p1)) {
	*u = (pdinfo->pd == 4)? 32 :
	    (pdinfo->pd == 12)? 96 : 8;
    } else {
	*u = state->bkbp_p1;
    }
}

int set_bkbp_periods (int l, int u)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (l > 0 && u > l) {
	state->bkbp_p0 = l;
	state->bkbp_p1 = u;
	return 0;
    } else {
	return 1;
    }
}

void unset_bkbp_periods (void)
{
    if (check_for_state()) {
	return;
    }

    state->bkbp_p0 = UNSET_INT;
    state->bkbp_p1 = UNSET_INT;
}

void set_gretl_echo (int e)
{
    if (check_for_state()) return;

    if (e) {
	state->flags |= STATE_ECHO_ON;
    } else {
	state->flags &= ~STATE_ECHO_ON;
    }
}

int gretl_echo_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, STATE_ECHO_ON);
}

void set_gretl_messages (int e)
{
    if (check_for_state()) return;

    if (e) {
	state->flags |= STATE_MSGS_ON;
    } else {
	state->flags &= ~STATE_MSGS_ON;
    }
}

int gretl_messages_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, STATE_MSGS_ON);
}

int gretl_warnings_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, STATE_WARN_ON);
}

int gretl_debugging_on (void)
{
    return gretl_debug;
}

int lists_protected (void)
{
    return protect_lists;
}

#define DEFAULT_MP_BITS 256
#define mp_bits_ok(b) (b >= 256 && b <= 8192)

void set_mp_bits (int b)
{
    if (mp_bits_ok(b)) {
	user_mp_bits = b;
    }
}

int get_mp_bits (void)
{
    if (user_mp_bits >= 256) {
	return user_mp_bits;
    } else {
	char *s = getenv("GRETL_MP_BITS");
	int b;

	if (s != NULL) {
	    b = atoi(s);
	    if (mp_bits_ok(b)) {
		return b;
	    }
	}
	return DEFAULT_MP_BITS;
    }
}

char get_csv_delim (const DATAINFO *pdinfo)
{
    check_for_state();
    if (state->delim > 0) {
	return state->delim;
    } else {
	return pdinfo->delim;
    }
}

const gretl_matrix *get_init_vals (void)
{
    check_for_state();
    return state->initvals;
}

void free_init_vals (void)
{
    if (state->initvals != NULL) {
	gretl_matrix_free(state->initvals);
	state->initvals = NULL;
    }
}

int n_init_vals (void)
{
    check_for_state();
    if (state->initvals != NULL) {
	return gretl_vector_get_length(state->initvals);
    } else {
	return 0;
    }
}

char *get_shelldir (void)
{
    check_for_state();

    if (state != NULL && *state->shelldir != '\0') {
	return state->shelldir;
    } else {
	return NULL;
    }
} 

int get_hac_lag (int T)
{
    int h = 0;

    check_for_state();

    /* Variants of Newey-West */

    if (state->ropts.user_lag >= 0 && state->ropts.user_lag < T - 2) {
	/* FIXME upper limit? */
	h = state->ropts.user_lag;
    } else if (state->ropts.auto_lag == AUTO_LAG_WOOLDRIDGE) {
	h = 4.0 * pow(T / 100.0, 2.0 / 9.0);
    } else {
	/* Stock-Watson default */
	h = 0.75 * pow(T, 1.0 / 3.0);
    }

    return h;
}

/* prewhitening implies nw3, but not vice versa */

int data_based_hac_bandwidth (void)
{
    if (is_unset(state->ropts.user_lag)) {
	if (state->ropts.auto_lag == AUTO_LAG_NEWEYWEST ||
	    (state->flags & STATE_PREWHITEN)) {
	    return 1;
	}
    }

    return 0;
}

static const char *hac_lag_string (void)
{
    check_for_state();

    if (state->ropts.user_lag >= 0 && state->ropts.user_lag < 1000) {
	static char lagstr[6];

	sprintf(lagstr, "%d", state->ropts.user_lag);
	return lagstr;
    } else if (state->ropts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return "nw1";
    } else {
	return "nw2";
    }
}

/* set max lag for HAC estimation */

static int parse_hac_lag_variant (const char *s)
{
    int err = E_DATA;

    if (!strcmp(s, "nw1")) {
	state->ropts.auto_lag = AUTO_LAG_STOCK_WATSON;
	state->ropts.user_lag = UNSET_INT;
	err = 0;
    } else if (!strcmp(s, "nw2")) {
	state->ropts.auto_lag = AUTO_LAG_WOOLDRIDGE;
	state->ropts.user_lag = UNSET_INT;
	err = 0;
    } else if (!strcmp(s, "nw3") ||
	       !strcmp(s, "auto")) {
	state->ropts.auto_lag = AUTO_LAG_NEWEYWEST;
	state->ropts.user_lag = UNSET_INT;
	err = 0;
    } else if (isdigit(*s)) {
	state->ropts.user_lag = atoi(s);
	err = 0;
    }

    return err;
}

static int 
libset_numeric_string (const char *s, int *pi, double *px, int *err)
{
    char *test;
    int ret = 1;

    if (s == NULL || *s == '\0' ||
	!strcmp(s, "inf") || !strcmp(s, "nan")) {
	return 0;
    }

    errno = 0;

    gretl_push_c_numeric_locale();

    if (px != NULL) {
	*px = strtod(s, &test);
	if (*test != '\0') {
	    ret = 0;
	} else if (errno == ERANGE) {
	    gretl_errmsg_set_from_errno(s);
	    *err = 1;
	}
    } else {
	long li = strtol(s, &test, 10);

	if (*test != '\0') {
	    ret = 0;
	} else if (errno == ERANGE) {
	    gretl_errmsg_set_from_errno(s);
	    *err = 1;
	} else {
	    *pi = (int) li;
	}
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

static int negval_invalid (const char *var)
{
    return (var == NULL || strcmp(var, BLAS_NMK_MIN)); 
}

static int libset_get_scalar (const char *var, const char *arg, 
			      int *pi, double *px)
{
    double x = NADBL;
    int err = 0;

    if (libset_numeric_string(arg, pi, px, &err)) {
	if (err) {
	    err = E_DATA;
	} else if (pi != NULL && negval_invalid(var) && *pi < 0) {
	    err = E_DATA;
	} else if (px != NULL && *px < 0.0) {
	    err = E_DATA;
	}
	return err;
    }

    if (gretl_is_scalar(arg)) {
	x = gretl_scalar_get_value(arg);
    } else {
	sprintf(gretl_errmsg, "'%s': not a scalar", arg);
	return E_UNKVAR;
    }

    if (negval_invalid(var) && x < 0.0) {
	return E_DATA;
    }

    if (px != NULL) {
	*px = x;
    } else if (pi != NULL) {
	if (na(x) || fabs(x) > (double) INT_MAX) {
	    err = E_DATA;
	} else {
	    *pi = (int) x;
	}
    }

    return err;
}

static int parse_hc_variant (const char *s)
{
    int i;

    check_for_state();

    if (!strncmp(s, "hc", 2)) {
	s += 2;
    }

    for (i=0; hc_version_strs[i] != NULL; i++) {
	if (!strcmp(s, hc_version_strs[i])) {
	    state->ropts.hc_version = i;
	    return 0;
	}
    }

    if (!strcmp(s, "4")) {
	state->ropts.hc_version = 4;
	return 0;
    }

    return 1;
}

static int parse_libset_int_code (const char *key, 
				  const char *val)
{
    int i, err = E_DATA;

    if (!strcmp(key, HC_VERSION)) {
	err = parse_hc_variant(val);
    } else if (!strcmp(key, HAC_LAG)) {
	err = parse_hac_lag_variant(val);
    } else if (!strcmp(key, GARCH_VCV)) {
	for (i=0; i<VCV_MAX; i++) {
	    if (!strcmp(val, garch_vcv_strs[i])) {
		state->garch_vcv = i;
		err = 0;
		break;
	    }
	}
    } else if (!strcmp(key, ARMA_VCV)) {
	if (!strcmp(val, "op")) {
	    state->arma_vcv = VCV_OP;
	    err = 0;
	} else if (!strcmp(val, "hessian")) {
	    state->arma_vcv = VCV_HESSIAN;
	    err = 0;
	}
    } else if (!strcmp(key, HAC_KERNEL)) {
	for (i=0; i<KERNEL_MAX; i++) {
	    if (!strcmp(val, hac_kernel_strs[i])) {
		state->ropts.hkern = i;
		err = 0;
		break;
	    }
	}
    } else if (!strcmp(key, VECM_NORM)) {
	for (i=0; i<NORM_MAX; i++) {
	    if (!strcmp(val, vecm_norm_strs[i])) {
		state->vecm_norm = i;
		err = 0;
		break;
	    }
	}
    }

    if (err) {
	gretl_errmsg_sprintf("%s: invalid value '%s'\n", key, val);
    }

    return err;
}	

void set_xsect_hccme (const char *s)
{
    char *scpy;

    if (check_for_state()) return;

    scpy = gretl_strdup(s);

    if (scpy != NULL) {
	lower(scpy);
	parse_hc_variant(scpy);
	free(scpy);
    }
}

void set_tseries_hccme (const char *s)
{
    char *scpy;

    if (check_for_state()) return;

    scpy = gretl_strdup(s);

    if (scpy != NULL) {
	lower(scpy);
	if (parse_hc_variant(scpy) == 0) {
	    libset_set_bool(FORCE_HC, 1);
	} else {
	    libset_set_bool(FORCE_HC, 0);
	}
	free(scpy);
    }
}

void set_panel_hccme (const char *s)
{
    if (check_for_state()) return;

    if (!strcmp(s, "Arellano")) {
	state->flags &= ~STATE_USE_PCSE;
    } else if (!strcmp(s, "PCSE")) {
	state->flags |= STATE_USE_PCSE;
    }
}

void set_garch_robust_vcv (const char *s)
{
    char *scpy;

    if (check_for_state()) return;

    scpy = gretl_strdup(s);

    if (scpy != NULL) {
	lower(scpy);
	if (!strcmp(s, "qml")) {
	    state->garch_robust_vcv = VCV_QML;
	} else if (!strcmp(s, "bw")) {
	    state->garch_robust_vcv = VCV_BW;
	}
	free(scpy);
    }
}

static int set_line_width (const char *s0, const char *s1,
			   DATAINFO *pdinfo, PRN *prn)
{
    int v, w, err = 0;

    if (!isdigit((unsigned char) *s1)) {
	return 1;
    }

    if (isdigit((unsigned char) *s0)) {
	v = atoi(s0);
    } else {
	v = current_series_index(pdinfo, s0);
    }

    if (v < 0) {
	return E_DATA;
    }

    w = atoi(s1);

    if (w < 0 || w > 32) {
	err = E_DATA;
    } else {
	var_set_linewidth(pdinfo, v, w);
	pprintf(prn, _("Line width for %s = %d\n"), 
		pdinfo->varname[v], w);
    }

    return err;
}

static int set_bkbp_limits (const char *s0, const char *s1,
			    PRN *prn)
{
    int p0, p1;
    int err = 0;

    err = libset_get_scalar(NULL, s0, &p0, NULL);
    if (!err) {
	err = libset_get_scalar(NULL, s1, &p1, NULL);
    }

    if (err) {
	return err;
    }

    if (p1 < p0) {
	/* 2nd entry should be bigger than 1st */
	int tmp = p1;

	p1 = p0;
	p0 = tmp;
    }

    pprintf(prn, _("Baxter-King band = %d-%d periods\n"), p0, p1);

    state->bkbp_p0 = p0;
    state->bkbp_p1 = p1;

    return 0;
}

static int set_initvals (const char *s, PRN *prn)
{
    gretl_matrix *m;
    char mname[VNAMELEN];
    int err = 0;

    /* skip past "set initvals" */
    s += 12;

    if (sscanf(s, "%15s", mname) != 1 || !strcmp(mname, "auto")) {
	gretl_matrix_free(state->initvals);
	state->initvals = NULL;
    } else {
	m = get_matrix_by_name(mname);
	if (m == NULL) {
	    pprintf(prn, _("'%s': no such matrix"), mname);
	    pputc(prn, '\n');
	    err = E_DATA;
	} else {
	    state->initvals = gretl_matrix_copy(m);
	    if (state->initvals == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

void shelldir_init (const char *s)
{
    if (s != NULL) {
	int n;

	*state->shelldir = '\0';
	strncat(state->shelldir, s, MAXLEN - 1);
	n = strlen(state->shelldir);
	if (n > 0 && (state->shelldir[n-1] == '\\' ||
		      state->shelldir[n-1] == '/')) {
	    state->shelldir[n-1] = '\0';
	}	
    } else {
	char *test = getcwd(state->shelldir, MAXLEN);

	if (test == NULL) {
	    *state->shelldir = '\0';
	} 
    }

    gretl_insert_builtin_string("shelldir", state->shelldir);
}

static int set_shelldir (const char *s)
{
    int len = 0, err = 0;

    /* skip past "set shelldir" and space */
    s += 12;
    s += strspn(s, " ");

    if (*s == '\0') {
	*state->shelldir = '\0';
	gretl_insert_builtin_string("shelldir", state->shelldir);
    } else if (*s == '"') {
	s++;
	len = haschar('"', s);
	if (len <= 0) {
	    err = E_PARSE;
	} 
    } else {
	len = strlen(s);
    }

    if (!err && len > 0) {
	char test[MAXLEN];
	char *home = NULL;
	int slen = len;

	if (*s == '~') {
	    home = getenv("HOME");
	    if (home != NULL) {
		s++;
		slen--;
		len = slen + strlen(home);
	    }
	} 

	*test = '\0';
    
	if (len >= MAXLEN) {
	    gretl_errmsg_set("shelldir: string is too long");
	    err = E_DATA;
	} else if (home != NULL) {
	    strcat(test, home);
	    strncat(test, s, slen);
	} else {
	    strncat(test, s, len);
	}

	if (!gretl_isdir(test)) {
	    gretl_errmsg_sprintf("shelldir: '%s' no such directory", test);
	    err = E_DATA;
	}

	if (!err) {
	    strcpy(state->shelldir, test);
	    gretl_insert_builtin_string("shelldir", state->shelldir);
	}
    }

    return err;
}

static int (*workdir_callback)();

void set_workdir_callback (int (*callback)())
{
    workdir_callback = callback;
}

static int set_workdir (const char *s)
{
    if (gretl_function_depth() > 0) {
	gretl_errmsg_set("set workdir: cannot be done inside a function");
	return 1;
    }

    if (workdir_callback == NULL) {
	return E_DATA;
    }

    /* skip past "set workdir" and space */
    s += 11;
    s += strspn(s, " ");

    if (*s != '\0') {
	char workdir[MAXLEN];
	int n = 0;

	if (*s == '"') {
	    n = sscanf(s+1, "%511[^\"]", workdir);
	} else {
	    n = sscanf(s, "%511s", workdir);
	}
	return (*workdir_callback)(workdir);
    } else {
	return E_DATA;
    }
}

static int parse_set_plotfile (const char *s)
{
    char *fname;
    int err = 0;

    while (isspace((unsigned char) *s)) {
	s++;
    }

    /* now skip two words, "set" and "plotfile" */
    s += strcspn(s, " ");
    s += strspn(s, " ");
    s += strcspn(s, " ");
    s += strspn(s, " ");

    fname = gretl_strdup(s);
    if (fname != NULL) {
	tailstrip(fname);
	set_gretl_plotfile(fname);
	free(fname);
    } else {
	err = E_ALLOC;
    }

    return err;
}

const char *csv_delims = ", \t;";

static char delim_from_arg (const char *s)
{
    int i;

    for (i=0; csv_delim_args[i] != NULL; i++) {
	if (!strcmp(s, csv_delim_args[i])) {
	    return csv_delims[i];
	}
    }

    return 0;
}

static const char *arg_from_delim (char c)
{
    int i;

    for (i=0; csv_delims[i] != '\0'; i++) {
	if (c == csv_delims[i]) {
	    return csv_delim_args[i];
	}
    }

    return "unset";
}

static void libset_print_bool (const char *s, PRN *prn)
{
    pprintf(prn, " %s = %d\n", s, libset_get_bool(s));
}

#define coded_intvar(s) (!strcmp(s, GARCH_VCV) || \
			 !strcmp(s, ARMA_VCV) || \
			 !strcmp(s, HAC_LAG) || \
			 !strcmp(s, HAC_KERNEL) || \
                         !strcmp(s, HC_VERSION) || \
			 !strcmp(s, VECM_NORM))

const char *intvar_code_string (const char *s)
{
    if (!strcmp(s, HAC_LAG)) {
	return hac_lag_string(); /* special */
    } else {
	return libset_option_string(s);
    }
}

static void libset_print_int (const char *s, PRN *prn)
{
    if (coded_intvar(s)) {
	pprintf(prn, " %s = %s\n", s, intvar_code_string(s));
    } else {
	int k = libset_get_int(s);

	if (is_unset(k)) {
	    pprintf(prn, " %s = auto\n", s);
	} else {
	    pprintf(prn, " %s = %d\n", s, k);
	}
    }
}

static void libset_print_double (const char *s, PRN *prn)
{
    double x = libset_get_double(s);

    if (na(x)) {
	pprintf(prn, " %s = auto\n", s);
    } else {
	pprintf(prn, " %s = %g\n", s, x);
    }
}

static void libset_header (char *s, PRN *prn) 
{
    pputs(prn, "\n --- ");
    pputs(prn, s);
    pputs(prn, " ---\n");
}

static int display_settings (PRN *prn)
{
    pputs(prn, _("Variables that can be set using \"set\""));
    pputs(prn, " (");
    pputs(prn, _("\"help set\" for details"));
    pputs(prn, "):\n");

    libset_header(_("Program interaction and behavior"), prn);

    pprintf(prn, " csv_delim = %s\n", arg_from_delim(state->delim));

    libset_print_bool(ECHO, prn);
    libset_print_bool(FORCE_DECP, prn);
    libset_print_bool(HALT_ON_ERR, prn);
    libset_print_int(LONGDIGITS, prn);
    libset_print_int(LOOP_MAXITER, prn);
    libset_print_bool(MAX_VERBOSE, prn);
    libset_print_bool(MESSAGES, prn);
    libset_print_bool(WARNINGS, prn);
    libset_print_int(GRETL_DEBUG, prn);
    libset_print_int(BLAS_NMK_MIN, prn);
    libset_print_bool(SHELL_OK, prn);

    if (*state->shelldir) {
	pprintf(prn, " shelldir = '%s'\n", state->shelldir);
    } else {
	pputs(prn, " shelldir = unset\n");
    }

    libset_print_bool(USE_CWD, prn);
    libset_print_bool(PROTECT_LISTS, prn);
    libset_print_bool(VERBOSE_INCLUDE, prn);
    libset_print_bool(SKIP_MISSING, prn);

    libset_header(_("Numerical methods"), prn);

    libset_print_int(BFGS_MAXITER, prn);
    libset_print_double(BFGS_TOLER, prn);
    libset_print_int(BHHH_MAXITER, prn);
    libset_print_double(BHHH_TOLER, prn);
    libset_print_int(RQ_MAXITER, prn);
    print_initvals(state->initvals, prn);
    libset_print_bool(USE_LBFGS, prn);
    libset_print_double(NLS_TOLER, prn);
    libset_print_bool(USE_SVD, prn);
    libset_print_bool(USE_FCP, prn);

    libset_header(_("Random number generation"), prn);

    pprintf(prn, " seed = %u\n", gretl_rand_get_seed());

    libset_header(_("Robust estimation"), prn);

    libset_print_int(BOOTREP, prn);
    libset_print_int(GARCH_VCV, prn);
    libset_print_int(ARMA_VCV, prn);
    libset_print_bool(FORCE_HC, prn);
    libset_print_int(HAC_LAG, prn);
    libset_print_int(HAC_KERNEL, prn);
    libset_print_bool(PREWHITEN, prn);
    libset_print_int(HC_VERSION, prn);
    libset_print_bool(PCSE, prn);
    libset_print_double(QS_BANDWIDTH, prn);

    libset_header(_("Filtering"), prn);

    if (is_unset(state->bkbp_p0) ||
	is_unset(state->bkbp_p1)) {
	pputs(prn, " bkbp_limits = auto\n");
    } else {
	pprintf(prn, " bkbp_limits = (%d, %d)\n", 
		state->bkbp_p0, 
		state->bkbp_p1);
    }

    libset_print_int(BKBP_K, prn);
    libset_print_double(HP_LAMBDA, prn);

    libset_header(_("Time series"), prn);

    libset_print_int(HORIZON, prn);
    libset_print_int(VECM_NORM, prn);

    pputc(prn, '\n');
    
    return 0;
}

static int 
libset_query_settings (const char *s, PRN *prn)
{
    int err = 0;

    if (libset_boolvar(s)) {
	pprintf(prn, "%s: boolean (on/off), currently %s\n", 
		s, libset_get_bool(s)? "on" : "off");
    } else if (coded_intvar(s)) {
	pprintf(prn, "%s: code, currently \"%s\"\n", s, intvar_code_string(s));
	coded_var_show_opts(s, prn);
    } else if (libset_int(s)) {
	int k = libset_get_int(s);

	if (is_unset(k)) {
	    pprintf(prn, "%s: positive integer, currently unset\n", s);
	} else {
	    pprintf(prn, "%s: positive integer, currently %d\n", s, k);
	}	    
    } else if (libset_double(s)) {
	double x = libset_get_double(s);

	if (na(x)) {
	    pprintf(prn, "%s: positive floating-point value, "
		    "currently automatic\n", s);
	} else {
	    pprintf(prn, "%s: positive floating-point value, "
		    "currently %g\n", s, x);
	}
    } else if (!strcmp(s, "initvals")) {
	if (state->initvals != NULL) {
	    pprintf(prn, "%s: matrix, currently\n", s);
	    gretl_matrix_print_to_prn(state->initvals, NULL, prn);
	} else {
	    pprintf(prn, "%s: matrix, currently null\n", s);
	}
    } else if (!strcmp(s, "seed")) {
	pprintf(prn, "%s: unsigned int, currently %u\n",
		s, state->seed);
    } else if (!strcmp(s, "csv_delim")) {
	pprintf(prn, "%s: named character, currently \"%s\"\n", s,
		arg_from_delim(state->delim));
	coded_var_show_opts(s, prn);
    } else if (!strcmp(s, "shelldir")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		state->shelldir);
    } else {
	err = 1;
    }

    return err;
}

#define default_ok(s) (!strcmp(s, BFGS_TOLER) || \
                       !strcmp(s, BHHH_TOLER) || \
                       !strcmp(s, HP_LAMBDA))

#define default_str(s) (!strcmp(s, "auto") || !strcmp(s, "default"))

#define boolean_on(s) (!strcmp(s, "on") || !strcmp(s, "1") || \
                       !strcmp(s, "true"))

#define boolean_off(s) (!strcmp(s, "off") || !strcmp(s, "0") || \
                        !strcmp(s, "false"))

int execute_set_line (const char *line, DATAINFO *pdinfo, PRN *prn)
{
    char setobj[32], setarg[32], setarg2[32];
    int k, nw, err = E_PARSE;
    double x;

    check_for_state();

    *setobj = *setarg = *setarg2 = '\0';

    nw = sscanf(line, "%*s %31s %31s %31s", setobj, setarg, setarg2);

    if (nw <= 0) {
	return display_settings(prn);
    }

    /* specials which need the whole line */
    if (nw > 1) {
	if (!strcmp(setobj, "plotfile")) {
	    return parse_set_plotfile(line);
	} else if (!strcmp(setobj, "initvals")) {
	    return set_initvals(line, prn);
	} else if (!strcmp(setobj, "shelldir")) {
	    return set_shelldir(line);
	} else if (!strcmp(setobj, "codevars")) {
	    return set_codevars(line);
	} else if (!strcmp(setobj, "workdir")) {
	    return set_workdir(line);
	}
    }

    if (nw == 1) {
	if (!strcmp(setobj, ECHO)) {
	    state->flags |= STATE_ECHO_ON;
	    err = 0;
	} else if (!strcmp(setobj, "stopwatch")) {
	    gretl_stopwatch();
	    err = 0;
	} else {
	    return libset_query_settings(setobj, prn);
	}
    } else if (nw == 2) {
	lower(setarg);

	if (libset_boolvar(setobj)) {
	    if (!strcmp(setobj, SHELL_OK)) {
		pprintf(prn, "You can only set this variable "
			"via the gretl GUI\n");
	    } else if (boolean_on(setarg)) {
		err = libset_set_bool(setobj, 1);
	    } else if (boolean_off(setarg)) {
		err = libset_set_bool(setobj, 0);
	    }
	} else if (libset_double(setobj)) {
	    if (default_ok(setobj) && default_str(setarg)) {
		libset_set_double(setobj, NADBL);
		err = 0;
	    } else {
		err = libset_get_scalar(NULL, setarg, NULL, &x);
		if (!err) {
		    err = libset_set_double(setobj, x);
		}
	    }
	} else if (!strcmp(setobj, "csv_delim")) {
	    char c = delim_from_arg(setarg);

	    if (c > 0) {
		state->delim = c;
		err = 0;
	    }
	} else if (!strcmp(setobj, "seed")) {
	    err = libset_get_scalar(NULL, setarg, &k, NULL);
	    if (!err) {
		gretl_rand_set_seed((unsigned int) k);
		if (gretl_messages_on() && !gretl_looping_quietly()) {
		    pprintf(prn, 
			    _("Pseudo-random number generator seeded with %d\n"), k);
		}
		state->seed = k;
	    }
	} else if (!strcmp(setobj, HORIZON)) {
	    /* horizon for VAR impulse responses */
	    if (!strcmp(setarg, "auto")) {
		state->horizon = UNSET_INT;
		err = 0;
	    } else {
		err = libset_get_scalar(NULL, setarg, &k, NULL);
		if (!err) {
		    state->horizon = k;
		} else {
		    state->horizon = UNSET_INT;
		}
	    }
	} else if (coded_intvar(setobj)) {
	    err = parse_libset_int_code(setobj, setarg);
	} else if (libset_int(setobj)) {
	    err = libset_get_scalar(setobj, setarg, &k, NULL);
	    if (!err) {
		err = libset_set_int(setobj, k);
	    }
	} else {
	    err = E_UNKVAR;
	}
    } else if (nw == 3) {
	if (!strcmp(setobj, "bkbp_limits")) {
	    err = set_bkbp_limits(setarg, setarg2, prn);
	} else if (!strcmp(setobj, "linewidth")) {
	    err = set_line_width(setarg, setarg2, pdinfo, prn);
	} else {
	    err = E_UNKVAR;
	}
    }
		    
    return err;
}

double libset_get_double (const char *key)
{
    if (check_for_state()) {
	return 1;
    }

    if (!strcmp(key, QS_BANDWIDTH)) {
	if (!na(state->ropts.qsband) && state->ropts.qsband > 0) {
	    return state->ropts.qsband;
	} else {
	    /* what's a sensible default here? */
	    return 2.0;
	}
    } else if (!strcmp(key, NLS_TOLER)) {
	if (na(state->nls_toler)) {
	    state->nls_toler = get_default_nls_toler();
	}
	return state->nls_toler;
    } else if (!strcmp(key, BHHH_TOLER)) {
	return state->bhhh_toler;
    } else if (!strcmp(key, BFGS_TOLER)) {
	if (na(state->bfgs_toler)) {
	    state->bfgs_toler = get_default_nls_toler();
	}
	return state->bfgs_toler;
    } else if (!strcmp(key, HP_LAMBDA)) {
	return state->hp_lambda;
    } else {
	fprintf(stderr, "libset_get_double: unrecognized "
		"variable '%s'\n", key);	
	return 0;
    }
}

int libset_set_double (const char *key, double val)
{
    int err = 0;

    if (check_for_state()) {
	return 1;
    }

    /* all the libset double vals must be positive */
    if (val <= 0.0) {
	return E_DATA;
    }

    if (!strcmp(key, QS_BANDWIDTH)) {
	state->ropts.qsband = val;
    } else if (!strcmp(key, NLS_TOLER)) {
	state->nls_toler = val;
    } else if (!strcmp(key, BHHH_TOLER)) {
	state->bhhh_toler = val;
    } else if (!strcmp(key, BFGS_TOLER)) {
	state->bfgs_toler = val;
    } else if (!strcmp(key, HP_LAMBDA)) {
	state->hp_lambda = val;
    } else {
	fprintf(stderr, "libset_set_double: unrecognized "
		"variable '%s'\n", key);	
	err = E_UNKVAR;
    }

    return err;
}

int libset_get_int (const char *key)
{
    if (check_for_state()) {
	return 0;
    }

    if (!strcmp(key, BFGS_MAXITER)) {
	return state->bfgs_maxiter;
    } else if (!strcmp(key, BHHH_MAXITER)) {
	return state->bhhh_maxiter;
    } else if (!strcmp(key, RQ_MAXITER)) {
	return state->rq_maxiter;
    } else if (!strcmp(key, BKBP_K)) {
	return state->bkbp_k;
    } else if (!strcmp(key, BOOTREP)) {
	return state->bootrep;
    } else if (!strcmp(key, GARCH_VCV)) {
	return state->garch_vcv;
    } else if (!strcmp(key, GARCH_ROBUST_VCV)) {
	return state->garch_robust_vcv;
    } else if (!strcmp(key, ARMA_VCV)) {
	return state->arma_vcv;
    } else if (!strcmp(key, HAC_KERNEL)) {
	return state->ropts.hkern;
    } else if (!strcmp(key, HC_VERSION)) {
	return state->ropts.hc_version;
    } else if (!strcmp(key, HORIZON)) {
	return state->horizon;
    } else if (!strcmp(key, LONGDIGITS)) {
	return state->longdigits;
    } else if (!strcmp(key, LOOP_MAXITER)) {
	return state->loop_maxiter;
    } else if (!strcmp(key, VECM_NORM)) {
	return state->vecm_norm;
    } else if (!strcmp(key, GRETL_DEBUG)) {
	return gretl_debug;
    } else if (!strcmp(key, BLAS_NMK_MIN)) {
	return get_blas_nmk_min();
    } else {
	fprintf(stderr, "libset_get_int: unrecognized "
		"variable '%s'\n", key);	
	return 0;
    }
}

static int intvar_min_max (const char *s, int *min, int *max,
			   int **var)
{
    *max = 100000;

    if (!strcmp(s, BFGS_MAXITER)) {
	*min = 1;
	*var = &state->bfgs_maxiter;
    } else if (!strcmp(s, BHHH_MAXITER)) {
	*min = 1;
	*var = &state->bhhh_maxiter;
    } else if (!strcmp(s, RQ_MAXITER)) {
	*min = 1;
	*var = &state->rq_maxiter;
    } else if (!strcmp(s, BKBP_K)) {
	*min = 1;
	*var = &state->bkbp_k;
    } else if (!strcmp(s, BOOTREP)) {
	*min = 1;
	*var = &state->bootrep;
    } else if (!strcmp(s, HAC_KERNEL)) {
	*min = 0;
	*max = KERNEL_MAX;
    } else if (!strcmp(s, HC_VERSION)) {
	*min = 0;
	*max = 4 + 1;
	*var = &state->ropts.hc_version;
    } else if (!strcmp(s, HORIZON)) {
	*min = 1;
	*var = &state->horizon;
    } else if (!strcmp(s, LONGDIGITS)) {
	*min = 1;
	*max = 21;
	*var = &state->longdigits;
    } else if (!strcmp(s, LOOP_MAXITER)) {
	*min = 1;
	*var = &state->loop_maxiter;
    } else if (!strcmp(s, VECM_NORM)) {
	*min = 0;
	*max = NORM_MAX;
	*var = &state->vecm_norm;
    } else if (!strcmp(s, GRETL_DEBUG)) {
	*min = 0;
	*max = 4;
	*var = &gretl_debug;
    } else {
	fprintf(stderr, "libset_set_int: unrecognized "
		"variable '%s'\n", s);	
	return E_UNKVAR;
    }

    return 0;
}

int libset_set_int (const char *key, int val)
{
    int min = 0, max = 0;
    int *ivar = NULL;
    int err = 0;

    if (check_for_state()) {
	return 1;
    }

    if (!strcmp(key, BLAS_NMK_MIN)) {
	set_blas_nmk_min(val);
	return 0;
    }

    err = intvar_min_max(key, &min, &max, &ivar);

    if (!err) {
	if (val < min || val >= max || ivar == NULL) {
	    err = E_DATA;
	} else {
	    *ivar = val;
	}
    }

    return err;
}

#ifndef WIN32
static int read_cli_shell_status (void)
{
    char shellstamp[FILENAME_MAX];
    FILE *fp;
    int ok = 0;

    sprintf(shellstamp, "%s.gretl_shell_stamp", gretl_dot_dir());
    fp = fopen(shellstamp, "r");
    if (fp != NULL) {
	ok = 1;
	fclose(fp);
    }

    return ok;
}
#endif

static int boolvar_get_flag (const char *s)
{
    if (!strcmp(s, ECHO)) {
	return STATE_ECHO_ON;
    } else if (!strcmp(s, MESSAGES)) {
	return STATE_MSGS_ON;
    } else if (!strcmp(s, WARNINGS)) {
	return STATE_WARN_ON;
    } else if (!strcmp(s, USE_SVD)) {
	return STATE_USE_SVD;
    } else if (!strcmp(s, USE_LBFGS)) {
	return STATE_USE_LBFGS;
    } else if (!strcmp(s, FORCE_DECP)) {
	return STATE_FORCE_DECPOINT;
    } else if (!strcmp(s, USE_CWD)) {
	return STATE_USE_CWD;
    } else if (!strcmp(s, USE_FCP)) {
	return STATE_USE_FCP;
    } else if (!strcmp(s, HALT_ON_ERR)) {
	return STATE_HALT_ON_ERR;
    } else if (!strcmp(s, MAX_VERBOSE)) {
	return STATE_MAX_VERBOSE;
    } else if (!strcmp(s, SHELL_OK)) {
	return STATE_SHELL_OK;
    } else if (!strcmp(s, FORCE_HC)) {
	return STATE_FORCE_HC;
    } else if (!strcmp(s, PREWHITEN)) {
	return STATE_PREWHITEN;
    } else if (!strcmp(s, PCSE)) {
	return STATE_USE_PCSE;
    } else if (!strcmp(s, VERBOSE_INCLUDE)) {
	return STATE_VERBOSE_INCLUDE;
    } else if (!strcmp(s, SKIP_MISSING)) {
	return STATE_SKIP_MISSING;
    } else {
	fprintf(stderr, "libset_get_bool: unrecognized "
		"variable '%s'\n", s);	
	return 0;
    }
}

static void set_flag_from_env (int flag, const char *s, int neg)
{
    char *e = getenv(s);
    int action = 0;

    if (e != NULL) {
	if (*e != '\0' && *e != '0') {
	    action = (neg)? -1 : 1;
	} else {
	    action = (neg)? 1 : -1;
	}
    }

    if (action > 0) {
	state->flags |= flag;
    } else if (action < 0) {
	state->flags &= ~flag;
    }
}

static void maybe_check_env (const char *s)
{
    if (!strcmp(s, USE_SVD)) {
	set_flag_from_env(STATE_USE_SVD, "GRETL_USE_SVD", 0);
    } else if (!strcmp(s, USE_LBFGS)) {
	set_flag_from_env(STATE_USE_LBFGS, "GRETL_USE_LBFGS", 0);
    } else if (!strcmp(s, HALT_ON_ERR)) {
	set_flag_from_env(STATE_HALT_ON_ERR, "GRETL_KEEP_GOING", 1);
    }

#ifndef WIN32
    if (!strcmp(s, SHELL_OK) && !gretl_in_gui_mode()) {
	if (read_cli_shell_status()) {
	    state->flags |= STATE_SHELL_OK;
	} else {
	    state->flags &= ~STATE_SHELL_OK;
	}
    }
#endif
}

int libset_get_bool (const char *key)
{
    int flag, ret = 0;

    /* global specials */

    if (!strcmp(key, PROTECT_LISTS)) {
	return protect_lists;
    } else if (!strcmp(key, R_FUNCTIONS)) {
	return R_functions;
    } else if (!strcmp(key, R_LIB)) {
	return R_lib;
    }

    if (!strcmp(key, MAX_VERBOSE) && gretl_debug > 1) {
	/* strong debugging turns on max_verbose */
	return 1;
    }

    if (check_for_state()) {
	return 0;
    }

    maybe_check_env(key);

    flag = boolvar_get_flag(key);
    if (flag == 0) {
	fprintf(stderr, "libset_get_bool: unrecognized "
		"variable '%s'\n", key);
	ret = 0;
    } else {
	ret = flag_to_bool(state, flag);
    }

    return ret;
}

static void libset_set_decpoint (int on)
{
#ifdef ENABLE_NLS
    static char num_locale[16];

    if (on) {
	char *orig = setlocale(LC_NUMERIC, "");

	*num_locale = '\0';
	strncat(num_locale, orig, 15);
	setlocale(LC_NUMERIC, "C");
    } else {
	setlocale(LC_NUMERIC, num_locale);
    }

    reset_local_decpoint();
#endif
}

static int check_R_setting (int *var, int val, const char *key)
{
    int err = 0;

#ifdef USE_RLIB
    *var = val;
#else
    if (val) {
	gretl_errmsg_sprintf("%s: not supported.", key);
	err = E_EXTERNAL;
    }
#endif

    return err;
}

int libset_set_bool (const char *key, int val)
{
    int flag, err = 0;

    if (check_for_state()) {
	return E_ALLOC;
    }

    /* global specials */

    if (!strcmp(key, PROTECT_LISTS)) {
	protect_lists = val;
	return 0;
    } else if (!strcmp(key, R_FUNCTIONS)) {
	return check_R_setting(&R_functions, val, key);
    } else if (!strcmp(key, R_LIB)) {
	return check_R_setting(&R_lib, val, key);
    }

    flag = boolvar_get_flag(key);

    if (flag == 0) {
	fprintf(stderr, "libset_set_bool: unrecognized "
		"variable '%s'\n", key);
	err = E_UNKVAR;
    } else if (val) {
	state->flags |= flag;
    } else {
	state->flags &= ~flag;
    }

    if (flag == STATE_FORCE_DECPOINT) {
	libset_set_decpoint(val);
    }

    return err;
}

/* Mechanism for pushing and popping program state for user-defined
   functions. push_program_state() is used when a function starts
   execution: the function gets a copy of the current program state,
   while that state is pushed onto the stack for restoration when the
   function exits.
*/

static int n_states;
static set_vars **state_stack;

int push_program_state (void)
{
    set_vars **sstack;
    set_vars *newstate;
    int ns = n_states;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "push_program_state: n_states = %d\n", ns);
#endif

    newstate = malloc(sizeof *newstate);

    if (newstate == NULL) {
	err = 1;
    } else {
	sstack = realloc(state_stack, (ns + 1) * sizeof *sstack);
	if (sstack == NULL) {
	    free(newstate);
	    err = 1;
	}
    }

    if (!err) {
	if (ns == 0) {
	    /* set all defaults */
	    state_vars_init(newstate);
	} else {
	    /* copy existing state */
	    state_vars_copy(newstate);
	}
	state_stack = sstack;
	state = state_stack[ns] = newstate;
	n_states++;
    }

#if PDEBUG
    if (!err) {
	fprintf(stderr, " state is now state_stack[%d]\n", ns);
    }
#endif

    return err;
}

static void free_state (set_vars *sv)
{
    if (sv->initvals != NULL) {
	gretl_matrix_free(sv->initvals);
    }

    free(sv);
}

/* Called when a user-defined function exits: restores the program
   state that was in force when the function started executing.
*/

int pop_program_state (void)
{
    set_vars **sstack;
    int ns = n_states;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "pop_program_state called: ns=%d\n", ns);
#endif

    if (ns < 2) {
	err = 1;
    } else {
	free_state(state_stack[ns - 1]);
	state_stack[ns - 1] = NULL;
	sstack = realloc(state_stack, (ns - 1) * sizeof *sstack);
	if (sstack == NULL) {
	    err = 1;
	}	
    }

    if (!err) {
	state_stack = sstack;
	state = state_stack[ns - 2];
	n_states--;
    }

#if PDEBUG
    fprintf(stderr, " state is now state_stack[%d]\n", ns - 2);
#endif

    return err;
}

/* initialization of all user-settable settings */

int libset_init (void)
{
    static int done;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "libset_init called, done=%d\n", done);
#endif

    if (!done) {
	err = push_program_state();
	done = 1;
    }

    return err;
}

void libset_cleanup (void)
{
    int i;

#if PDEBUG
    fprintf(stderr, "libset_cleanup called\n");
#endif

    for (i=0; i<n_states; i++) {
	free_state(state_stack[i]);
    }

    free(state_stack);
    state_stack = NULL;
    n_states = 0;
}

/* switches for looping and batch mode: output: these depend on the
   state of the program calling libgretl, but are not user-settable
*/

static int batch_mode;
static int gui_mode;

void set_loop_on (int quiet)
{
    state->flags |= STATE_LOOPING;
    if (quiet) {
	state->flags |= STATE_LOOP_QUIET;
    }
}

void set_loop_off (void)
{
    state->flags &= ~STATE_LOOPING;
    
    /* If we're not currently governed by "loop quietness" at
       caller level, turn such quietness off too 
    */
    if (state->flags & STATE_LOOP_QUIET) {
	int i = n_states - 1;

	if (i <= 0 || !(state_stack[i-1]->flags & STATE_LOOP_QUIET)) {
	    state->flags ^= STATE_LOOP_QUIET;
	}
    }
}

/* returns 1 if there's a loop going on anywhere in the "caller
   ancestry" of the current execution level, else 0.
*/

int gretl_looping (void)
{
    int i, ns = n_states;

    for (i=0; i<ns; i++) {
	if (state_stack[i]->flags & STATE_LOOPING) {
	    return 1;
	}
    }

    return 0;
}

/* returns 1 if there's a loop going on at the current execution
   stack level, else 0.
*/

int gretl_looping_currently (void)
{
    return (state->flags & STATE_LOOPING)? 1 : 0;
}

int gretl_looping_quietly (void)
{
    return (state->flags & STATE_LOOP_QUIET)? 1 : 0;
}

void gretl_set_batch_mode (int b)
{
    batch_mode = b;
}

int gretl_in_batch_mode (void)
{
    return batch_mode;
}

void gretl_set_gui_mode (int g)
{
    gui_mode = g;
}

int gretl_in_gui_mode (void)
{
    return gui_mode;
}

static ITER_PRINT_FUNC ifunc;

void set_iter_print_func (ITER_PRINT_FUNC func)
{
    ifunc = func;
}

int iter_print_func_installed (void)
{
    return ifunc != NULL;
}

int iter_print_callback (int i, PRN *prn)
{
    int ret = 0;

    if (ifunc != NULL) {
	ret = (*ifunc)(i, prn);
    }

    return ret;
}

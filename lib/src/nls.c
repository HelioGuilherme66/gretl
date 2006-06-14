/*
 *  Copyright (c) 2003-2005 by Allin Cottrell
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

/* Nonlinear least squares for libgretl, using minpack */

#include "libgretl.h" 
#include "libset.h"

#include "f2c.h"
#include "../../minpack/minpack.h"  

#define NLS_DEBUG 0

enum {
    NUMERIC_DERIVS,
    ANALYTIC_DERIVS
} nls_modes;

typedef struct _nls_param nls_param;

struct _nls_param {
    char name[VNAMELEN]; /* name of parameter (scalar variable) */
    char *deriv;         /* string representation of derivative of regression
			    function with respect to param (or NULL) */
    int varnum;          /* ID number of the scalar variable in dataset */
    int dernum;          /* ID number of the variable holding the derivative */
};

struct _nls_spec {
    int ci;             /* NLS or MLE */
    int mode;           /* derivatives: numeric or analytic */
    gretlopt opt;       /* can include OPT_V for verbose output */
    int depvar;         /* ID number of dependent variable */
    int uhatnum;        /* ID number of variable holding residuals */
    char *nlfunc;       /* string representation of nonlinear function,
			   expressed in terms of the residuals */
    int nparam;         /* number of parameters to be estimated */
    int naux;           /* number of auxiliary commands */
    int ngenrs;         /* number of variable-generating formulae */
    int iters;          /* number of iterations performed */
    int fncount;        /* number of function evaluations (ML) */
    int grcount;        /* number of gradient evaluations (ML) */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int nobs;           /* number of observations used */
    double ess;         /* error sum of squares */
    double ll;          /* log likelihood */
    double tol;         /* tolerance for stopping iteration */
    nls_param *params;  /* array of information on function parameters
			   (see the _nls_param struct above) */
    doublereal *coeff;  /* coefficient estimates */
    char **aux;         /* auxiliary commands */
    GENERATOR **genrs;  /* variable-generation pointers */
};

/* file-scope global variables: we need to access these variables in
   the context of the callback functions that we register with the
   minpack routines, and there is no provision in minpack for passing
   additional pointers into the callbacks.
*/

static double ***nZ;
static DATAINFO *ndinfo;
static PRN *nprn;

static nls_spec private_spec;
static nls_spec *pspec;
static integer one = 1;
static int genr_err;

static void update_nls_param_values (const double *x);

static void destroy_genrs_array (GENERATOR **genrs, int n)
{
    int i;

    for (i=0; i<n; i++) {
	destroy_genr(genrs[i]);
    }

    free(genrs);
}

/* new-style: "compile" the required equations first, so we can
   subsequently execute the compiled versions for greater
   efficiency */

static int nls_genr_setup (void)
{
    GENERATOR **genrs;
    char formula[MAXLINE];
    int i, j, n_gen, nparam;
    int v, err = 0;

    pspec->ngenrs = 0;

    nparam = (pspec->mode == ANALYTIC_DERIVS)? pspec->nparam : 0;

    n_gen = 1 + pspec->naux + nparam;

#if NLS_DEBUG
    fprintf(stderr, "nls_genr_setup: current v = %d, n_gen = %d\n", 
	    ndinfo->v, n_gen);
#endif

    genrs = malloc(n_gen * sizeof *genrs);
    if (genrs == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n_gen; i++) {
	genrs[i] = NULL;
    }

    j = 0;

    for (i=0; i<n_gen && !err; i++) {
	if (i < pspec->naux) {
	    /* auxiliary variables */
	    strcpy(formula, pspec->aux[i]);
	} else if (i == pspec->naux) {
	    /* residual/likelihood function */
	    sprintf(formula, "$nl_y = %s", pspec->nlfunc); 
	} else {
	    /* derivatives/gradients */
	    sprintf(formula, "$nl_x%d = %s", i, pspec->params[j++].deriv);
	}
	
	genrs[i] = genr_compile(formula, nZ, ndinfo, OPT_P, NULL);
	err = genr_get_err(genrs[i]);
#if NLS_DEBUG
	fprintf(stderr, "genrs[%d] = %p, err = %d\n", i, (void *) genrs[i], err);
#endif

	/* first pass calculating residual: ensure all parameters are
	   non-zero so we can get a valid reading on possible missing
	   values */
	if (i == pspec->naux) {
	    int k;

	    for (k=0; k<pspec->nparam; k++) {
		v = pspec->params[k].varnum;
		if ((*nZ)[v][0] == 0.0) {
		    (*nZ)[v][0] = 0.0001;
		}
	    }	    
	}

	if (!err) {
	    err = execute_genr(genrs[i], ndinfo->v);
	}

	if (!err) {
	    v = genr_get_varnum(genrs[i]);
	    if (i == pspec->naux) {
		pspec->uhatnum = v;
	    } else if (j > 0) {
		pspec->params[j-1].dernum = v;
	    }
#if NLS_DEBUG
	    fprintf(stderr, " genr->varnum = %d\n", v);
	    fprintf(stderr, " formula '%s'\n", formula);
#endif
	} else {
	    fprintf(stderr, "execute_genr: formula '%s', error = %d\n", formula, err);
	}
    }

    if (err) {
	destroy_genrs_array(genrs, n_gen);
    } else {
	pspec->ngenrs = n_gen;
	pspec->genrs = genrs;
    }

    return err;
}

static int nls_auto_genr (int i)
{
    int j;

    if (pspec->genrs == NULL) {
	genr_err = nls_genr_setup();
	if (genr_err) {
	    fprintf(stderr, "nls_genr_setup failed\n");
	}
	return genr_err;
    }

    for (j=0; j<pspec->naux; j++) {
#if NLS_DEBUG
	fprintf(stderr, "nls_auto_genr: generating aux var:\n %s\n", pspec->aux[j]);
#endif
	genr_err = execute_genr(pspec->genrs[j], ndinfo->v);
    }

    j = pspec->naux + i;
#if NLS_DEBUG
    fprintf(stderr, "nls_auto_genr: executing genr[%d] at %p, with oldv = %d\n",
	    j, (void *) pspec->genrs[j], ndinfo->v);
#endif
    genr_err = execute_genr(pspec->genrs[j], ndinfo->v);


#if NLS_DEBUG
    if (genr_err) {
	fprintf(stderr, " varnum = %d, err = %d\n", genr_get_varnum(pspec->genrs[j]), 
		genr_err);
	errmsg(genr_err, nprn);
    } 
#endif

    return genr_err;    
}

/* wrappers for the above to enhance comprehensibility below */

static int nls_calculate_fvec (void)
{
    return nls_auto_genr(0);
}

static int nls_calculate_deriv (int i)
{
    return nls_auto_genr(i + 1);
}

static int nlspec_allocate_param (nls_spec *spec)
{
    nls_param *params;
    double *coeff;
    int nt = spec->nparam + 1;

    params = realloc(spec->params, nt * sizeof *spec->params);
    if (params == NULL) {
	return 1;
    }

    spec->params = params;

    coeff = realloc(spec->coeff, nt * sizeof *spec->coeff);
    if (coeff == NULL) {
	free(params);
	return 1;
    }

    spec->coeff = coeff;

    spec->nparam += 1;

    return 0;
}

/* allocate space for an additional regression parameter in the
   nls_spec struct and add its info */

static int 
real_add_param_to_spec (const char *vname, int vnum, double initval,
			nls_spec *spec)
{
    int i;

#if NLS_DEBUG
    fprintf(stderr, "real_add_param: adding '%s'\n", vname);
#endif

    if (nlspec_allocate_param(spec)) {
	return E_ALLOC;
    }

    i = spec->nparam - 1;
    
    spec->params[i].varnum = vnum;
    spec->params[i].dernum = 0;
    strcpy(spec->params[i].name, vname);
    spec->params[i].deriv = NULL;

    spec->coeff[i] = initval;

    return 0;
}

/* scrutinize word and see if it's a new scalar that should
   be added to the NLS specification; if so, add it */

static int 
maybe_add_param_to_spec (nls_spec *spec, const char *word, 
			 const double **Z, const DATAINFO *pdinfo)
{
    int i, v;

#if NLS_DEBUG
    fprintf(stderr, "maybe_add_param: looking at '%s'\n", word);
#endif

    /* if word represents a math function or constant, skip it */
    if (genr_function_from_string(word) || !strcmp(word, "pi")) {
	return 0;
    }

    /* try looking up word as the name of a variable */
    v = varindex(pdinfo, word);

    if (v < pdinfo->v) {
	/* existing variable */
	if (var_is_series(pdinfo, v)) {
	    /* if term is not a scalar, skip it: only scalars can figure
	       as regression parameters */
	    return 0;
	}
    } else {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), word);
	return E_UNKVAR;
    }

    /* if this term is already present in the specification, skip it */
    for (i=0; i<spec->nparam; i++) {
	if (strcmp(word, spec->params[i].name) == 0) {
	    return 0;
	}
    }

#if NLS_DEBUG
    fprintf(stderr, "maybe_add_param: adding '%s'\n", word);
#endif

    /* else: add this param to the NLS specification */
    return real_add_param_to_spec(word, v, Z[v][0], spec);
}

/* Parse NLS function specification string to find names of variables
   that may figure as parameters of the regression function.  We need
   to do this only if analytical derivatives have not been supplied.
*/

static int 
get_params_from_nlfunc (nls_spec *spec, const double **Z,
			const DATAINFO *pdinfo)
{
    const char *s = spec->nlfunc;
    const char *p;
    char vname[VNAMELEN];
    int n, np = 0;
    int err = 0;

#if NLS_DEBUG
    fprintf(stderr, "get_params: looking at '%s'\n", s);
#endif

    while (*s && !err) {
	p = s;
	if (isalpha(*s) && *(s + 1)) { 
	    /* find a variable name */
	    n = gretl_varchar_spn(s);
	    p += n;
	    if (n > VNAMELEN - 1) {
		/* variable name is too long */
		return 1;
	    }
	    *vname = 0;
	    strncat(vname, s, n);
	    if (np > 0) {
		/* right-hand side term */
		err = maybe_add_param_to_spec(spec, vname, Z, pdinfo);
	    }
	    np++;
	} else {
	    p++;
	}
	s = p;
    }

    return err;
}

static int nls_spec_allocate_params (nls_spec *spec, int np)
{
    spec->params = malloc(np * sizeof *spec->params);
    if (spec->params == NULL) {
	return E_ALLOC;
    }

    spec->coeff = malloc(np * sizeof *spec->coeff);
    if (spec->coeff == NULL) {
	free(spec->params);
	spec->params = NULL;
	return E_ALLOC;
    }    

    return 0;
}

/* For case where analytical derivatives are not given, the user
   may supply a line like:

     params b0 b1 b2 b3

   specifying the parameters to be estimated.  Here we parse
   such a list and add the parameter info to spec.  The terms
   in the list must be pre-existing scalar variables.
*/

static int 
nls_spec_add_params_from_line (nls_spec *spec, const char *s,
			       const double **Z, const DATAINFO *pdinfo)
{
    int i, nf = count_fields(s);
    const char *p = s;
    int err = 0;

    if (spec->params != NULL || nf == 0) {
	return E_DATA;
    }

    err = nls_spec_allocate_params(spec, nf);
    if (err) {
	return err;
    }

    for (i=0; i<nf && !err; i++) {
	char *pname = gretl_word_strdup(p, &p);
	int v;

	if (pname != NULL) {
	    if (strlen(pname) > VNAMELEN - 1) {
		pname[VNAMELEN - 1] = '\0';
	    }
	    v = varindex(pdinfo, pname);
	    if (v >= pdinfo->v || var_is_series(pdinfo, v)) {
		err = E_DATA;
	    } else {
		spec->params[i].varnum = v;
		spec->params[i].dernum = 0;
		strcpy(spec->params[i].name, pname);
		spec->params[i].deriv = NULL;
		spec->coeff[i] = Z[v][0];
	    }
	    free(pname);
	} else {
	    err = E_ALLOC;
	}
    }

    if (err) {
	free(spec->params);
	spec->params = NULL;
	free(spec->coeff);
	spec->coeff = NULL;
	spec->nparam = 0;
    } else {
	spec->nparam = nf;
    }

    return err;
}

/**
 * nls_spec_add_param_list:
 * @spec: nls specification.
 * @list: list of variables by ID number.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Adds to @spec a list of (scalar) parameters to be estimated, as
 * given in @list.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nls_spec_add_param_list (nls_spec *spec, const int *list,
			     const double **Z, const DATAINFO *pdinfo)
{
    int i, np = list[0];
    int err = 0;

    if (spec->params != NULL || np == 0) {
	return E_DATA;
    }

    err = nls_spec_allocate_params(spec, np);
    if (err) {
	return err;
    }

    for (i=0; i<np && !err; i++) {
	int v = list[i+1];

	if (v >= pdinfo->v || var_is_series(pdinfo, v)) {
	    err = E_DATA;
	} else {
	    spec->params[i].varnum = v;
	    spec->params[i].dernum = 0;
	    strcpy(spec->params[i].name, pdinfo->varname[v]);
	    spec->params[i].deriv = NULL;
	    spec->coeff[i] = Z[v][0];
	}
    }

    if (err) {
	free(spec->params);
	spec->params = NULL;
	free(spec->coeff);
	spec->coeff = NULL;
	spec->nparam = 0;
    } else {
	spec->nparam = np;
    }

    return err;
}

/* Adjust starting and ending points of sample if need be, to avoid
   missing values; abort if there are missing values within the
   (possibly reduced) sample range.  For this purpose we generate
   the nls residual variable.
*/

static int nls_missval_check (nls_spec *spec)
{
    int t, v, miss = 0;
    int t1 = spec->t1, t2 = spec->t2;
    int err = 0;

#if NLS_DEBUG
    fprintf(stderr, "nls_missval_check: calling nls_calculate_fvec\n");
#endif

    /* generate the nls residual variable */
    err = nls_calculate_fvec();
    if (err) {
	return err;
    }
	
    /* ID number of LHS variable */
    v = pspec->uhatnum;

#if NLS_DEBUG
    fprintf(stderr, "nls_missval_check: checking var %d (%s)\n",
	    v, ndinfo->varname[v]);
    fprintf(stderr, " before trimming: spec->t1 = %d, spec->t2 = %d\n",
	    spec->t1, spec->t2);
#endif

    for (t=spec->t1; t<=spec->t2; t++) {
	if (na((*nZ)[v][t])) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=spec->t2; t>=spec->t1; t--) {
	if (na((*nZ)[v][t])) {
	    t2--;
	} else {
	    break;
	}
    }

    if (t2 - t1 + 1 < spec->nparam) {
	return E_DF;
    }

    for (t=t1; t<=t2; t++) {
	if (na((*nZ)[v][t])) {
	    fprintf(stderr, " nls_missval_check: after setting t1=%d, t2=%d, "
		    "got NA for var %d at obs %d\n", t1, t2, v, t);
	    miss = 1;
	    break;
	}
    }  

    if (miss) {
	strcpy(gretl_errmsg, _("There were missing data values"));
	return 1;
    }

    spec->t1 = t1;
    spec->t2 = t2;
    spec->nobs = t2 - t1 + 1;

#if NLS_DEBUG
    fprintf(stderr, " after: spec->t1 = %d, spec->t2 = %d, spec->nobs = %d\n\n",
	    spec->t1, spec->t2, spec->nobs);
#endif

    return 0;
}

/* this function is used in the context of BFGS */

static double get_mle_ll (const double *b, void *unused)
{
    int t, v = pspec->uhatnum;

    update_nls_param_values(b);

    /* calculate log-likelihood given current parameter estimates */
    if (nls_calculate_fvec()) {
	return NADBL;
    }

    pspec->ll = 0.0;
    for (t=pspec->t1; t<=pspec->t2; t++) {
	pspec->ll -= (*nZ)[v][t];
    }

    return pspec->ll;
}

/* analytical derivatives, used in the context of BFGS */

static int get_mle_gradient (double *b, double *g, int n, 
			     BFGS_LL_FUNC llfunc,
			     void *unused)
{
    int i, t, v;
    int err = 0;

    update_nls_param_values(b);

    for (i=0; i<pspec->nparam; i++) {
	if (nls_calculate_deriv(i)) {
	    fprintf(stderr, "error calculating deriv\n");
	    return 1;
	}

	v = pspec->params[i].dernum;

	g[i] = 0.0;
	if (var_is_series(ndinfo, v)) {
	    /* derivative may be vector or scalar */
	    for (t=pspec->t1; t<=pspec->t2; t++) {
		if (na((*nZ)[v][t])) {
		    fprintf(stderr, "NA in gradient calculation\n");
		    err = 1;
		} else {
		    g[i] += (*nZ)[v][t];
		}
	    }
	} else {
	    g[i] += (*nZ)[v][0];
	}
#if NLS_DEBUG > 1
	fprintf(stderr, "g[%d] = %g, based on nZ[%d]\n", i, g[i], v);
#endif
    }

    return err;
}

/* default numerical calculation of gradient in context of BFGS */

static int BFGS_numeric_gradient (double *b, double *g, int n,
				  BFGS_LL_FUNC llfunc,
				  void *p)
{
    const double h = 1.0e-8;
    double b0, f1, f2;
    int i;

    /* consider Richardson extrapolation here? */

    for (i=0; i<n; i++) {
	b0 = b[i];
	b[i] = b0 - h;
	f1 = llfunc(b, p);
	if (na(f1)) {
	    b[i] = b0;
	    return 1;
	}
	b[i] = b0 + h;
	f2 = llfunc(b, p);
	if (na(f2)) {
	    b[i] = b0;
	    return 1;
	}
	g[i] = (f2 - f1) / (2.0 * h);
    }

    return 0;
}

static void 
print_mle_iter_stats (double ll, int nparam, const double *b, const double *g, 
		      int iter, double sl, PRN *prn)
{
    int i;

    if (na(ll)) {
	pprintf(prn, "Iteration %d: log likelihood = NA", iter);	
    } else {
	pprintf(prn, "Iteration %d: log likelihood = %.8g", iter, ll);
    }
    if (iter > 1) {
	pprintf(prn, " (steplength = %.8g)", sl);
    }	
    pputc(prn, '\n');
	
    pputs(prn, "Parameters: ");
    for (i=0; i<nparam; i++) {
	pprintf(prn, "%#12.5g", b[i]);
    }
    pputc(prn, '\n');

    pputs(prn, "Gradients:  ");
    for (i=0; i<nparam; i++) {
	pprintf(prn, "%#12.5g", -g[i]);
    }
    pputs(prn, "\n\n");
}

/* this function is used in the context of the minpack callback, and
   also for checking derivatives in the MLE case
*/

static int get_nls_fvec (double *fvec)
{
    int j, t, v = pspec->uhatnum;

#if NLS_DEBUG > 1
    fprintf(stderr, "*** get_nls_fvec called\n");
#endif

    /* calculate residual given current parameter estimates */
    if (nls_calculate_fvec()) {
	return 1;
    }

    pspec->ess = 0.0;
    pspec->ll = 0.0;

    j = 0;

    /* transcribe from dataset to fvec array */
    for (t=pspec->t1; t<=pspec->t2; t++) {
	if (na((*nZ)[v][t])) {
	    fprintf(stderr, "nls_calculate_fvec: produced NA at obs %d\n", t);
	    return 1;
	}
	fvec[j] = (*nZ)[v][t];
#if NLS_DEBUG > 1
	fprintf(stderr, "fvec[%d] = nZ[%d][%d] = %g\n", j, v, t, fvec[j]);
#endif
	if (pspec->ci == MLE) {
	    pspec->ll -= fvec[j];
	} else {
	    pspec->ess += fvec[j] * fvec[j];
	}
	j++;
    }

    pspec->iters += 1;

    if (pspec->ci == NLS && (pspec->opt & OPT_V)) {
	pprintf(nprn, "iteration %2d: SSR = %.8g\n", pspec->iters, pspec->ess);
    }

    return 0;
}

/* this function is used in the context of the minpack callback */

static int get_nls_deriv (int i, double *deriv)
{
    int j, t, vec, v = pspec->params[i].dernum;

#if NLS_DEBUG
    fprintf(stderr, "get_nls_deriv: getting deriv %d\n", i);
#endif

    /* calculate value of deriv with respect to param */
    if (nls_calculate_deriv(i)) {
	return 1;
    }

    if (v == 0) { /* FIXME */
	v = pspec->params[i].dernum = ndinfo->v - 1;
    }

    /* derivative may be vector or scalar */
    vec = var_is_series(ndinfo, v);

#if NLS_DEBUG
    fprintf(stderr, " v = %d, vec = %d\n", v, vec);
#endif

    j = 0;
    /* transcribe from dataset to deriv array */
    for (t=pspec->t1; t<=pspec->t2; t++) {
	if (vec) {
	    deriv[j] = - (*nZ)[v][t];
	} else {
	    deriv[j] = - (*nZ)[v][0];
	}
#if NLS_DEBUG > 1
	fprintf(stderr, " set deriv[%d] = nZ[%d][%d] = %g\n", j, v, 
		(vec)? t : 0, deriv[j]);
#endif
	j++;
    }

    return 0;
}

/* compute auxiliary statistics and add them to the NLS 
   model struct */

static void add_stats_to_model (MODEL *pmod, nls_spec *spec,
				const double **Z)
{
    int dv = spec->depvar;
    double d, tss;
    int t;

    pmod->ess = spec->ess;
    pmod->sigma = sqrt(spec->ess / (pmod->nobs - spec->nparam));
    
    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, Z[dv]);

    tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	d = Z[dv][t] - pmod->ybar;
	tss += d * d;
    }  

    if (tss == 0.0) {
	pmod->rsq = NADBL;
    } else {
	pmod->rsq = 1.0 - pmod->ess / tss;
    }

    pmod->adjrsq = NADBL;
}

/* MLE: add variance matrix elements and standard errors based
   on OPG: (GG')^{-1} */

static int add_OPG_vcv (MODEL *pmod, nls_spec *spec)
{
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    double x = 0.0;
    int k = spec->nparam;
    int T = spec->nobs;
    int i, v, t, err = 0;

    G = gretl_matrix_alloc(k, T);
    V = gretl_matrix_alloc(k, k);
    if (G == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (pspec->mode == NUMERIC_DERIVS) {
	const double eps = 1e-8;
	double bi, x0, x1;

	/* construct approximation to G matrix based
	   on finite differences */

	for (i=0; i<pspec->nparam; i++) {
	    bi = spec->coeff[i];
	    spec->coeff[i] = bi - eps;
	    update_nls_param_values(spec->coeff);
	    nls_calculate_fvec();
	    for (t=0; t<T; t++) {
		x0 = (*nZ)[pspec->uhatnum][t + spec->t1];
		gretl_matrix_set(G, i, t, x0);
	    }
	    spec->coeff[i] = bi + eps;
	    update_nls_param_values(spec->coeff);
	    nls_calculate_fvec();
	    for (t=0; t<T; t++) {
		x1 = (*nZ)[pspec->uhatnum][t + spec->t1];
		x0 = gretl_matrix_get(G, i, t);
		gretl_matrix_set(G, i, t, (x1 - x0) / (2.0 * eps));
	    }
	    spec->coeff[i] = bi;
	}
    } else {
	for (i=0; i<k; i++) {
	    v = spec->params[i].dernum;
	    if (var_is_scalar(ndinfo, v)) {
		x = (*nZ)[v][0];
	    }
	    for (t=0; t<T; t++) {
		if (var_is_series(ndinfo, v)) {
		    x = (*nZ)[v][t + spec->t1];
		}
		gretl_matrix_set(G, i, t, x);
	    }
	}
    }

    gretl_matrix_multiply_mod(G, GRETL_MOD_NONE,
			      G, GRETL_MOD_TRANSPOSE,
			      V);

    err = gretl_invert_symmetric_matrix(V);

    if (!err) {
	int j;

	for (i=0; i<k; i++) {
	    for (j=0; j<=i; j++) {
		x = gretl_matrix_get(V, i, j);
		pmod->vcv[ijton(i, j, k)] = x;
		if (i == j) {
		    pmod->sderr[i] = sqrt(x);
		} 
	    }
	}
    }

 bailout:

    gretl_matrix_free(G);
    gretl_matrix_free(V);

    return err;
}

/* NLS: add coefficient covariance matrix and standard errors 
   based on GNR */

static int add_nls_std_errs_to_model (MODEL *pmod)
{
    int i, k;

    if (pmod->vcv == NULL && makevcv(pmod)) {
	return E_ALLOC;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	k = ijton(i, i, pmod->ncoeff);
	if (pmod->vcv[k] == 0.0) {
	    pmod->sderr[i] = 0.0;
	} else if (pmod->vcv[k] > 0.0) {
	    pmod->sderr[i] = sqrt(pmod->vcv[k]);
	} else {
	    pmod->sderr[i] = NADBL;
	}
    }

    return 0;
}

/* transcribe coefficient estimates into model struct */

static void add_coeffs_to_model (MODEL *pmod, double *coeff)
{
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = coeff[i];
    }
}

static int 
add_param_names_to_model (MODEL *pmod, nls_spec *spec, const DATAINFO *pdinfo)
{
    int i, np = pmod->ncoeff + 1;

    pmod->params = malloc(np * sizeof *pmod->params);
    if (pmod->params == NULL) {
	return 1;
    }
    pmod->nparams = np;

    if (spec->ci == MLE) {
	int n = strlen(spec->nlfunc);

	pmod->params[0] = malloc(n + 3);
	if (pmod->params[0] != NULL) {
	    sprintf(pmod->params[0], "l = %s", spec->nlfunc + 2);
	    n = strlen(pmod->params[0]);
	    pmod->params[0][n-1] = '\0';
	} 
    } else {
	pmod->params[0] = gretl_strdup(pdinfo->varname[spec->depvar]);
    } 

    if (pmod->params[0] == NULL) {
	free(pmod->params);
	return 1;
    }    

    for (i=1; i<=pmod->ncoeff; i++) {
	pmod->params[i] = gretl_strdup(spec->params[i-1].name);
	if (pmod->params[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) {
		free(pmod->params[j]);
	    }
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    return 1;
	}
    }

    return 0;
}

static void 
add_fit_resid_to_model (MODEL *pmod, nls_spec *spec, double *uhat, 
			const double **Z, int perfect)
{
    int t, j = 0;

    if (perfect) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] = 0.0;
	    pmod->yhat[t] = Z[spec->depvar][t];
	}
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] = uhat[j];
	    pmod->yhat[t] = Z[spec->depvar][t] - uhat[j];
	    j++;
	}
    }
}

/* this may be used later for generating out-of-sample forecasts --
   see nls_forecast() in forecast.c
*/

static int transcribe_nls_function (MODEL *pmod, const char *s)
{
    char *formula;
    int err = 0;

    /* skip "depvar - " */
    s += strcspn(s, " ") + 3;

    formula = gretl_strdup(s);
    if (s != NULL) {
	gretl_model_set_string_as_data(pmod, "nl_regfunc", formula); 
    } else {
	err = E_ALLOC;
    }

    return err;
}

#if NLS_DEBUG > 1
static void 
print_GNR_dataset (const int *list, double **gZ, DATAINFO *gdinfo)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
    int t1 = gdinfo->t1;

    fprintf(stderr, "gdinfo->t1 = %d, gdinfo->t2 = %d\n",
	    gdinfo->t1, gdinfo->t2);
    gdinfo->t1 = 0;
    printdata(list, (const double **) gZ, gdinfo, OPT_O, prn);
    gdinfo->t1 = t1;
    gretl_print_destroy(prn);
}
#endif

/* Gauss-Newton regression to calculate standard errors for the NLS
   parameters (see Davidson and MacKinnon).  This model is also taken
   as the basis for the model struct returned by the nls function.
*/

static MODEL GNR (double *uhat, double *jac, nls_spec *spec,
		  const double **Z, const DATAINFO *pdinfo, 
		  PRN *prn)
{
    double **gZ = NULL;
    DATAINFO *gdinfo;
    int *glist;
    MODEL gnr;
    gretlopt lsqopt;
    int i, j, t;
    int T = spec->nobs;
    int iters = spec->iters;
    int perfect = 0;
    int err = 0;

    if (gretl_iszero(0, spec->nobs - 1, uhat)) {
	pputs(prn, _("Perfect fit achieved\n"));
	perfect = 1;
	for (t=0; t<spec->nobs; t++) {
	    uhat[t] = 1.0;
	}
	spec->ess = 0.0;
    }

    /* number of variables = 1 (const) + 1 (depvar) + spec->nparam
       (derivatives) */
    gdinfo = create_new_dataset(&gZ, spec->nparam + 2, pdinfo->n, 0);
    if (gdinfo == NULL) {
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    /* transcribe sample info */
    gdinfo->t1 = spec->t1;
    gdinfo->t2 = spec->t2;

#if 0
    fprintf(stderr, "pdinfo->n = %d, gdinfo->t1 = %d, gdinfo->t2 = %d\n",
	    pdinfo->n, gdinfo->t1, gdinfo->t2);
#endif
    
    glist = gretl_list_new(spec->nparam + 1);

    if (glist == NULL) {
	destroy_dataset(gZ, gdinfo);
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    j = 0;

    /* dependent variable (NLS residual) */
    glist[1] = 1;
    strcpy(gdinfo->varname[1], "gnr_y");
    for (t=0; t<gdinfo->n; t++) {
	if (t < gdinfo->t1 || t > gdinfo->t2) {
	    gZ[1][t] = NADBL;
	} else {
	    gZ[1][t] = uhat[j++];
	}
    }

    for (i=0; i<spec->nparam; i++) {
	/* independent vars: derivatives wrt NLS params */
	int v = i + 2;

	glist[v] = v;
	sprintf(gdinfo->varname[v], "gnr_x%d", i + 1);
	if (spec->mode == ANALYTIC_DERIVS) {
	    for (t=0; t<gdinfo->t1; t++) {
		gZ[v][t] = NADBL;
	    }
	    for (t=gdinfo->t2; t<gdinfo->n; t++) {
		gZ[v][t] = NADBL;
	    }
	    get_nls_deriv(i, gZ[v] + gdinfo->t1);
	} else {
	    j = T * i; /* calculate offset into jac */
	    for (t=0; t<gdinfo->n; t++) {
		if (t < gdinfo->t1 || t > gdinfo->t2) {
		    gZ[v][t] = NADBL;
		} else {
		    gZ[v][t] = jac[j++];
		}
	    }
	}
    }

#if NLS_DEBUG > 1
    print_GNR_dataset(glist, gZ, gdinfo);
#endif

    lsqopt = OPT_A;
    if (spec->opt & OPT_R) {
	/* robust variance matrix, if wanted */
	lsqopt |= OPT_R;
    }

    gnr = lsq(glist, &gZ, gdinfo, OLS, lsqopt);

#if NLS_DEBUG
    gnr.name = gretl_strdup("GNR for NLS");
    printmodel(&gnr, gdinfo, OPT_NONE, prn);
    free(gnr.name);
    gnr.name = NULL;
#endif

    if (gnr.errcode) {
	pputs(prn, _("In Gauss-Newton Regression:\n"));
	errmsg(gnr.errcode, prn);
	err = 1;
    } 

    if (gnr.list[0] != glist[0]) {
	strcpy(gretl_errmsg, _("Failed to calculate Jacobian"));
	gnr.errcode = E_DATA;
    }

    if (gnr.errcode == 0) {
	gnr.ci = spec->ci;
	add_stats_to_model(&gnr, spec, Z);
	if (add_nls_std_errs_to_model(&gnr)) {
	    gnr.errcode = E_ALLOC;
	}
    }

    if (gnr.errcode == 0) {
	ls_criteria(&gnr);
	add_coeffs_to_model(&gnr, spec->coeff);
	add_param_names_to_model(&gnr, spec, pdinfo);
	add_fit_resid_to_model(&gnr, spec, uhat, Z, perfect);
	gnr.list[1] = spec->depvar;

	/* set relevant data on model to be shipped out */
	gretl_model_set_int(&gnr, "iters", iters);
	gretl_model_set_double(&gnr, "tol", spec->tol);
	transcribe_nls_function(&gnr, pspec->nlfunc);
    }

    destroy_dataset(gZ, gdinfo);
    free(glist);

    return gnr;
}

/* allocate space to copy info into MLE model struct */

static int mle_model_allocate (MODEL *pmod, nls_spec *spec)
{
    int k = spec->nparam;
    int nvc = (k * k + k) / 2;

    pmod->coeff = malloc(k * sizeof *pmod->coeff);
    pmod->sderr = malloc(k * sizeof *pmod->sderr);
    pmod->vcv = malloc(nvc * sizeof *pmod->vcv);

    if (pmod->coeff == NULL || pmod->sderr == NULL || pmod->vcv == NULL) {
	pmod->errcode = E_ALLOC;
    } else {
	pmod->ncoeff = k;
    }

    return pmod->errcode;
}

/* work up the results of ML estimation into the form of a gretl
   MODEL */

static int make_mle_model (MODEL *pmod, nls_spec *spec, 
			   const DATAINFO *pdinfo)
{
    mle_model_allocate(pmod, spec);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    pmod->t1 = pspec->t1;
    pmod->t2 = pspec->t2;
    pmod->nobs = pspec->nobs;
    
    /* hmm */
    pmod->dfn = pmod->ncoeff;
    pmod->dfd = pmod->nobs - pmod->ncoeff;

    pmod->ci = MLE;
    pmod->lnL = pspec->ll;
    mle_criteria(pmod, 0);

    add_coeffs_to_model(pmod, spec->coeff);

    pmod->errcode = add_param_names_to_model(pmod, spec, pdinfo);

    if (!pmod->errcode) {
	pmod->errcode = add_OPG_vcv(pmod, spec);
    }

    if (!pmod->errcode) {
	gretl_model_set_int(pmod, "fncount", pspec->fncount);
	gretl_model_set_int(pmod, "grcount", pspec->grcount);
	gretl_model_set_double(pmod, "tol", spec->tol);
    }

    return pmod->errcode;
}

static int add_nls_coeffs (MODEL *pmod, nls_spec *spec)
{
    pmod->ncoeff = spec->nparam;
    pmod->full_n = 0;

    pmod->errcode = 
	gretl_model_allocate_storage(pmod);

    if (!pmod->errcode) {
	add_coeffs_to_model(pmod, spec->coeff);
    }

    return pmod->errcode;
}

/* free up resources associated with the nlspec struct */

static void clear_nls_spec (nls_spec *spec)
{
    int i;

    if (spec == NULL) {
	return;
    }

    if (spec->params != NULL) {
	for (i=0; i<spec->nparam; i++) {
	    free(spec->params[i].deriv);
	}
	free(spec->params);
	spec->params = NULL;
    }

    if (spec->aux != NULL) {
	for (i=0; i<spec->naux; i++) {
	    free(spec->aux[i]);
	}
	free(spec->aux);
	spec->aux = NULL;
    }

    if (spec->genrs != NULL) {
	for (i=0; i<spec->ngenrs; i++) {
	    destroy_genr(spec->genrs[i]);
	}
	free(spec->genrs);
	spec->genrs = NULL;
    }    

    free(spec->nlfunc);
    spec->nlfunc = NULL;

    free(spec->coeff);
    spec->coeff = NULL;

    spec->ci = NLS;
    spec->mode = NUMERIC_DERIVS;
    spec->opt = OPT_NONE;

    spec->nparam = 0;
    spec->naux = 0;
    spec->ngenrs = 0;

    spec->depvar = 0;
    spec->uhatnum = 0;

    spec->iters = 0;
    spec->fncount = 0;
    spec->grcount = 0;

    spec->t1 = spec->t2 = 0;
    spec->nobs = 0;
}

/* 
   Next block: functions that interface with minpack.

   The details below may be obscure, but here's the basic idea: The
   minpack functions are passed an array ("fvec") that holds the
   calculated values of the function to be minimized, at a given value
   of the parameters, and also (in the case of analytic derivatives)
   an array holding the Jacobian ("jac").  Minpack is also passed a
   callback function that will recompute the values in these arrays,
   given a revised vector of parameter estimates.

   As minpack does its iterative thing, at each step it invokes the
   callback function, supplying its updated parameter estimates and
   saying via a flag variable ("iflag") whether it wants the function
   itself or the Jacobian re-evaluated.

   The libgretl strategy involves holding all the relevant values (nls
   residual, nls derivatives, and nls parameters) as variables in a
   gretl dataset (Z-array and datainfo-struct pair).  The callback
   function that we supply to minpack first transcribes the revised
   parameter estimates into the dataset, then invokes genr() to
   recalculate the residual and derivatives, then transcribes the
   results back into the fvec and jac arrays.
*/

static void update_nls_param_values (const double *x)
{
    int i, v;

    /* write the values produced by minpack into the dataset */
    for (i=0; i<pspec->nparam; i++) {
	v = pspec->params[i].varnum;
	(*nZ)[v][0] = x[i];
#if NLS_DEBUG
	fprintf(stderr, "revised param[%d] = %g\n", i, x[i]);
#endif
    }
}

/* callback for lm_calculate (below) to be used by minpack */

static int nls_calc (integer *m, integer *n, double *x, double *fvec, 
		     double *jac, integer *ldjac, integer *iflag)
{
    int T = *m;
    int i;

#if NLS_DEBUG
    fprintf(stderr, "nls_calc called by minpack with iflag = %d\n", 
	    (int) *iflag);
#endif

    /* write current parameter values into dataset Z */
    update_nls_param_values(x);

    if (*iflag == 1) {
	/* calculate function at x, results into fvec */
	if (get_nls_fvec(fvec)) {
	    *iflag = -1;
	}
    } else if (*iflag == 2) {
	/* calculate jacobian at x, results into jac */
	for (i=0; i<*n; i++) {
	    if (get_nls_deriv(i, &jac[i*T])) {
		*iflag = -1; 
	    }
	}	
    }

    return 0;
}

/* in case the user supplied analytical derivatives for the
   parameters, check them for sanity */

static int check_derivatives (integer m, integer n, double *x,
			      double *fvec, double *jac,
			      integer ldjac, PRN *prn)
{
#if NLS_DEBUG > 1
    int T = pspec->nobs * pspec->nparam;
#endif
    integer mode, iflag;
    doublereal *xp = NULL;
    doublereal *err = NULL;
    doublereal *fvecp = NULL;
    int i, badcount = 0, zerocount = 0;

    xp = malloc(n * sizeof *xp);
    err = malloc(m * sizeof *err);
    fvecp = malloc(m * sizeof *fvecp);

    if (xp == NULL || err == NULL || fvecp == NULL) {
	free(err);
	free(xp);
	free(fvecp);
	return 1;
    }

#if NLS_DEBUG > 1
    fprintf(stderr, "\nchkder, starting: m=%d, n=%d, ldjac=%d\n",
	    (int) m, (int) n, (int) ldjac);
    for (i=0; i<pspec->nparam; i++) {
	fprintf(stderr, "x[%d] = %g\n", i, x[i]);
    }    
    for (i=0; i<pspec->nobs; i++) {
	fprintf(stderr, "fvec[%d] = %g\n", i, fvec[i]);
    }
#endif

    /* mode 1: x contains the point of evaluation of the function; on
       output xp is set to a neighboring point. */
    mode = 1;
    chkder_(&m, &n, x, fvec, jac, &ldjac, xp, fvecp, &mode, err);

    /* calculate gradient */
    iflag = 2;
    nls_calc(&m, &n, x, fvec, jac, &ldjac, &iflag);
    if (iflag == -1) goto chkderiv_abort;

#if NLS_DEBUG > 1
    fprintf(stderr, "\nchkder, calculated gradient\n");
    for (i=0; i<T; i++) {
	fprintf(stderr, "jac[%d] = %g\n", i, jac[i]);
    }
#endif

    /* calculate function, at neighboring point xp */
    iflag = 1;
    nls_calc(&m, &n, xp, fvecp, jac, &ldjac, &iflag);
    if (iflag == -1) goto chkderiv_abort; 

    /* mode 2: on input, fvec must contain the functions, the rows of
       fjac must contain the gradients evaluated at x, and fvecp must
       contain the functions evaluated at xp.  On output, err contains
       measures of correctness of the respective gradients.
    */
    mode = 2;
    chkder_(&m, &n, x, fvec, jac, &ldjac, xp, fvecp, &mode, err);

#if NLS_DEBUG > 1
    fprintf(stderr, "\nchkder, done mode 2:\n");
    for (i=0; i<m; i++) {
	fprintf(stderr, "%d: fvec = %.12g, fvecp = %.12g, err = %g\n", i, 
		fvec[i], fvecp[i], err[i]);
    }
#endif

    /* examine "err" vector */
    for (i=0; i<m; i++) {
	if (err[i] == 0.0) {
	    zerocount++;
	} else if (err[i] < 0.35) {
	    badcount++;
	}
    }

    if (zerocount > 0) {
	strcpy(gretl_errmsg, 
	       _("NLS: The supplied derivatives seem to be incorrect"));
	fprintf(stderr, "%d out of %d tests gave zero\n", zerocount, (int) m);
    } else if (badcount > 0) {
	pputs(prn, _("Warning: The supplied derivatives may be incorrect, or perhaps\n"
		     "the data are ill-conditioned for this function.\n"));
	pprintf(prn, _("%d out of %d gradients looked suspicious.\n\n"),
		badcount, (int) m);
    }

 chkderiv_abort:
    free(xp);
    free(err);
    free(fvecp);

    return (zerocount > m/4);
}

/* driver for BFGS code below */

static int mle_calculate (nls_spec *spec, double *fvec, double *jac, PRN *prn)
{
    integer n = spec->nparam;
    int maxit = 200; /* arbitrary? */
    int err = 0;

    if (spec->mode == ANALYTIC_DERIVS) {
	integer m = spec->nobs;
	integer ldjac = m; 

	err = check_derivatives(m, n, spec->coeff, fvec, jac, ldjac, prn);
    }

    if (!err) {
	BFGS_GRAD_FUNC gradfun = (spec->mode == ANALYTIC_DERIVS)?
	    get_mle_gradient : NULL;

	err = BFGS_max(n, spec->coeff, maxit, pspec->tol, 
		       &spec->fncount, &spec->grcount, NULL,
		       get_mle_ll, gradfun, NULL,
		       pspec->opt, nprn);
    }

    return err;    
}

/* driver for minpack levenberg-marquandt code for use when analytical
   derivatives have been supplied */

static int lm_calculate (nls_spec *spec, double *fvec, double *jac, PRN *prn)
{
    integer info, lwa;
    integer m, n, ldjac;
    integer *ipvt;
    doublereal *wa;
    int err = 0;

    m = spec->nobs;              /* number of observations */
    n = spec->nparam;            /* number of parameters */
    lwa = 5 * n + m;             /* work array size */
    ldjac = m;                   /* leading dimension of jac array */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);

    if (wa == NULL || ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    err = check_derivatives(m, n, spec->coeff, fvec, jac, ldjac, prn);
    if (err) {
	goto nls_cleanup; 
    }

    /* call minpack */
    lmder1_(nls_calc, &m, &n, spec->coeff, fvec, jac, &ldjac, &spec->tol, 
	    &info, ipvt, wa, &lwa);

    switch ((int) info) {
    case -1: 
	err = 1;
	break;
    case 0:
	strcpy(gretl_errmsg, _("Invalid NLS specification"));
	err = 1;
	break;
    case 1:
    case 2:
    case 3:
    case 4: /* is this right? */
	pprintf(prn, _("Convergence achieved after %d iterations\n"),
		spec->iters);
	break;
    case 5:
    case 6:
    case 7:
	sprintf(gretl_errmsg, 
		_("NLS: failed to converge after %d iterations"),
		spec->iters);
	err = 1;
	break;
    default:
	break;
    }

 nls_cleanup:

    free(wa);
    free(ipvt);

    return err;    
}

/* callback for lm_approximate (below) to be used by minpack */

static int 
nls_calc_approx (integer *m, integer *n, double *x, double *fvec,
		 integer *iflag)
{
    /* write current parameter values into dataset Z */
    update_nls_param_values(x);

    /* calculate function at x, results into fvec */    
    if (get_nls_fvec(fvec)) {
	*iflag = -1;
    }

    return 0;
}

/* driver for minpack levenberg-marquandt code for use when the
   Jacobian must be approximated numerically */

static int 
lm_approximate (nls_spec *spec, double *fvec, double *jac, PRN *prn)
{
    integer info, m, n, ldjac;
    integer maxfev, mode = 1, nprint = 0, nfev = 0;
    integer iflag = 0;
    integer *ipvt;
    doublereal gtol = 0.0;
    doublereal epsfcn = 0.0, factor = 100.;
    doublereal *diag, *qtf;
    doublereal *wa1, *wa2, *wa3, *wa4;
    int err = 0;
    
    m = spec->nobs;              /* number of observations */
    n = spec->nparam;            /* number of parameters */
    ldjac = m;                   /* leading dimension of jac array */

    maxfev = 200 * (n + 1);

    diag = malloc(n * sizeof *diag);
    qtf = malloc(n * sizeof *qtf);
    wa1 = malloc(n * sizeof *wa1);
    wa2 = malloc(n * sizeof *wa2);
    wa3 = malloc(n * sizeof *wa3);
    wa4 = malloc(m * sizeof *wa4);
    ipvt = malloc(n * sizeof *ipvt);

    if (diag == NULL || qtf == NULL ||
	wa1 == NULL || wa2 == NULL || wa3 == NULL || wa4 == NULL ||
	ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    /* call minpack */
    lmdif_(nls_calc_approx, &m, &n, spec->coeff, fvec, 
	   &spec->tol, &spec->tol, &gtol, &maxfev, &epsfcn, diag, &mode, &factor,
	   &nprint, &info, &nfev, jac, &ldjac, 
	   ipvt, qtf, wa1, wa2, wa3, wa4);

    spec->iters = nfev;

    switch ((int) info) {
    case -1: 
	err = 1;
	break;
    case 0:
	strcpy(gretl_errmsg, _("Invalid NLS specification"));
	err = 1;
	break;
    case 1:
    case 2:
    case 3:
    case 4:
	pprintf(prn, _("Convergence achieved after %d iterations\n"),
		spec->iters);
	break;
    case 5:
    case 6:
    case 7:
    case 8:
	sprintf(gretl_errmsg, 
		_("NLS: failed to converge after %d iterations"),
		spec->iters);
	err = 1;
	break;
    default:
	break;
    }

    if (!err) {
	double ess = spec->ess;
	int iters = spec->iters;
	gretlopt opt = spec->opt;

	spec->opt = OPT_NONE;

	/* call minpack again */
	fdjac2_(nls_calc_approx, &m, &n, spec->coeff, fvec, jac, 
		&ldjac, &iflag, &epsfcn, wa4);
	spec->ess = ess;
	spec->iters = iters;
	spec->opt = opt;
    }

 nls_cleanup:

    free(diag);
    free(qtf);
    free(wa1);
    free(wa2);
    free(wa3);
    free(wa4);
    free(ipvt);

    return err;    
}

/* below: public functions */

/**
 * nls_spec_add_param_with_deriv:
 * @spec: pointer to nls specification.
 * @dstr: string specifying a derivative with respect to a
 *   parameter of the regression function.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Adds an analytical derivative to @spec.  This pointer must
 * have previously been obtained by a call to #nls_spec_new.
 * The required format for @dstr is "%varname = %formula", where
 * %varname is the name of the (scalar) variable holding the parameter
 * in question, and %formula is an expression, of the sort that
 * is fed to gretl's %genr command, giving the derivative of the
 * regression function in @spec with respect to the parameter.
 * The variable holding the parameter must be already present in
 * the dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int 
nls_spec_add_param_with_deriv (nls_spec *spec, const char *dstr,
			       const double **Z, const DATAINFO *pdinfo)
{
    nls_param *param = NULL;
    const char *p = dstr;
    char *vname = NULL;
    int i, v, err = 0;

    if (nlspec_allocate_param(spec)) {
	return E_ALLOC;
    }

    i = spec->nparam - 1;

    param = &spec->params[i];

    if (!strncmp(p, "deriv ", 6)) {
	/* make starting with "deriv" optional */
	p += 6;
    }

    err = equation_get_lhs_and_rhs(p, &vname, &param->deriv);
    if (err) {
	fprintf(stderr, "parse error in deriv string: '%s'\n", dstr);
	return E_PARSE;
    }

    *param->name = '\0';
    strncat(param->name, vname, VNAMELEN - 1);
    free(vname);

    v = varindex(pdinfo, param->name);
    if (v < pdinfo->v) {
	param->varnum = v;
	param->dernum = 0;
	spec->coeff[i] = Z[v][0];
    } else {
	free(param->deriv);
	param->deriv = NULL;
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), param->name);
	err = E_UNKVAR;
    }

    if (!err) {
	spec->mode = ANALYTIC_DERIVS;
    }

#if NLS_DEBUG
    if (param->deriv != NULL) {
	fprintf(stderr, "add_param_with_deriv: '%s'\n"
		" set varnum = %d, initial value = %g\n", dstr, 
		param->varnum, spec->coeff[i]);
    }
#endif

    return err;
}

/**
 * nls_spec_add_aux:
 * @spec: pointer to nls specification.
 * @s: string specifying generation of an auxiliary variable
 * (for use in calculating function or derivatives).
 *
 * Adds the specification of an auxiliary variable to @spec, 
 * which pointer must have previously been obtained by a call 
 * to #nls_spec_new.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nls_spec_add_aux (nls_spec *spec, const char *s)
{
    char **aux;
    char *this;
    int nx = spec->naux + 1;
    int err = 0;

    this = gretl_strdup(s);
    if (this == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	aux = realloc(spec->aux, nx * sizeof *spec->aux);
	if (aux == NULL) {
	    free(this);
	    err = E_ALLOC;
	}
    }

    if (!err) {
	spec->aux = aux;
	spec->aux[nx - 1] = this;
	spec->naux += 1;
    }

    return err;
}

/**
 * nls_spec_set_regression_function:
 * @spec: pointer to nls specification.
 * @fnstr: string specifying nonlinear regression function.
 * @pdinfo: information on dataset.
 *
 * Adds the regression function to @spec.  This pointer must
 * have previously been obtained by a call to #nls_spec_new.
 * The required format for @fnstr is "%varname = %formula", where
 * %varname is the name of the dependent variable and %formula
 * is an expression of the sort that is fed to gretl's %genr command.
 * The dependent variable must be already present in the
 * dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int 
nls_spec_set_regression_function (nls_spec *spec, const char *fnstr, 
				  const DATAINFO *pdinfo)
{    
    const char *p = fnstr;
    char *vname = NULL;
    char *rhs = NULL;
    int flen, err = 0;

    if (spec->nlfunc != NULL) {
	free(spec->nlfunc);
	spec->nlfunc = NULL;
    }

    if (!strncmp(p, "nls ", 4) || !strncmp(p, "mle ", 4)) {
	/* starting with "nls" / "mle" is optional here */
	p += 4;
    }    

    if (equation_get_lhs_and_rhs(p, &vname, &rhs)) { 
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), fnstr);
	err =  E_PARSE;
    } else {
	spec->depvar = varindex(pdinfo, vname);
	if (spec->depvar == pdinfo->v) {
	    if (spec->ci == NLS) {
		sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
		err = E_UNKVAR;
	    } else {
		/* MLE: don't need depvar */
		spec->depvar = 0;
	    }
	}
    }

    if (!err) {
	if (spec->ci == MLE) {
	    flen = strlen(rhs) + 4;
	} else {
	    flen = strlen(vname) + strlen(rhs) + 6;
	}

	spec->nlfunc = malloc(flen);

	if (spec->nlfunc == NULL) {
	    err = E_ALLOC;
	} else {
	    if (spec->ci == MLE) {
		sprintf(spec->nlfunc, "-(%s)", rhs);
	    } else {
		sprintf(spec->nlfunc, "%s - (%s)", vname, rhs);
	    }
	}
    }

    free(vname);
    free(rhs);

    return err;
}

/**
 * nls_spec_set_t1_t2:
 * @spec: pointer to nls specification.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Sets the sample range for estimation of @spec.  This pointer must
 * have previously been obtained by a call to #nls_spec_new.
 */

void nls_spec_set_t1_t2 (nls_spec *spec, int t1, int t2)
{
    if (spec != NULL) {
	spec->t1 = t1;
	spec->t2 = t2;
    }
}

#define genr_line(s) (!strncmp(s, "series", 6) || \
                      !strncmp(s, "scalar", 6) || \
                      !strncmp(s, "genr", 4))

#define param_line(s) (!strncmp(s, "deriv", 5) || \
                       !strncmp(s, "params", 6))

/**
 * nls_parse_line:
 * @ci: either %NLS or %MLE (docs not finished on this)
 * @line: specification of regression function or derivative
 *        of this function with respect to a parameter.
 * @Z: data array.
 * @pdinfo: information on dataset.
 * @prn: gretl printing struct (for warning messages).
 *
 * This function is used to create the specification of a
 * nonlinear regression function, to be estimated via #nls.
 * It should first be called with a @line containing a
 * string specification of the regression function.  Optionally,
 * it can then be called one or more times to specify 
 * analytical derivatives for the parameters of the regression
 * function.  
 *
 * The format of @line should be that used in the #genr function.
 * When specifying the regression function, the formula may
 * optionally be preceded by the string %nls.  When specifying
 * a derivative, @line must start with the string %deriv.  See
 * the gretl manual for details.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int nls_parse_line (int ci, const char *line, const double **Z,
		    const DATAINFO *pdinfo, PRN *prn)
{
    int err = 0;

    pspec = &private_spec;
    pspec->ci = ci;

    if (genr_line(line)) {
	err = nls_spec_add_aux(pspec, line);
    } else if (param_line(line)) {
	if (pspec->nlfunc == NULL) {
	    strcpy(gretl_errmsg, _("No regression function has been specified"));
	    err = E_PARSE;
	} else {
	    if (*line == 'd') {
		/* "deriv" */
		if (pspec->mode != ANALYTIC_DERIVS && pspec->params != NULL) {
		    strcpy(gretl_errmsg, "You cannot supply both a \"params\" "
			   "line and analytical derivatives");
		    err = E_PARSE;
		} else {
		    err = nls_spec_add_param_with_deriv(pspec, line, Z, pdinfo);
		}
	    } else {
		/* "params" */
		if (pspec->mode != ANALYTIC_DERIVS) {
		    err = nls_spec_add_params_from_line(pspec, line + 6, Z, pdinfo);
		} else {
		    pprintf(prn, _("Analytical derivatives supplied: "
				   "\"params\" line will be ignored"));
		    pputc(prn, '\n');
		}
	    }
	}
    } else {
	/* do we already have an nls specification under way? */
	if (pspec->nlfunc != NULL) {
	    clear_nls_spec(pspec);
	}
	err = nls_spec_set_regression_function(pspec, line, pdinfo);
	if (!err) {
	    nls_spec_set_t1_t2(pspec, pdinfo->t1, pdinfo->t2);
	}
    }

    return err;
}

static double default_nls_toler;

/**
 * get_default_nls_toler:
 *
 * Returns: the default value used in the convergence criterion
 * for estimation of models using nonlinear least squares.
 */

double get_default_nls_toler (void)
{
    if (default_nls_toler == 0.0) {
	default_nls_toler = pow(dpmpar_(&one), .75);
    }

    return default_nls_toler;
}

static void 
save_likelihood_vector (nls_spec *spec, double ***pZ, DATAINFO *pdinfo)
{
    int t;

    for (t=0; t<pdinfo->n; t++) {
	if (t < spec->t1 || t > spec->t2) {
	    (*pZ)[spec->depvar][t] = NADBL;
	} else {
	    (*pZ)[spec->depvar][t] = - (*pZ)[spec->uhatnum][t];
	}
    }
}

/* static function providing the real content for the two public
   wrapper functions below */

static MODEL real_nls (nls_spec *spec, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    MODEL nlsmod;
    double *fvec = NULL;
    double *jac = NULL;
    int origv = pdinfo->v;
    int i, t, err = 0;

    genr_err = 0;
    gretl_model_init(&nlsmod);
    gretl_model_smpl_init(&nlsmod, pdinfo);

    if (spec != NULL) {
	/* the caller supplied an nls specification directly */
	pspec = spec;
    } else {
	/* we use the static spec composed via nls_parse_line() */
	pspec = &private_spec;
    }

    if (pspec->nlfunc == NULL) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	nlsmod.errcode = E_PARSE;
	goto bailout;
    } 

    pspec->opt = opt;

    /* publish pZ, pdinfo and prn */
    nZ = pZ;
    ndinfo = pdinfo;
    nprn = prn;

    if (pspec->mode == NUMERIC_DERIVS && pspec->nparam == 0) {
	err = get_params_from_nlfunc(pspec, (const double **) *pZ, pdinfo);
	if (err) {
	    if (err == 1) {
		nlsmod.errcode = E_PARSE;
	    } else {
		nlsmod.errcode = err;
	    }
	    goto bailout;
	}
    }

    if (pspec->nparam == 0) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	nlsmod.errcode = E_PARSE;
	goto bailout;
    } 

    err = nls_missval_check(pspec);
    if (err) {
	nlsmod.errcode = err;
	*pZ = *nZ;
	goto bailout;
    }

    /* allocate arrays to be passed to minpack */
    fvec = malloc(pspec->nobs * sizeof *fvec);
    jac = malloc(pspec->nobs * pspec->nparam * sizeof *jac);

    if (fvec == NULL || jac == NULL) {
	nlsmod.errcode = E_ALLOC;
	*pZ = *nZ;
	goto bailout;
    }

    i = 0;
    for (t=pspec->t1; t<=pspec->t2; t++) {
	fvec[i++] = (*nZ)[pspec->uhatnum][t];
    }

    /* get tolerance from user setting or default */
    pspec->tol = get_nls_toler();

    pputs(prn, (pspec->mode == NUMERIC_DERIVS)?
	  _("Using numerical derivatives\n") :
	  _("Using analytical derivatives\n"));

    if (pspec->ci == MLE) {
	err = mle_calculate(pspec, fvec, jac, prn);
    } else {
	/* invoke appropriate minpack driver function */
	if (pspec->mode == NUMERIC_DERIVS) {
	    err = lm_approximate(pspec, fvec, jac, prn);
	} else {
	    err = lm_calculate(pspec, fvec, jac, prn);
	    if (err) {
		fprintf(stderr, "lm_calculate returned %d\n", err);
	    }
	}
    }

    /* re-attach data array pointer: may have moved! */
    *pZ = *nZ;

    pprintf(prn, _("Tolerance = %g\n"), pspec->tol);

    if (!err) {
	if (pspec->ci == NLS) {
	    if (pspec->opt & OPT_C) {
		/* coefficients only: don't bother with GNR */
		add_nls_coeffs(&nlsmod, pspec);
	    } else {
		/* Use Gauss-Newton Regression for covariance matrix,
		   standard errors */
		nlsmod = GNR(fvec, jac, pspec, (const double **) *pZ, 
			     pdinfo, prn);
	    }
	} else {
	    make_mle_model(&nlsmod, pspec, pdinfo);
	    if (pspec->depvar > 0) {
		save_likelihood_vector(pspec, pZ, pdinfo);
	    }
	}
    } else {
	if (nlsmod.errcode == 0) { 
	    if (genr_err != 0) {
		nlsmod.errcode = genr_err;
	    } else {
		nlsmod.errcode = E_NOCONV;
	    }
	}
    }

 bailout:

    free(fvec);
    free(jac);

    if (spec == NULL) {
	clear_nls_spec(pspec);
    }

    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    if (nlsmod.errcode == 0 && !(opt & OPT_A)) {
	set_model_id(&nlsmod);
    }

    return nlsmod;
}

/**
 * nls:
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @opt: may include %OPT_V for verbose output, %OPT_R
 * for robust covariance matrix.
 * @prn: printing struct.
 *
 * Computes estimates of a model via nonlinear least squares.
 * The model must have been specified previously, via calls to
 * the function #nls_parse_line.  
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL nls (double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    return real_nls(NULL, pZ, pdinfo, opt, prn);
}

/**
 * model_from_nls_spec:
 * @spec: nls specification.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @opt: may include %OPT_V for verbose output, %OPT_A to
 * treat as an auxiliary model, %OPT_C to produce coefficient
 * estimates only (don't bother with GNR to produce standard
 * errors).
 * @prn: printing struct.
 *
 * Computes estimates of the model specified in @spec, via nonlinear 
 * least squares. The @spec must first be obtained using #nls_spec_new, and
 * initialized using #nls_spec_set_regression_function.  If analytical
 * derivatives are to be used (which is optional but recommended)
 * these are set using #nls_spec_add_param_with_deriv.
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL model_from_nls_spec (nls_spec *spec, double ***pZ, DATAINFO *pdinfo, 
			   gretlopt opt, PRN *prn)
{
    return real_nls(spec, pZ, pdinfo, opt, prn);
}

/**
 * nls_spec_new:
 * @ci: either %NLS or %MLE.
 * @pdinfo: information on dataset.
 *
 * Returns: a pointer to a newly allocated nls specification,
 * or %NULL on failure.
 */

nls_spec *nls_spec_new (int ci, const DATAINFO *pdinfo)
{
    nls_spec *spec;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->nlfunc = NULL;

    spec->params = NULL;
    spec->nparam = 0;

    spec->aux = NULL;
    spec->naux = 0;
    
    spec->genrs = NULL;
    spec->ngenrs = 0;
    
    spec->coeff = NULL;

    spec->ci = ci;
    spec->mode = NUMERIC_DERIVS;
    spec->opt = OPT_NONE;

    spec->nparam = 0;
    spec->depvar = 0;

    spec->iters = 0;
    spec->fncount = 0;
    spec->grcount = 0;

    spec->t1 = pdinfo->t1;
    spec->t2 = pdinfo->t2;
    spec->nobs = 0;

    return spec;
}

/**
 * nls_spec_destroy:
 * @spec: pointer to nls specification.
 *
 * Frees all resources associated with @spec, and frees the
 * pointer itself.
 */

void nls_spec_destroy (nls_spec *spec)
{
    clear_nls_spec(spec);
    free(spec);
}

static double **triangular_array_new (int n)
{
    double **m;
    int i;

    m = malloc(n * sizeof *m);

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    m[i] = malloc((i + 1) * sizeof **m);
	}
    }

    return m;
}

static void free_triangular_array (double **m, int n)
{
    int i;

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    free(m[i]);
	}
	free(m);
    }
}

#define BFGS_DEBUG 0

#define stepfrac	0.2
#define acctol		0.0001 
#define reltest		10.0

/* FIXME: it would be nice to trash this, but that depends on
   switching all the signs correctly below, which (to me!) is
   more difficult than it looks -- AC */

static void reverse_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
	g[i] = -g[i];
    }
}

/* apparatus for constructing numerical approximation to
   the Hessian */

static void hess_h_init (double *h, double *h0, int n)
{
    int i;

    for (i=0; i<n; i++) {
	h[i] = h0[i];
    }
}

static void hess_h_reduce (double *h, double v, int n)
{
    int i;

    for (i=0; i<n; i++) {
	h[i] /= v;
    }
}

static void hess_b_adjust_i (double *c, double *b, double *h, int n, 
			     int i, double sgn)
{
    int k;

    for (k=0; k<n; k++) {
	c[k] = b[k] + (k == i) * sgn * h[i];
    }
}

static void hess_b_adjust_ij (double *c, double *b, double *h, int n, 
			      int i, int j, double sgn)
{
    int k;

    for (k=0; k<n; k++) {
	c[k] = b[k] + (k == i) * sgn * h[i] +
	    (k == j) * sgn * h[j];
    }
}

/* The algorithm below implements the method of Richardson
   Extrapolation.  It is derived from code in the gnu R package
   "numDeriv" by Paul Gilbert, which was in turn derived from C code
   by Xinqiao Liu.  Turned back into C and modified for gretl by
   Allin Cottrell, June 2006.
*/

static double *numerical_hessian (double *b, double *c, int n,
				  BFGS_LL_FUNC func, void *data)
{
    double *D = NULL;
    double *h0 = NULL;
    double *h = NULL;
    double *Dx = NULL;
    double *Hx = NULL;
    double *Hd = NULL;

    gretl_matrix *V = NULL;
    double *vcv = NULL;

    /* numerical parameters */
    double eps = 1.0e-4;
    double d = 0.0001;
    double v = 2.0;      /* reduction factor for h */
    int r = 4;           /* number of Richardson steps */

    double f0, f1, f2;
    double p4m;

    int vn = (n * (n + 1)) / 2;
    int dn = vn + n;
    int i, j, k, m, u;
    int err = 0;

    h0 = malloc(n * sizeof *h0);
    h  = malloc(n * sizeof *h);
    Dx = malloc(r * sizeof *Dx);
    Hx = malloc(r * sizeof *Hx);
    Hd = malloc(n * sizeof *Hd);
    D  = malloc(dn * sizeof *D);

    if (h0 == NULL || h == NULL || Dx == NULL || 
	Hx == NULL || Hd == NULL || D == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* vech form of variance matrix */
    V = gretl_column_vector_alloc(vn);
    if (V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }	

    for (i=0; i<n; i++) {
	h0[i] = d * b[i] + eps * (b[i] == 0.0);
    }

    f0 = func(b, data);

    /* first derivatives and Hessian diagonal */

    for (i=0; i<n; i++) {
	hess_h_init(h, h0, n);
	for (k=0; k<r; k++) {
	    hess_b_adjust_i(c, b, h, n, i, 1.0);
	    f1 = func(c, data);
	    hess_b_adjust_i(c, b, h, n, i, -1.0);
	    f2 = func(c, data);
	    /* F'(i) */
	    Dx[k] = (f1 - f2) / (2.0 * h[i]); 
	    /* F''(i) */
	    Hx[k] = (f1 - 2.0 * f0 + f2) / (h[i] * h[i]);
	    hess_h_reduce(h, v, n);
	}
	p4m = 4.0;
	for (m=0; m<r-1; m++) {
	    for (k=0; k<r-m; k++) {
		Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1.0);
		Hx[k] = (Hx[k+1] * p4m - Hx[k]) / (p4m - 1.0);
	    }
	    p4m *= 4.0;
	}
	D[i] = Dx[0];
	Hd[i] = Hx[0];
    }

    /* second derivatives: lower half of Hessian only */

    u = n;
    for (i=0; i<n; i++) {
	for (j=0; j<=i; j++) {
	    if (i == j) {
		D[u] = Hd[i];
	    } else {
		hess_h_init(h, h0, n);
		for (k=0; k<r; k++) {
		    hess_b_adjust_ij(c, b, h, n, i, j, 1.0);
		    f1 = func(c, data);
		    hess_b_adjust_ij(c, b, h, n, i, j, -1.0);
		    f2 = func(c, data);
		    /* cross-partial */
		    Dx[k] = (f1 - 2.0 * f0 + f2 - Hd[i] * h[i] * h[i]
			     - Hd[j] * h[j] * h[j]) / (2.0 * h[i] * h[j]);
		    hess_h_reduce(h, v, n);
		}
		p4m = 4.0;
		for (m=0; m<r-1; m++) {
		    for (k=0; k<r-m; k++) {
			Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1.0);
		    }
		    p4m *= 4.0;
		}
		D[u] = Dx[0];
	    }
	    u++;
	}
    }

    /* transcribe the negative of the Hessian */
    u = n;
    for (i=0; i<n; i++) {
	for (j=0; j<=i; j++) {
	    k = ijton(i, j, n);
	    V->val[k] = -D[u++];
	}
    }

    err = gretl_invert_packed_symmetric_matrix(V);
    if (!err) {
	vcv = gretl_matrix_steal_data(V);
    }  

    gretl_matrix_free(V);

 bailout:

    free(D);
    free(h0);
    free(h);
    free(Dx);
    free(Hx);
    free(Hd);

    return vcv;
}

/**
 * BFGS_max:
 * @n: number elements in array @b.
 * @b: array of adjustable coefficients.
 * @maxit: the maximum number of iterations to allow.
 * @reltol: relative tolerance for terminating iteration.
 * @fncount: location to receive count of function evaluations.
 * @grcount: location to receive count of gradient evaluations.
 * @hessvcv: location to receive packed triangular representation
 * of covariance matrix, based on numerical approximation to the
 * Hessian at the last iteration, or %NULL.
 * @get_ll: pointer to function used to calculate log
 * likelihood.
 * @get_gradient: pointer to function used to calculate the 
 * gradient, or %NULL for default numerical calculation.
 * @callback_data: pointer that will be passed as the last
 * parameter to the callback functions @get_ll and @get_gradient.
 * @opt: may contain %OPT_V for verbose operation.
 * @prn: printing struct (or %NULL).
 *
 * Obtains the set of values for @b which jointly maximize the
 * log-likelihood as calculated by @get_ll.  Uses the BFGS
 * variable-metric method.  Based on Pascal code in J. C. Nash,
 * "Compact Numerical Methods for Computers," 2nd edition, converted
 * by p2c then re-crafted by B. D. Ripley for gnu R.  Revised for 
 * gretl by Allin Cottrell.
 * 
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int BFGS_max (int n, double *b, int maxit, double reltol,
	      int *fncount, int *grcount, double **hessvcv,
	      BFGS_LL_FUNC get_ll, BFGS_GRAD_FUNC get_gradient, 
	      void *callback_data, gretlopt opt, PRN *prn)
{
    int ll_ok, done;
    double *g = NULL, *t = NULL, *X = NULL, *c = NULL, **H = NULL;
    int ndelta, fcount, gcount;
    double fmax, f, sumgrad;
    int i, j, ilast, iter;
    double s, steplen = 0.0;
    double D1, D2;
    int err = 0;

    if (get_gradient == NULL) {
	get_gradient = BFGS_numeric_gradient;
    }

    g = malloc(n * sizeof *g);
    t = malloc(n * sizeof *t);
    X = malloc(n * sizeof *X);
    c = malloc(n * sizeof *c);
    H = triangular_array_new(n);

    if (g == NULL || t == NULL || X == NULL || c == NULL || H == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    f = get_ll(b, callback_data);

    if (na(f)) {
	fprintf(stderr, "initial value of f is not finite\n");
	err = E_DATA;
	goto bailout;
    }

    fmax = f;
    iter = ilast = fcount = gcount = 1;
    get_gradient(b, g, n, get_ll, callback_data);
    reverse_gradient(g, n);

    do {
	if (opt & OPT_V) {
	    print_mle_iter_stats(f, n, b, g, iter, steplen, prn);
	}
	if (ilast == gcount) {
	    /* (re-)start: initialize curvature matrix */
	    for (i=0; i<n; i++) {
		for (j=0; j<i; j++) {
		    H[i][j] = 0.0;
		}
		H[i][i] = 1.0;
	    }
	}
	for (i=0; i<n; i++) {
	    /* copy coefficients to X, gradient to c */
	    X[i] = b[i];
	    c[i] = g[i];
	}
	sumgrad = 0.0;
	for (i=0; i<n; i++) {
	    s = 0.0;
	    for (j=0; j<=i; j++) {
		s -= H[i][j] * g[j];
	    }
	    for (j=i+1; j<n; j++) {
		s -= H[j][i] * g[j];
	    }
	    t[i] = s;
	    sumgrad += s * g[i];
	}

	if (sumgrad < 0.0) {	
	    /* search direction is uphill, actually */
	    steplen = 1.0;
	    ll_ok = 0;
	    do {
		/* loop so long as (a) we haven't achieved a definite
		   improvement in the log-likelihood and (b) there is
		   still some prospect of doing so */
		ndelta = n;
		for (i=0; i<n; i++) {
		    b[i] = X[i] + steplen * t[i];
		    if (reltest + X[i] == reltest + b[i]) {
			/* no change in coefficient */
			ndelta--;
		    }
		}
		if (ndelta > 0) {
		    f = get_ll(b, callback_data);
		    fcount++;
		    ll_ok = !na(f) &&
			(f >= fmax + sumgrad * steplen * acctol);
		    if (!ll_ok) {
			/* loglik not good: try smaller step */
			steplen *= stepfrac;
		    }
		}
	    } while (ndelta != 0 && !ll_ok);

	    done = fabs(fmax - f) <= reltol * (fabs(fmax) + reltol);

#if BFGS_DEBUG
	    fprintf(stderr, "LHS=%g, RHS=%g; done = %d\n",
		    fabs(fmax - f), reltol * (fabs(fmax) + reltol),
		    done);
#endif

	    /* prepare to stop if relative change is small enough */
	    if (done) {
		ndelta = 0;
		fmax = f;
	    }

	    if (ndelta > 0) {
		/* making progress */
		fmax = f;
		get_gradient(b, g, n, get_ll, callback_data);
		reverse_gradient(g, n);
		gcount++;
		iter++;
		D1 = 0.0;
		for (i=0; i<n; i++) {
		    t[i] = steplen * t[i];
		    c[i] = g[i] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0.0) {
		    D2 = 0.0;
		    for (i=0; i<n; i++) {
			s = 0.0;
			for (j=0; j<=i; j++) {
			    s += H[i][j] * c[j];
			}
			for (j=i+1; j<n; j++) {
			    s += H[j][i] * c[j];
			}
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i=0; i<n; i++) {
			for (j=0; j<=i; j++) {
			    H[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
			}
		    }
		} else {
		    /* D1 <= 0.0 */
		    ilast = gcount;
		}
	    } else if (ilast < gcount) {
		ndelta = n;
		ilast = gcount;
	    }
	} else {
	    /* downhill search */
	    if (ilast == gcount) {
		/* we just reset: don't reset again; set ndelta = 0 so
		   that we exit the main loop
		*/
		ndelta = 0;
	    } else {
		/* reset for another attempt */
		ilast = gcount;
		ndelta = n;
	    }
	}

	if (iter >= maxit) {
	    break;
	}

	if (gcount - ilast > 2 * n) {
	    /* periodic restart of curvature computation */
	    ilast = gcount;
	}

    } while (ndelta > 0 || ilast < gcount);

#if BFGS_DEBUG
    fprintf(stderr, "terminated: ndelta=%d, ilast=%d, gcount=%d\n",
	    ndelta, ilast, gcount);
#endif

    if (iter >= maxit) {
	fprintf(stderr, "stopped after %d iterations\n", iter);
	err = E_NOCONV;
    }

    *fncount = fcount;
    *grcount = gcount;

    if (!err && hessvcv != NULL) {
	*hessvcv = numerical_hessian(b, c, n, get_ll, callback_data);
    }

    if (opt & OPT_V) {
	pputs(nprn, "\n--- FINAL VALUES: \n");	
	print_mle_iter_stats(f, n, b, g, iter, steplen, prn);
	pputs(nprn, "\n\n");	
    }

 bailout:

    free(g);
    free(t);
    free(X);
    free(c);
    free_triangular_array(H, n);

    return err;
}

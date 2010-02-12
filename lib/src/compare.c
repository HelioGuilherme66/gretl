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

/*
    compare.c - gretl model comparison procedures
*/

#include "libgretl.h"
#include "libset.h"
#include "gretl_panel.h"
#include "var.h"
#include "system.h"
#include "missing_private.h"
#include "matrix_extra.h"
#include "plotspec.h"
#include "tsls.h"

#define WDEBUG 0

#define print_add_omit_model(m,o) (!(opt & OPT_Q) && !(opt & OPT_I))

enum {
    CHISQ_FORM,
    F_FORM
};

enum {
    COMPARE_ADD,
    COMPARE_OMIT,
    COMPARE_OMIT_WALD,
    COMPARE_ADD_WALD
};

struct COMPARE {
    int cmd;       /* ADD or OMIT */
    int m1;        /* ID for first model */
    int m2;        /* ID for second model */
    int ci;        /* estimator code for the first model */
    int dfn;       /* numerator degrees of freedom */
    int dfd;       /* denominator degrees of freedom */ 
    double F;      /* F test statistic */
    double wald;   /* Wald chi-square test statistic */
    double trsq;   /* T*R^2 test statistic */
    double LR;     /* likelihood ratio test statistic */
    int score;     /* number of info stats showing improvement */
    int robust;    /* = 1 when robust vcv is in use, else 0 */
    int err;       /* error code */
};

/* Critical values for Quandt likelihood ratio (break) test:
   columns contain the 10%, 5% and 1% values respectively,
   for the given number of restrictions, q, and 15% trimming.  
   Taken from Stock and Watson, Introduction to Econometrics
   (Addison-Wesley, 2003), based on Andrews.
*/

#define QLR_QMAX 20

static double QLR_critvals[QLR_QMAX][3] = {
    { 7.12, 8.68, 12.16 }, /*  1 */
    { 5.00, 5.86,  7.78 }, /*  2 */
    { 4.09, 4.71,  6.02 }, /*  3 */
    { 3.59, 4.09,  5.12 }, /*  4 */
    { 3.26, 3.66,  4.53 }, /*  5 */
    { 3.02, 3.37,  4.12 }, /*  6 */
    { 2.84, 3.15,  3.82 }, /*  7 */
    { 2.69, 2.98,  3.57 }, /*  8 */
    { 2.58, 2.84,  3.38 }, /*  9 */
    { 2.48, 2.71,  3.23 }, /* 10 */
    { 2.40, 2.62,  3.09 }, /* 11 */
    { 2.33, 2.54,  2.97 }, /* 12 */
    { 2.27, 2.46,  2.87 }, /* 13 */
    { 2.21, 2.40,  2.78 }, /* 14 */
    { 2.16, 2.34,  2.71 }, /* 15 */
    { 2.12, 2.29,  2.64 }, /* 16 */
    { 2.08, 2.25,  2.58 }, /* 17 */
    { 2.05, 2.20,  2.53 }, /* 18 */
    { 2.01, 2.17,  2.48 }, /* 19 */
    { 1.99, 2.13,  2.43 }  /* 20 */
};

/* Given a list of variables, check them against the independent
   variables included in a model, and construct a mask with 1s in
   positions where there is a match, 0s otherwise.  If the test list
   is NULL, match all variables except the constant.
*/

static char *
mask_from_test_list (const int *list, const MODEL *pmod, int *err)
{
    char *mask;
    int off1 = 2, off2 = 0;
    int nmask = 0;
    int i, j;

    mask = calloc(pmod->ncoeff, 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (pmod->ci == ARBOND) {
	/* find correct offset into independent vars in list */
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == LISTSEP) {
		off1 = i + 2;
	    }
	}
	off2 = pmod->list[1];
    }

    for (i=0; i<pmod->ncoeff; i++) {
	if (list != NULL) {
	    for (j=1; j<=list[0]; j++) {
		if (pmod->list[i + off1] == list[j]) {
#if WDEBUG
		    fprintf(stderr, "matched var %d at pmod->list[%d]: set mask[%d] = 1\n", 
			    list[j], i + off1, i + off2);
#endif
		    mask[i + off2] = 1;
		    nmask++;
		}
	    }
	} else if (pmod->list[i + off1] != 0) {
	    mask[i + off2] = 1;
	}
    }

    if (list != NULL && nmask != list[0]) {
	fprintf(stderr, "mask from list: list[0] = %d but nmask = %d\n",
		list[0], nmask);
	*err = E_DATA;
    }

    return mask;
}

/* Wald (chi-square and/or F) test for a set of zero restrictions on
   the parameters of a given model, based on the covariance matrix of
   the unrestricted model. Suitable for use where the original model
   is estimated by FGLS or IV.  Note that if list is NULL, we do an
   automatic test, for the significance of all vars but the constant.
*/

static int 
wald_test (const int *list, MODEL *pmod, double *chisq, double *F)
{
    char *mask = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    double wX = NADBL;
    double wF = NADBL;
    int err = 0;

    mask = mask_from_test_list(list, pmod, &err);
    if (err) {
	return err;
    }

    if (!err) {
	C = gretl_vcv_matrix_from_model(pmod, mask, &err);
    }

    if (!err) {
	b = gretl_coeff_vector_from_model(pmod, mask, &err);
    }  

    if (!err) {
#if WDEBUG
	gretl_matrix_print(C, "Wald VCV matrix");
	gretl_matrix_print(b, "Wald coeff vector");
#endif
	err = gretl_invert_symmetric_matrix(C);
#if WDEBUG
	gretl_matrix_print(C, "VCV-inverse");
#endif
    }

    if (!err) {
	wX = gretl_scalar_qform(b, C, &err);
    }

#if WDEBUG
    fprintf(stderr, "wX (quad form) = %g\n", wX);
#endif

    if (!err) {
	wF = wX / gretl_vector_get_length(b);
    }

    if (!err) {
	if (chisq != NULL) *chisq = wX;
	if (F != NULL) *F = wF;
    }

#if WDEBUG
    fprintf(stderr, "Wald test: F = %g, Chi^2 = %g\n", wF, wX);
#endif

    free(mask);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return err;
}

/**
 * wald_omit_F:
 * @list: list of variables to omit (or NULL).
 * @pmod: model to be tested.
 *
 * Simple form of Wald F-test for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated F-value, or #NADBL on failure.
 */

double wald_omit_F (const int *list, MODEL *pmod)
{
    double F = NADBL;

    wald_test(list, pmod, NULL, &F);
    return F;
}

/**
 * wald_omit_chisq:
 * @list: list of variables to omit (or NULL).
 * @pmod: model to be tested.
 *
 * Simple form of Wald chi-square for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated chi-square value, or #NADBL on failure.
 */

double wald_omit_chisq (const int *list, MODEL *pmod)
{
    double X = NADBL;

    wald_test(list, pmod, &X, NULL);
    return X;
}

/* ----------------------------------------------------- */

static int 
add_diffvars_to_test (ModelTest *test, const int *list, 
		      const DATAINFO *pdinfo)
{
    char *vnames;
    int i, len = 0;
    int err = 0;

    for (i=1; i<=list[0]; i++) {
	len += strlen(pdinfo->varname[list[i]]) + 1;
    }

    vnames = malloc(len);

    if (vnames == NULL) {
	err = 1;
    } else {
	*vnames = '\0';
	for (i=1; i<=list[0]; i++) {
	    strcat(vnames, pdinfo->varname[list[i]]);
	    if (i < list[0]) {
		strcat(vnames, " ");
	    }
	}
	model_test_set_allocated_param(test, vnames);
    }

    return err;
}

static void
gretl_make_compare (const struct COMPARE *cmp, const int *diffvars, 
		    MODEL *newmod, const DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn)
{
    ModelTest *test = NULL;
    const char *tstr = NULL;
    int tcode = 0, statcode = 0;
    double testval = NADBL;
    double pval = NADBL;
    int verbosity = 2;
    int i, stat_ok = 1;

    if (cmp->err) {
	return;
    }

    if (opt & OPT_I) {
	verbosity = 0;
    } else if (opt & OPT_Q) {
	verbosity = 1;
    }

    if (!na(cmp->wald)) {
	testval = cmp->wald;
	statcode = GRETL_STAT_WALD_CHISQ;
    } else if (!na(cmp->LR)) {
	testval = cmp->LR;
	statcode = GRETL_STAT_LR;
    } else if (!na(cmp->F)) {
	testval = cmp->F;
	statcode = GRETL_STAT_F;
    } else {
	stat_ok = 0;
    }

    if (!stat_ok && (!verbosity || cmp->score < 0)) {
	/* nothing to show */
	return;
    }

    if (cmp->cmd == OMIT) {
	tcode = GRETL_TEST_OMIT;
	tstr = N_("omit");
    } else {
	tcode = GRETL_TEST_ADD;
	tstr = N_("add");
    }

    if (stat_ok && (opt & OPT_S)) {
	test = model_test_new(tcode);
    }

    /* preamble */

    if (verbosity > 1 && cmp->m1 >= 0 && cmp->m2 >= 0) {
	pprintf(prn, _("Comparison of Model %d and Model %d:\n"), 
		cmp->m1, cmp->m2);
    } 

    if (verbosity > 0 && stat_ok) {
	if (diffvars[0] == 1) {
	    pputs(prn, "\n  ");
	    pprintf(prn, _("Null hypothesis: the regression parameter is zero for %s"), 
		    pdinfo->varname[diffvars[1]]);
	    pputc(prn, '\n');
	} else {
	    const char *vname;
	    int nc = 0;

	    pputs(prn, _("\n  Null hypothesis: the regression parameters are "
			 "zero for the variables\n"));
	    pputs(prn, "    ");
	    for (i=1; i<=diffvars[0]; i++) {
		vname = pdinfo->varname[diffvars[i]];
		nc += strlen(vname) + 2;
		pprintf(prn, "%s", vname);
		if (i < diffvars[0]) {
		    if (nc > 60) {
			pputs(prn, ",\n    ");
			nc = 0;
		    } else {
			pputs(prn, ", ");
		    }
		}
	    }
	    pputc(prn, '\n');
	}
    }

    /* computation of p-value */

    if (statcode == GRETL_STAT_WALD_CHISQ || statcode == GRETL_STAT_LR) {
	pval = chisq_cdf_comp(cmp->dfn, testval);
    } else if (statcode == GRETL_STAT_F) {
	pval = snedecor_cdf_comp(cmp->dfn, cmp->dfd, cmp->F);
    }

    /* printout, if wanted */

    if (verbosity > 0) {
	if (diffvars[0] > 1 || statcode == GRETL_STAT_WALD_CHISQ) {
	    pputc(prn, '\n');
	}
	if (statcode == GRETL_STAT_WALD_CHISQ) {
	    pprintf(prn, "  %s:\n    %s(%d) = %g, ",  _("Asymptotic test statistic"),
		    _("Wald chi-square"), cmp->dfn, testval);
	    pprintf(prn, _("with p-value = %g\n"), pval);
	    if (verbosity <= 1) {
		pputc(prn, '\n');
	    }
	    if (!na(cmp->F)) {
		/* alternate F-form for Wald */
		pval = snedecor_cdf_comp(cmp->dfn, cmp->dfd, cmp->F);
		pputs(prn, "    ");
		pprintf(prn, "%s: F(%d, %d) = %g, ", _("F-form"), 
			cmp->dfn, cmp->dfd, cmp->F);
		pprintf(prn, _("with p-value = %g\n"), pval);
	    }
	} else if (statcode == GRETL_STAT_LR) {
	    pprintf(prn, "  %s:\n    %s(%d) = %g, ",  _("Likelihood ratio test"),
		    _("Chi-square"), cmp->dfn, cmp->LR);
	    pprintf(prn, _("with p-value = %g\n"), pval);
	    if (verbosity <= 1) {
		pputc(prn, '\n');
	    }
	} else if (statcode == GRETL_STAT_F) {
	    pprintf(prn, "  %s: %s(%d, %d) = %g, ", _("Test statistic"), 
		    (cmp->robust)? _("Robust F") : "F",
		    cmp->dfn, cmp->dfd, cmp->F);
	    pprintf(prn, _("with p-value = %g\n"), pval);
	    if (verbosity <= 1) {
		pputc(prn, '\n');
	    }
	}
    }

    /* record test result */

    if (stat_ok) {
	record_test_result(testval, pval, _(tstr));
	if (test != NULL) {
	    model_test_set_teststat(test, statcode);
	    model_test_set_dfn(test, cmp->dfn);
	    if (statcode == GRETL_STAT_F) {
		model_test_set_dfd(test, cmp->dfd);
	    }
	    model_test_set_value(test, testval);
	    model_test_set_pvalue(test, pval);
	}
    }

    if (verbosity > 0 && statcode == GRETL_STAT_WALD_CHISQ && !na(cmp->LR)) {
	/* print auxiliary LR test if available */
	pprintf(prn, "\n  %s:\n    %s(%d) = %g, ",  _("Likelihood ratio test"),
		_("Chi-square"), cmp->dfn, cmp->LR);
	pprintf(prn, _("with p-value = %g\n"), chisq_cdf_comp(cmp->dfn, cmp->LR));
	if (verbosity <= 1) {
	    pputc(prn, '\n');
	}
    }	

    if (test != NULL) {
	add_diffvars_to_test(test, diffvars, pdinfo);
	maybe_add_test_to_model(newmod, test);
    }

    if (verbosity > 1) {
	if (cmp->score >= 0) {
	    pprintf(prn, _("  Of the %d model selection statistics, %d "), 
		    C_MAX, cmp->score);
	    if (cmp->score == 1) {
		pputs(prn, _("has improved.\n"));
		pputc(prn, '\n');
	    } else {
		pputs(prn, _("have improved.\n\n"));
	    }
	} else {
	    pputc(prn, '\n');
	}
    }
}

static struct COMPARE 
add_or_omit_compare (MODEL *pmodA, MODEL *pmodB, int flag,
		     const int *testvars)
{
    struct COMPARE cmp;
    MODEL *umod, *rmod;
    int i;	

    cmp.err = 0;

    if (flag == COMPARE_ADD || flag == COMPARE_ADD_WALD) {
	umod = pmodB;
	rmod = pmodA;
	cmp.cmd = ADD;
    } else {
	umod = pmodA;
	rmod = pmodB;
	cmp.cmd = OMIT;
    }

    cmp.F = cmp.wald = cmp.trsq = cmp.LR = NADBL;
    cmp.score = -1;
    cmp.robust = 0;
    cmp.dfn = testvars[0];

    if (pmodA != NULL) {
	/* may be NULL for Wald add test */
	cmp.m1 = pmodA->ID;
	cmp.ci = pmodA->ci;
    } else {
	cmp.m1 = -1;
	cmp.ci = pmodB->ci;
    }

    if (pmodB != NULL) {
	/* may be NULL for Wald omit test */
	cmp.m2 = pmodB->ID;
	if (pmodA != NULL) {
	    cmp.dfn = umod->ncoeff - rmod->ncoeff;
	}
    } else {
	cmp.m2 = -1;
    }

    cmp.dfd = umod->dfd;

    if (flag == COMPARE_OMIT_WALD || flag == COMPARE_ADD_WALD) {
	cmp.err = wald_test(testvars, umod, &cmp.wald, &cmp.F);
    } else if (pmodA != NULL && (pmodA->opt & OPT_R)) {
	cmp.F = wald_omit_F(testvars, umod);
	cmp.robust = 1;
    } else if (LIMDEP(cmp.ci) && umod != NULL && rmod != NULL) {
	cmp.LR = 2.0 * (umod->lnL - rmod->lnL);
    } else if (cmp.ci == OLS) {
	cmp.F = ((rmod->ess - umod->ess) / umod->ess) * cmp.dfd / cmp.dfn;
    } else if (cmp.dfn >= 1) {
	cmp.err = wald_test(testvars, umod, &cmp.wald, NULL);
    }

    if (pmodA != NULL && pmodB != NULL && flag != COMPARE_ADD_WALD) {
	int miss = 0;

	cmp.score = 0;
	for (i=0; i<C_MAX; i++) { 
	    if (na(pmodB->criterion[i]) || na(pmodA->criterion[i])) {
		miss++;
		continue;
	    }
	    if (pmodB->criterion[i] < pmodA->criterion[i]) {
		cmp.score++;
	    }
	}
	if (miss == C_MAX) {
	    cmp.score = -1;
	}
    }

    if (na(cmp.LR) && na(cmp.F) && 
	umod != NULL && rmod != NULL && 
	!na(umod->lnL) && !na(rmod->lnL) &&
	umod->lnL > rmod->lnL) {
	/* if we're doing a Wald test and a likelihood ratio test is
	   also supported, add the LR test */
	cmp.LR = 2.0 * (umod->lnL - rmod->lnL);
    }

    return cmp;
}

/* reconstitute full varlist for WLS, AR and count models */

static int *
full_model_list (const MODEL *pmod, const int *inlist)
{
    int i, len = 0, pos = 0;
    int *flist = NULL;

    if (pmod->ci == WLS) { 
	len = inlist[0] + 2;
    } else if (COUNT_MODEL(pmod->ci)) {
	len = inlist[0] + 3;
    } else if (pmod->ci == AR) {
	pos = pmod->arinfo->arlist[0] + 1;
	len = pos + inlist[0] + 2;
    }

    flist = malloc(len * sizeof *flist);
    if (flist == NULL) {
	return NULL;
    }

    if (pmod->ci == WLS) { 
	flist[0] = len - 1;
	flist[1] = pmod->nwt;
	for (i=1; i<=inlist[0]; i++) {
	    flist[i+1] = inlist[i];
	}
    } else if (COUNT_MODEL(pmod->ci)) {
	int offvar = gretl_model_get_int(pmod, "offset_var");

	flist[0] = len - 1;
	for (i=1; i<=inlist[0]; i++) {
	    flist[i] = inlist[i];
	}
	flist[flist[0] - 1] = LISTSEP;
	flist[flist[0]] = offvar;
    } else if (pmod->ci == AR) {
	flist[0] = len - 2;
	for (i=1; i<pos; i++) {
	    flist[i] = pmod->arinfo->arlist[i];
	}
	flist[pos] = LISTSEP;
	for (i=1; i<=inlist[0]; i++) {
	    flist[pos+i] = inlist[i];
	}
    }

    return flist;
}

static gretlopt retrieve_arbond_opts (const MODEL *pmod)
{
    gretlopt opt = OPT_NONE;

    if (pmod->opt & OPT_D) {
	opt |= OPT_D;
    }

    if (gretl_model_get_int(pmod, "asy")) {
	opt |= OPT_A;
    }

    if (gretl_model_get_int(pmod, "step") == 2) {
	opt |=OPT_T;
    }    

    return opt;
}

static int obs_diff_ok (const MODEL *m_old, const MODEL *m_new)
{
    int tdiff, ndiff = m_new->nobs - m_old->nobs;

    if (m_old->ci == AR1) {
	return 0;
    }

    if (ndiff > 0) {
	tdiff = (m_new->t2 - m_new->t1) - (m_old->t2 - m_old->t1);
	if (ndiff == tdiff) {
	    return 1;
	}
    }

    return 0;
}

static int has_offset_var (const MODEL *pmod)
{
    return (COUNT_MODEL(pmod->ci) &&
	    gretl_model_get_int(pmod, "offset_var") > 0);   
}

#define be_quiet(o) ((o & OPT_A) || (o & OPT_Q))

#define SMPL_DEBUG 0

static MODEL replicate_estimator (const MODEL *orig, int **plist,
				  double ***pZ, DATAINFO *pdinfo,
				  gretlopt myopt, PRN *prn)
{
    MODEL rep;
    const char *param = NULL;
    char altparm[32] = {0};
    double rho = 0.0;
    int *list = *plist;
    int mc = get_model_count();
    int repci = orig->ci;
    int order = 0;
    int first = 1;

    gretl_model_init(&rep);

    /* recreate options and auxiliary vars, if required */

    transcribe_option_flags(&myopt, orig->opt,
			    OPT_D | OPT_J | OPT_R);

    if (orig->ci == AR1) {
	if (orig->opt & OPT_H) {
	    myopt |= OPT_H;
	    if (gretl_model_get_int(orig, "no-corc")) {
		myopt |= OPT_B;
	    }
	} else if (orig->opt & OPT_P) {
	    myopt |= OPT_P;
	}
	rho = estimate_rho(list, pZ, pdinfo, myopt, prn, &rep.errcode);
    } else if (orig->ci == WLS || orig->ci == AR || has_offset_var(orig)) {
	int *full_list = full_model_list(orig, list);

	free(list);
	if (full_list == NULL) {
	    rep.errcode = E_ALLOC;
	} else {
	    list = *plist = full_list;
	}
    } else if (orig->ci == ARBOND) {
	param = gretl_model_get_data(orig, "istr");
	myopt |= retrieve_arbond_opts(orig);
    } else if (orig->ci == ARCH) {
	order = gretl_model_get_int(orig, "arch_order");
    } else if (orig->ci == LOGIT || orig->ci == PROBIT) {
	if (gretl_model_get_int(orig, "ordered")) {
	    myopt |= OPT_D;
	} else if (gretl_model_get_int(orig, "multinom")) {
	    myopt |= OPT_M;
	}
    } else if (orig->ci == PANEL) {
	if (gretl_model_get_int(orig, "pooled")) {
	    myopt |= OPT_P;
	} else if (orig->opt & OPT_U) {
	    myopt |= OPT_U;
	} else if (orig->opt & OPT_W) {
	    myopt |= OPT_W;
	    if (gretl_model_get_int(orig, "iters")) {
		myopt |= OPT_I;
	    }
	}
    } else if (orig->ci == LAD && gretl_model_get_int(orig, "rq")) {
	double x;

	x = gretl_model_get_double(orig, "tau");
	sprintf(altparm, "%g", x);
	if (gretl_model_get_int(orig, "rq_nid")) {
	    myopt |= OPT_R;
	}
	x = gretl_model_get_double(orig, "rq_alpha");
	if (!na(x)) {
	    myopt |= OPT_I;
	    set_optval_double(QUANTREG, OPT_I, x);
	}
    } else if (orig->ci == IVREG) {
	transcribe_option_flags(&myopt, orig->opt, OPT_L | OPT_G);
    }

    if (rep.errcode) {
	return rep;
    }

 try_again:

    switch (orig->ci) {

    case AR:
	rep = ar_func(list, pZ, pdinfo, myopt, prn);
	break;
    case ARBOND:
	rep = arbond_model(list, param, (const double **) *pZ, 
			   pdinfo, myopt, prn);
	break;
    case ARCH:
	rep = arch_model(list, order, pZ, pdinfo, myopt, prn);
	break;
    case LOGIT:
    case PROBIT:
	rep = logit_probit(list, pZ, pdinfo, orig->ci, myopt, NULL);
	break;
    case TOBIT:
	rep = tobit_model(list, pZ, pdinfo, NULL);
	break;
    case LAD:
	if (gretl_model_get_int(orig, "rq")) {
	    rep = quantreg(altparm, list, pZ, pdinfo, myopt, NULL);
	} else {
	    rep = lad(list, pZ, pdinfo);
	}
	break;
    case POISSON:
    case NEGBIN:
	rep = count_model(list, orig->ci, pZ, pdinfo, myopt, NULL);
	break;
    case HECKIT:
	rep = heckit_model(list, pZ, pdinfo, myopt, NULL);
	break;
    case IVREG:
	rep = ivreg(list, pZ, pdinfo, myopt);
	break;
    case LOGISTIC: 
	{
	    char lmaxstr[32];
	    double lmax;

	    lmax = gretl_model_get_double(orig, "lmax");
	    sprintf(lmaxstr, "lmax=%g", lmax);
	    rep = logistic_model(list, pZ, pdinfo, lmaxstr);
	}
	break;
    case PANEL:
	rep = panel_model(list, pZ, pdinfo, myopt, prn);
	break;
    default:
	/* handles OLS, AR1, WLS, HSK, etc. */
	if (rho != 0.0) {
	    rep = ar1_lsq(list, pZ, pdinfo, repci, myopt, rho);
	} else if (gretl_model_get_int(orig, "pooled")) {
	    myopt |= OPT_P;
	    rep = panel_model(list, pZ, pdinfo, myopt, prn);
	} else {
	    rep = lsq(list, pZ, pdinfo, repci, myopt);
	}
	break;
    }

#if SMPL_DEBUG
    fprintf(stderr, "replicate_estimator:\n"
	    " orig: t1=%d, t2=%d, nobs = %d\n"
	    " rep:  t1=%d, t2=%d, nobs = %d\n",
	    orig->t1, orig->t2, orig->nobs,
	    rep.t1, rep.t2, rep.nobs);
#endif

    /* check that we got the same sample as the original */
    if (!rep.errcode && rep.nobs != orig->nobs) {
	if (first && obs_diff_ok(orig, &rep)) {
	    pdinfo->t1 = orig->t1;
	    pdinfo->t2 = orig->t2;
	    clear_model(&rep);
	    first = 0;
	    goto try_again;
	} else {
	    fprintf(stderr, "Original obs = %d but new = %d\n", orig->nobs, rep.nobs);
	    rep.errcode = E_DATA;
	}
    } 

    /* if the model count went up for an aux regression,
       bring it back down */
    if (be_quiet(myopt) && get_model_count() > mc) {
	model_count_minus();
    }

    return rep;
}

static void nonlin_test_header (int code, PRN *prn)
{
    pputc(prn, '\n');
    if (code == AUX_SQ) {
	pputs(prn, _("Non-linearity test (squared terms)"));
    } else {
	pputs(prn, _("Non-linearity test (log terms)"));
    }
    pputs(prn, "\n\n");
}

static int
real_nonlinearity_test (MODEL *pmod, int *list,
			double ***pZ, DATAINFO *pdinfo,
			int aux_code, gretlopt opt, 
			PRN *prn)
{
    MODEL aux;
    int t, err = 0;

    /* grow data set to accommodate new dependent var */
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[pdinfo->v - 1][t] = pmod->uhat[t];
    }

    /* replace the dependent var */
    list[1] = pdinfo->v - 1;

    aux = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (aux.errcode) {
	err = aux.errcode;
	fprintf(stderr, "auxiliary regression failed\n");
    } else {
	double trsq = aux.rsq * aux.nobs;
	int df = list[0] - pmod->list[0];
	double pval = chisq_cdf_comp(df, trsq);

	aux.aux = aux_code;

	if (opt & OPT_Q) {
	    nonlin_test_header(aux_code, prn);
	} else {
	    printmodel(&aux, pdinfo, opt, prn);
	    pputc(prn, '\n');
	}

	pprintf(prn, "  %s: TR^2 = %g,\n  ", _("Test statistic"), trsq);
	pprintf(prn, _("with p-value = prob(Chi-square(%d) > %g) = %g\n\n"), 
		df, trsq, pval);

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = model_test_new((aux_code == AUX_SQ)?
				  GRETL_TEST_SQUARES : GRETL_TEST_LOGS);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LM);
		model_test_set_dfn(test, df);
		model_test_set_value(test, trsq);
		model_test_set_pvalue(test, chisq_cdf_comp(df, trsq));
		maybe_add_test_to_model(pmod, test);
	    }
	}

	record_test_result(trsq, pval, _("non-linearity"));
    } 

    clear_model(&aux);

    return err;
}

/**
 * nonlinearity_test:
 * @pmod: pointer to original model.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @aux_code: %AUX_SQ for squares or %AUX_LOG for logs
 * @opt: if contains %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Run an auxiliary regression to test @pmod for non-linearity,
 * via the addition of either squares or logs of the original
 * indepdendent variables.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int nonlinearity_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		       int aux_code, gretlopt opt, PRN *prn) 
{
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int *tmplist = NULL;
    const int orig_nvar = pdinfo->v; 
    int err = 0;

    if (!command_ok_for_model(ADD, 0, pmod->ci)) {
	return E_NOTIMP;
    }

    if (pmod->ci == LOGISTIC || pmod->ci == LAD) {
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod->list, pdinfo, pmod->ID);
    if (err) {
	return err;
    }

    /* re-impose the sample that was in force when the original model
       was estimated */
    impose_model_smpl(pmod, pdinfo);

    /* add squares or logs */
    tmplist = augment_regression_list(pmod->list, aux_code, 
				      pZ, pdinfo);
    if (tmplist == NULL) {
	return E_ALLOC;
    } else if (tmplist[0] == pmod->list[0]) {
	/* no vars were added */
	if (aux_code == AUX_SQ) {
	    fprintf(stderr, "gretl: generation of squares failed\n");
	    err = E_SQUARES;
	} else if (aux_code == AUX_LOG) {
	    fprintf(stderr, "gretl: generation of logs failed\n");
	    err = E_LOGS;
	}
    }

    if (!err) {
	err = real_nonlinearity_test(pmod, tmplist, pZ, pdinfo, aux_code, 
				     opt, prn);
    }
	
    /* trash any extra variables generated (squares, logs) */
    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    free(tmplist);

    return err;
}

static void remove_special_flags (gretlopt *popt)
{
    gretlopt rem = OPT_S | OPT_T | OPT_B | OPT_P | OPT_A;

    delete_option_flags(popt, rem);
}

static int add_vars_missing (const MODEL *pmod, const int *list,
			     const double **Z)
{
    int i, t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!model_missing(pmod, t)) {
	    for (i=1; i<=list[0]; i++) {
		if (na(Z[list[i]][t])) {
		    return E_MISSDATA;
		}
	    }
	}
    }

    return 0;
}

/**
 * add_test:
 * @addvars: list of variables to add to original model.
 * @orig: pointer to original model.
 * @pmod: pointer to receive new model, with vars added.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: can contain %OPT_Q (quiet) to suppress printing
 * of the new model, %OPT_O to print covariance matrix,
 * %OPT_I for silent operation.  Additional options that
 * are applicable only for 2SLS models: %OPT_B (add to
 * both list of regressors and list of instruments), %OPT_T
 * (add to list of instruments only).
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after adding the specified
 * variables.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int add_test (const int *addvars, MODEL *orig, MODEL *pmod, 
	      double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn)
{
    gretlopt est_opt = opt;
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int *tmplist = NULL;
    const int orig_nvar = pdinfo->v; 
    int err = 0;

    if (orig == NULL || orig->list == NULL || addvars == NULL) {
	return 1;
    }

    if (!command_ok_for_model(ADD, 0, orig->ci)) {
	return E_NOTIMP;
    }

    if (exact_fit_check(orig, prn)) {
	return 0;
    }

    /* check for changes in original list members */
    err = list_members_replaced(orig->list, pdinfo, orig->ID);
    if (err) {
	return err;
    }

    /* check for NAs in add list relative to model */
    err = add_vars_missing(orig, addvars, (const double **) *pZ);
    if (err) {
	return err;
    }

    /* create augmented regression list */
    if (orig->ci == IVREG) {
	tmplist = ivreg_list_add(orig->list, addvars, opt, &err);
    } else if (orig->ci == PANEL || orig->ci == ARBOND) {
	tmplist = panel_list_add(orig, addvars, &err);
    } else {
	tmplist = gretl_list_add(orig->list, addvars, &err);
    }

    if (err) {
	return err;
    }

    /* impose as sample range the sample range in force
       when the original model was estimated */
    impose_model_smpl(orig, pdinfo);

    /* don't pass special opts to replicate_estimator() */
    remove_special_flags(&est_opt);

    /* Run augmented regression, matching the original estimation
       method; use OPT_Z to suppress the elimination of perfectly
       collinear variables.
    */
    *pmod = replicate_estimator(orig, &tmplist, pZ, pdinfo, 
				(est_opt | OPT_Z), prn);

    if (pmod->errcode) {
	err = pmod->errcode;
	errmsg(err, prn);
    }

    if (!err) {
	int flag = COMPARE_ADD_WALD;
	struct COMPARE cmp;
	MODEL *origmod = orig;
	int *addlist;

	pmod->aux = AUX_ADD;

	if (print_add_omit_model(pmod, opt)) {
	    printmodel(pmod, pdinfo, est_opt, prn);
	}

	if (orig->ci == OLS && pmod->nobs == orig->nobs) {
	    flag = COMPARE_ADD;
	} else if (opt & OPT_B) {
	    /* ivreg: adding as both regressor and instrument */
	    origmod = NULL;
	}

	addlist = gretl_list_diff_new(pmod->list, orig->list, 2);
	if (addlist != NULL) {
	    cmp = add_or_omit_compare(origmod, pmod, flag, addlist);
	    gretl_make_compare(&cmp, addlist, orig, pdinfo, opt, prn);
	    free(addlist);
	}
    }

    /* trash any extra variables generated (squares, logs) */
    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    free(tmplist);

    return err;
}

static int wald_omit_test (const int *list, MODEL *pmod, 
			   const DATAINFO *pdinfo, gretlopt opt,
			   PRN *prn)
{
    struct COMPARE cmp;
    int *test = NULL;
    int err = 0;

    /* test validity of omissions */

    if (pmod->ci == IVREG) {
	test = ivreg_list_omit(pmod->list, list, opt, &err);
    } else if (pmod->ci == PANEL || pmod->ci == ARBOND) {
	test = panel_list_omit(pmod, list, &err);
    } else {
	test = gretl_list_omit(pmod->list, list, 2, &err);
    }

    if (err) {
	return err;
    }

    free(test);

    cmp = add_or_omit_compare(pmod, NULL, COMPARE_OMIT_WALD, list);
    if (cmp.err) {
	return cmp.err;
    }

    gretl_make_compare(&cmp, list, pmod, pdinfo, opt, prn);

    return err;
}

/* Check whether coefficient @i corresponds to a variable
   that is removable from the model: this is the case if either
   (a) the @cands list is NULL, or (b) coefficient @i is the
   coefficient on one of the variables in @cands.
*/

static int coeff_is_removable (const int *cands, const MODEL *pmod,
			       DATAINFO *pdinfo, int i)
{
    int ret = 1;

    if (cands != NULL) {
	const char *vname;
	int j, pj;

	ret = 0; /* reverse the presumption */

	for (j=1; j<=cands[0]; j++) {
	    vname = pdinfo->varname[cands[j]];
	    pj = gretl_model_get_param_number(pmod, pdinfo, vname);
	    if (pj == i) {
		ret = 1;
		break;
	    }
	}
    }

    return ret;
}

/* Determine if a model contains a variable with p-value
   greater than some cutoff alpha_max; and if so, remove
   this variable from the regression list.  If the list
   @cands is non-NULL then confine the search to candidate
   variables in that list.
*/

static int auto_drop_var (const int *cands, const MODEL *pmod, 
			  int *list, DATAINFO *pdinfo, 
			  double alpha_max, int d0, 
			  PRN *prn, int *err)
{
    double tstat, pv = 0.0, tmin = 4.0;
    int i, k = -1;
    int ret = 0;

    if (pmod->ncoeff == 1) {
	return 0;
    }
    
    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
	if (coeff_is_removable(cands, pmod, pdinfo, i)) {
	    tstat = fabs(pmod->coeff[i] / pmod->sderr[i]);
	    if (tstat < tmin) {
		tmin = tstat;
		k = i;
	    }
	}
    }

    if (k >= 0) {
	pv = coeff_pval(pmod->ci, tmin, pmod->dfd);
    }

    if (pv > alpha_max) {
	char pname[VNAMELEN];

	if (d0) {
	    pputc(prn, '\n');
	    pprintf(prn, _("Sequential elimination using two-sided alpha = %.2f"),
		    alpha_max);
	    pputs(prn, "\n\n");
	}

	gretl_model_get_param_name(pmod, pdinfo, k, pname);
	pprintf(prn, _(" Dropping %-16s (p-value %.3f)\n"), pname, pv);
	*err = gretl_list_delete_at_pos(list, k + 2);
	ret = 1;
    }

    return ret;
}

static void list_copy_values (int *targ, const int *src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }
}

/* run a loop in which the least significant variable is dropped
   from the regression list, provided its p-value exceeds some
   specified cutoff.  FIXME this probably still needs work for 
   estimators other than OLS.  If @omitlist is non-NULL the
   routine is confined to members of the list.
*/

static int auto_omit (const int *omitlist, 
		      MODEL *orig, MODEL *pmod, 
		      double ***pZ, DATAINFO *pdinfo, 
		      gretlopt est_opt, gretlopt opt,
		      PRN *prn)
{
    double amax;
    int *tmplist = NULL;
    int i, err = 0;

    tmplist = gretl_list_copy(orig->list);
    if (tmplist == NULL) {
	return E_ALLOC;
    }

    amax = get_optval_double(OMIT, OPT_A);
    if (na(amax) || amax <= 0.0 || amax >= 1.0) {
	amax = 0.10;
    }

    if (!auto_drop_var(omitlist, orig, tmplist, pdinfo, amax, 1, prn, &err)) {
	free(tmplist);
	return (err)? err : E_NOOMIT;
    }    

    for (i=0; !err; i++) {
	if (i > 0) {
	    set_reference_missmask_from_model(orig);
	}
	*pmod = replicate_estimator(orig, &tmplist, pZ, pdinfo, 
				    est_opt, prn);
	if (pmod->errcode) {
	    err = pmod->errcode;
	    fprintf(stderr, "auto_omit: error %d from replicate_estimator\n", err);
	} else {
	    list_copy_values(tmplist, pmod->list);
	    if (auto_drop_var(omitlist, pmod, tmplist, pdinfo, 
			      amax, 0, prn, &err)) {
		model_count_minus();
		clear_model(pmod);
	    } else {
		break;
	    }
	}
    }

    free(tmplist);

    return err;
}

/* create reduced list for "omit" test on model, based on
   the list of variables to be dropped, omitvars
*/

static int make_short_list (MODEL *orig, const int *omitvars,
			    gretlopt opt, int *omitlast,
			    int **plist)
{
    int *list = NULL;
    int err = 0;

    if (omitvars == NULL || omitvars[0] == 0) {
	if (orig->ci == IVREG) {
	    return E_PARSE;
	} else {
	    *omitlast = 1;
	}
    }

    if (orig->ci == IVREG) {
	list = ivreg_list_omit(orig->list, omitvars, opt, &err);
    } else if (orig->ci == PANEL || orig->ci == ARBOND) {
	list = panel_list_omit(orig, omitvars, &err);
    } else if (*omitlast) {
	/* special: just drop the last variable */
	list = gretl_list_omit_last(orig->list, &err);
    } else {
	list = gretl_list_omit(orig->list, omitvars, 2, &err);
    }

    if (list != NULL && list[0] == 1) {
	err = E_NOVARS;
    }

    *plist = list;

    return err;
}

static int omit_options_inconsistent (gretlopt opt)
{
    if ((opt & OPT_T) || (opt & OPT_B)) {
	/* 2sls: omitting variable as instrument */
	if (opt & OPT_W) {
	    /* can't use Wald method on original VCV */
	    return 1;
	}
    }

    if ((opt & OPT_A) && (opt & OPT_W)) {
	/* auto and Wald options incompatible */
	return 1;
    }

    return 0;
}

/**
 * omit_test:
 * @omitvars: list of variables to omit from original model.
 * @orig: pointer to original model.
 * @pmod: pointer to receive new model, with vars omitted.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt: can contain %OPT_Q (quiet) to suppress printing
 * of the new model, %OPT_O to print covariance matrix,
 * %OPT_I for silent operation; for %OPT_A, see below.
 * Additional options that are applicable only for IV models: 
 * %OPT_B (omit from both list of regressors and list of 
 * instruments), %OPT_T (omit from list of instruments only).
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after removing the variables
 * specified in @omitvars.  Or if %OPT_A is given, proceed 
 * sequentially, at each step dropping the least significant 
 * variable provided its p-value is above a certain threshold 
 * (currently 0.10, two-sided).  Or if @omitvars is %NULL 
 * and @orig was not estimated using two-stage least squares,
 * drop the last independent variable in @orig.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test (const int *omitvars, MODEL *orig, MODEL *pmod, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    gretlopt est_opt = opt;
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int omitlast = 0;
    int *tmplist = NULL;
    int err = 0;

    if (orig == NULL || orig->list == NULL) {
	err = E_DATA;
    } else if (!command_ok_for_model(OMIT, 0, orig->ci)) {
	err = E_NOTIMP;
    } else if (omit_options_inconsistent(opt)) {
	err = E_BADOPT;
    }

    if (err) {
	return err;
    }

    if (opt & OPT_W) {
	return wald_omit_test(omitvars, orig, pdinfo, opt, prn);
    }

    /* check that vars to omit have not been redefined */
    if ((err = list_members_replaced(orig->list, pdinfo, orig->ID))) {
	return err;
    }

    if (!(opt & OPT_A)) {
	err = make_short_list(orig, omitvars, opt, &omitlast, &tmplist);
	if (err) {
	    return err;
	}
    }

    /* impose the sample range used for the original model */ 
    impose_model_smpl(orig, pdinfo);

    /* set the mask for missing obs within the sample range, based
       on the original model */
    set_reference_missmask_from_model(orig);

    /* extract option flags that should not be passed to estimator
       functions */
    remove_special_flags(&est_opt);

    if (opt & OPT_A) {
	if (omitvars != NULL && omitvars[0] > 0) {
	    err = auto_omit(omitvars, orig, pmod, pZ, pdinfo, est_opt, opt, prn);
	} else {
	    err = auto_omit(NULL, orig, pmod, pZ, pdinfo, est_opt, opt, prn);
	}
    } else {
	*pmod = replicate_estimator(orig, &tmplist, pZ, pdinfo, est_opt, prn);
	err = pmod->errcode;
    }

    if (err) {
	errmsg(err, prn);
    } else {
	if (orig->ci == LOGIT || orig->ci == PROBIT) {
	    pmod->aux = AUX_OMIT;
	}

	if (print_add_omit_model(orig, opt)) {
	    printmodel(pmod, pdinfo, est_opt, prn); 
	}	

	if (!omitlast) {
	    int flag = COMPARE_OMIT;
	    MODEL *newmod = pmod;
	    struct COMPARE cmp;
	    int *omitlist;

	    if (opt & OPT_B) {
		/* ivreg: omitting both as regressor and instrument */
		newmod = NULL;
		flag = COMPARE_OMIT_WALD;
	    } 

	    omitlist = gretl_list_diff_new(orig->list, pmod->list, 2);
	    if (omitlist != NULL) {
		cmp = add_or_omit_compare(orig, newmod, flag, omitlist);
		gretl_make_compare(&cmp, omitlist, orig, pdinfo, opt, prn); 
		free(omitlist);
	    }
	}

	if (orig->ci == LOGIT || orig->ci == PROBIT) {
	    pmod->aux = AUX_NONE;
	}
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    free(tmplist);

    return err;
}

double get_dw_pvalue (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		      int *err)
{
    MODEL dwmod;
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int *list = NULL;
    double pv = NADBL;

    if (pmod == NULL || pmod->list == NULL) {
	*err = E_DATA;
    } else if ((pmod->ci != OLS && pmod->ci != PANEL) || 
	       pmod->missmask != NULL || na(pmod->dw)) {
	*err = E_BADSTAT;
    } else {
	/* check that relevant vars have not been redefined */
	*err = list_members_replaced(pmod->list, pdinfo, pmod->ID);
    }

    if (!*err) {
	list = gretl_list_copy(pmod->list);
	if (list == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	return NADBL;
    }

    gretl_model_init(&dwmod);

    /* impose the sample range used for the original model */ 
    impose_model_smpl(pmod, pdinfo);

    dwmod = replicate_estimator(pmod, &list, pZ, pdinfo, OPT_A | OPT_I, NULL);
    *err = dwmod.errcode;

    if (!*err) {
	pv = gretl_model_get_double(&dwmod, "dw_pval");
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;
    
    clear_model(&dwmod);
    free(list);

    return pv;
}

/**
 * reset_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if contains %OPT_S, save test results to model. %OPT_Q
 * suppresses the printout of the auxiliary regression. %OPT_R and
 * %OPT_C stand for "squares only" and "cubes only", respectively.
 * @prn: gretl printing struct.
 *
 * Carries out Ramsey's RESET test for model specification.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int reset_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn)
{
    int *newlist = NULL;
    MODEL aux;
    double RF;
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int i, t, v = pdinfo->v; 
    int addcols;
    const char *mode;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    err = incompatible_options(opt, OPT_C | OPT_R);

    if (err) {
	return err;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    gretl_model_init(&aux);

    if (opt & OPT_C) {
	addcols = 1;
	mode = N_("squares only");
    } else if (opt & OPT_R) {
	addcols = 1;
	mode = N_("cubes only");
    } else {
	addcols = 2;
	mode = N_("squares and cubes");
    }

    impose_model_smpl(pmod, pdinfo);

    if (pmod->ncoeff + addcols >= pdinfo->t2 - pdinfo->t1) {
	err = E_DF;
    }

    if (!err) {
	newlist = malloc((pmod->list[0] + addcols + 1) * sizeof *newlist);
	if (newlist == NULL) {
	    err = E_ALLOC;
	}
    }

    if (newlist != NULL) {
	newlist[0] = pmod->list[0] + addcols;
	for (i=1; i<=pmod->list[0]; i++) {
	    newlist[i] = pmod->list[i];
	}
	if (dataset_add_series(addcols, pZ, pdinfo)) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add yhat^2 and/or yhat^3 to data set */
	int sqcol = v;
	int cubecol = (opt & OPT_C) ? v : v + 1;

	for (t = pmod->t1; t<=pmod->t2; t++) {
	    double xx = pmod->yhat[t];

	    if (!(opt & OPT_C)) {
		(*pZ)[sqcol][t] = xx * xx;
	    }
	    if (!(opt & OPT_R)) {
		(*pZ)[cubecol][t] = xx * xx * xx;
	    }
	}

	if (!(opt & OPT_C)) {
	    strcpy(pdinfo->varname[sqcol], "yhat^2");
	    newlist[pmod->list[0] + 1] = sqcol;
	}

	if (!(opt & OPT_R)) {
	    strcpy(pdinfo->varname[cubecol], "yhat^3");
	    newlist[newlist[0]] = cubecol;
	}
    }

    if (!err) {
	aux = lsq(newlist, pZ, pdinfo, OLS, OPT_A);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    } 

    if (!err) {
	double pval;

	aux.aux = AUX_RESET;

	if (!(opt & OPT_Q)) {
	    printmodel(&aux, pdinfo, OPT_NONE, prn);
	} else {
	    if (!(opt & OPT_G)) {
		pputc(prn, '\n');
	    }
	    pputs(prn, _("RESET test for specification"));
	    pprintf(prn, " (%s)\n", _(mode));
	}

	RF = ((pmod->ess - aux.ess) / addcols) / (aux.ess / aux.dfd);
	pval = snedecor_cdf_comp(addcols, aux.dfd, RF);

	pprintf(prn, "%s: F = %f,\n", _("Test statistic"), RF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		addcols, aux.dfd, RF, pval);
	pputc(prn, '\n');

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_RESET);
	    gretlopt topt = OPT_NONE;

	    if (test != NULL) {
		if (opt & OPT_R) {
		    topt = OPT_R;
		} else if (opt & OPT_C) {
		    topt = OPT_C;
		}
		model_test_set_teststat(test, GRETL_STAT_RESET);
		model_test_set_dfn(test, addcols);
		model_test_set_dfd(test, aux.dfd);
		model_test_set_value(test, RF);
		model_test_set_pvalue(test, pval);
		model_test_set_opt(test, topt);
		maybe_add_test_to_model(pmod, test);
	    }	    
	}

	record_test_result(RF, pval, "RESET");
    }

    free(newlist);
    dataset_drop_last_variables(addcols, pZ, pdinfo); 
    clear_model(&aux); 

    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    return err;
}

static void bg_test_header (int order, PRN *prn, int ivreg)
{
    if (ivreg) {
	pprintf(prn, "\n%s ", _("Godfrey (1994) test for"));
    } else {
	pprintf(prn, "\n%s ", _("Breusch-Godfrey test for")); 
    }

    if (order > 1) {
	pprintf(prn, "%s %d\n", _("autocorrelation up to order"),
		order);
    } else {
	pprintf(prn, "%s\n", _("first-order autocorrelation"));
    }

    pputc(prn, '\n');
}

static double ivreg_autocorr_wald_stat (MODEL *aux, int order, int *err)
{
    gretl_vector *b = gretl_vector_alloc(order);
    gretl_matrix *V1 = gretl_matrix_alloc(order, order);
    gretl_vector *WT = gretl_vector_alloc(1);
    gretl_matrix *V0 = NULL;
    double x = NADBL;
    int i, j, ki, kj;

    if (b == NULL || V1 == NULL || WT == NULL) {
	*err = E_ALLOC;
    } else {
	V0 = gretl_model_get_matrix(aux, M_VCV, err);
    }

    if (!*err) {
	ki = aux->ncoeff - order;

	for (i=0; i<order; i++) {
	    x = aux->coeff[ki];
	    gretl_vector_set(b, i, x);
	    kj = ki;
	    for (j=i; j<order; j++) {
		x = gretl_matrix_get(V0, ki, kj); 
		gretl_matrix_set(V1, i, j, x);
		gretl_matrix_set(V1, j, i, x);
		kj++;
	    }
	    ki++;
	}
	*err = gretl_invert_symmetric_matrix(V1);
    }

    if (!*err) {
	gretl_matrix_qform(b, GRETL_MOD_NONE, V1, WT, GRETL_MOD_NONE);
	x = gretl_vector_get(WT, 0) / order;
    }
	
    gretl_vector_free(WT);
    gretl_vector_free(b);
    gretl_matrix_free(V1);
    gretl_matrix_free(V0);

    return x;
}

/**
 * ivreg_autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save test results to model;
 * if %OPT_Q, be less verbose.
 * @prn: gretl printing struct.
 *
 * Tests the given IV model for autocorrelation of order equal
 * to the specified value, or equal to the frequency of the data if
 * the supplied @order is zero, as per Godfrey (1999), "Testing for
 * Serial Correlation by Variable Addition in Dynamic Models Estimated
 * by Instrumental Variables", RES. Note that none of the
 * asymptotically equivalent tests given on page 553 is used
 * here. Instead, we estimate the model augmented with lags and then
 * perform a Wald-type test. The resulting chi-square statistic is
 * divided by its degrees of freedom as a finite-sample adjustment and
 * compared to an F distribution.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

static int ivreg_autocorr_test (MODEL *pmod, int order, 
				double ***pZ, DATAINFO *pdinfo, 
				gretlopt opt, PRN *prn)
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int n = pdinfo->n;
    int v = pdinfo->v;
    int *addlist = NULL;
    int *testlist = NULL;
    double x, pval = 1.0;
    MODEL aux;
    int i, t;
    int err = 0;

    if (dataset_is_panel(pdinfo)) { 
	return E_NOTIMP;
    }

    if (pmod->missmask != NULL) {
	return E_DATA;
    }

    /* impose original sample range */
    impose_model_smpl(pmod, pdinfo);

    gretl_model_init(&aux);

    if (order <= 0) {
	order = pdinfo->pd;
    }

    if (pmod->ncoeff + order >= pdinfo->t2 - pdinfo->t1) {
	return E_DF;
    }

    addlist = gretl_list_new(order);

    if (addlist == NULL) {
	err = E_ALLOC;
    } else {
	err = dataset_add_series(1, pZ, pdinfo);
    }

    if (!err) {
	/* add uhat to data set */
	for (t=0; t<n; t++) {
	    (*pZ)[v][t] = pmod->uhat[t];
	}
	strcpy(pdinfo->varname[v], "uhat");
	strcpy(VARLABEL(pdinfo, v), _("residual"));
	/* then lags of same */
	for (i=1; i<=order && !err; i++) {
	    int lnum;

	    lnum = laggenr(v, i, pZ, pdinfo);

	    if (lnum < 0) {
		gretl_errmsg_set(_("lagging uhat failed"));
		err = E_LAGS;
	    } else {
		/* set the first entries to 0 for compatibility with Godfrey (1994)
		   and PcGive (not perfect) */
		for (t=smpl_t1; t<smpl_t1+i; t++) {
		    (*pZ)[lnum][t] = 0;
		}
		addlist[i] = lnum;
	    }
	}
	if (!err) {
	    /* compose augmented regression list */
	    testlist = ivreg_list_add(pmod->list, addlist, OPT_B, &err);
	}
    } 

    if (!err) {
	gretlopt ivopt = OPT_A;

	transcribe_option_flags(&ivopt, pmod->opt,
				OPT_L | OPT_G | OPT_R);
	aux = ivreg(testlist, pZ, pdinfo, ivopt);
	err = aux.errcode;
    }

    if (!err) {
	x = ivreg_autocorr_wald_stat(&aux, order, &err);
    }

    if (!err) {
	aux.aux = AUX_AR;
	gretl_model_set_int(&aux, "BG_order", order);
	pval = snedecor_cdf_comp(order, aux.nobs - pmod->ncoeff - order, x);

	if (opt & OPT_Q) {
	    bg_test_header(order, prn, 1);
	} else {
	    printmodel(&aux, pdinfo, OPT_S, prn);
	} 

	pputc(prn, '\n');
	pprintf(prn, "%s: Pseudo-LMF = %f,\n", _("Test statistic"), x);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, aux.nobs - pmod->ncoeff, x, pval);
	pputc(prn, '\n');
	record_test_result(x / order, pval, _("autocorrelation"));

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_AUTOCORR);
	    
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LMF);
		model_test_set_dfn(test, order);
		model_test_set_dfd(test, aux.nobs - pmod->ncoeff);
		model_test_set_order(test, order);
		model_test_set_value(test, x);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	    
	}
    }

    free(addlist);
    free(testlist);

    dataset_drop_last_variables(pdinfo->v - v, pZ, pdinfo); 
    clear_model(&aux); 

    /* reset sample as it was */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

    return err;
}

/**
 * autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save test results to model;
 * if %OPT_Q, be less verbose.
 * @prn: gretl printing struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero. Gives TR^2 and LMF test statistics.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int autocorr_test (MODEL *pmod, int order, 
		   double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn)
{
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int *newlist = NULL;
    MODEL aux;
    double RSSxe, RSSx = pmod->ess;
    int i, t, n = pdinfo->n, v = pdinfo->v; 
    double trsq, LMF, lb, pval = 1.0;
    int err = 0;

    if (pmod->ci == IVREG) {
	return ivreg_autocorr_test(pmod, order, pZ, pdinfo, opt, prn);
    }

    if (pmod->ci != OLS && pmod->ci != VAR) { 
	return E_NOTIMP;
    }

    if (pmod->missmask != NULL) {
	return E_DATA;
    }

    if (dataset_is_panel(pdinfo)) {
#if 1 /* FIXME */
	return E_NOTIMP;
#else
	return panel_autocorr_test(pmod, order, *pZ, pdinfo, opt, prn);
#endif
    }

    gretl_model_init(&aux);

    if (order <= 0) {
	order = pdinfo->pd;
    }

    if (pmod->ncoeff + order >= pmod->t2 - pmod->t1) {
	return E_DF;
    }

    newlist = gretl_list_new(pmod->list[0] + order);

    if (newlist == NULL) {
	err = E_ALLOC;
    } else {
	newlist[0] = pmod->list[0] + order;
	for (i=2; i<=pmod->list[0]; i++) {
	    newlist[i] = pmod->list[i];
	}
	if (dataset_add_series(1 + order, pZ, pdinfo)) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add uhat to data set: substitute zeros for 
	   pre-sample values */
	for (t=0; t<n; t++) {
	    if (t < pmod->t1) {
		(*pZ)[v][t] = 0.0;
	    } else {
		(*pZ)[v][t] = pmod->uhat[t];
	    }
	}
	strcpy(pdinfo->varname[v], "uhat");
	strcpy(VARLABEL(pdinfo, v), _("residual"));
	/* then order lags of same */
	for (i=1; i<=order; i++) {
	    int s, lv = v + i;
	    double ul;

	    sprintf(pdinfo->varname[lv], "uhat_%d", i);
	    newlist[pmod->list[0] + i] = lv;
	    for (t=0; t<pdinfo->n; t++) {
		s = t - i;
		if (s < 0) {
		    (*pZ)[lv][t] = 0.0;
		} else {
		    ul = (*pZ)[v][s];
		    (*pZ)[lv][t] = (na(ul))? 0.0 : ul;
		}
	    }
	} 
    }

    /* LMF apparatus: see Kiviet, Review of Economic Studies,
       53/2, 1986, equation (5), p. 245.
    */

    if (!err) {
	/* regression on [X~E], using original sample */
	impose_model_smpl(pmod, pdinfo);
	newlist[1] = v;
	aux = lsq(newlist, pZ, pdinfo, OLS, OPT_A);
	err = aux.errcode;
	if (err) {
	   errmsg(err, prn);
	} else { 
	    RSSxe = aux.ess;
	}
    }

    if (!err) {
	int dfd = aux.nobs - pmod->ncoeff - order;
	int lberr;

	aux.aux = AUX_AR;
	gretl_model_set_int(&aux, "BG_order", order);
	trsq = aux.rsq * aux.nobs;
	LMF = ((RSSx - RSSxe) / RSSxe) * dfd / order;
	pval = snedecor_cdf_comp(order, dfd, LMF);

	if (pmod->aux != AUX_VAR) {
	    if (opt & OPT_Q) {
		bg_test_header(order, prn, 0);
	    } else {
		printmodel(&aux, pdinfo, OPT_NONE, prn);
		pputc(prn, '\n');
	    } 
	    pprintf(prn, "%s: LMF = %f,\n", _("Test statistic"), LMF);
	    pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		    order, aux.nobs - pmod->ncoeff - order, LMF, pval);

	    pprintf(prn, "\n%s: TR^2 = %f,\n", 
		    _("Alternative statistic"), trsq);
	    pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", _("with p-value"), 
		    _("Chi-square"), order, trsq, chisq_cdf_comp(order, trsq));

	    lb = ljung_box(order, pmod->t1, pmod->t2, (*pZ)[v], &lberr);
	    if (!na(lb)) {
		pprintf(prn, "Ljung-Box Q' = %g,\n", lb);
		pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n", _("with p-value"), 
			_("Chi-square"), order, lb, chisq_cdf_comp(order, lb));
	    }

	    pputc(prn, '\n');
	    record_test_result(LMF, pval, _("autocorrelation"));
	}

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_AUTOCORR);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LMF);
		model_test_set_dfn(test, order);
		model_test_set_dfd(test, aux.nobs - pmod->ncoeff - order);
		model_test_set_order(test, order);
		model_test_set_value(test, LMF);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	    
	}
    }

    free(newlist);
    dataset_drop_last_variables(pdinfo->v - v, pZ, pdinfo); 
    clear_model(&aux); 

    /* reset sample as it was */
    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    return err;
}

static int chow_active (int split, const double *x, int t)
{
    if (x != NULL) {
	return x[t] == 1.0;
    } else {
	return (t >= split);
    } 
}

/* compose list of variables to be used for the Chow test and add
   them to the data set */

static int *
make_chow_list (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		int split, int dumv, int *err)
{
    int *chowlist = NULL;
    int l0 = pmod->list[0];
    int ninter = pmod->ncoeff - pmod->ifc; /* number of interaction terms */
    int havedum = (dumv > 0);
    int newvars = ninter + 1 - havedum;
    int i, t, v = pdinfo->v;

    if (havedum && in_gretl_list(pmod->list, dumv)) {
	gretl_errmsg_sprintf(_("The model already contains %s"), 
			     pdinfo->varname[dumv]);
	*err = E_DATA;
	return NULL;
    }

    if (dataset_add_series(newvars, pZ, pdinfo)) {
	*err = E_ALLOC;
    } else {
	chowlist = gretl_list_new(pmod->list[0] + ninter + 1);
	if (chowlist == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	const double *cdum = NULL;

	for (i=1; i<=l0; i++) { 
	    chowlist[i] = pmod->list[i];
	}

	if (dumv > 0) {
	    /* we have a user-supplied dummy var */
	    cdum = (*pZ)[dumv];
	} else {
	    /* generate the split variable */
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[v][t] = (double) (t >= split); 
	    }
	    strcpy(pdinfo->varname[v], "splitdum");
	    strcpy(VARLABEL(pdinfo, v), _("dummy variable for Chow test"));
	}

	chowlist[l0 + 1] = (dumv > 0)? dumv : v;

	/* and the interaction terms */
	for (i=0; i<ninter; i++) {
	    int pv = pmod->list[i + 2 + pmod->ifc];
	    int sv = v + i + 1 - havedum;

	    for (t=0; t<pdinfo->n; t++) {
		if (model_missing(pmod, t)) {
		    (*pZ)[sv][t] = NADBL;
		} else if (chow_active(split, cdum, t)) {
		    (*pZ)[sv][t] = (*pZ)[pv][t];
		} else {
		    (*pZ)[sv][t] = 0.0;
		}
	    }

	    if (havedum) {
		sprintf(pdinfo->varname[sv], "%.2s_", pdinfo->varname[dumv]);
	    } else {
		strcpy(pdinfo->varname[sv], "sd_");
	    }
	    strncat(pdinfo->varname[sv], pdinfo->varname[pv], VNAMELEN - 4);
	    chowlist[l0 + 2 + i] = sv;
	}
    }

    return chowlist;
}

static double QLR_get_critval (double Fmax, int dfn, int *a, int *approx)
{
    double crit = 0.0;
    int i, j = dfn - 1;

    *approx = 0;

    if (j >= QLR_QMAX) {
	j = QLR_QMAX - 1;
	*approx = 1;
    } 

    for (i=2; i>=0; i--) {
	if (Fmax > QLR_critvals[j][i]) {
	    *a = (i == 2)? 1 : (i == 1)? 5 : 10;
	    crit = QLR_critvals[j][i];
	    break;
	}
    }

    if (crit == 0.0) {
	crit = QLR_critvals[j][0];
    }

    return crit;
}

static int QLR_graph (const double *Ft, int t1, int t2, 
		      int tmax, int dfn, const DATAINFO *pdinfo)
{
    const double *x = gretl_plotx(pdinfo);
    FILE *fp;
    int t, err = 0;

    fp = get_plot_input_stream(PLOT_REGULAR, &err);
    if (err) {
	return err;
    }

    print_keypos_string(GP_KEY_LEFT_TOP, fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, "plot \\\n"
	    "'-' using 1:2 title '%s' w lines\n",
	    _("Chow F-test for break"));
    for (t=t1; t<=t2; t++) {
	fprintf(fp, "%g %g\n", x[t], Ft[t-t1]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();
    
    fclose(fp);

    return gnuplot_make_graph();
}

static void save_QLR_test (MODEL *pmod, const char *datestr,
			   double Fmax, double crit, double alpha,
			   int dfn, int dfd)
{
    ModelTest *test = model_test_new(GRETL_TEST_QLR);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_SUP_WALD);
	model_test_set_param(test, datestr);
	model_test_set_value(test, Fmax);
	model_test_set_crit_and_alpha(test, crit, alpha);
	model_test_set_dfn(test, dfn);
	model_test_set_dfd(test, dfd);
	maybe_add_test_to_model(pmod, test);
    }	  
}

static void QLR_print_result (MODEL *pmod,
			      double Fmax, int tmax, int dfn, int dfd,
			      const DATAINFO *pdinfo, gretlopt opt,
			      PRN *prn)
{
    char datestr[OBSLEN];
    double crit = 0.0;
    int a = 0, approx = 0;

    ntodate(datestr, tmax, pdinfo);

    pputs(prn, _("Quandt likelihood ratio test for structural break at an "
		 "unknown point,\nwith 15 percent trimming"));
    pputs(prn, ":\n\n");
    pprintf(prn, _("The maximum F(%d, %d) = %g occurs "
		   "at observation %s"), dfn, dfd, Fmax, datestr);
    pputc(prn, '\n');

    crit = QLR_get_critval(Fmax, dfn, &a, &approx);

    if (a > 0) {
	pprintf(prn, _("Significant at the %d percent level "), a);
	pprintf(prn, _("(%d%% critical value %s %.2f)"), a, 
		(approx)? "<" : "=", crit);
    } else {
	if (approx) {
	    pprintf(prn, _("10%% critical value for q = 20 is %.2f"),
		    crit);
	} else {
	    pputs(prn, _("Not significant at the 10 percent level "));
	    pprintf(prn, _("(10%% value = %.2f)"), crit);
	}
	a = 10;
    }

    pputs(prn, "\n\n");
    pputs(prn, _("This statistic does not follow the standard "
	  "F distribution;\ncritical values are from Stock and Watson "
	  "(2003)."));
    pputs(prn, "\n\n");
	  

    if (opt & OPT_S) {
	save_QLR_test(pmod, datestr, Fmax, crit, a / 100.0,
		      dfn, dfd);
    }
}

static double robust_chow_test (MODEL *pmod, const int *list,
				int *err)
{
    double test = NADBL;
    int *tlist;

    tlist = gretl_list_diff_new(pmod->list, list, 2);

    if (tlist == NULL) {
	*err = E_ALLOC;
    } else {
	test = wald_omit_F(tlist, pmod);
	if (!na(test)) {
	    test *= tlist[0]; /* chi-square form */
	}
	free(tlist);
    }

    return test;
}

static void save_chow_test (MODEL *pmod, char *chowparm,
			    int dumv, double test, double pval,
			    int dfn, int dfd)
{
    int ttype = (dumv > 0)? GRETL_TEST_CHOWDUM : GRETL_TEST_CHOW;
    ModelTest *mt = model_test_new(ttype);

    if (mt != NULL) {
	if (dfd == 0) {
	    model_test_set_teststat(mt, GRETL_STAT_WALD_CHISQ);
	} else {
	    model_test_set_teststat(mt, GRETL_STAT_F);
	}
	model_test_set_param(mt, chowparm);
	model_test_set_value(mt, test);
	model_test_set_pvalue(mt, pval);
	model_test_set_dfn(mt, dfn);
	model_test_set_dfd(mt, dfd);
	maybe_add_test_to_model(pmod, mt);
    }	  
}

static int get_chow_dummy (const char *s, const double **Z,
			   const DATAINFO *pdinfo, int *err)
{
    int v = current_series_index(pdinfo, s);

    if (v < 0) {
	*err = E_UNKVAR;
	return 0;
    } else if (!gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[v])) {
	*err = E_DATA;
	return 0;
    }

    return v;
}

/**
 * chow_test:
 * @line: command line for parsing.
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save test results to model;
 * if %OPT_D included, do the Chow test based on a given dummy
 * variable; if %OPT_T included, do the QLR test.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability (Chow test).
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test (const char *line, MODEL *pmod, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int *chowlist = NULL;
    int origv = pdinfo->v;
    char chowparm[VNAMELEN];
    MODEL chow_mod;
    int QLR = (opt & OPT_T);
    int dumv = 0, split = 0, smax = 0;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    impose_model_smpl(pmod, pdinfo);

    gretl_model_init(&chow_mod);

    if (QLR) {
	/* 15 percent trimming */
	split = pmod->t1 + 0.15 * pmod->nobs;
	smax = pmod->t1 + 0.85 * pmod->nobs;
    } else {
	if (sscanf(line, "%*s %15s", chowparm) != 1) {
	    err = E_PARSE;
	} else if (opt & OPT_D) {
	    dumv = get_chow_dummy(chowparm, (const double **) *pZ, pdinfo, &err);
	} else {
	    split = dateton(chowparm, pdinfo);
	    if (split <= 0 || split >= pdinfo->n) { 
		gretl_errmsg_set(_("Invalid sample split for Chow test"));
		err = E_DATA;
	    }
	    smax = split;
	}
    }

    if (!err) {
	chowlist = make_chow_list(pmod, pZ, pdinfo, split, dumv, &err);
    }

    if (err) {
	goto bailout;
    }

    if (QLR) {
	/* Quandt likelihood ratio */
	double F, Fmax = 0.0;
	double *Ft = NULL;
	int dfn = 0, dfd = 0;
	int tmax = 0;
	int i, t;

	if (gretl_in_gui_mode()) {
	    Ft = malloc((smax - split + 1) * sizeof *Ft);
	}
	
	for (t=split; t<=smax; t++) {
	    chow_mod = lsq(chowlist, pZ, pdinfo, OLS, OPT_A);
	    if (chow_mod.errcode) {
		err = chow_mod.errcode;
		errmsg(err, prn);
		break;
	    }
	    dfn = chow_mod.ncoeff - pmod->ncoeff;
	    dfd = chow_mod.dfd;
	    F = (pmod->ess - chow_mod.ess) * dfd / (chow_mod.ess * dfn);
	    if (F > Fmax) {
		tmax = t;
		Fmax = F;
	    }
	    if (Ft != NULL) {
		Ft[t - split] = F;
	    }
#if 0
	    fprintf(stderr, "split at t=%d: F(%d,%d)=%g\n", t, 
		    dfn, dfd, F);
	    fprintf(stderr, " pmod->ess = %g, chow_mod.ess = %g\n", 
		    pmod->ess, chow_mod.ess);
#endif
	    clear_model(&chow_mod);
	    for (i=0; i<pmod->ncoeff; i++) {
		(*pZ)[origv+i][t] = 0.0;
	    }
	}

	if (!err) {
	    QLR_print_result(pmod, Fmax, tmax, dfn, dfd, pdinfo, opt, prn);
	    record_test_result(Fmax, NADBL, "QLR");
	    if (Ft != NULL) {
		QLR_graph(Ft, split, smax, tmax, dfn, pdinfo);
	    }
	}

	if (Ft != NULL) {
	    free(Ft);
	}
    } else {
	/* regular (or robust) Chow test */
	int robust = (pmod->opt & OPT_R);
	gretlopt lsqopt = OPT_A;

	if (robust) {
	    lsqopt |= OPT_R;
	}

	chow_mod = lsq(chowlist, pZ, pdinfo, OLS, lsqopt);

	if (chow_mod.errcode) {
	    err = chow_mod.errcode;
	    errmsg(err, prn);
	} else if (chow_mod.ncoeff <= pmod->ncoeff) {
	    err = chow_mod.errcode = E_DATA;
	    errmsg(err, prn);
	} else {
	    int dfd = (robust)? 0 : chow_mod.dfd;
	    int dfn = chow_mod.ncoeff - pmod->ncoeff;
	    double test, pval = NADBL;

	    if (!(opt & OPT_Q)) {
		chow_mod.aux = AUX_CHOW;
		printmodel(&chow_mod, pdinfo, OPT_NONE, prn);
	    }

	    if (robust) {
		test = robust_chow_test(&chow_mod, pmod->list, &err);
		if (!na(test)) {
		    pval = chisq_cdf_comp(dfn, test);
		}
	    } else {
		test = (pmod->ess - chow_mod.ess) * dfd / 
		    (chow_mod.ess * dfn);
		if (!na(test)) {
		    pval = snedecor_cdf_comp(dfn, dfd, test);
		}
	    }

	    if (!na(test) && !na(pval)) {
		if (opt & OPT_Q) {
		    pputc(prn, '\n');
		}
		if (opt & OPT_D) {
		    pprintf(prn, _("Chow test for structural difference with respect to %s"),
			    chowparm);
		} else {
		    pprintf(prn, _("Chow test for structural break at observation %s"),
			    chowparm);
		} 
		pputc(prn, '\n');

		if (robust) {
		    pprintf(prn, "  %s(%d) = %g %s %.4f\n", _("Chi-square"),
			    dfn, test, _("with p-value"), pval);
		    pprintf(prn, "  %s: F(%d, %d) = %g %s %.4f\n\n", _("F-form"), 
			    dfn, chow_mod.dfd, test / dfn, _("with p-value"),
			    snedecor_cdf_comp(dfn, chow_mod.dfd, test / dfn));
		} else {
		    pprintf(prn, "  F(%d, %d) = %g %s %.4f\n\n", 
			    dfn, dfd, test, _("with p-value"), pval);
		}

		if (opt & OPT_S) {
		    save_chow_test(pmod, chowparm, dumv, test, pval, dfn, dfd);
		}

		record_test_result(test, pval, "Chow");
	    } 
	}
	clear_model(&chow_mod);
    }

 bailout:

    /* clean up extra variables */
    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);
    free(chowlist);

    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    return err;
}

/* compute v'Mv, for symmetric M */

static double vprime_M_v (double *v, double *M, int n)
{
    int i, j, jmin, jmax, k;
    double xx, val = 0.0;

    k = jmin = 0;
    for (i=0; i<n; i++) {
	xx = 0.0;
	for (j=jmin; j<n; j++) {
	    xx += v[j] * M[k++];
	}
	val += xx * v[i];
	jmin++;
    }

    jmax = 1;
    for (i=1; i<n; i++) {
	k = i;
	xx = 0.0;
	for (j=0; j<jmax; j++) {
	    xx += v[j] * M[k];
	    k += n - j - 1;
	}
	val += xx * v[i];
	jmax++;
    }

    return val;
}

static int cusum_compute (MODEL *pmod, double *cresid, int T, int k,
			  double *wbar, double ***pZ, DATAINFO *pdinfo) 
{
    MODEL cmod;
    gretlopt opt = OPT_X | OPT_A;
    double *xvec;
    double xx;
    int n = T - k;
    int i, j, t;
    int err = 0;

    xvec = malloc(k * sizeof *xvec);
    if (xvec == NULL) {
	return E_ALLOC;
    }

    for (j=0; j<n && !err; j++) {
	cmod = lsq(pmod->list, pZ, pdinfo, OLS, opt);
	if (cmod.errcode) {
	    err = cmod.errcode;
	} else {
	    /* compute ex post prediction error */
	    t = pdinfo->t2 + 1;
	    xx = 0.0;
	    for (i=0; i<cmod.ncoeff; i++) {
		xvec[i] = (*pZ)[cmod.list[i+2]][t];
		xx += cmod.coeff[i] * xvec[i];
	    }
	    cresid[j] = (*pZ)[pmod->list[1]][t] - xx;
	    cmod.ci = CUSUM;
	    err = makevcv(&cmod, 1.0);
	    xx = vprime_M_v(xvec, cmod.vcv, cmod.ncoeff);
	    cresid[j] /= sqrt(1.0 + xx);
	    *wbar += cresid[j];
	    pdinfo->t2 += 1;
	}
	clear_model(&cmod); 
    }

    free(xvec);

    return err;
}

#define okfreq(p) (p == 1 || p == 4 || p == 12 || p == 24 || p == 52)

static int cusum_do_graph (double a, double b, const double *W, 
			   int t1, int k, int m, 
			   DATAINFO *pdinfo, gretlopt opt)
{
    FILE *fp = NULL;
    const double *obs = NULL;
    double frac = 1.0;
    double x0 = 0.0;
    int j, t, err = 0;

    fp = get_plot_input_stream(PLOT_CUSUM, &err);
    if (err) {
	return err;
    }

    if (dataset_is_time_series(pdinfo) && okfreq(pdinfo->pd)) {
	b *= pdinfo->pd;
	frac /= pdinfo->pd;
        obs = gretl_plotx(pdinfo);
	if (obs != NULL) {
	    x0 = obs[t1 + k];
	}
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xlabel '%s'\n", _("Observation"));
    fputs("set nokey\n", fp);

    if (opt & OPT_R) {
	fprintf(fp, "set title '%s'\n",
		/* xgettext:no-c-format */
		_("CUSUMSQ plot with 95% confidence band"));
	fprintf(fp, "plot \\\n%g*(x-%g) title '' w dots lt 2, \\\n", b, x0 - frac);
	fprintf(fp, "%g+%g*(x-%g) title '' w lines lt 2, \\\n", -a, b, x0 - frac);
	fprintf(fp, "%g+%g*(x-%g) title '' w lines lt 2, \\\n", a, b, x0 - frac);
    } else {
	fputs("set xzeroaxis\n", fp);
	fprintf(fp, "set title '%s'\n",
		/* xgettext:no-c-format */
		_("CUSUM plot with 95% confidence band"));
	fprintf(fp, "plot \\\n%g+%g*(x-%g) title '' w lines lt 2, \\\n", a, b, x0);
	fprintf(fp, "%g-%g*(x-%g) title '' w lines lt 2, \\\n", -a, b, x0);
    }	

    fputs("'-' using 1:2 w linespoints lt 1\n", fp);

    for (j=0; j<m; j++) { 
	t = t1 + k + j;
	if (obs != NULL) {
	    fprintf(fp, "%g %g\n", obs[t], W[j]);
	} else {
	    fprintf(fp, "%d %g\n", t, W[j]);
	}
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    err = gnuplot_make_graph();

    return err;
}

static void cusum_harvey_collier (double wbar, double sigma, int m,
				  MODEL *pmod, gretlopt opt,
				  PRN *prn)
{
    double hct, pval;

    hct = (sqrt((double) m) * wbar) / sigma;
    pval = student_pvalue_2(m - 1, hct);
    pprintf(prn, _("\nHarvey-Collier t(%d) = %g with p-value %.4g\n\n"), 
	    m - 1, hct, pval);

    if (opt & OPT_S) {
	ModelTest *test = model_test_new(GRETL_TEST_CUSUM);

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_HARVEY_COLLIER);
	    model_test_set_dfn(test, m - 1);
	    model_test_set_value(test, hct);
	    model_test_set_pvalue(test, pval);
	    maybe_add_test_to_model(pmod, test);
	}
    }

    record_test_result(hct, pval, "Harvey-Collier");
}

/**
 * cusum_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save results of test to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for parameter stability via the CUSUM test,
 * or if @opt includes %OPT_R, via the CUSUMSQ test; %OPT_Q makes
 * the test quiet.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int cusum_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn) 
{
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    char cumdate[OBSLEN];
    double wbar = 0.0;
    double *cresid = NULL, *W = NULL;
    int quiet = opt & OPT_Q;
    int m, i, j;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    if (has_missing_obs(pmod)) {
	return E_DATA;
    }

    /* number of forecasts */
    m = T - k; 

    /* set sample based on model to be tested */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t1 + k - 1;    

    cresid = malloc(m * sizeof *cresid);
    W = malloc(m * sizeof *W);

    if (cresid == NULL || W == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = cusum_compute(pmod, cresid, T, k, &wbar, pZ, pdinfo);
	if (err) {
	    errmsg(err, prn);
	}
    }

    if (!err) {
	double a, b, den = 0.0, sigma = 0.0;
	int sig;

	if (opt & OPT_R) {
	    double a1, a2, a3;
	    double n = 0.5 * m - 1;

	    pprintf(prn, "\n%s\n\n", _("CUSUMSQ test for stability of parameters"));

	    for (j=0; j<m; j++) {
		den += cresid[j] * cresid[j];
	    }

	    /* see Edgerton and Wells, Oxford Bulletin of Economics and
	       Statistics, 56, 1994, pp. 355-365 */
	    a1 = 1.358015 / sqrt(n);
	    a2 = -0.6701218 / n;
	    a3 = -0.8858694 / pow(n, 1.5);

	    /* 0.5 * total height of band */
	    a = a1 + a2 + a3;
	    /* slope of expectation wrt time */
	    b = 1.0 / m;
	    if (!quiet) {
		pputs(prn, _("Cumulated sum of squared residuals"));
	    }
	} else {
	    wbar /= T - k;
	    pprintf(prn, "\n%s\n\n", _("CUSUM test for stability of parameters"));
	    pprintf(prn, _("mean of scaled residuals = %g\n"), wbar);

	    for (j=0; j<m; j++) {
		sigma += (cresid[j] - wbar) * (cresid[j] - wbar);
	    }
	    sigma /= T - k - 1;
	    sigma = sqrt(sigma);
	    pprintf(prn, _("sigmahat                 = %g\n\n"), sigma);

	    /* height of confidence band for first prediction */
	    a = 0.948 * sqrt((double) m);
	    /* slope of confidence band limit wrt time */
	    b = 2.0 * a / m;
	    if (!quiet) {
		pputs(prn, _("Cumulated sum of scaled residuals"));
	    }
	}

	pputc(prn, '\n');
	pputs(prn, /* xgettext:no-c-format */
	      _("('*' indicates a value outside of 95% confidence band)"));
	pputs(prn, "\n\n");
    
	for (j=0; j<m; j++) {
	    W[j] = 0.0;
	    if (opt & OPT_R) {
		for (i=0; i<=j; i++) {
		    W[j] += cresid[i] * cresid[i] / den;
		}
		sig = fabs(W[j] - (j+1) / (double) m) > a;
	    } else {
		for (i=0; i<=j; i++) {
		    W[j] += cresid[i];
		}
		W[j] /= sigma;
		sig = fabs(W[j]) > a + j * b;
	    }
	    if (!quiet) {
		ntodate(cumdate, pmod->t1 + k + j, pdinfo);
		pprintf(prn, " %s %9.3f %s\n", cumdate, W[j], sig? "*" : "");
	    }
	}

	if (!(opt & OPT_R)) {
	    cusum_harvey_collier(wbar, sigma, m, pmod, opt, prn);
	}

	/* plot with 95% confidence bands if not in batch mode */
	if (!gretl_in_batch_mode()) {
	    err = cusum_do_graph(a, b, W, pmod->t1, k, m, pdinfo, opt);
	}
    }

    /* restore original sample range */
    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;
    
    free(cresid);
    free(W);

    return err;
}

/**
 * comfac_test:
 * @pmod: pointer to original model.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt: if contains %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * If @pmod was estimated via an AR(1) estimator, run an
 * auxiliary regression to test the implied common-factor
 * restriction.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int comfac_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    MODEL cmod;
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int v = pdinfo->v;
    int *biglist = NULL;
    int clearit = 0;
    int nadd, i, k, t;
    int err;

    if (pmod->ci != AR1 || (pmod->opt & OPT_P)) {
	/* can't handle Prais-Winsten? */
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod->list, pdinfo, pmod->ID);
    if (err) {
	return err;
    }

    biglist = gretl_list_copy(pmod->list);
    if (biglist == NULL) {
	return E_ALLOC;
    }

    nadd = 1 + pmod->ncoeff - pmod->ifc;

    err = dataset_add_series(nadd, pZ, pdinfo);
    if (err) {
	free(biglist);
	return err;
    }

    /* add lags of the dependent variable and all regressors: some of
       these may be redundant but we'll just let the redundat terms be
       eliminated automatically via gretl's collinearity checking.
    */

    k = v;
    for (i=1; i<=pmod->list[0]; i++) {
	int src, lag, parent;

	src = pmod->list[i];
	if (src == 0) {
	    continue;
	}
	for (t=0; t<pdinfo->n; t++) {
	    if (t == 0 || na((*pZ)[src][t-1])) {
		(*pZ)[k][t] = NADBL;
	    } else {
		(*pZ)[k][t] = (*pZ)[src][t-1];
	    }
	}
	biglist = gretl_list_append_term(&biglist, k);
	if (biglist == NULL) {
	    err = E_ALLOC;
	    break;
	}
	lag = is_standard_lag(src, pdinfo, &parent);
	if (lag && parent) {
	    char tmp[8];

	    sprintf(tmp, "_%d", lag + 1);
	    strcpy(pdinfo->varname[k], pdinfo->varname[parent]);
	    gretl_trunc(pdinfo->varname[k], 15 - strlen(tmp));
	    strcat(pdinfo->varname[k], tmp);
	} else {
	    strcpy(pdinfo->varname[k], pdinfo->varname[src]);
	    gretl_trunc(pdinfo->varname[k], 13);
	    strcat(pdinfo->varname[k], "_1");
	}
	k++;
    }

    if (!err) {
	/* re-impose the sample that was in force when the original model
	   was estimated */
	impose_model_smpl(pmod, pdinfo);
	cmod = lsq(biglist, pZ, pdinfo, OLS, OPT_A);
	clearit = 1;
	err = cmod.errcode;
    }

    if (!err) {
	if (cmod.nobs != pmod->nobs || cmod.ess > pmod->ess || cmod.dfd >= pmod->dfd) {
	    /* something has gone wrong */
	    err = E_DATA;
	}
    }

    if (!err) {
	/* construct an F-test based on the SSR from the original
	   AR(1) model and the SSR from the unrestricted model, cmod.
	*/
	int dfd = cmod.dfd;
	int dfn = pmod->dfd - dfd - 1; /* account for rho */
	double SSRr = pmod->ess;
	double SSRu = cmod.ess;
	double Ftest = ((SSRr - SSRu)/dfn) / (SSRu/dfd);
	double pval = snedecor_cdf_comp(dfn, dfd, Ftest);

	if (!(opt & OPT_Q)) {
	    cmod.aux = AUX_COMFAC;
	    printmodel(&cmod, pdinfo, OPT_S, prn);
	    pputc(prn, '\n');
	}

	pputs(prn, _("Test of common factor restriction"));
	pputs(prn, "\n\n");

	pprintf(prn, "  %s: %s(%d, %d) = %g, ", _("Test statistic"), 
		"F", dfn, dfd, Ftest);
	pprintf(prn, _("with p-value = %g\n"), pval);
	pputc(prn, '\n');

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = model_test_new(GRETL_TEST_COMFAC);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_F);
		model_test_set_dfn(test, dfn);
		model_test_set_dfd(test, dfd);
		model_test_set_value(test, Ftest);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }
	}

	record_test_result(Ftest, pval, _("common factor restriction"));
    }

    if (clearit) {
	clear_model(&cmod);
    }

    /* delete the added variables and restore the original
       sample range */

    dataset_drop_last_variables(nadd, pZ, pdinfo);   
    free(biglist);

    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    return err;
}

/**
 * panel_hausman_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Tests the given pooled model for fixed and random effects.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int panel_hausman_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn) 
{
    if (pmod->ci != OLS || !dataset_is_panel(pdinfo)) {
	pputs(prn, _("This test is only relevant for pooled models\n"));
	return 1;
    }

    if (pmod->ifc == 0) {
	pputs(prn, _("This test requires that the model contains a constant\n"));
	return 1;
    }

    if (!balanced_panel(pdinfo)) { /* ?? */
	pputs(prn, _("Sorry, can't do this test on an unbalanced panel.\n"
		"You need to have the same number of observations\n"
		"for each cross-sectional unit"));
	return 1;
    } else {
	panel_diagnostics(pmod, pZ, pdinfo, opt, prn);
    }

    return 0;
}

/**
 * add_leverage_values_to_dataset:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @m: matrix containing leverage values.
 * @flags: option flags: combination of %SAVE_LEVERAGE, %SAVE_INFLUENCE,
 * and %SAVE_DFFITS.
 *
 * Adds to the working dataset one or more series calculated by
 * the gretl test for leverage/influence of data points.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int add_leverage_values_to_dataset (double ***pZ, DATAINFO *pdinfo,
				    gretl_matrix *m, unsigned char flags)
{
    int t1, t2;
    int addvars = 0;

    if (flags & SAVE_LEVERAGE) addvars++;
    if (flags & SAVE_INFLUENCE) addvars++;
    if (flags & SAVE_DFFITS) addvars++;

    if (dataset_add_series(addvars, pZ, pdinfo)) {
	return E_ALLOC;
    }

    t1 = gretl_matrix_get_t1(m);
    t2 = t1 + gretl_matrix_rows(m);

    /* add leverage? */
    if (flags & SAVE_LEVERAGE) {
	int t, v = pdinfo->v - addvars;
	int j = 0;

	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t >= t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = gretl_matrix_get(m, j++, 0);
	    }
	}
	strcpy(pdinfo->varname[v], "lever");
	make_varname_unique(pdinfo->varname[v], v, pdinfo);
	strcpy(VARLABEL(pdinfo, v), "leverage values");
    }

    /* add influence? */
    if (flags & SAVE_INFLUENCE) {
	int t, v = pdinfo->v - (addvars - 1);
	int j = 0;

	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t >= t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = gretl_matrix_get(m, j++, 1);
	    }
	}	
	strcpy(pdinfo->varname[v], "influ");
	make_varname_unique(pdinfo->varname[v], v, pdinfo);
	strcpy(VARLABEL(pdinfo, v), "influence values");
    }

    /* add DFFITS? */
    if (flags & SAVE_DFFITS) {
	int t, v = pdinfo->v - (addvars - 2);
	int j = 0;

	for (t=0; t<pdinfo->n; t++) {
	    double s, h;

	    if (t < t1 || t >= t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		/* s = studentized residuals */
		h = gretl_matrix_get(m, j, 0);
		s = gretl_matrix_get(m, j, 2);
		if (na(h) || na(s)) {
		    (*pZ)[v][t] = NADBL;
		} else {
		    (*pZ)[v][t] = s * sqrt(h / (1.0 - h));
		}
		j++;
	    }
	}	
	strcpy(pdinfo->varname[v], "dffits");
	make_varname_unique(pdinfo->varname[v], v, pdinfo);
	strcpy(VARLABEL(pdinfo, v), "DFFITS values");
    }

    return 0;
}

/**
 * leverage_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if %OPT_S, add calculated series to data set.
 * @prn: gretl printing struct.
 *
 * Tests the data used in the given model for points with
 * high leverage and influence on the estimates.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int leverage_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn)
{
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, double ***, 
				     const DATAINFO *, gretlopt,
				     PRN *, int *);
    gretl_matrix *m;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    model_leverage = get_plugin_function("model_leverage", &handle);
    if (model_leverage == NULL) {
	return 1;
    }

    m = (*model_leverage)(pmod, pZ, pdinfo, OPT_NONE, prn, &err);

    if (!err && (opt & OPT_S)) {
	err = add_leverage_values_to_dataset(pZ, pdinfo, m, 
					     SAVE_LEVERAGE |
					     SAVE_INFLUENCE| 
					     SAVE_DFFITS);
    }

    gretl_matrix_free(m);

    close_plugin(handle);

    return err;
}

/**
 * vif_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Calculates and displays the Variance Inflation Factors for
 * the independent variables in the given model.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int vif_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    void *handle;
    int (*print_vifs) (MODEL *, double ***, DATAINFO *, PRN *);
    int err;

    gretl_error_clear();

    print_vifs = get_plugin_function("print_vifs", &handle);
    if (print_vifs == NULL) {
	return 1;
    }

    err = (*print_vifs)(pmod, pZ, pdinfo, prn);

    close_plugin(handle);

    return err;
}

/**
 * model_test_driver:
 * @param: auxiliary parameter for some uses.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: controls which test(s) will be performed; %OPT_Q
 * gives less verbose results.
 * @prn: gretl printing struct.
 * 
 * Performs some subset of gretl's "modtest" tests on the
 * model last estimated, and prints the results to @prn.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int model_test_driver (const char *param, 
		       double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    GretlObjType type;
    gretlopt testopt;
    void *ptr;
    int k = 0;
    int err = 0;

    if (opt == OPT_NONE || opt == OPT_Q) {
	pprintf(prn, "modtest: no options selected\n");
	return 0;
    }

    err = incompatible_options(opt, OPT_A | OPT_H | OPT_L | OPT_S |
			       OPT_N | OPT_P | OPT_W | OPT_X);
    if (err) {
	return err;
    }

    ptr = get_last_model(&type);  
    if (ptr == NULL) {
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN && exact_fit_check(ptr, prn)) {
	return 0;
    }

    if (opt & (OPT_A | OPT_H)) {
	/* autocorrelation and arch: lag order */
	k = atoi(param);
	if (k == 0) {
	    k = pdinfo->pd;
	}
    }

    testopt = (opt & OPT_Q)? OPT_Q : OPT_NONE;

    /* non-linearity (squares) */
    if (!err && (opt & OPT_S)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, pZ, pdinfo, 
				    AUX_SQ, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* non-linearity (logs) */
    if (!err && (opt & OPT_L)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, pZ, pdinfo, 
				    AUX_LOG, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* heteroskedasticity (White or Breusch-Pagan) */
    if (!err && (opt & (OPT_W | OPT_X | OPT_B))) {
	if (type == GRETL_OBJ_EQN) {
	    transcribe_option_flags(&testopt, opt, OPT_B | OPT_X);
	    if ((opt & OPT_B) && (opt & OPT_R)) {
		testopt |= OPT_R;
	    }
	    err = whites_test(ptr, pZ, pdinfo, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* autocorrelation */
    if (!err && (opt & OPT_A)) {
	if (type == GRETL_OBJ_EQN) {
	    err = autocorr_test(ptr, k, pZ, pdinfo, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = gretl_VAR_autocorrelation_test(ptr, k, pZ, pdinfo, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_autocorrelation_test(ptr, k, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* ARCH */
    if (!err && (opt & OPT_H)) {
	if (type == GRETL_OBJ_EQN) {
	    err = arch_test(ptr, k, pdinfo, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = gretl_VAR_arch_test(ptr, k, pdinfo, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_arch_test(ptr, k, prn);
	} else {
	    err = E_NOTIMP;
	}
    }    

    /* normality of residual */
    if (!err && (opt & OPT_N)) {
	err = last_model_test_uhat(pZ, pdinfo, testopt, prn);
    }

    /* groupwise heteroskedasticity */
    if (!err && (opt & OPT_P)) {
	if (type == GRETL_OBJ_EQN) {
	    err = groupwise_hetero_test(ptr, pdinfo, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* common factor restriction */
    if (!err && (opt & OPT_C)) {
	if (type == GRETL_OBJ_EQN) {
	    err = comfac_test(ptr, pZ, pdinfo, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }    

    return err;
}

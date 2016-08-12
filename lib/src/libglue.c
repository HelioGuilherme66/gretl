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

/* Miscellaneous functions to bridge between gretl commands and
   the corresponding libgretl functions. This glue allows a cleaner
   interface for the latter, hiding some parsing of strings
   coming from the command line.
*/

#include "libgretl.h"
#include "var.h"
#include "system.h"
#include "gretl_panel.h"
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "boxplots.h"
#include "libglue.h"

/*
 * model_test_driver:
 * @order: lag order for --autocorr and --arch.
 * @dset: dataset struct.
 * @opt: controls which test(s) will be performed; OPT_Q
 * gives less verbose results, OPT_I gives silent operation.
 * @prn: gretl printing struct.
 * 
 * Performs some subset of gretl's "modtest" tests on the
 * model last estimated, and prints the results to @prn.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int model_test_driver (int order, DATASET *dset, 
		       gretlopt opt, PRN *prn)
{
    GretlObjType type;
    gretlopt testopt = OPT_NONE;
    void *ptr;
    int k = 0;
    int err = 0;

    if (opt == OPT_NONE || opt == OPT_Q || opt == OPT_I) {
	/* note: OPT_Q and OPT_I are just quiet and silent respectively */
	pprintf(prn, "modtest: no options selected\n");
	return 0;
    }

    err = incompatible_options(opt, OPT_A | OPT_H | OPT_L | OPT_S |
			       OPT_N | OPT_P | OPT_W | OPT_X | OPT_D);
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
	k = order > 0 ? order : dset->pd;
    }

    /* transcribe the quietness flags */
    if (opt & OPT_I) {
	testopt = OPT_I | OPT_Q;
    } else if (opt & OPT_Q) {
	testopt = OPT_Q;
    }

    /* non-linearity (squares) */
    if (!err && (opt & OPT_S)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, dset, AUX_SQ, 
				    testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* non-linearity (logs) */
    if (!err && (opt & OPT_L)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, dset, AUX_LOG, 
				    testopt, prn);
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
	    err = whites_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* autocorrelation */
    if (!err && (opt & OPT_A)) {
	if (type == GRETL_OBJ_EQN) {
	    err = autocorr_test(ptr, k, dset, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = gretl_VAR_autocorrelation_test(ptr, k, dset, 
						 testopt, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_autocorrelation_test(ptr, k, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* ARCH */
    if (!err && (opt & OPT_H)) {
	if (type == GRETL_OBJ_EQN) {
	    err = arch_test(ptr, k, dset, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = gretl_VAR_arch_test(ptr, k, dset, 
				      testopt, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_arch_test(ptr, k, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }    

    /* normality of residual */
    if (!err && (opt & OPT_N)) {
	err = last_model_test_uhat(dset, testopt, prn);
    }

    /* groupwise heteroskedasticity */
    if (!err && (opt & OPT_P)) {
	if (type == GRETL_OBJ_EQN) {
	    err = groupwise_hetero_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* common factor restriction */
    if (!err && (opt & OPT_C)) {
	if (type == GRETL_OBJ_EQN) {
	    err = comfac_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* cross-sectional dependence */
    if (!err && (opt & OPT_D)) {
	if (type == GRETL_OBJ_EQN) {
	    err = panel_xdepend_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }    

    return err;
}

static int get_chow_dummy (const char *s, const DATASET *dset, 
			   int *err)
{
    int v = current_series_index(dset, s);

    if (v < 0) {
	*err = E_UNKVAR;
    } else if (!gretl_isdummy(dset->t1, dset->t2, dset->Z[v])) {
	*err = E_DATA;
    }

    return v;
}

/*
 * chow_test_driver:
 * @param: parameter (observation or name of dummy)
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_D included, do the Chow test based on a given dummy
 * variable.
 * @prn: gretl printing struct.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test_driver (const char *param, MODEL *pmod, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    int chowparm = 0;
    int err = 0;

    if (param == NULL || *param == '\0') {
	return E_DATA;
    }

    if (opt & OPT_D) {
	chowparm = get_chow_dummy(param, dset, &err);
    } else {
	chowparm = dateton(param, dset);
    }

    if (!err) {
	if (opt & OPT_D) {
	    err = chow_test_from_dummy(chowparm, pmod, dset, opt, prn);
	} else {
	    err = chow_test(chowparm, pmod, dset, opt, prn);
	}
    }

    return err;
}

/* The @param here may contain a scalar or a matrix: in either case,
   convert to a list of lag-orders before handing off to the real
   Levin-Lin-Chu code.
*/

int llc_test_driver (const char *param, const int *list, 
		     DATASET *dset, gretlopt opt, PRN *prn)
{
    gretl_matrix *m = NULL;
    int *plist = NULL;
    int p0 = -1;
    int err = 0;

    if (param == NULL) {
	err = E_DATA;
    } else if (*param == '{') {
	m = generate_matrix(param, dset, &err);
	if (!err) {
	    plist = gretl_list_from_vector(m, &err);
	}
	gretl_matrix_free(m);
    } else if (gretl_is_matrix(param)) {
	m = get_matrix_by_name(param);
	plist = gretl_list_from_vector(m, &err);
    } else if (integer_string(param)) {
	p0 = atoi(param);
    } else if (gretl_is_scalar(param)) {
	p0 = gretl_scalar_get_value(param, NULL);
    } else {
	err = E_DATA;
    }

    if (!err) {
	if (plist != NULL) {
	    err = levin_lin_test(list[1], plist, dset, opt, prn);
	    free(plist);
	} else {
	    int tmplist[2] = {1, p0};

	    err = levin_lin_test(list[1], tmplist, dset, opt, prn);
	}
    }

    return err;
}

/* parse the tau vector out of @param before calling the
   "real" quantreg function
*/

MODEL quantreg_driver (const char *param, const int *list, 
		       DATASET *dset, gretlopt opt, PRN *prn)
{
    gretl_vector *tau;
    MODEL mod;
    int err = 0;

    tau = generate_matrix(param, dset, &err);

    if (!err && gretl_vector_get_length(tau) == 0) {
	err = E_DATA;
    }

    if (err) {
	gretl_model_init(&mod, dset);
	mod.errcode = err;
    } else {
	mod = quantreg(tau, list, dset, opt, prn);
    }

    gretl_matrix_free(tau);

    return mod;
}

/* wrapper for the various sorts of logit and probit models
   that gretl supports
*/

MODEL logit_probit (int *list, DATASET *dset, int ci, 
		    gretlopt opt, PRN *prn)
{
    int yv = list[1];

    if (ci == LOGIT && (opt & OPT_M)) {
	return multinomial_logit(list, dset, opt, prn);
    } else if (ci == PROBIT && (opt & OPT_E)) {
	return reprobit_model(list, dset, opt, prn);
    } else if (gretl_isdummy(dset->t1, dset->t2, dset->Z[yv])) {
	if (ci == LOGIT) {
	    return binary_logit(list, dset, opt, prn);
	} else {
	    return binary_probit(list, dset, opt, prn);
	}
    } else {
	if (ci == LOGIT) {
	    return ordered_logit(list, dset, opt, prn);
	} else {
	    return ordered_probit(list, dset, opt, prn);
	}
    } 
}

/* parse out optional "ymax=..." parameter before calling the real
   logistic model function 
*/

MODEL logistic_driver (const int *list, DATASET *dset,
		       gretlopt opt) 
{
    double lmax;
    int err = 0;

    lmax = get_optval_double(LOGISTIC, OPT_M, &err);
    
    if (err) {
	MODEL mdl;

	gretl_model_init(&mdl, dset);
	mdl.errcode = err;
	return mdl;
    }

    return logistic_model(list, lmax, dset);
}

/* assemble the left and right limits for tobit using gretl's
   option apparatus before calling the real tobit function
*/

MODEL tobit_driver (const int *list, DATASET *dset, 
		    gretlopt opt, PRN *prn)
{
    MODEL model;
    double llim = -1.0e300;
    double rlim = NADBL;
    int err = 0;

    if (opt & OPT_L) {
	/* we should have an explicit lower limit */
	llim = get_optval_double(TOBIT, OPT_L, &err);
	if (!err && na(llim)) {
	    err = E_INVARG;
	} 
    }

    if (!err && (opt & OPT_M)) {
	/* we should have an explicit upper limit */
	rlim = get_optval_double(TOBIT, OPT_M, &err);
	if (!err && (na(rlim) || rlim <= llim)) {
	    err = E_INVARG; 
	}	
    }

    if (err) {
	gretl_model_init(&model, dset);
	model.errcode = err;
	return model;
    }

    if (!(opt & (OPT_L | OPT_M))) {
	/* the default: left-censoring at zero */
	llim = 0;
    }

    return tobit_model(list, llim, rlim, dset, opt, prn);
}

/*
 * do_modprint:
 * @line: command line.
 * @opt: may contain %OPT_O for specifying output, and if
 * TeX output is called for then %OPT_C calls for
 * a complete LaTeX document.
 * @prn: gretl printer.
 *
 * Prints to @prn the coefficient table and optional additional statistics
 * for a model estimated "by hand". Mainly useful for user-written functions.
 * 
 * The string @line must contain, in order: (1) the name of a k x 2 matrix
 * containing k coefficients and k associated standard errors and (2) the
 * name of a string variable containing at least k comma- or space-
 * separated names for the coefficients (or a string literal on that 
 * pattern). 
 *
 * Optionally, @line may contain a third element, the name of a vector 
 * containing p additional statistics.  In that case element (2) should 
 * contain k + p names, the additional p names to be associated with the 
 * additional statistics. 
 *
 * Returns: 0 on success, non-zero on failure.
 */

int do_modprint (const char *mname, const char *names, 
		 gretlopt opt, PRN *prn)
{
    gretl_matrix *coef_se = NULL;
    gretl_matrix *addstats = NULL;
    const char *parnames = NULL;
    int err = 0;

    if (mname == NULL || names == NULL) {
	return E_ARGS;
    }

    /* first: name of k x 2 matrix */
    coef_se = get_matrix_by_name(mname);
    if (coef_se == NULL) {
	err = E_UNKVAR;
    } else if (gretl_matrix_cols(coef_se) != 2) {
	gretl_errmsg_set(_("modprint: the first matrix argument must have 2 columns"));
	err = E_DATA;
    }
 
    if (!err) {
	/* second: string containing names */
	if (opt & OPT_L) {
	    /* treat as string _L_iteral */
	    parnames = names;
	} else {
	    /* FIXME accept array of strings */
	    parnames = get_string_by_name(names);
	    if (parnames == NULL) {
		err = E_PARSE;
	    }
	}
    }

    if (!err && (opt & OPT_A)) {
	/* optional third field: extra matrix */
	const char *aname = get_optval_string(MODPRINT, OPT_A);

	if (aname != NULL) {
	    addstats = get_matrix_by_name(aname);
	    if (addstats == NULL) {
		err = E_UNKVAR;
	    }	    
	}
    }

    if (!err) {
	PrnFormat fmt = GRETL_FORMAT_TXT;
	char fname[FILENAME_MAX];

	*fname = '\0';

	if (opt & OPT_O) {
	    /* try for --output=filename, and if found let
	       the suffix determine the output type
	    */
	    const char *s = get_optval_string(MODPRINT, OPT_O);

	    if (s != NULL && *s != '\0') {
		strcpy(fname, s);
		if (has_suffix(fname, ".tex")) {
		    fmt = GRETL_FORMAT_TEX;
		    if (opt & OPT_C) {
			fmt |= GRETL_FORMAT_DOC;
		    }		    
		} else if (has_suffix(fname, ".rtf")) {
		    fmt = GRETL_FORMAT_RTF;
		} else if (has_suffix(fname, ".csv")) {
		    fmt = GRETL_FORMAT_CSV;
		}
	    }
	}	

	if (*fname != '\0') {
	    PRN *myprn;
	    
	    gretl_maybe_switch_dir(fname);
	    myprn = gretl_print_new_with_filename(fname, &err);
	    if (!err) {
		gretl_print_set_format(myprn, fmt);
		err = print_model_from_matrices(coef_se, addstats,
						parnames, myprn);
		gretl_print_destroy(myprn);
	    }
	} else {
	    gretl_print_set_format(prn, fmt);
	    err = print_model_from_matrices(coef_se, addstats, parnames, prn);
	}
    }

    return err;
}

int *matrix_bandplot_biglist (int ci,
			      const gretl_matrix *m,
			      const int *list,
			      int *err)
{
    const char *s = get_optval_string(ci, OPT_N);
    gchar **S = NULL;
    int *biglist = NULL;
    int ccol = 0, wcol = 0;
    int c, i;

    if (s == NULL) {
	*err = E_INVARG;
	return NULL;
    }

    S = g_strsplit(s, ",", -1);

    for (i=0; i<2 && !*err; i++) {
	c = 0;
	if (S[i] == NULL) {
	    *err = E_DATA;
	} else if (integer_string(S[i])) {
	    c = atoi(S[i]);
	} else {
	    c = get_scalar_value_by_name(S[i], err);
	}
	if (!*err && c >= 1 && c <= m->cols) {
	    if (i == 0) {
		ccol = c;
	    } else {
		wcol = c;
	    }
	} else {
	    c = 0;
	}	
	if (!*err && c == 0) {
	    *err = invalid_field_error(S[i]);
	}
    }

    g_strfreev(S);

    if (!*err) {
	biglist = gretl_list_copy(list);
	gretl_list_append_term(&biglist, ccol);
	gretl_list_append_term(&biglist, wcol);
    }

    return biglist;
}

int matrix_command_driver (int ci, 
			   const int *list, 
			   const char *param,
			   const DATASET *dset, 
			   gretlopt opt,
			   PRN *prn)
{
    gretl_matrix *m = NULL;
    DATASET *mdset = NULL;
    int *collist = NULL;
    const char *mname;
    int cmax = 0;
    int err = 0;

    mname = get_optval_string(ci, OPT_X);

    if (mname != NULL) {
	m = get_matrix_by_name(mname);
    }

    if (gretl_is_null_matrix(m)) {
	return E_DATA;
    }

    if (ci == GNUPLOT && (opt & OPT_N)) {
	/* --band=... */
	int *biglist = matrix_bandplot_biglist(ci, m, list, &err);

	if (!err) {
	    mdset = gretl_dataset_from_matrix(m, biglist, OPT_B, &err);
	    cmax = mdset->v - 3;
	    free(biglist);
	} 
    } else if (ci == SCATTERS) {
	/* note: this is a special case, for now */
	return matrix_scatters(m, list, dset, opt);
    } else if (list != NULL && list[0] == 0) {
	/* use all columns of the matrix */
	mdset = gretl_dataset_from_matrix(m, NULL, OPT_B, &err);
    } else if (list != NULL && list[0] == 1 && ci == SUMMARY) {
	/* summary stats for a single specified column */
	mdset = gretl_dataset_from_matrix(m, list, OPT_B | OPT_N, &err);
    } else {
	/* note that a NULL list is OK here */
	mdset = gretl_dataset_from_matrix(m, list, OPT_B, &err);
    }

    if (!err) {
	if (cmax == 0) {
	    cmax = mdset->v - 1;
	}
	dataset_set_matrix_name(mdset, mname);
	collist = gretl_consecutive_list_new(1, cmax);
	if (collist == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	if (ci != GNUPLOT) {
	    opt &= ~OPT_X;
	}
	if (ci == BXPLOT) {
	    err = boxplots(collist, param, mdset, opt);
	} else if (ci == GNUPLOT) {
	    err = gnuplot(collist, param, mdset, opt);
	} else if (ci == SUMMARY) {
	    err = list_summary(collist, 0, mdset, opt, prn);
	} else {
	    err = E_DATA;
	}
    }

    destroy_dataset(mdset);   
    free(collist);

    return err;
}

int matrix_freq_driver (const int *list,
			gretlopt opt,
			PRN *prn)
{
    gretl_matrix *m = NULL;
    DATASET *mdset = NULL;
    const char *mname;
    int err = 0;

    if (list != NULL && list[0] != 1) {
	return E_DATA;
    }

    mname = get_optval_string(FREQ, OPT_X);

    if (mname != NULL) {
	m = get_matrix_by_name(mname);
    }

    if (gretl_is_null_matrix(m)) {
	err = E_DATA;
    } else {
	if (list == NULL) {
	    /* this is OK if m is a column vector */
	    if (m->cols == 1) {
		int mlist[2] = {1, 1};
		
		mdset = gretl_dataset_from_matrix(m, mlist, OPT_B, &err);
	    } else {
		err = E_ARGS;
	    }
	} else {
	    mdset = gretl_dataset_from_matrix(m, list, OPT_B, &err);
	}
    }

    if (!err) {
	err = freqdist(1, mdset, opt, prn);
    }

    destroy_dataset(mdset);   

    return err;
}

int list_summary_driver (const int *list, const DATASET *dset, 
			 gretlopt opt, PRN *prn)
{
    int wtvar = 0;
    int err = 0;

    if (opt & OPT_W) {
	const char *wname = get_optval_string(SUMMARY, OPT_W);

	if (wname == NULL) {
	    err = E_DATA;
	} else {
	    wtvar = current_series_index(dset, wname);
	    if (wtvar < 0) {
		err = E_UNKVAR;
	    }
	}
    }

    if (!err) {
	err = list_summary(list, wtvar, dset, opt, prn);
    }

    return err;
}



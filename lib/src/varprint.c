/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#include "libgretl.h" 
#include "var.h"  
#include "varprint.h"
#include "libset.h"

/**
 * gretl_VAR_print_VCV:
 * @var:
 * @prn:
 *
 *
 * Returns:
 */

int gretl_VAR_print_VCV (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->S == NULL) {
	err = 1;
    } else {
	print_contemp_covariance_matrix(var->S, var->ldet, prn);
    }

    return err;
}

static void tex_print_double (double x, PRN *prn)
{
    char number[16];

    x = screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);

    if (x < 0.) {
	pprintf(prn, "$-$%s", number + 1);
    } else {
	pputs(prn, number);
    }
}

/* printing of impulse responses and variance decompositions */

#define IRF_ROW_MAX 4
#define VDC_ROW_MAX 5

enum {
    IRF,
    VDC
};

static void VAR_RTF_row_spec (int ncols, PRN *prn)
{
    int lcol = 800, colwid = 1600;
    int i, cellx = lcol;

    pputs(prn, "{\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");
    for (i=0; i<ncols; i++) {
	pprintf(prn, "\\cellx%d", cellx);
	cellx += colwid;
    }
    pputc(prn, '\n');
}

static void VAR_info_header_block (int code, int v, int block, 
				   const DATAINFO *pdinfo, 
				   PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[16];

    if (tex) {
	pputs(prn, "\\vspace{1em}\n\n");
	if (code == IRF) {
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    tex_escape(vname, pdinfo->varname[v]));
	} else {
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    tex_escape(vname, pdinfo->varname[v]));
	}
	if (block == 0) {
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, " (%s)\n\n", I_("continued"));
	}
	pprintf(prn, "\\vspace{1em}\n\n\\begin{longtable}{%s}\n",
		(code == IRF)? "rrrrr" : "rrrrrr");
    } else if (rtf) {
	pputs(prn, "\\par\n\n");
	if (code == IRF) {
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[v]);
	} else {
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    pdinfo->varname[v]);
	}
	if (block == 0) {
	    pputs(prn, "\\par\n\n");
	} else {
	    pprintf(prn, " (%s)\\par\n\n", I_("continued"));
	}
	/* FIXME */
	VAR_RTF_row_spec((code == IRF)? IRF_ROW_MAX : VDC_ROW_MAX, prn);
    } else {
	if (code == IRF) {	
	    pprintf(prn, _("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[v]);
	} else {
	    pprintf(prn, _("Decomposition of variance for %s"), 
		    pdinfo->varname[v]);
	}
	if (block == 0) {
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, " (%s)\n\n", _("continued"));
	}
    }

    /* first column: period number header */
    if (tex) {
	pprintf(prn, "%s & ", I_("period"));
    } else if (rtf) {
	pprintf(prn, "\\intbl \\qc %s\\cell ", I_("period"));
    } else {
	pputs(prn, _("period"));
    }
}

static void VAR_info_print_vname (int code, int i, int v, int endrow, 
				  const DATAINFO *pdinfo, PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[16];

    if (tex) {
	pprintf(prn, " %s ", tex_escape(vname, pdinfo->varname[v]));
	if (endrow) {
	   pputs(prn, "\\\\");
	} else { 
	    pputs(prn, "& ");
	}
    } else if (rtf) {
	pprintf(prn, "\\qc %s\\cell", pdinfo->varname[v]);
	if (endrow) {
	    pputs(prn, " \\intbl \\row");
	} 
    } else {
	int w = (code == IRF)? 13 : 11;

	if (code == IRF && i == 0) w--; /* FIXME width for i18n */
	pprintf(prn, "%*s", w, pdinfo->varname[v]);
    }
}

static void VAR_info_print_period (int t, PRN *prn)
{
    if (tex_format(prn)) {
	pprintf(prn, "%d & ", t);
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\intbl \\qc %d\\cell ", t);
    } else {
	pprintf(prn, " %3d  ", t);
    }
}

static void VAR_info_end_row (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\intbl \\row\n");
    } else {
	pputc(prn, '\n');
    }
}

static void VAR_info_end_table (PRN *prn)	
{
    if (tex_format(prn)) {
	pputs(prn, "\\end{longtable}\n\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "}\n");
    } else {
	pputc(prn, '\n');
    }
}

#define IRF_ROW_MAX 4
#define VDC_ROW_MAX 5

/**
 * gretl_VAR_print_impulse_response:
 * @var: pointer to VAR struct.
 * @shock:
 * @periods: number of periods over which to print response.
 * @pdinfo: dataset information.
 * @pause: if non-zero, pause between sections of output.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  int pause, PRN *prn)
{
    int i, t;
    int vsrc;
    int rows = var->neqns * (var->order + var->ecm);
    gretl_matrix *rtmp, *ctmp;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return 1;
    }  

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) {
	return E_ALLOC;
    }

    ctmp = gretl_matrix_alloc(rows, var->neqns);
    if (ctmp == NULL) {
	gretl_matrix_free(rtmp);
	return E_ALLOC;
    }

    if (var->ci == VECM) {
	vsrc = var->jinfo->list[shock + 1];
    } else {
	vsrc = (var->models[shock])->list[1];
    }

    blockmax = var->neqns / IRF_ROW_MAX;
    if (var->neqns % IRF_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax && !err; block++) {
	int k, vtarg, endrow;
	double r;

	VAR_info_header_block(IRF, vsrc, block, pdinfo, prn);

	for (i=0; i<IRF_ROW_MAX; i++) {
	    k = IRF_ROW_MAX * block + i;
	    if (k >= var->neqns) {
		break;
	    }
	    if (var->ci == VECM) {
		vtarg = var->jinfo->list[k + 1];
	    } else {
		vtarg = (var->models[k])->list[1];
	    }
	    endrow = !(i < IRF_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(IRF, i, vtarg, endrow, pdinfo, prn);
	}

	if (tex || rtf) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}

	for (t=0; t<periods && !err; t++) {
	    VAR_info_print_period(t + 1, prn);
	    if (t == 0) {
		/* calculate initial estimated responses */
		err = gretl_matrix_copy_values(rtmp, var->C);
	    } else {
		/* calculate further estimated responses */
		err = gretl_matrix_multiply(var->A, rtmp, ctmp);
		gretl_matrix_copy_values(rtmp, ctmp);
	    }

	    if (err) break;

	    /* matrix rtmp holds the responses */

	    for (i=0; i<IRF_ROW_MAX; i++) {
		k = IRF_ROW_MAX * block + i;
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(rtmp, k, shock);
		if (tex) {
		    tex_print_double(r, prn);
		    if (i < IRF_ROW_MAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.5g\\cell ", r);
		} else {
		    pprintf(prn, "%#12.5g ", r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (rtmp != NULL) gretl_matrix_free(rtmp);
    if (ctmp != NULL) gretl_matrix_free(ctmp);

    return err;
}

int gretl_VAR_print_all_impulse_responses (GRETL_VAR *var, const DATAINFO *pdinfo, 
					   int horizon, PRN *prn)
{
    int i, pause = 0, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (plain_format(prn)) {
	pause = gretl_get_text_pause();
    } else if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_impulse_response(var, i, horizon, pdinfo, 
					       pause, prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }   

    return err;
}

#define VD_ROW_MAX 5

/**
 * gretl_VAR_print_fcast_decomp:
 * @var: pointer to VAR struct.
 * @targ:
 * @periods: number of periods over which to print decomposition.
 * @pdinfo: dataset information.
 * @pause: if non-zero, pause between sections of output.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      int pause, PRN *prn)
{
    int i, t;
    int vtarg;
    gretl_matrix *vd = NULL;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return 1;
    } 

    vd = gretl_VAR_get_fcast_decomp(var, targ, periods);
    if (vd == NULL) {
	return E_ALLOC;
    }

    if (var->ci == VECM) {
	vtarg = var->jinfo->list[targ + 1];
    } else {
	vtarg = (var->models[targ])->list[1];
    }

    blockmax = (var->neqns + 1) / VDC_ROW_MAX;
    if ((var->neqns + 1) % VDC_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax; block++) {
	int k, vsrc, endrow;
	double r;

	VAR_info_header_block(VDC, vtarg, block, pdinfo, prn);

	for (i=0; i<VDC_ROW_MAX; i++) {
	    k = VDC_ROW_MAX * block + i - 1;
	    if (k < 0) {
		if (tex) {
		    pprintf(prn, " %s & ", I_("std. error"));
		} else if (rtf) {
		    pprintf(prn, " \\qc %s\\cell ", I_("std. error"));
		} else {
		    pprintf(prn, " %12s ", _("std. error"));
		}
		continue;
	    }
	    if (k >= var->neqns) {
		break;
	    }
	    if (var->ci == VECM) {
		vsrc = var->jinfo->list[k + 1];
	    } else {
		vsrc = (var->models[k])->list[1];
	    }
	    endrow = !(i < VDC_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(VDC, i, vsrc, endrow, pdinfo, prn);
	}

	if (tex || rtf) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}

	for (t=0; t<periods && !err; t++) {
	    VAR_info_print_period(t + 1, prn);
	    for (i=0; i<VDC_ROW_MAX; i++) {
		k = VDC_ROW_MAX * block + i - 1;
		if (k < 0) {
		    r = gretl_matrix_get(vd, t, var->neqns);
		    if (tex) {
			pprintf(prn, "%g & ", r);
		    } else if (rtf) {
			pprintf(prn, "\\qc %g\\cell", r);
		    } else {
			pprintf(prn, " %14g ", r);
		    }
		    continue;
		}
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(vd, t, k);
		if (tex) {
		    pprintf(prn, "$%.4f$", r);
		    if (i < VDC_ROW_MAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.4f\\cell", r);
		} else {
		    pprintf(prn, "%10.4f ", r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (vd != NULL) {
	gretl_matrix_free(vd);
    }

    return err;
}

int gretl_VAR_print_all_fcast_decomps (GRETL_VAR *var, const DATAINFO *pdinfo, 
				       int horizon, PRN *prn)
{
    int i, pause = 0, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (plain_format(prn)) {
	pause = gretl_get_text_pause();
    } else if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_fcast_decomp(var, i, horizon, pdinfo, 
					   pause, prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }

    return err;
}

void print_Johansen_test_case (JohansenCode jcode, PRN *prn)
{
    const char *jcase[] = {
	N_("Case 1: No constant"),
	N_("Case 2: Restricted constant"),
	N_("Case 3: Unrestricted constant"),
	N_("Case 4: Restricted trend, unrestricted constant"),
	N_("Case 5: Unrestricted trend and constant")
    };

    if (jcode <= J_UNREST_TREND) {
	if (plain_format(prn)) {
	    pputs(prn, _(jcase[jcode]));
	} else {
	    pputs(prn, I_(jcase[jcode]));
	}
    }
}

static void 
print_VECM_coint_eqns (JohansenInfo *jv, const DATAINFO *pdinfo, PRN *prn)
{
    int rtf = rtf_format(prn);
    char s[16];
    int rows = gretl_matrix_rows(jv->Beta);
    int i, j;
    double x;

    pputs(prn, _("Cointegrating vectors"));
    if (jv->Bse != NULL) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    } 
    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    for (i=0; i<rows; i++) {
	char vname[16];

	if (i < jv->list[0]) {
	    sprintf(vname, "%s(-1)", pdinfo->varname[jv->list[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    strcpy(vname, "const");
	} else if (jv->code == J_REST_TREND) {
	    strcpy(vname, "trend");
	}
	if (rtf) {
	    pputs(prn, vname);
	} else {
	    pprintf(prn, "%-12s", vname);
	}

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Beta, i, j);
	    if (jv->Bse == NULL) {
		x /= gretl_matrix_get(jv->Beta, j, j);
	    }
	    if (rtf) {
		pprintf(prn, "\t%#.5g ", x);
	    } else {
		pprintf(prn, "%#12.5g ", x);
	    }
	}
	gretl_prn_newline(prn);

	if (jv->Bse != NULL) {
	    /* standard errors */
	    if (rtf) {
		pputs(prn, "\t");
	    } else {
		pprintf(prn, "%13s", " ");
	    }
	    for (j=0; j<jv->rank; j++) {
		if (i < jv->rank) {
		    x = 0.0;
		} else {
		    x = gretl_matrix_get(jv->Bse, i - jv->rank, j);
		}
		sprintf(s, "(%#.5g)", x);
		if (rtf) {
		    pprintf(prn, "\t%s", s);
		} else {
		    pprintf(prn, "%12s ", s);
		}
	    }
	    gretl_prn_newline(prn);
	}
    }

    gretl_prn_newline(prn);
}

static void print_VECM_omega (GRETL_VAR *jvar, const DATAINFO *pdinfo, PRN *prn)
{
    int rtf = rtf_format(prn);
    int *list = jvar->jinfo->list;
    char s[32];
    int i, j;

    pprintf(prn, "%s\n", _("Cross-equation covariance matrix"));
    gretl_prn_newline(prn);

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[list[i+1]]);
	if (i == 0) {
	    if (rtf) {
		pprintf(prn, "\t\t%s", s);
	    } else {
		pprintf(prn, "%25s", s);
	    }
	} else {
	    if (rtf) {
		pprintf(prn, "\t%s", s);
	    } else {
		pprintf(prn, "%13s", s);
	    }
	}
    }
    gretl_prn_newline(prn);

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[list[i+1]]);
	if (rtf) {
	    pputs(prn, s);
	    if (strlen(s) < 8) {
		pputc(prn, '\t');
	    }	    
	} else {
	    pprintf(prn, "%-13s", s);
	}
	for (j=0; j<jvar->neqns; j++) {
	    if (rtf) {
		pprintf(prn, "\t%#.5g", gretl_matrix_get(jvar->S, i, j));
	    } else {
		pprintf(prn, "%#12.5g ", gretl_matrix_get(jvar->S, i, j));
	    }
	}
	gretl_prn_newline(prn);
    }

    gretl_prn_newline(prn);

    pprintf(prn, "%s = %g", _("determinant"), exp(jvar->ldet));
    gretl_prn_newline(prn);
}

static void 
print_vecm_header_info (GRETL_VAR *vecm, PRN *prn)
{
    gretl_prn_newline(prn);
    pprintf(prn, "%s = %d", 
	    (plain_format(prn))? _("Cointegration rank") : I_("Cointegration rank"),
	    jrank(vecm));
    gretl_prn_newline(prn);
    print_Johansen_test_case(jcode(vecm), prn); 
}

/**
 * gretl_VAR_print:
 * @var: pointer to VAR struct.
 * @pdinfo: dataset information.
 * @opt: if includes %OPT_I, include impulse responses; if
 * includes %OPT_F, include forecast variance decompositions;
 * if includes %OPT_Q, don't print individual regressions.
 * @prn: pointer to printing struct.
 *
 * Prints the models in @var, along with relevant F-tests and
 * possibly impulse responses and variance decompositions.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print (GRETL_VAR *var, const DATAINFO *pdinfo, gretlopt opt, 
		     PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    char Vstr[72];
    int vecm = var->ci == VECM;
    int dfd = var->models[0]->dfd;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int quiet = (opt & OPT_Q);
    int pause = 0;
    int i, j, k, v;

    if (prn == NULL) {
	return 0;
    }

    ntodate(startdate, var->t1, pdinfo);
    ntodate(enddate, var->t2, pdinfo);

    if (rtf) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    if (vecm) {
	if (tex || rtf) {
	    sprintf(Vstr, I_("VECM system, lag order %d"), var->order + 1);
	} else {
	    sprintf(Vstr, _("VECM system, lag order %d"), var->order + 1);
	}
    } else {
	if (tex || rtf) {
	    sprintf(Vstr, I_("VAR system, lag order %d"), var->order);
	} else {
	    sprintf(Vstr, _("VAR system, lag order %d"), var->order);
	}
    }

    if (tex) {
	pputs(prn, "\\begin{center}");
	pprintf(prn, "\n%s\\\\\n", Vstr);
	pprintf(prn, I_("%s estimates, observations %s--%s ($T=%d$)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, prn);
	}
	pputs(prn, "\n\\end{center}\n");
    } else if (rtf) {
	gretl_print_toggle_doc_flag(prn);
	pprintf(prn, "\n%s\\par\n", Vstr);
	pprintf(prn, I_("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, prn);
	}	
	pputs(prn, "\\par\n\n");
    } else {
	pause = gretl_get_text_pause();
	pprintf(prn, "\n%s\n", Vstr);
	pprintf(prn, _("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? ("Maximum likelihood") : _("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, prn);
	}
	pputs(prn, "\n\n");
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_coint_eqns(var, pdinfo, prn);
	} else {
	    print_VECM_coint_eqns(var->jinfo, pdinfo, prn);
	}
    }

    if (tex) {
	tex_print_VAR_ll_stats(var, prn);
    } else if (rtf) {
	pprintf(prn, "%s = %#g\\par\n", I_("Log-likelihood"), var->ll);
	pprintf(prn, "%s = %#g\\par\n", I_("Determinant of covariance matrix"), 
		exp(var->ldet));
	pprintf(prn, "%s = %.4f\\par\n", I_("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\\par\n", I_("BIC"), var->BIC);
    } else {
	pprintf(prn, "%s = %#g\n", _("Log-likelihood"), var->ll);
	pprintf(prn, "%s = %#g\n", _("Determinant of covariance matrix"), exp(var->ldet));
	pprintf(prn, "%s = %.4f\n", _("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\n", _("BIC"), var->BIC);
    }

    if (vecm) {
	pputc(prn, '\n');
    }

    k = 0;

    for (i=0; i<var->neqns; i++) {
	char Fstr[24];

	if(!quiet) {
	    printmodel(var->models[i], pdinfo, OPT_NONE, prn);
	} else {
	    if (!var->ecm) {
	        v = var->models[i]->list[1];
		if (tex) {
		    pputs(prn, "\n\\begin{center}\n");
		    pprintf(prn, "%s\\\\[1em]\n", I_("Equation for "));
		    pprintf(prn, "%s\\\n", pdinfo->varname[v]);
		    pputs(prn, "\n\\end{center}\n");
		} else if (rtf) {
		    pprintf(prn, "\\par\n%s", I_("Equation for "));
		    pprintf(prn, "%s:\\par\n\n", pdinfo->varname[v]);
		} else {
		    pprintf(prn, "\n%s", I_("Equation for "));
		    pprintf(prn, "%s:\n", pdinfo->varname[v]);
		}
	    }
	}

	if (pause) {
	    scroll_pause();
	}

	if (vecm) {
	    continue;
	}

	if (tex) {
	    pputs(prn, "\n\\begin{center}\n");
	    pprintf(prn, "%s\\\\[1em]\n", I_("F-tests of zero restrictions"));
	    pputs(prn, "\\begin{tabular}{lll}\n");
	} else if (rtf) {
	    pprintf(prn, "%s:\\par\n\n", I_("F-tests of zero restrictions"));
	} else {
	    pprintf(prn, "  %s:\n", _("F-tests of zero restrictions"));
	}

	for (j=0; j<var->neqns; j++) {
	    v = (var->models[j])->list[1];
	    if (tex) {
		pprintf(prn, I_("All lags of %-8s "), pdinfo->varname[v]);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\\\\n", I_("p-value"), 
			fdist(var->Fvals[k], var->order, dfd));
	    } else if (rtf) {
		pprintf(prn, I_("All lags of %-8s "), pdinfo->varname[v]);
		pprintf(prn, "F(%d, %d) = %10g, ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\par\n", I_("p-value"), 
			fdist(var->Fvals[k], var->order, dfd));
	    } else {
		pputs(prn, "    ");
		pprintf(prn, _("All lags of %-8s "), pdinfo->varname[v]);
		sprintf(Fstr, "F(%d, %d)", var->order, dfd);
		pprintf(prn, "%12s = %#10.5g, ", Fstr, var->Fvals[k]);
		pprintf(prn, "%s %.4f\n", _("p-value"), 
			fdist(var->Fvals[k], var->order, dfd));
	    }
	    k++;
	}

	if (var->order > 1) {
	    if (tex) {
		pprintf(prn, I_("All vars, lag %-6d "), var->order);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\\\\n", I_("p-value"), 
			fdist(var->Fvals[k], var->neqns, dfd));
	    } else if (rtf) {
		pprintf(prn, I_("All vars, lag %-6d "), var->order);
		pprintf(prn, "F(%d, %d) = %10g, ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\par\n", I_("p-value"), 
			fdist(var->Fvals[k], var->neqns, dfd));
	    } else {
		pputs(prn, "    ");
		pprintf(prn, _("All vars, lag %-6d "), var->order);
		sprintf(Fstr, "F(%d, %d)", var->neqns, dfd);
		pprintf(prn, "%12s = %#10.5g, ", Fstr, var->Fvals[k]);
		pprintf(prn, "%s %.4f\n", _("p-value"), 
			fdist(var->Fvals[k], var->neqns, dfd));
	    } 
	    k++;
	}

	if (tex) {
	    pputs(prn, "\\end{tabular}\n"
		  "\\end{center}\n\n"
		  "\\clearpage\n\n");
	} else if (rtf) {
	    pputs(prn, "\\par\\n\n");
	} else if (pause) {
	    scroll_pause();
	}
    }

    pputc(prn, '\n');

    /* global LR test on max lag */
    if (!na(var->LR)) {
	char h0str[64];
	char h1str[64];
	int df = var->neqns * var->neqns;

	if (tex || rtf) {
	    sprintf(h0str, I_("the longest lag is %d"), var->order - 1);
	    sprintf(h1str, I_("the longest lag is %d"), var->order);
	} else {
	    sprintf(h0str, _("the longest lag is %d"), var->order - 1);
	    sprintf(h1str, _("the longest lag is %d"), var->order);
	}	    

	if (tex) {
	    pprintf(prn, "\\noindent %s ---\\par\n", I_("For the system as a whole"));
	    pprintf(prn, "%s: %s\\par\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "%s: %s\\par\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f (%s %f)\\par\n",
		    I_("Likelihood ratio test"), 
		    df, var->LR, I_("p-value"), chisq(var->LR, df));
	} else if (rtf) {
	    pprintf(prn, "\\par %s\n", I_("For the system as a whole"));
	    pprintf(prn, "\\par %s: %s\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "\\par %s: %s\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "\\par %s: %s(%d) = %g (%s %f)\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    df, var->LR, I_("p-value"), chisq(var->LR, df));
	} else {
	    pprintf(prn, "%s\n", _("For the system as a whole"));
	    pprintf(prn, "  %s: %s\n", _("Null hypothesis"), h0str);
	    pprintf(prn, "  %s: %s\n", _("Alternative hypothesis"), h1str);
	    pprintf(prn, "  %s: %s(%d) = %g (%s %f)\n",
		    _("Likelihood ratio test"), _("Chi-square"), 
		    df, var->LR, _("p-value"), chisq(var->LR, df));
	}
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_omega(var, pdinfo, prn);
	} else {
	    print_VECM_omega(var, pdinfo, prn);
	    pputc(prn, '\n');
	}
    } else {
	pputc(prn, '\n');
    }

    if (opt & OPT_I) {
	gretl_VAR_print_all_impulse_responses(var, pdinfo, 0, prn);
    }

    if (opt & OPT_F) {
	gretl_VAR_print_all_fcast_decomps(var, pdinfo, 0, prn);
    }

    if (rtf) {
	pputs(prn, "}\n");
    }

    return 0;
}

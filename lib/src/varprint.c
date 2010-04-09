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

#include "libgretl.h" 
#include "var.h"  
#include "vartest.h" 
#include "johansen.h"
#include "varprint.h"
#include "libset.h"
#include "texprint.h"

#include <glib.h>

/**
 * gretl_VAR_print_sigma:
 * @var: pointer to gretl VAR structure.
 * @prn: printing object.
 *
 * Prints to @prn the contemporaneous (cross-equation) variance
 * matrix for @var.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print_sigma (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->S == NULL) {
	err = 1;
    } else {
	print_contemp_covariance_matrix(var->S, var->ldet, prn);
    }

    return err;
}

/* printing of impulse responses and variance decompositions */

#define IRF_ROW_MAX 4
#define IRF_WIDTH  13

#define VDC_ROW_MAX 5
#define VDC_WIDTH  11

#define VECM_WIDTH 13

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
    char vname[32];

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

static void VAR_info_print_vname (int i, int v, int endrow, int width,
				  const DATAINFO *pdinfo, PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[32];

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
	pprintf(prn, "%*s", width, pdinfo->varname[v]);
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

static int varprint_namelen (const GRETL_VAR *var, const DATAINFO *pdinfo,
			     int rmax, int block)
{
    int len, maxlen = 0;
    int i, k, v;

    for (i=0; i<rmax; i++) {
	k = rmax * block + i - 1;
	if (k < 0) {
	    continue;
	}
	if (k >= var->neqns) {
	    break;
	}
	v = var->ylist[k+1];
	len = strlen(pdinfo->varname[v]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    return maxlen;
}

/**
 * gretl_VAR_print_impulse_response:
 * @var: pointer to VAR struct.
 * @shock: index number of the "shock" variable.
 * @periods: number of periods over which to print response.
 * @pdinfo: dataset information.
 * @prn: gretl printing object.
 *
 * Prints to @prn the estimated responses of the endogenous
 * variables in @var to a one-standard deviation shock in
 * the specified variable: @shock is a zero-based index into
 * the equations of the VAR so for example if @shock = 1,
 * the responses are to a shock in the second endogenous
 * variable in the VAR specification.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  PRN *prn)
{
    gretl_matrix *rtmp = NULL, *ctmp = NULL;
    gretl_matrix *C = var->C;
    int rows = var->neqns * effective_order(var);
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int i, t, vsrc, err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return 1;
    }  

    if (var->ord != NULL) {
	C = reorder_responses(var, &err);
	if (err) {
	    return err;
	}
    } 

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    ctmp = gretl_matrix_alloc(rows, var->neqns);

    if (rtmp == NULL || ctmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    vsrc = var->ylist[shock + 1];

    blockmax = var->neqns / IRF_ROW_MAX;
    if (var->neqns % IRF_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax && !err; block++) {
	int k, vtarg, endrow;
	int namelen, width;
	double r;

	VAR_info_header_block(IRF, vsrc, block, pdinfo, prn);

	namelen = varprint_namelen(var, pdinfo, IRF_ROW_MAX, block);
	width = (namelen < IRF_WIDTH - 1)? IRF_WIDTH : namelen + 1;

	for (i=0; i<IRF_ROW_MAX; i++) {
	    k = IRF_ROW_MAX * block + i;
	    if (k >= var->neqns) {
		break;
	    }
	    vtarg = var->ylist[k+1];
	    endrow = !(i < IRF_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(i, vtarg, endrow, width, pdinfo, prn);
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
		err = gretl_matrix_copy_values(rtmp, C);
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
		    if (i == 0) pputc(prn, ' ');
		    pprintf(prn, "%#*.5g ", width - 1, r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);
    }

 bailout:

    if (rtmp != NULL) gretl_matrix_free(rtmp);
    if (ctmp != NULL) gretl_matrix_free(ctmp);

    if (C != var->C) {
	gretl_matrix_free(C);
    }

    return err;
}

int gretl_VAR_print_all_impulse_responses (GRETL_VAR *var, const DATAINFO *pdinfo, 
					   int horizon, PRN *prn)
{
    int i, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_impulse_response(var, i, horizon, pdinfo, 
					       prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }   

    return err;
}

/**
 * gretl_VAR_print_fcast_decomp:
 * @var: pointer to VAR struct.
 * @targ:
 * @periods: number of periods over which to print decomposition.
 * @pdinfo: dataset information.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      PRN *prn)
{
    int vtarg;
    gretl_matrix *vd = NULL;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int i, t, err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return 1;
    } 

    vd = gretl_VAR_get_fcast_decomp(var, targ, periods, &err);
    if (err) {
	pprintf(prn, "Forecast decomposition failed\n");
	return err;
    }

    vtarg = var->ylist[targ + 1];

    blockmax = (var->neqns + 1) / VDC_ROW_MAX;
    if ((var->neqns + 1) % VDC_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax; block++) {
	int k, vsrc, endrow;
	int namelen, width;
	double r;

	VAR_info_header_block(VDC, vtarg, block, pdinfo, prn);

	namelen = varprint_namelen(var, pdinfo, VDC_ROW_MAX, block);
	width = (namelen < VDC_WIDTH - 1)? VDC_WIDTH : namelen + 1;

	for (i=0; i<VDC_ROW_MAX; i++) {
	    k = VDC_ROW_MAX * block + i - 1;
	    if (k < 0) {
		if (tex) {
		    pprintf(prn, " %s & ", I_("std. error"));
		} else if (rtf) {
		    pprintf(prn, " \\qc %s\\cell ", I_("std. error"));
		} else {
		    pprintf(prn, " %14s", _("std. error"));
		}
		continue;
	    }
	    if (k >= var->neqns) {
		break;
	    }
	    vsrc = var->ylist[k+1];
	    endrow = !(i < VDC_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(i, vsrc, endrow, width, pdinfo, prn);
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
		    /* standard error column */
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
		    pprintf(prn, "%*.4f ", width - 1, r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);
    }

    if (vd != NULL) {
	gretl_matrix_free(vd);
    }

    return err;
}

int gretl_VAR_print_all_fcast_decomps (GRETL_VAR *var, const DATAINFO *pdinfo, 
				       int horizon, PRN *prn)
{
    int i, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_fcast_decomp(var, i, horizon, pdinfo, 
					   prn);
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

static char *make_beta_vname (char *vname,
			      const GRETL_VAR *v,
			      const DATAINFO *pdinfo,
			      int i)
{
    if (i < v->neqns) {
	strcpy(vname, pdinfo->varname[v->ylist[i+1]]);
    } else if (auto_restr(v) && i == v->neqns) {
	strcpy(vname, (jcode(v) == J_REST_CONST)? "const" : "trend");
    } else if (v->rlist != NULL) {
	int k = i - v->ylist[0] - auto_restr(v) + 1;

	strcpy(vname, pdinfo->varname[v->rlist[k]]);
    } 

    return vname;
}

static int max_beta_namelen (GRETL_VAR *v, 
			     const DATAINFO *pdinfo)
{
    int r = gretl_matrix_rows(v->jinfo->Beta);
    char s[32];
    int i, ni, n = 0;

    for (i=0; i<r; i++) {    
	make_beta_vname(s, v, pdinfo, i);
	ni = strlen(s);
	if (ni > n) {
	    n = ni;
	}
    }

    return n;
}

static void 
print_VECM_coint_eqns (GRETL_VAR *jvar, 
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    int utf = plain_format(prn);
    int rtf = rtf_format(prn);
    char namefmt[8];
    char s[16], vname[32];
    int rows = gretl_matrix_rows(jv->Beta);
    int nwid;
    int i, j;
    double x;

    pprintf(prn, "beta (%s", (utf)? _("cointegrating vectors") :
	    I_("cointegrating vectors"));
    if (jv->Bse != NULL) {
	pprintf(prn, ", %s)", (utf)? _("standard errors in parentheses") :
		I_("standard errors in parentheses"));
    } else {
	pputc(prn, ')');
    }

    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    nwid = max_beta_namelen(jvar, pdinfo) + 1;
    sprintf(namefmt, "%%-%ds", nwid);

    for (i=0; i<rows; i++) {
	make_beta_vname(vname, jvar, pdinfo, i);
	if (rtf) {
	    pputs(prn, vname);
	} else {
	    pprintf(prn, namefmt, vname);
	}

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Beta, i, j);
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
		bufspace(nwid + 1, prn);
	    }
	    for (j=0; j<jv->rank; j++) {
		x = gretl_matrix_get(jv->Bse, i, j);
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

    rows = gretl_matrix_rows(jv->Alpha);

    pprintf(prn, "alpha (%s", (utf)? _("adjustment vectors") :
	    I_("adjustment vectors"));
    if (jv->Ase != NULL) {
	pprintf(prn, ", %s)", (utf)? _("standard errors in parentheses") :
		I_("standard errors in parentheses"));
    } else {
	pputc(prn, ')');
    }

    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    for (i=0; i<rows; i++) {
	sprintf(vname, "%s", pdinfo->varname[jvar->ylist[i+1]]);
	if (rtf) {
	    pputs(prn, vname);
	} else {
	    pprintf(prn, namefmt, vname);
	}

	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Alpha, i, j);
	    if (rtf) {
		pprintf(prn, "\t%#.5g ", x);
	    } else {
		pprintf(prn, "%#12.5g ", x);
	    }
	}
	gretl_prn_newline(prn);

	if (jv->Ase != NULL) {
	    if (rtf) {
		pputs(prn, "\t");
	    } else {
		bufspace(nwid + 1, prn);
	    }
	    for (j=0; j<jv->rank; j++) {
		x = gretl_matrix_get(jv->Ase, i, j);
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

static int max_vlen (const int *list, const DATAINFO *pdinfo)
{
    int i, ni, n = 0;

    for (i=1; i<=list[0]; i++) {
	ni = strlen(pdinfo->varname[list[i]]);
	if (ni > n) {
	    n = ni;
	}
    }

    return n;
}

static void print_VECM_omega (GRETL_VAR *jvar, const DATAINFO *pdinfo, PRN *prn)
{
    int utf = plain_format(prn);
    int rtf = rtf_format(prn);
    int *list = jvar->ylist;
    char s[32];
    int w0, wi = 12;
    int i, j;

    pprintf(prn, "%s:\n", (utf)? _("Cross-equation covariance matrix") :
	    I_("Cross-equation covariance matrix"));
    gretl_prn_newline(prn);

    w0 = max_vlen(list, pdinfo) + 1;

    /* top row: names of Y variables */

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "%s", pdinfo->varname[list[i+1]]);
	if (i == 0) {
	    if (rtf) {
		pprintf(prn, "\t\t%s", s);
	    } else {
		wi = strlen(s);
		if (wi < VECM_WIDTH) wi = VECM_WIDTH;
		pprintf(prn, "%*s", w0 + wi, s);
	    }
	} else {
	    if (rtf) {
		pprintf(prn, "\t%s", s);
	    } else {
		wi = strlen(s) + 1;
		if (wi < VECM_WIDTH) wi = VECM_WIDTH;
		pprintf(prn, "%*s", wi, s);
	    }
	}
    }

    gretl_prn_newline(prn);

    /* subsequent rows: Y name plus values */

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "%s", pdinfo->varname[list[i+1]]);
	if (rtf) {
	    pputs(prn, s);
	    if (strlen(s) < 8) {
		pputc(prn, '\t');
	    }	    
	} else {
	    pprintf(prn, "%-*s ", w0, s);
	}
	for (j=0; j<jvar->neqns; j++) {
	    if (rtf) {
		pprintf(prn, "\t%#.5g", gretl_matrix_get(jvar->S, i, j));
	    } else {
		wi = strlen(pdinfo->varname[list[j+1]]);
		if (wi < VECM_WIDTH - 1) wi = VECM_WIDTH - 1;
		pprintf(prn, "%#*.5g ", wi, gretl_matrix_get(jvar->S, i, j));
	    }
	}
	gretl_prn_newline(prn);
    }

    gretl_prn_newline(prn);

    pprintf(prn, "%s = %g", (utf)? _("determinant") : I_("determinant"),
	    exp(jvar->ldet));

    gretl_prn_newline(prn);
}

/* FIXME TeX and RTF */

static void vecm_print_restrictions (GRETL_VAR *vecm, PRN *prn)
{
    if (vecm->jinfo->R != NULL) {
	pputs(prn, "\n\n");
	pputs(prn, _("Restrictions on beta:"));
	pputc(prn, '\n');
	print_restriction_from_matrices(vecm->jinfo->R, vecm->jinfo->q, 
					'b', gretl_VECM_n_beta(vecm), 
					prn);
	pputc(prn, '\n');
    }

    if (vecm->jinfo->Ra != NULL) {
	if (vecm->jinfo->R == NULL) {
	    pputs(prn, "\n\n");
	}
	pputs(prn, _("Restrictions on alpha:"));
	pputc(prn, '\n');
	print_restriction_from_matrices(vecm->jinfo->Ra, vecm->jinfo->qa, 
					'a', gretl_VECM_n_alpha(vecm), 
					prn);
	pputc(prn, '\n');
    }	  
}

static void print_LR_stat (double x, int df, PRN *prn)
{
    double pv = chisq_cdf_comp(df, x);

    if (na(pv)) {
	return;
    }

    if (tex_format(prn)) {
	pprintf(prn, "$2 (l_u - l_r) = %g$", x);
	gretl_prn_newline(prn);
	pprintf(prn, "$P(\\chi^2_{%d} > %g) = %g$", df, x, pv);
    } else if (rtf_format(prn)) {
	pprintf(prn, "2 * (lu - lr) = %g", x);
	gretl_prn_newline(prn);
	pprintf(prn, I_("P(Chi-Square(%d) > %g) = %g"), df, x, pv);
    } else {
	pprintf(prn, "2 * (lu - lr) = %g", x);
	gretl_prn_newline(prn);
	pprintf(prn, _("P(Chi-Square(%d) > %g) = %g"), df, x, pv);
    }
}

enum {
    LR_TOTAL,
    LR_RELATIVE
};

static void 
vecm_print_LR_test (GRETL_VAR *vecm, PRN *prn, int code)
{
    double ll0, x;
    int df;

    if (code == LR_RELATIVE) {
	ll0 = vecm->jinfo->prior_ll;
	df = vecm->jinfo->lrdf - vecm->jinfo->prior_df;
	gretl_prn_newline(prn);
	if (tex_format(prn) || rtf_format(prn)) {
	    pputs(prn, I_("Relative to prior restriction"));
	} else {
	    pputs(prn, _("Relative to prior restriction"));
	} 
	pputc(prn, ':');
	gretl_prn_newline(prn);
    } else {
	ll0 = vecm->jinfo->ll0;
	df = vecm->jinfo->lrdf;
    }

    x = 2.0 * (ll0 - vecm->ll);

    if (tex_format(prn)) {
	pprintf(prn, I_("Unrestricted loglikelihood $(l_u) = %.8g$"), ll0);
	gretl_prn_newline(prn);
	pprintf(prn, I_("Restricted loglikelihood $(l_r) = %.8g$"), vecm->ll);
    } else if (rtf_format(prn)) {
	pprintf(prn, I_("Unrestricted loglikelihood (lu) = %.8g"), ll0);
	gretl_prn_newline(prn);
	pprintf(prn, I_("Restricted loglikelihood (lr) = %.8g"), vecm->ll);
    } else {
	pprintf(prn, _("Unrestricted loglikelihood (lu) = %.8g"), ll0);
	gretl_prn_newline(prn);
	pprintf(prn, _("Restricted loglikelihood (lr) = %.8g"), vecm->ll);
    }

    gretl_prn_newline(prn);
    print_LR_stat(x, df, prn);
    gretl_prn_newline(prn);
}

static void 
print_vecm_header_info (GRETL_VAR *vecm, int *lldone, PRN *prn)
{
    JohansenInfo *J = vecm->jinfo;
    gretl_prn_newline(prn);

    if (vecm->jinfo == NULL) {
	return;
    }
    
    pprintf(prn, "%s = %d", 
	    (plain_format(prn))? _("Cointegration rank") : I_("Cointegration rank"),
	    jrank(vecm));
    gretl_prn_newline(prn);
    print_Johansen_test_case(jcode(vecm), prn); 

    if (J->R != NULL || J->Ra != NULL) {
	vecm_print_restrictions(vecm, prn);
	if (!na(J->ll0) && J->lrdf > 0) {
	    vecm_print_LR_test(vecm, prn, LR_TOTAL);
	    *lldone = 1;
	} 
	if (!na(J->prior_ll) && J->prior_df > 0) {
	    vecm_print_LR_test(vecm, prn, LR_RELATIVE);
	    *lldone = 1;
	} 	
	if (!*lldone) {
	    pputc(prn, '\n');
	}
    } else {
	pputc(prn, '\n');
    }
}

static void VAR_print_LB_stat (const GRETL_VAR *var, PRN *prn)
{
    int k = var->order + (var->ci == VECM);
    int df = var->neqns * var->neqns * (var->LBs - k);
    double pv = chisq_cdf_comp(df, var->LB);

    if (tex_format(prn)) {
	pprintf(prn, "\\noindent\n%s: LB(%d) = %g, %s = %d [%.4f]\\par\n",
		I_("Portmanteau test"), var->LBs, var->LB,
		I_("df"), df, pv);	    
    } else if (rtf_format(prn)) {
	pprintf(prn, "%s: LB(%d) = %g, %s = %d [%.4f]\\par\n",
		I_("Portmanteau test"), var->LBs, var->LB,
		I_("df"), df, pv);
    } else {
	pprintf(prn, "%s: LB(%d) = %g, %s = %d [%.4f]\n",
		_("Portmanteau test"), var->LBs, var->LB,
		_("df"), df, pv);
    }
}

static int Ivals_ok (const GRETL_VAR *var)
{
    return var->Ivals != NULL && !na(var->Ivals[0])
	&& !na(var->Ivals[1]) && !na(var->Ivals[2]);
}

static int max_Ftest_label_len (GRETL_VAR *var, const DATAINFO *pdinfo,
				char *s, int maxlag)
{
    int len, len1, len2;
    int maxnamelen = 0;
    int i, v;

    for (i=0; i<var->neqns; i++) {
	v = (var->models[i])->list[1];
	len = strlen(pdinfo->varname[v]);
	if (len > maxnamelen) {
	    maxnamelen = len;
	}
    }

    sprintf(s, _("All lags of %s"), "x");
    len1 = g_utf8_strlen(s, -1) + maxnamelen - 1;

    sprintf(s, _("All vars, lag %d"), maxlag);
    len2 = g_utf8_strlen(s, -1);

    len = (len1 > len2)? len1 : len2;

    return (len > 25)? (len + 1) : 25;
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
    char label[72];
    int vecm = (var->ci == VECM);
    int dfd = var->models[0]->dfd;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int quiet = (opt & OPT_Q);
    int nlags, maxlag, nextlag;
    int llen, fwidth = 0;
    int lldone = 0;
    double pv;
    int i, j, k, v;

    if (prn == NULL) {
	return 0;
    }

    nlags = var_n_lags(var);
    maxlag = var_max_lag(var);
    
    if (var->lags != NULL) {
	nextlag = var->lags[nlags - 1];
    } else {
	nextlag = maxlag - 1;
    }

    ntodate(startdate, var->t1, pdinfo);
    ntodate(enddate, var->t2, pdinfo);

    if (rtf) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    if (vecm) {
	if (tex || rtf) {
	    sprintf(label, I_("VECM system, lag order %d"), var->order + 1);
	} else {
	    sprintf(label, _("VECM system, lag order %d"), var->order + 1);
	}
    } else {
	if (tex || rtf) {
	    sprintf(label, I_("VAR system, lag order %d"), var->order);
	} else {
	    sprintf(label, _("VAR system, lag order %d"), var->order);
	}
    }

    if (tex) {
	pputs(prn, "\\begin{center}");
	pprintf(prn, "\n%s\\\\\n", label);
	pprintf(prn, I_("%s estimates, observations %s--%s ($T=%d$)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, &lldone, prn);
	}
	pputs(prn, "\n\\end{center}\n");
    } else if (rtf) {
	gretl_print_toggle_doc_flag(prn);
	pprintf(prn, "\n%s\\par\n", label);
	pprintf(prn, I_("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, &lldone, prn);
	}	
	pputs(prn, "\\par\n\n");
    } else {
	pprintf(prn, "\n%s\n", label);
	pprintf(prn, _("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? _("Maximum likelihood") : _("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, &lldone, prn);
	}
	pputc(prn, '\n');
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_coint_eqns(var, pdinfo, prn);
	} else {
	    print_VECM_coint_eqns(var, pdinfo, prn);
	}
    }

    if (tex) {
	tex_print_VAR_ll_stats(var, prn);
    } else if (rtf) {
	if (!lldone) {
	    pprintf(prn, "%s = %.8g\\par\n", I_("Log-likelihood"), var->ll);
	}
	pprintf(prn, "%s = %.8g\\par\n", I_("Determinant of covariance matrix"), 
		exp(var->ldet));
	pprintf(prn, "%s = %.4f\\par\n", I_("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\\par\n", I_("BIC"), var->BIC);
	pprintf(prn, "%s = %.4f\\par\n", I_("HQC"), var->HQC);
    } else {
	if (!lldone) {
	    pprintf(prn, "%s = %.8g\n", _("Log-likelihood"), var->ll);
	}
	pprintf(prn, "%s = %.8g\n", _("Determinant of covariance matrix"), exp(var->ldet));
	pprintf(prn, "%s = %.4f\n", _("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\n", _("BIC"), var->BIC);
	pprintf(prn, "%s = %.4f\n", _("HQC"), var->HQC);
    }

    if (var->LBs > 0 && !na(var->LB)) {
	VAR_print_LB_stat(var, prn);
    }

    if (vecm && !quiet) {
	pputc(prn, '\n');
    }

    k = 0;

    for (i=0; i<var->neqns; i++) {
	char Fstr[24];

	if (!quiet) {
	    printmodel(var->models[i], pdinfo, OPT_NONE, prn);
	} else {
	    if (var->ci != VECM) {
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
	    pprintf(prn, "%s:\n\n", _("F-tests of zero restrictions"));
	}

	if (!tex) {
	    fwidth = max_Ftest_label_len(var, pdinfo, label, maxlag);
	}

	for (j=0; j<var->neqns; j++) {
	    pv = snedecor_cdf_comp(nlags, dfd, var->Fvals[k]);
	    v = (var->models[j])->list[1];
	    if (tex) {
		pprintf(prn, I_("All lags of %s"), pdinfo->varname[v]);
		pputs(prn, " & ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", nlags, dfd, var->Fvals[k]);
		pprintf(prn, "[%.4f]\\\\\n", pv);
	    } else if (rtf) {
		sprintf(label, I_("All lags of %s"), pdinfo->varname[v]);
		llen = strlen(label);
		pputs(prn, label);
		bufspace(fwidth - llen, prn);
		pprintf(prn, "F(%d, %d) = %8.5g ", nlags, dfd, var->Fvals[k]);
		pprintf(prn, "[%.4f]\\par\n", pv);
	    } else {
		sprintf(label, _("All lags of %s"), pdinfo->varname[v]);
		llen = g_utf8_strlen(label, -1);
		pputs(prn, label);
		bufspace(fwidth - llen, prn);
		sprintf(Fstr, "F(%d, %d)", nlags, dfd);
		pprintf(prn, "%12s = %#8.5g [%.4f]\n", Fstr, var->Fvals[k], pv);
	    }
	    k++;
	}

	if (var->order > 1) {
	    pv = snedecor_cdf_comp(var->neqns, dfd, var->Fvals[k]);
	    if (tex) {
		pprintf(prn, I_("All vars, lag %d"), maxlag);
		pputs(prn, " & ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "[%.4f]\\\\\n", pv);
	    } else if (rtf) {
		sprintf(label, I_("All vars, lag %d"), maxlag);
		llen = strlen(label);
		pputs(prn, label);
		bufspace(fwidth - llen, prn);
		pprintf(prn, "F(%d, %d) = %8.5g ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "[%.4f]\\par\n", pv);
	    } else {
		sprintf(label, _("All vars, lag %d"), maxlag);
		llen = g_utf8_strlen(label, -1);
		pputs(prn, label);
		bufspace(fwidth - llen, prn);
		sprintf(Fstr, "F(%d, %d)", var->neqns, dfd);
		pprintf(prn, "%12s = %#8.5g [%.4f]\n", Fstr, var->Fvals[k], pv);
	    } 
	    k++;
	}

	if (tex) {
	    pputs(prn, "\\end{tabular}\n"
		  "\\end{center}\n\n"
		  "\\clearpage\n\n");
	} else if (rtf) {
	    pputs(prn, "\\par\\n\n");
	} 
    }

    /* global LR test on max lag */
    if (!na(var->LR)) {
	char h0str[64];
	char h1str[64];
	int df = var->neqns * var->neqns;

	pputc(prn, '\n');

	if (rtf || tex) {
	    sprintf(h0str, I_("the longest lag is %d"), nextlag);
	    sprintf(h1str, I_("the longest lag is %d"), maxlag);
	} else {
	    sprintf(h0str, _("the longest lag is %d"), nextlag);
	    sprintf(h1str, _("the longest lag is %d"), maxlag);
	}

	pv = chisq_cdf_comp(df, var->LR);

	if (tex) {
	    pprintf(prn, "\\noindent %s ---\\par\n", I_("For the system as a whole"));
	    pprintf(prn, "%s: %s\\par\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "%s: %s\\par\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f [%.4f]\\par\n",
		    I_("Likelihood ratio test"), df, var->LR, pv);
	} else if (rtf) {
	    pprintf(prn, "\\par %s\n", I_("For the system as a whole"));
	    pprintf(prn, "\\par %s: %s\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "\\par %s: %s\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "\\par %s: %s(%d) = %g [%.4f]\n", I_("Likelihood ratio test"), 
		    I_("Chi-square"), df, var->LR, pv);
	} else {
	    int ordlen = (var->order > 10)? 2 : 1;

	    pprintf(prn, "%s:\n\n", _("For the system as a whole"));
	    pprintf(prn, "  %s: %s\n", _("Null hypothesis"), h0str);
	    pprintf(prn, "  %s: %s\n", _("Alternative hypothesis"), h1str);
	    pprintf(prn, "  %s: %s(%d) = %g [%.4f]\n", _("Likelihood ratio test"), 
		    _("Chi-square"), df, var->LR, pv);
	    if (Ivals_ok(var)) {
		/* Info criteria comparison */
		pprintf(prn, "\n  %s:\n", _("Comparison of information criteria"));
		pputs(prn, "  ");
		pprintf(prn, _("Lag order %*d"), ordlen, var->order);
		pprintf(prn, ": AIC = %#.6g, BIC = %#.6g, HQC = %#.6g\n", 
			var->AIC, var->BIC, var->HQC);
		pputs(prn, "  ");
		pprintf(prn, _("Lag order %*d"), ordlen, var->order - 1);
		pprintf(prn, ": AIC = %#.6g, BIC = %#.6g, HQC = %#.6g\n", 
			var->Ivals[0], var->Ivals[1], var->Ivals[2]);
	    }
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

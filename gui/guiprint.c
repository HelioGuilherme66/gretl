/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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

/*  guiprint.c - RTF and LaTeX generation for gretl, plus native
    printing */ 

#include "gretl.h"

static int r_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
                         PRN *prn);

#if defined(USE_GNOME)

#include <libgnomeprint/gnome-print.h>
#include <libgnomeprint/gnome-printer-dialog.h>

static void time_string (char *s)
{
    time_t prntime = time(NULL);
    
    sprintf(s, _("gretl output %s"), ctime(&prntime));
    s[strlen(s)-1] = '\0';
}

void winprint (char *fullbuf, char *selbuf)
{
    GnomePrinter *printer;
    GnomePrintContext *pc;    
    GnomeFont *font;
    char *p, linebuf[90], hdrstart[48], hdr[70];
    int page_lines = 47;
    int x, y, line, page;
    size_t len;

    printer = gnome_printer_dialog_new_modal();

    if (!printer) {
	free(fullbuf);
	free(selbuf);
	return;
    }

    pc = gnome_print_context_new_with_paper_size(printer, "US-Letter");

    gnome_print_beginpage (pc, _("gretl output"));

    /* could use GNOME_FONT_MEDIUM below */
    /* font = gnome_font_new_closest("Courier", GNOME_FONT_BOOK, FALSE, 10); */
    font = gnome_font_new("Courier", 10.0);
    gnome_print_setfont(pc, font);
    gnome_print_setrgbcolor(pc, 0, 0, 0);

    time_string(hdrstart);
    if (selbuf != NULL) p = selbuf;
    else p = fullbuf;
    page = 1;
    x = 72;
    time_string(hdrstart);
    while (*p) { /* pages loop */
	line = 0;
	y = 756;
	if (page > 1) 
	    gnome_print_beginpage (pc, _("gretl output"));
	sprintf(hdr, _("%s, page %d"), hdrstart, page++);
	gnome_print_moveto(pc, x, y);
	gnome_print_show(pc, hdr);
	y = 720;
	while (*p && line < page_lines) { /* lines loop */
	    len = strcspn(p, "\n");
	    *linebuf = '\0';
	    strncat(linebuf, p, len);
	    gnome_print_moveto(pc, x, y);
	    gnome_print_show(pc, linebuf);
	    p += len + 1;
	    y -= 14; /* line spacing */
	    line++;
	}
	gnome_print_showpage(pc);
    }

    /* clean up */
    gnome_print_context_close(pc);
    gtk_object_unref(GTK_OBJECT(font));
    gtk_object_unref(GTK_OBJECT(printer));
    free(fullbuf);
    if (selbuf) 
	free(selbuf);
}

void gnome_print_graph (const char *fname)
{
    GnomePrinter *printer;
    GnomePrintContext *pc; 
    GdkPixbuf *pbuf;
    char tmp[MAXLEN];
    int image_left_x = 530, image_bottom_y = 50;
    int width, height;

    printer = gnome_printer_dialog_new_modal();

    if (!printer) return;

    /* run gnuplot on the plotfile to generate pngtmp */
    sprintf(tmp, "\"%s\" \"%s\"", paths.gnuplot, fname);
    if (system(tmp)) {
	errbox("Failed to generate graph");
	gtk_object_unref(GTK_OBJECT(printer));
	return;
    }

    build_path(paths.userdir, "gretltmp.png", tmp, NULL);
    pbuf = gdk_pixbuf_new_from_file(tmp);
    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);
    remove(tmp);

    pc = gnome_print_context_new_with_paper_size(printer, "US-Letter");

    gnome_print_beginpage(pc, _("gretl output"));

    gnome_print_gsave(pc);
    gnome_print_translate(pc, image_left_x, image_bottom_y);
    gnome_print_rotate(pc, 90);
    gnome_print_scale(pc, width, height);
    gnome_print_pixbuf(pc, pbuf);
    gnome_print_grestore(pc);
    gnome_print_showpage(pc);

    /* clean up */
    gnome_print_context_close(pc);
    gtk_object_unref(GTK_OBJECT(printer));
}

#endif /* USE_GNOME */

void model_to_rtf (MODEL *pmod)
{
    PRN *prn;

    if (bufopen(&prn)) return;
    
    r_printmodel(pmod, datainfo, prn);

    prn_to_clipboard(prn);

    gretl_print_destroy(prn);
}

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

/* ......................................................... */ 

static int r_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    prn->format = GRETL_PRINT_FORMAT_RTF;
    
    return printmodel (pmod, pdinfo, prn);
}

/* ......................................................... */ 

static void printfrtf (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\qr %s\\cell\\intbl \\row\n",
		    I_("undefined"));
	} else {
	    pprintf(prn, "\\qr %s\\cell", I_("undefined"));
	}
	return;
    }

    if (endrow) {
	pprintf(prn, "\\qr %#*g\\cell\\intbl \\row\n", GRETL_DIGITS, zz);
    } else {
	pprintf(prn, "\\qr %#*g\\cell", GRETL_DIGITS, zz);
    }
}

#define SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx1600\\cellx3200\\cellx4800\\cellx6400" \
                   "\\cellx8000\n"

#define VAR_SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2000\\cellx4000\\cellx6000\\cellx8000\n"

/* ............................................................. */

void rtfprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
{
    char date1[9], date2[9], tmp[128];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s - %s"),
	    date1, date2);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
    
    if (lo == 1) {
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "%s\\par\n\n", tmp);
	pprintf(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	strcpy(tmp, I_("(missing values denoted by -999 will be skipped)"));
	pprintf(prn, "%s\\par\n\n", tmp);
	pprintf(prn, "{" SUMM_ROW
		"\\intbl \\qc %s\\cell", I_("Variable"));
    }

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1)
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	printfrtf(xbar, prn, 0);
	printfrtf(summ->xmedian[v], prn, 0);
	printfrtf(summ->xpx[v], prn, 0);
	printfrtf(summ->xpy[v], prn, 1);
    }

    if (lo > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
			I_("Variable"));
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Std. Dev."), I_("C.V."), I_("Skewness"), I_("Ex. kurtosis"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	if (lo > 1)
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = -999;
	printfrtf(std, prn, 0);
	printfrtf(xcv, prn, 0);
	printfrtf(summ->xskew[v], prn, 0);
	printfrtf(summ->xkurt[v], prn, 1);
    }

    pprintf(prn, "}}\n");
}

/* ............................................................. */

static void printftex (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow)
	    pprintf(prn, "\\multicolumn{1}{c}{%s}\\\\", I_("undefined"));
	else
	    pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
    } else {
	char s[32];

	tex_dcolumn_double(zz, s);

	if (endrow) 
	    pprintf(prn, "%s\\\\", s);
	else
	    pprintf(prn, "%s & ", s);
    }	
}

/* ............................................................. */

void texprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
{
    char date1[9], date2[9], vname[16], tmp[128];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;
    char pt = get_local_decpoint();

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s--%s"),
	    date1, date2);

    pprintf(prn, "\\begin{center}\n%s\\\\\n", tmp);
    
    if (lo == 1) {
	tex_escape(vname, pdinfo->varname[summ->list[1]]);
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		vname, summ->n);
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pprintf(prn, "\\begin{tabular}{rrrr}\n");
    } else {
	strcpy(tmp, I_("(missing values denoted by $-999$ will be "
		"skipped)"));
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pprintf(prn, "\\begin{tabular}{lD{%c}{%c}{-1}"
		"D{%c}{%c}{-1}D{%c}{%c}{-1}D{%c}{%c}{-1}}\n", 
		pt, pt, pt, pt, pt, pt, pt, pt);
	pprintf(prn, "%s &", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1) {
	    tex_escape(vname, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(xbar, prn, 0);
	printftex(summ->xmedian[v], prn, 0);
	printftex(summ->xpx[v], prn, 0);
	printftex(summ->xpy[v], prn, 1);
	if (v == lo) pprintf(prn, "[10pt]\n\n");
	else pprintf(prn, "\n");
    }

    if (lo > 1) pprintf(prn, "%s & ", I_("Variable"));

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Std.\\ Dev."), I_("C.V."), I_("Skewness"), I_("Ex.\\ kurtosis"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	if (lo > 1) {
	    tex_escape(vname, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", vname);
	}
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = -999;
	printftex(std, prn, 0);
	printftex(xcv, prn, 0);
	printftex(summ->xskew[v], prn, 0);
	printftex(summ->xkurt[v], prn, 1);
	pprintf(prn, "\n");
    }

    pprintf(prn, "\\end{tabular}\n\\end{center}\n");
    
}

/* ......................................................... */ 

static void outxx (const double xx, PRN *prn)
{
    if (na(xx)) pprintf(prn, "%s & ", I_("undefined"));
    else pprintf(prn, "$%.4f$ & ", xx);
}

/* ......................................................... */ 

static void rtf_outxx (const double xx, PRN *prn)
{
    if (na(xx)) pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    else pprintf(prn, "\\qc %.4f\\cell ", xx);
}

/* ......................................................... */ 

static void rtf_corr_row (int lo, PRN *prn)
{
    pprintf(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    if (lo == 2) {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx3500\n");
    }
    else if (lo == 3) {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx5000\n");
    }
    else if (lo == 4) {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx6500\n");
    }
    else {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx7500\\cellx8000\n");
    }

    pprintf(prn, "\\intbl ");
}

/* ......................................................... */ 

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pprintf(prn, "\\cell ");
}

/* ......................................................... */ 

void rtfprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char date1[9], date2[9], tmp[128];
    enum { FIELDS = 5 };

    ntodate(date1, corr->t1, pdinfo);
    ntodate(date2, corr->t2, pdinfo);

    sprintf(tmp, I_("Correlation coefficients, using the observations "
		    "%s - %s"), date1, date2);
    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n(%s)\\par\n",
	    tmp, I_("skipping any missing values"));

    sprintf(tmp, I_("5%% critical value (two-tailed) = %.4f for n = %d"), 
	    rhocrit95(corr->n), corr->n);
    pprintf(prn, "%s\\par\n\\par\n{", tmp);
    
    lo = corr->list[0];

    for (i=0; i<=lo/FIELDS; i++) {
	int pad;

	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pad = (lo > FIELDS)? FIELDS - p : lo - p;

	rtf_corr_row(lo, prn);

	if (pad) rtf_table_pad(pad, prn);

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = corr->list[j + nf];
	    pprintf(prn, "%d) %s\\cell %s", ljnf, pdinfo->varname[ljnf],
		    (j == p)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    pprintf(prn, "\\intbl "); 
	    if (pad) rtf_table_pad(pad, prn);
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		rtf_outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", corr->list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    pprintf(prn, "\\intbl "); 
	    rtf_table_pad(pad + j - 1, prn);
	    ij2 = nf + j;
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		rtf_outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", corr->list[ij2]);
	}
    }
    pprintf(prn, "}}\n");
}

/* ......................................................... */ 

void texprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char date1[9], date2[9], vname[16], tmp[128];
    enum { FIELDS = 5 };

    ntodate(date1, corr->t1, pdinfo);
    ntodate(date2, corr->t2, pdinfo);

    lo = corr->list[0];

    sprintf(tmp, I_("Correlation coefficients, using the observations "
		    "%s--%s"), date1, date2);
    pprintf(prn, "\\begin{center}\n%s\\\\\n(%s)\\\\\n", 
	    tmp, I_("skipping any missing values"));

    sprintf(tmp, I_("5\\%% critical value (two-tailed) = %.4f for n = %d"), 
	    rhocrit95(corr->n), corr->n);
    pprintf(prn, "%s\\\\\n", tmp);

    pprintf(prn, "\\vspace{8pt}\n\\begin{tabular}{rrr%s}\n",
	    (lo == 3)? "r" : (lo == 4)? "rr" : "rrr");

    for (i=0; i<=lo/FIELDS; i++) {
	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = corr->list[j + nf];
	    tex_escape(vname, pdinfo->varname[ljnf]);
	    pprintf(prn, "%d) %s%s", ljnf, vname,
		    (j == p)? " &\\\\" : " & ");
	}
	pprintf(prn, "\n");
	
	/* insert spacers */
	for (j=1; j<=p; ++j) 
	    pprintf(prn, "\\rule{13ex}{0pt} & ");
	pprintf(prn, "\\\\\[-6pt]\n");    

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "(%d\\\\\n", corr->list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j-1; k++) pprintf(prn, " & ");
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "(%d\\\\\n", corr->list[ij2]);
	}
	pprintf(prn, "\\\\\n");
    }
    pprintf(prn, "\\end{tabular}\n\\end{center}\n");
}

/* ........................................................... */

static 
void tex_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[9], date2[9]; 

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pprintf(prn, "\\begin{raggedright}\n");
    pprintf(prn, I_("Full data range: %s--%s ($n$ = %d)\\\\\n"),
	    pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    pprintf(prn, I_("Model estimation range: %s--%s"), date1, date2);

    if (fr->nobs == pdinfo->n) pprintf(prn, "\\\\\n");
    else pprintf(prn, " ($n$ = %d)\\\\\n", fr->nobs); 

    pprintf(prn, I_("Standard error of residuals = %f\n"), fr->sigma);
    pprintf(prn, "\\end{raggedright}\n");
}

/* ........................................................... */

void 
texprint_fit_resid (const FITRESID *fr, const DATAINFO *pdinfo, PRN *prn)
{
    int t, anyast = 0;
    int n = pdinfo->n;
    double xx;

    tex_fit_resid_head(fr, pdinfo, prn); 

    pprintf(prn, "\n\\begin{center}\n"
	    "\\begin{tabular}{rrrrl}\n"
	    "\\multicolumn{1}{c}{%s} & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    I_("Obs"), fr->depvar,
	    I_("fitted"), I_("residuals"));

    for (t=0; t<n; t++) {
	if (t == fr->t1 && t) pprintf(prn, "\\\\\n");
	if (t == fr->t2 + 1) pprintf(prn, "\\\\\n");

	print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & ");
	
	if (na(fr->actual[t]) || na(fr->fitted[t])) { 
	    pprintf(prn, "\\\\\n");
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    pprintf(prn, "%13.*f & %13.*f & %13.*f & %s \\\\\n", 
		    fr->pmax, fr->actual[t],
		    fr->pmax, fr->fitted[t], fr->pmax, xx,
		    (ast)? " *" : "");
	}
    }

    pprintf(prn, "\\end{tabular}\n"
	    "\\end{center}\n\n");

    if (anyast) pprintf(prn, I_("\\textit{Note}: * denotes a residual "
				"in excess of 2.5 standard errors\n\n"));
}

/* .................................................................. */

void texprint_fcast_with_errs (const FITRESID *fr, 
			       const DATAINFO *pdinfo, 
			       PRN *prn)
{
    int t;
    double *maxerr;

    maxerr = mymalloc(fr->nobs * sizeof *maxerr);
    if (maxerr == NULL) return;

    pprintf(prn, I_(" For 95%% confidence intervals, t(%d, .025) = %.3f\n"), 
	    fr->pmax, fr->tval);
    pprintf(prn, "\n     Obs ");
    pprintf(prn, "%12s", fr->depvar);
    pprintf(prn, "%*s", UTF_WIDTH(I_("prediction"), 14), _("prediction"));
    pprintf(prn, "%*s", UTF_WIDTH(_(" std. error"), 14), _(" std. error"));
    pprintf(prn, _("   95%% confidence interval\n"));
    pprintf(prn, "\n");

    for (t=0; t<fr->nobs; t++) {
	print_obs_marker(t + fr->t1, pdinfo, prn);
	_printxs(fr->actual[t], 15, PRINT, prn);
	_printxs(fr->fitted[t], 15, PRINT, prn);
	_printxs(fr->sderr[t], 15, PRINT, prn);
	maxerr[t] = fr->tval * fr->sderr[t];
	_printxs(fr->fitted[t] - maxerr[t], 15, PRINT, prn);
	pprintf(prn, " -");
	_printxs(fr->fitted[t] + maxerr[t], 10, PRINT, prn);
	pprintf(prn, "\n");
    }
}

/* .................................................................. */

void augment_copy_menu (windata_t *vwin)
{
    GtkItemFactoryEntry item;
    const char *itempaths[] = {
	N_("/Edit/Copy _all"),
	N_("/Edit/Copy all/as plain _text"),
	N_("/Edit/Copy all/as _LaTeX"),
	N_("/Edit/Copy all/as _RTF")
    };

    item.path = NULL;

    if (gtk_item_factory_get_item(vwin->ifac, "/Edit/Copy all")) {
	gtk_item_factory_delete_item(vwin->ifac, "/Edit/Copy all");
    }

    item.path = mymalloc(64);
    if (item.path == NULL) return;

    /* menu branch */
    sprintf(item.path, itempaths[0]);
    item.callback = NULL;
    item.callback_action = 0;
    item.item_type = "<Branch>";
    item.accelerator = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);

    /* common for menu items */
    item.item_type = NULL;    
    item.accelerator = NULL;
    
    /* plain text option */
    sprintf(item.path, itempaths[1]);
    item.callback = text_copy;
    item.callback_action = COPY_TEXT;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);    

    /* LaTeX option */
    sprintf(item.path, itempaths[2]);
    item.callback = text_copy;
    item.callback_action = COPY_LATEX;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1); 

    /* RTF option */
    sprintf(item.path, itempaths[3]);
    item.callback = text_copy;
    item.callback_action = COPY_RTF;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1); 

    free(item.path);
} 



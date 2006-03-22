#ifdef USE_GNOME

static GdkPixbuf *png_mono_pixbuf (const char *fname)
{
    FILE *fsrc, *ftmp;
    char cmd[MAXLEN], temp[MAXLEN], fline[MAXLEN];
    GdkPixbuf *pbuf = NULL;

    sprintf(temp, "%sgpttmp.XXXXXX", paths.userdir);
    if (mktemp(temp) == NULL) {
	return NULL;
    }

    ftmp = gretl_fopen(temp, "w");
    if (ftmp == NULL) {
	return NULL;
    }

    fsrc = gretl_fopen(fname, "r");
    if (fsrc == NULL) {
	fclose(ftmp);
	remove(temp);
	return NULL;
    }

    fprintf(ftmp, "set term pbm mono\n"
	    "set output '%s%s'\n", 
	    paths.userdir, GRETL_PBM_TMP);

    while (fgets(fline, MAXLEN-1, fsrc)) {
	if (strncmp(fline, "set term", 8) && 
	    strncmp(fline, "set output", 10)) {
	    fputs(fline, ftmp);
	}
    }

    fclose(fsrc);
    fclose(ftmp);

    /* run gnuplot on the temp plotfile */
    sprintf(cmd, "\"%s\" \"%s\"", paths.gnuplot, temp);
    if (system(cmd)) {
	remove(temp);
	return NULL;
    }

    remove(temp);

    build_path(paths.userdir, GRETL_PBM_TMP, temp, NULL);
#if GTK_MAJOR_VERSION >= 2
    pbuf = gdk_pixbuf_new_from_file(temp, NULL);
#else
    pbuf = gdk_pixbuf_new_from_file(temp);
#endif
    remove(temp);

    return pbuf;
}

#endif /* USE_GNOME */

void rtf_print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn)
{
    const char *obs;

    if (pdinfo->markers) { 
	obs = pdinfo->S[t];
    } else {
	char tmp[OBSLEN]; 

	ntodate(tmp, t, pdinfo);
	obs = tmp;
    }

    pprintf(prn, "\\intbl \\ql %s\\cell", obs);
}

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

static void printfrtf (double zz, PRN *prn, int endrow)
{
    /* was using "qr", for right alignment */

    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\qc %s\\cell\\intbl \\row\n",
		    I_("undefined"));
	} else {
	    pprintf(prn, "\\qc %s\\cell", I_("undefined"));
	}
	return;
    }

    if (endrow) {
	pprintf(prn, "\\qc %#.*g\\cell\\intbl \\row\n", GRETL_DIGITS, zz);
    } else {
	pprintf(prn, "\\qc %#.*g\\cell", GRETL_DIGITS, zz);
    }
}

#define SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx1600\\cellx3200\\cellx4800\\cellx6400" \
                   "\\cellx8000\n"

#define VAR_SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                      "\\cellx2000\\cellx4000\\cellx6000\\cellx8000\n"

static void 
rtfprint_summary (const Summary *summ, const DATAINFO *pdinfo, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN], tmp[128];
    int i, vi;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s - %s"),
	    date1, date2);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
    
    if (summ->list[0] == 1) {
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "%s\\par\n\n", tmp);
	pputs(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	if (summ->missing) {
	    strcpy(tmp, I_("(missing values were skipped)"));
	    pprintf(prn, "%s\\par\n\n", tmp); /* FIXME */
	}
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

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printfrtf(summ->mean[i], prn, 0);
	printfrtf(summ->median[i], prn, 0);
	printfrtf(summ->low[i], prn, 0);
	printfrtf(summ->high[i], prn, 1);
    }

    if (summ->list[0] > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
				   I_("Variable"));

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Std. Dev."), I_("C.V."), I_("Skewness"), I_("Ex. kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printfrtf(summ->sd[i], prn, 0);
	printfrtf(summ->cv[i], prn, 0);
	printfrtf(summ->skew[i], prn, 0);
	printfrtf(summ->xkurt[i], prn, 1);
    }

    pputs(prn, "}}\n");
}

static void printftex (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\multicolumn{1}{c}{%s}\\\\", I_("undefined"));
	} else {
	    pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
	}
    } else {
	char s[32];

	tex_dcolumn_double(zz, s);

	if (endrow) {
	    pprintf(prn, "$%s$\\\\", s);
	} else {
	    pprintf(prn, "$%s$ & ", s);
	}
    }	
}

static void 
texprint_summary (const Summary *summ, const DATAINFO *pdinfo, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN], vname[16], tmp[128];
    int i, vi;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s--%s"),
	    date1, date2);

    pprintf(prn, "\\begin{center}\n%s\\\\\n", tmp);
    
    if (summ->list[0] == 1) {
	tex_escape(vname, pdinfo->varname[summ->list[1]]);
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		vname, summ->n);
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pputs(prn, "\\begin{tabular}{rrrr}\n");
    } else {
	if (summ->missing) {
	    pprintf(prn, "%s\\\\[8pt]\n\n", I_("(missing values were skipped)"));
	} else {
	    pputs(prn, "\n\\vspace{8pt}\n\n");
	}
	pputs(prn, "\\begin{tabular}{lrrrr}\n");
	pprintf(prn, "%s &", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(summ->mean[i], prn, 0);
	printftex(summ->median[i], prn, 0);
	printftex(summ->low[i], prn, 0);
	printftex(summ->high[i], prn, 1);
	if (i == summ->list[0] - 1) {
	    pputs(prn, "[10pt]\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    if (summ->list[0] > 1) {
	pprintf(prn, "%s & ", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Std.\\ Dev."), I_("C.V."), I_("Skewness"), I_("Ex.\\ kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(summ->sd[i], prn, 0);
	printftex(summ->cv[i], prn, 0);
	printftex(summ->skew[i], prn, 0);
	printftex(summ->xkurt[i], prn, 1);
	pputc(prn, '\n');
    }

    pputs(prn, "\\end{tabular}\n\\end{center}\n");
}

void special_print_summary (const Summary *summ, const DATAINFO *pdinfo,
			    PRN *prn)
{
    if (tex_format(prn)) {
	texprint_summary(summ, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_summary(summ, pdinfo, prn);
    }
}

static void tex_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "%s & ", I_("undefined"));
    } else {
	pprintf(prn, "$%.4f$ & ", xx);
    }
}

static void rtf_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc %.4f\\cell ", xx);	
    }
}

static void rtf_vmat_row (int lo, PRN *prn)
{
    int i, w = 1400;
    int cmax = (lo + 1 > 6)? 6 : lo + 1;

    pputs(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    for (i=1; i<=cmax; i++) {
	pprintf(prn, "\\cellx%d", w * i);
    }

    pputs(prn, "\n\\intbl ");
}

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pputs(prn, "\\cell ");
}

static void rtf_vmat_blank_row (int lo, int n, PRN *prn)
{
    rtf_vmat_row(lo, prn);
    while (n--) pputs(prn, "\\cell ");
    pputs(prn, "\\intbl \\row\n");
}

#define FIELDS 5

static void
rtfprint_vmatrix (const VMatrix *vmat, const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int n = vmat->t2 - vmat->t1 + 1;
    int blockmax = vmat->dim / FIELDS;
    int nf, li2, p, k, idx, ij2;
    char tmp[128];

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, vmat->t1, pdinfo);
	ntodate(date2, vmat->t2, pdinfo);

	sprintf(tmp, I_("Correlation coefficients, using the observations "
			"%s - %s"), date1, date2);
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
	if (vmat->missing) {
	    pprintf(prn, "%s\\par\n", I_("(missing values were skipped)"));
	}
	sprintf(tmp, I_("5%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pprintf(prn, "%s\\par\n\\par\n{", tmp);
    } else {
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n{",
		I_("Coefficient covariance matrix"));
    }
    
    for (i=0; i<=blockmax; i++) {
	int pad;

	nf = i * FIELDS;
	li2 = vmat->dim - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pad = (vmat->dim > FIELDS)? FIELDS - p : vmat->dim - p;

	rtf_vmat_row(vmat->dim, prn);

	if (pad) rtf_table_pad(pad, prn);

	/* print the varname headings */
	for (j=0; j<p; j++)  {
	    pprintf(prn, "%s\\cell %s", vmat->names[j + nf],
		    (j == p - 1)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    pputs(prn, "\\intbl "); 
	    if (pad) {
		rtf_table_pad(pad, prn);
	    }
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, vmat->dim);
		if (vmat->ci == CORR) {
		    rtf_outxx(vmat->vec[idx], prn);
		} else {
		    printfrtf(vmat->vec[idx], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", vmat->names[j]);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    pputs(prn, "\\intbl "); 
	    rtf_table_pad(pad + j, prn);
	    ij2 = nf + j;
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, vmat->dim);
		if (vmat->ci == CORR) {
		    rtf_outxx(vmat->vec[idx], prn);
		} else {
		    printfrtf(vmat->vec[idx], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", vmat->names[ij2]);
	}

	if (i < blockmax) {
	    rtf_vmat_blank_row(vmat->dim, pad + p + 1, prn);
	}
    }

    pputs(prn, "}}\n");
}

static void
texprint_vmatrix (const VMatrix *vmat, const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int n = vmat->t2 - vmat->t1 + 1;
    int lo, nf, li2, p, k, idx, ij2;
    char vname[16];
    int fields = 5;

    lo = vmat->dim;

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, vmat->t1, pdinfo);
	ntodate(date2, vmat->t2, pdinfo);

	pputs(prn, "\\begin{center}\n");
	pprintf(prn, I_("Correlation coefficients, using the observations "
			"%s--%s"), date1, date2);
	pputs(prn, "\\\\\n");
	if (vmat->missing) {
	    pputs(prn, I_("(missing values were skipped)"));
	    pputs(prn, "\\\\\n");
	}
	pprintf(prn, I_("5\\%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pputs(prn, "\\\\\n");
    } else {
	pprintf(prn, "\\begin{center}\n%s\\\\\n", 
		I_("Coefficient covariance matrix"));
    }

    pputs(prn, "\\vspace{8pt}\n");

    for (i=0; i<=lo/fields; i++) {
	nf = i * fields;
	li2 = lo - nf;
	/* p = number of cols we'll print */
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	pputs(prn, "\\begin{tabular}{");
	for (j=0; j<p; j++) {
	    pputc(prn, 'r');
	}
	pputs(prn, "l}\n");

	/* print the varname headings */
	for (j=0; j<p; j++)  {
	    tex_escape(vname, vmat->names[j + nf]);
	    if (vmat->ci == CORR) {
		pprintf(prn, "%s%s", vname,
			(j == p - 1)? " &\\\\\n" : " & ");
	    } else {
		pprintf(prn, "\\multicolumn{1}{c}{%s}%s", vname,
			(j == p - 1)? " &\\\\\n" : " &\n");
	    }
	}
	
	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, lo);
		if (vmat->ci == CORR) {
		    tex_outxx(vmat->vec[idx], prn);
		} else {
		    printftex(vmat->vec[idx], prn, 0);
		}
	    }
	    tex_escape(vname, vmat->names[j]);
	    pprintf(prn, "%s\\\\\n", vname);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j; k++) {
		pputs(prn, " & ");
	    }
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, lo);
		if (vmat->ci == CORR) {
		    tex_outxx(vmat->vec[idx], prn);
		} else {
		    printftex(vmat->vec[idx], prn, 0);
		}
	    }
	    tex_escape(vname, vmat->names[ij2]);
	    pprintf(prn, "%s\\\\\n", vname);
	}

	pputs(prn, "\\end{tabular}\n\n");
    }

    pputs(prn, "\\end{center}\n");
}

void special_print_vmatrix (const VMatrix *vmat, const DATAINFO *pdinfo, 
			    PRN *prn)
{
    if (tex_format(prn)) {
	texprint_vmatrix(vmat, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_vmatrix(vmat, pdinfo, prn);
    }
}

static 
void tex_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN]; 

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "\\begin{raggedright}\n");
    pputs(prn, I_("Model estimation range:"));
    pprintf(prn, " %s--%s", date1, date2);

    pprintf(prn, I_("Standard error of residuals = %g"), fr->sigma);
    pputs(prn, "\n\\end{raggedright}\n");
}

static 
void rtf_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN]; 
    char tmp[128];

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pputs(prn, I_("Model estimation range:")); 
    pprintf(prn, " %s - %s\\par\n", date1, date2);

    sprintf(tmp, I_("Standard error of residuals = %g"), 
	    fr->sigma);
    pprintf(prn, "\\qc %s\\par\n\\par\n", tmp);
}

static void tex_print_x (double x, int pmax, PRN *prn)
{
    if (x < 0) {
	pputs(prn, "$-$");
    } 

    x = fabs(x);

    if (pmax != PMAX_NOT_AVAILABLE) {
	pprintf(prn, "%.*f", pmax, x);
    } else {
	pprintf(prn, "%g", x);
    }

    pputs(prn, " & ");
}

static void texprint_fit_resid (const FITRESID *fr, 
				const DATAINFO *pdinfo, 
				PRN *prn)
{
    int t, anyast = 0;
    double xx;
    char vname[16];

    tex_fit_resid_head(fr, pdinfo, prn); 

    tex_escape(vname, fr->depvar);

    pprintf(prn, "\n\\begin{center}\n"
	    "\\begin{longtable}{rrrrl}\n"
	    " & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    vname, I_("fitted"), I_("residual"));

    for (t=fr->t1; t<=fr->t2; t++) {
	tex_print_obs_marker(t, pdinfo, prn);
	pputs(prn, " & ");

	if (na(fr->actual[t])) {
	    ;
	} else if (na(fr->fitted[t])) {
	    tex_print_x(fr->actual[t], fr->pmax, prn);
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    tex_print_x(fr->actual[t], fr->pmax, prn);
	    tex_print_x(fr->fitted[t], fr->pmax, prn);
	    tex_print_x(xx, fr->pmax, prn);
	    if (ast) {
		pputs(prn, " *");
	    }
	}
	pputs(prn, " \\\\\n");
    }

    pputs(prn, "\\end{longtable}\n\\end{center}\n\n");

    if (anyast) {
	pputs(prn, I_("\\textit{Note}: * denotes a residual "
		      "in excess of 2.5 standard errors\n\n"));
    }
}

#define FR_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2400\\cellx4000\\cellx5600" \
                "\\cellx6100\n"

static void rtfprint_fit_resid (const FITRESID *fr, 
				const DATAINFO *pdinfo, 
				PRN *prn)
{
    double xx;
    int anyast = 0;
    int t;

    rtf_fit_resid_head(fr, pdinfo, prn);

    pputs(prn, "{" FR_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc \\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\ql \\cell"
	    " \\intbl \\row\n",
	    fr->depvar, I_("fitted"), I_("residual"));

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	if (na(fr->actual[t])) {
	    pputs(prn, "\\qc \\cell \\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else if (na(fr->fitted[t])) {	 
	    printfrtf(fr->actual[t], prn, 0);
	    pputs(prn, "\\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) {
		anyast = 1;
	    }
	    printfrtf(fr->actual[t], prn, 0);
	    printfrtf(fr->fitted[t], prn, 0);
	    printfrtf(xx, prn, 0);
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", 
		    (ast)? "*" : "");
	}
    }

    pputs(prn, "}\n");
    if (anyast) {
	pprintf(prn, "\\par\n\\qc %s \\par\n",
		I_("Note: * denotes a residual in excess of 2.5 standard errors"));
    }
    pputs(prn, "}\n");
}

void special_print_fit_resid (const FITRESID *fr, 
			      const DATAINFO *pdinfo, 
			      PRN *prn)
{
    if (tex_format(prn)) {
	texprint_fit_resid(fr, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_fit_resid(fr, pdinfo, prn);
    }
}

/* .................................................................. */

static void texprint_fcast_x (double x, int places, char *str)
{
    if (places != PMAX_NOT_AVAILABLE && !na(x)) {
	sprintf(str, "%.*f", places, x);
    } else {
	tex_dcolumn_double(x, str);
    }
}

static void texprint_fcast_without_errs (const FITRESID *fr, 
					 const DATAINFO *pdinfo, 
					 PRN *prn)
{
    char actual[32], fitted[32];
    char vname[16];
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "%% The table below needs the \"dcolumn\" and "
	  "\"longtable\" packages\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{longtable}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    D{%c}{%c}{-1}}%% col 3: fitted\n",
	    pt, pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s} \\\\ [4pt] \n",
	    I_("Obs"), vname, I_("prediction"));

    for (t=fr->t1; t<=fr->t2; t++) {
	texprint_fcast_x(fr->actual[t], fr->pmax, actual);
	texprint_fcast_x(fr->fitted[t], fr->pmax, fitted);
	tex_print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & %s & %s \\\\\n",
		actual, fitted);
    }

    pputs(prn, "\\end{longtable}\n\\end{center}\n\n");
}

static void texprint_fcast_with_errs (const FITRESID *fr, 
				      const DATAINFO *pdinfo, 
				      PRN *prn)
{
    double maxerr;
    int pmax = fr->pmax;
    int errpmax = fr->pmax;
    char actual[32], fitted[32], sderr[32], lo[32], hi[32];
    char vname[16];
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "\\begin{center}\n");
    if (fr->model_ci == ARMA) {
	pprintf(prn, _("For 95\\%% confidence intervals, $z(.025) = %.2f$\n\n"), 
		1.96);
    } else {
	pprintf(prn, I_("For 95\\%% confidence intervals, $t(%d, .025) = %.3f$\n\n"), 
		fr->df, fr->tval);
    }
    pputs(prn, "\\end{center}\n");

    pputs(prn, "%% The table below needs the \"dcolumn\" and "
	  "\"longtable\" packages\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{longtable}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    D{%c}{%c}{-1}%% col 3: fitted\n"
	    "      D{%c}{%c}{-1}%% col 4: std error\n"
	    "        D{%c}{%c}{-1}%% col 5: conf int lo\n"
	    "         D{%c}{%c}{-1}}%% col 5: conf int hi\n",
	    pt, pt, pt, pt, pt, pt, pt, pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s}\n"
	    " & \\multicolumn{1}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s} \\\\\n",
	    I_("Obs"), vname,
	    I_("prediction"), I_("std. error"),
	    /* xgettext:no-c-format */
	    I_("95\\% confidence interval"));

    pputs(prn, "& & & & \\multicolumn{1}{c}{low} & "
	  "\\multicolumn{1}{c}{high} \\\\\n");

    if (pmax < 4) {
	errpmax = pmax + 1;
    }

    for (t=fr->t1; t<=fr->t2; t++) {
	double xlo, xhi;

	if (na(fr->sderr[t])) {
	    xlo = xhi = NADBL;
	} else {
	    maxerr = fr->tval * fr->sderr[t];
	    xlo = fr->fitted[t] - maxerr;
	    xhi = fr->fitted[t] + maxerr;
	}
	texprint_fcast_x(fr->actual[t], pmax, actual);
	texprint_fcast_x(fr->fitted[t], pmax, fitted);
	texprint_fcast_x(fr->sderr[t], errpmax, sderr);
	texprint_fcast_x(xlo, pmax, lo);
	texprint_fcast_x(xhi, pmax, hi);
	tex_print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & %s & %s & %s & %s & %s \\\\\n",
		actual, fitted, sderr, lo, hi);
    }

    pputs(prn, "\\end{longtable}\n\\end{center}\n\n");
}

#define FC_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2200\\cellx3600\n"

#define FCE_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                 "\\cellx800\\cellx2200\\cellx3600\\cellx5000" \
                 "\\cellx7800\n"

static void rtfprint_fcast_without_errs (const FITRESID *fr, 
					 const DATAINFO *pdinfo, 
					 PRN *prn)
{
    int t;

    pputs(prn, "{\\rtf1\\par\n\n");

    pputs(prn, "{" FC_ROW "\\intbl ");

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction")); 

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
    }

    pputs(prn, "}}\n");
}

static void rtfprint_fcast_with_errs (const FITRESID *fr, 
				      const DATAINFO *pdinfo, 
				      PRN *prn)
{
    int t;
    double maxerr;
    char tmp[128];

    if (fr->model_ci == ARMA) {
	sprintf(tmp, I_("For 95%% confidence intervals, z(.025) = %.2f"), 
		1.96);
    } else {
	sprintf(tmp, I_("For 95%% confidence intervals, t(%d, .025) = %.3f"), 
		fr->df, fr->tval);
    }
    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n", tmp);

    pputs(prn, "{" FCE_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction"), 
	    I_("std. error"),
	    /* xgettext:no-c-format */
	    I_("95% confidence interval"));

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	maxerr = fr->tval * fr->sderr[t];
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
	printfrtf(fr->sderr[t], prn, 0);
	if (na(fr->sderr[t])) {
	    pputs(prn, "\\qc \\cell \\intbl \\row\n");
	} else {
	    maxerr = fr->tval * fr->sderr[t];
	    pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell \\intbl \\row\n", 
		    GRETL_DIGITS, fr->fitted[t] - maxerr, 
		    GRETL_DIGITS, fr->fitted[t] + maxerr);
	}
    }

    pputs(prn, "}}\n");
}

void special_print_forecast (const FITRESID *fr, 
			     const DATAINFO *pdinfo, 
			     PRN *prn)
{
    if (tex_format(prn)) {
	if (fr->sderr != NULL) {
	    texprint_fcast_with_errs(fr, pdinfo, prn);
	} else {
	    texprint_fcast_without_errs(fr, pdinfo, prn);
	}
    } else if (rtf_format(prn)) {
	if (fr->sderr != NULL) {
	    rtfprint_fcast_with_errs(fr, pdinfo, prn);
	} else {
	    rtfprint_fcast_without_errs(fr, pdinfo, prn);
	}
    }
}

static void 
texprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    char vname[16];

    tex_escape(vname, cf->names[i]);
    pprintf(prn, " %s & ", vname);

    if (isnan(cf->coeff[i])) {
	pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
    } else {
	char coeff[32];

	tex_dcolumn_double(cf->coeff[i], coeff);
	pprintf(prn, "%s & ", coeff);
    }

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "\\multicolumn{2}{c}{%s}", I_("undefined"));
    } else {
	char lo[32], hi[32];

	tex_dcolumn_double(cf->coeff[i] - cf->maxerr[i], lo);
	tex_dcolumn_double(cf->coeff[i] + cf->maxerr[i], hi);
	pprintf(prn, "%s & %s", lo, hi);
    }
    pputs(prn, "\\\\\n");
}

static void texprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    char pt = get_local_decpoint();
    int i;

    pprintf(prn, "$t(%d, .025) = %.3f$\n\n", cf->df, tcrit95(cf->df));

    pputs(prn, "%% The table below needs the \"dcolumn\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{rD{%c}{%c}{-1}D{%c}{%c}{-1}D{%c}{%c}{-1}}\n",
	    pt, pt, pt, pt, pt, pt);

    pprintf(prn, " %s%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{2}{c}{%s}\\\\\n",
	    I_("Variable"), I_("Coefficient"),
	    /* xgettext:no-c-format */
	    I_("95\\% confidence interval"));

    pprintf(prn, " & & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}\\\\\n",
	    I_("low"), I_("high"));

    for (i=0; i<cf->ncoeff; i++) {
	texprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n");
}

static void 
rtfprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    pprintf(prn, "\\qc %s\\cell", cf->names[i]);

    printfrtf(cf->coeff[i], prn, 0);

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell ", 
		GRETL_DIGITS, cf->coeff[i] - cf->maxerr[i], 
		GRETL_DIGITS, cf->coeff[i] + cf->maxerr[i]);
    }
    pputs(prn, " \\intbl \\row\n");
}

#define CF_ROW  "\\trowd \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx2400\\cellx4000\\cellx7200\n" 

static void rtfprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    int i;

    pprintf(prn, "{\\rtf1\\par\n\\qc t(%d, .025) = %.3f\\par\n\\par\n", 
	    cf->df, tcrit95(cf->df));

    pputs(prn, "{" CF_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Variable"), I_("Coefficient"), 
	    /* xgettext:no-c-format */
	    I_("95% confidence interval"));

    for (i=0; i<cf->ncoeff; i++) {
	rtfprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "}}\n");
}

void special_print_confints (const CoeffIntervals *cf, PRN *prn)
{
    if (tex_format(prn)) {
	texprint_confints(cf, prn);
    } else if (rtf_format(prn)) {
	rtfprint_confints(cf, prn);
    }
}

/* copy data to buffer in CSV format and place on clipboard */

#define SCALAR_DIGITS 12

static int data_to_buf_as_csv (const int *list, PRN *prn)
{
    int i, t, l0 = list[0];
    int *pmax = NULL;
    int tsamp = datainfo->t2 - datainfo->t1 + 1;
    double xx;

    if (l0 == 0) return 1;

    if (datainfo->delim == ',' && ',' == datainfo->decpoint) {
	errbox(_("You can't use the same character for "
		 "the column delimiter and the decimal point"));
	return 1;
    }

    pmax = mymalloc(l0 * sizeof *pmax);
    if (pmax == NULL) {
	return 1;
    }

    for (i=1; i<=l0; i++) {
	if (datainfo->vector[list[i]]) {
	    pmax[i-1] = get_precision(&Z[list[i]][datainfo->t1], 
				      tsamp, 8);
	} else {
	    pmax[i-1] = SCALAR_DIGITS;
	}
    }	

    if (datainfo->decpoint != ',') {
	gretl_push_c_numeric_locale();
    }

    /* obs column heading? */
    if (datainfo->S != NULL || datainfo->structure != CROSS_SECTION) {
	pprintf(prn, "obs%c", datainfo->delim);
    }

    /* variable names */
    for (i=1; i<=l0; i++) {
	pprintf(prn, "%s", datainfo->varname[list[i]]);
	pputc(prn, (i < l0)? datainfo->delim : '\n');
    }

    /* actual data values */
    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	csv_obs_to_prn(t, datainfo, prn);
	for (i=1; i<=l0; i++) { 
	    xx = (datainfo->vector[list[i]])? 
		Z[list[i]][t] : Z[list[i]][0];
	    if (na(xx)) {
		pputs(prn, "NA");
	    } else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%.10g", xx);
	    } else {
		pprintf(prn, "%.*f", pmax[i-1], xx);
	    }
	    pputc(prn, (i < l0)? datainfo->delim : '\n');
	}
    }

    if (datainfo->decpoint != ',') {
	gretl_pop_c_numeric_locale();
    }

    if (pmax != NULL) {
	free(pmax);
    }

    return 0;
}

static int real_csv_to_clipboard (const char *liststr)
{
    char line[MAXLINE];
    int *list = NULL;
    PRN *prn = NULL;
    int err = 0;

    sprintf(line, "store csv %s", liststr);
    list = command_list_from_string(line);

    if (list != NULL) {
	err = bufopen(&prn);
	if (!err) {
	    err = data_to_buf_as_csv(list, prn);
	}
	if (!err) {
	    err = prn_to_clipboard(prn, GRETL_FORMAT_CSV);
	}
    }

    free(list);
    gretl_print_destroy(prn);

    return err;
}

int csv_to_clipboard (void)
{
    gretlopt opt = OPT_NONE;
    int err = 0;

    delimiter_dialog(&opt);
    data_save_selection_wrapper(COPY_CSV, GINT_TO_POINTER(opt));

    if (storelist != NULL && *storelist != 0) {
	err = real_csv_to_clipboard(storelist);
	free(storelist);
	storelist = NULL;
    }

    return err;
}

int csv_selected_to_clipboard (void)
{
    char *liststr;
    int err = 0;

    liststr = main_window_selection_as_string();

    if (liststr != NULL) {
	delimiter_dialog(NULL);
	err = real_csv_to_clipboard(liststr);
	free(liststr);
    }

    return err;
}

#include "series_view.h"
#include "fileselect.h"

int csv_copy_listed_vars (windata_t *vwin, int fmt, int action)
{
    int *list = series_view_get_list(vwin);
    PRN *prn = NULL;
    char save_delim = datainfo->delim;
    char save_decpoint = datainfo->decpoint;
    int i, err = 0;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] >= datainfo->v) {
		gui_errmsg(E_DATA);
		return E_DATA;
	    }
	}

	if (fmt == GRETL_FORMAT_CSV) {
	    datainfo->delim = ',';
	    datainfo->decpoint = '.';
	} else {
	    datainfo->delim = '\t';
	}

	if (series_view_is_sorted(vwin)) {
	    prn = vwin_print_sorted_as_csv(vwin);
	    if (prn == NULL) {
		err = 1;
	    }
	} else {
	    err = bufopen(&prn);
	    if (!err) {
		err = data_to_buf_as_csv(list, prn);
	    }
	}
	if (!err) {
	    if (action == W_COPY) {
		err = prn_to_clipboard(prn, fmt);
	    } else {
		file_selector(_("Save data"), EXPORT_CSV, FSEL_DATA_PRN, prn);
	    }
	}
	gretl_print_destroy(prn);
	free(list);
    }

    datainfo->delim = save_delim;
    datainfo->decpoint = save_decpoint;

    return err;
}

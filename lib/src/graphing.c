/*
 *  Copyright (c) by Allin Cottrell and Riccardo Lucchetti
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

/* graphing.c for gretl */

#include "libgretl.h"
#include "var.h"
#include "libset.h"
#include "../../cephes/libprob.h"

#include <unistd.h>

#define GP_DEBUG 0

#ifdef _WIN32
# include <windows.h>
#else
# ifdef USE_GTK2
#  define USE_GSPAWN
#  include <glib.h>
#  include <signal.h>
#  if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#  endif
#  ifndef WEXITSTATUS
#   define WEXITSTATUS(stat_val) ((unsigned)(stat_val) >> 8)
#  endif
#  ifndef WIFEXITED
#   define WIFEXITED(stat_val) (((stat_val) & 255) == 0)
#  endif
# endif /* GTK2 */
#endif /* ! _WIN32 */


static char gnuplot_path[MAXLEN];
static const char *auto_ols_string = "# plot includes automatic OLS line\n";
static int gp_small_font_size;

struct gnuplot_info {
    GnuplotFlags flags;
    int ts_plot;
    int yscale;
    int impulses;
    int lo;
    int ols_ok;
    int t1;
    int t2;
    int npoints;
    int miss;
    int oddman;
    int toomany;
    double xrange;
    const char *yformula;
    double *yvar1;
    double *yvar2;
};

struct plot_type_info {
    PlotType ptype;
    const char *pstr;
};

struct plot_type_info ptinfo[] = {
    { PLOT_REGULAR,        NULL },
    { PLOT_CORRELOGRAM,    "correlogram" },
    { PLOT_CUSUM,          "CUSUM test" },
    { PLOT_FORECAST,       "forecasts with 95 pc conf. interval" },
    { PLOT_FREQ_SIMPLE,    "frequency plot (simple)" },
    { PLOT_FREQ_NORMAL,    "frequency plot (against normal)" },
    { PLOT_FREQ_GAMMA,     "frequency plot (against gamma)" },
    { PLOT_GARCH,          "GARCH residual plot" },
    { PLOT_HURST,          "rescaled range plot" },
    { PLOT_IRFBOOT,        "impulse response plot with quantiles" },
    { PLOT_KERNEL,         "kernel density plot" },
    { PLOT_LEVERAGE,       "leverage/influence plot" },
    { PLOT_MULTI_SCATTER,  "multiple scatterplots" },
    { PLOT_PERIODOGRAM,    "periodogram" },
    { PLOT_RANGE_MEAN,     "range-mean plot" },
    { PLOT_SAMPLING_DIST,  "sampling distribution" },
    { PLOT_TRI_GRAPH,      "TRAMO / X12A tri-graph" },
    { PLOT_VAR_ROOTS,      "VAR inverse roots plot" },
    { PLOT_ELLIPSE,        "confidence ellipse plot" },
    { PLOT_MULTI_IRF,      "multiple impulse responses" },
    { PLOT_PANEL,          "multiple panel plots" },
    { PLOT_TYPE_MAX,       NULL }
};
    
#ifdef USE_GSPAWN

# define SPAWN_DEBUG 0  /* define != 0 for debugging statements */

/**
 * gnuplot_test_command:
 * @cmd: gnuplot command string.
 * 
 * See if the installed version of gnuplot will accept a given
 * command.
 *
 * Returns: 0 if gnuplot successfully handles the given command,
 * 1 on error.
 */

int gnuplot_test_command (const char *cmd)
{
    int ok, ret = 1;
    char errbuf[32];
    int child_pid = 0, sinp = 0, serr = 0;
    GError *error = NULL;
    gchar *argv[] = {
	(*gnuplot_path == 0)? "gnuplot" : gnuplot_path,
	NULL
    };

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_async_with_pipes (NULL,
				   argv,
				   NULL,
				   G_SPAWN_SEARCH_PATH |
				   G_SPAWN_STDOUT_TO_DEV_NULL |
				   G_SPAWN_DO_NOT_REAP_CHILD,
				   NULL,
				   NULL,
				   &child_pid,
				   &sinp,
				   NULL,
				   &serr,
				   &error);

# if SPAWN_DEBUG
    fprintf(stderr, "Testing gnuplot command '%s'\n", cmd);
    fprintf(stderr, "ok=%d, child_pid=%d, sinp=%d\n",
	    ok, child_pid, sinp);
# endif

    if (ok) {
	int test, status;
	int errbytes;

	write(sinp, cmd, strlen(cmd));
	write(sinp, "\n", 1);
	close(sinp);
	test = waitpid(child_pid, &status, 0);
# if SPAWN_DEBUG
	fprintf(stderr, "waitpid returned %d, WIFEXITED %d, "
		"WEXITSTATUS %d\n", test, WIFEXITED(status),
		WEXITSTATUS(status));
# endif
	if (test == child_pid && WIFEXITED(status)) {
	    ret = WEXITSTATUS(status);
	}
	errbytes = read(serr, errbuf, sizeof errbuf - 1);
	if (errbytes > 0) {
	    errbuf[errbytes] = '\0';
	    if (strstr(errbuf, "not find/open font")) {
		ret = 1;
	    }
	} 
	close(serr);
    } else {
	fprintf(stderr, "error: '%s'\n", error->message);
	g_error_free(error);
    }

# if SPAWN_DEBUG
    fprintf(stderr, "gnuplot test: ret = %d\n", ret);
# endif

    return ret;
}

#elif !defined(_WIN32)

#include <signal.h>

int gnuplot_test_command (const char *cmd)
{
    char fullcmd[512];
    int err;

    signal(SIGCHLD, SIG_DFL);

    sprintf(fullcmd, "echo \"%s\" | %s 2>/dev/null", cmd,
	    (*gnuplot_path == 0)? "gnuplot" : gnuplot_path);
    err =  system(fullcmd);

    return err;
}

#endif /* ! USE_GSPAWN, ! WIN32 */

static int printvars (FILE *fp, int t, const int *list, const double **Z,
		      const char *label, double offset)
{
    int i, miss = 0;
    double xx;

    for (i=1; i<=list[0]; i++) {
	xx = Z[list[i]][t];
	if (na(xx)) {
	    fputs("? ", fp);
	    miss = 1;
	} else {
	    if (i == 1) { 
		/* the x variable */
		xx += offset;
	    }
	    fprintf(fp, "%.8g ", xx);
	}
    }

    if (label != NULL) {
	fprintf(fp, "# %s", label);
    }

    fputc('\n', fp);
    
    return miss;
}

static int factorized_vars (struct gnuplot_info *gpinfo,
			    const double **Z, int ynum, int dum)
{
    int fn, t, i = 0;
    double xx;

    fn = gpinfo->t2 - gpinfo->t1 + 1;

    gpinfo->yvar1 = malloc(fn * sizeof *gpinfo->yvar1);
    if (gpinfo->yvar1 == NULL) {
	return 1;
    }

    gpinfo->yvar2 = malloc(fn * sizeof *gpinfo->yvar2);
    if (gpinfo->yvar2 == NULL) {
	free(gpinfo->yvar1);
	return 1;
    }

    for (t=gpinfo->t1; t<=gpinfo->t2; t++) {
	if (na(Z[ynum][t])) {
	    gpinfo->yvar1[i] = NADBL;
	    gpinfo->yvar2[i] = NADBL;
	} else {
	    xx = Z[dum][t];
	    if (floateq(xx, 1.)) {
		gpinfo->yvar1[i] = Z[ynum][t];
		gpinfo->yvar2[i] = NADBL;
	    } else {
		gpinfo->yvar1[i] = NADBL;
		gpinfo->yvar2[i] = Z[ynum][t];
	    }
	}
	i++;
    }

    return 0;
}

#ifdef WIN32

int gnuplot_has_ttf (void)
{
    /* we know the gnuplot supplied with gretl for win32
       does TrueType fonts */
    return 1;
}

int gnuplot_has_pdf (void)
{
    /* ... and we know it does PDF output */
    return 1;
}

int gnuplot_has_specified_colors (void)
{
    /* ... and we know it does specified colors */
    return 1;
}

static int gnuplot_has_specified_emf_colors (void)
{
    /* ... and we know it does specified emf colors */
    return 1;
}

static int gnuplot_has_style_fill (void)
{
    /* ... and that it does style fill */
    return 1;
}

static char *label_front (void)
{
    /* ... and that it handles "front" for labels */
    return " front";
}

#else

int gnuplot_has_ttf (void)
{
    static int err = -1; 

    /* try a range of ttf fonts that might plausibly be installed
       with X11 */

    if (err == -1) {
	err = gnuplot_test_command("set term png font luxisr 8");
	if (err) {
	    err = gnuplot_test_command("set term png font Vera 8");
	}
	if (err) {
	    err = gnuplot_test_command("set term png font arial 8");
	}
    }

    return !err;
}

int gnuplot_has_pdf (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term pdf");
    }

    return !err;
}

int gnuplot_has_specified_colors (void)
{
    static int err = -1; 

    if (err == -1) {
	/* try the old-style command: 
	   if it fails, we have the new driver, we hope */
	err = gnuplot_test_command("set term png color");
    }

    return err;
}

static int gnuplot_has_specified_emf_colors (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term emf color xff0000");
    }

    return !err;
}

static int gnuplot_has_style_fill (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set style fill solid");
    }

    return !err;
}

static char *label_front (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set label 'foo' at 0,0 front");
    }

    if (err) {
	return "";
    } else {
	return " front";
    }
}

#endif /* WIN32 */

static void 
write_gnuplot_font_string (char *fstr, const char *grfont, PlotType ptype)
{
    int shrink = 0;

    if ((ptype == PLOT_MULTI_IRF || ptype == PLOT_MULTI_SCATTER) 
	&& gp_small_font_size > 0) {
	char fname[64];
	int fsize;

	if (sscanf(grfont, "%s %d", fname, &fsize) == 2) {
	    sprintf(fstr, " font %s %d", fname, gp_small_font_size);
	    shrink = 1;
	}
    }

    if (!shrink) {
	sprintf(fstr, " font %s", grfont);
    }
}

/**
 * get_gretl_png_term_line:
 * @ptype: indication of the sort of plot to be made, which
 * may made a difference to the color palette chosen.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the PNG "terminal".  Checks the environment for setting of 
 * %GRETL_PNG_GRAPH_FONT.  Also appends a color-specification string 
 * if the gnuplot PNG driver supports this.
 *
 * Returns: the terminal string, "set term png ..."
 */

const char *get_gretl_png_term_line (PlotType ptype)
{
    static char png_term_line[256];
    char font_string[128];
    char size_string[16];
    char color_string[64];
    int gpcolors = 1, gpttf = 1;
    const char *grfont = NULL;

    *font_string = 0;
    *size_string = 0;
    *color_string = 0;

#ifndef WIN32
    gpcolors = gnuplot_has_specified_colors();
    gpttf = gnuplot_has_ttf();
#endif

    /* plot font setup */
    if (gpttf) {
	grfont = gretl_png_font();
	if (*grfont == 0) {
	    grfont = getenv("GRETL_PNG_GRAPH_FONT");
	}
	if (grfont != NULL && *grfont != 0) {
	    write_gnuplot_font_string(font_string, grfont, ptype);
	}
    }

    /* plot color setup */
    if (gpcolors) {
	int i;

	strcpy(color_string, " xffffff x000000 x202020");
	for (i=0; i<3; i++) {
	    strcat(color_string, " ");
	    strcat(color_string, get_gnuplot_pallette(i, ptype));
	}
    } else {
	strcpy(color_string, " color"); /* old PNG driver */
    }

    if (ptype == PLOT_VAR_ROOTS) {
	strcpy(size_string, " size 480,480");
    }

    sprintf(png_term_line, "set term png%s%s%s",
	    font_string, size_string, color_string);

#if GP_DEBUG
    fprintf(stderr, "png term line:\n'%s'\n", png_term_line);
#endif

    return png_term_line;
}

static void png_font_to_emf (const char *pngfont, char *emfline)
{
    char name[32];
    int pt;

    if (sscanf(pngfont, "%31s %d", name, &pt) == 2) {
	char ptstr[8];

	if (pt <= 8) {
	    pt = 12;
	} else {
	    pt = 16;
	}

	strcat(emfline, "'");
	strcat(emfline, name);
	strcat(emfline, "' ");
	sprintf(ptstr, "%d ", pt);
	strcat(emfline, ptstr);
    }
}

/**
 * get_gretl_emf_term_line:
 * @ptype: indication of the sort of plot to be made.
 * @color: 1 if graph is to be in color, else 0.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the EMF "terminal".  
 *
 * Returns: the term string, "set term emf ..."
 */

const char *get_gretl_emf_term_line (PlotType ptype, int color)
{
    static char emf_term_line[256];
    const char *grfont = NULL;
    
    strcpy(emf_term_line, "set term emf ");

    if (color) {
	strcat(emf_term_line, "color ");
    } else {
	strcat(emf_term_line, "mono dash ");
    }

    /* font spec */
    grfont = gretl_png_font();
    if (grfont != NULL && *grfont != 0) {
	png_font_to_emf(grfont, emf_term_line);
    }

    if (color && gnuplot_has_specified_emf_colors()) {
	if (frequency_plot_code(ptype)) {
	    strcat(emf_term_line, get_gnuplot_pallette(0, PLOT_FREQ_SIMPLE));
	} else {
	    int i;

	    for (i=0; i<3; i++) {
		const char *colstr = get_gnuplot_pallette(i, 0);

		if (*colstr != '\0') {
		    strcat(emf_term_line, colstr);
		    strcat(emf_term_line, " ");
		}
	    }
	}
    } 

    return emf_term_line;
}

/**
 * plot_type_from_string:
 * @str: initial comment line from plot file.
 *
 * Returns: the special plot code corresponding to the initial
 * comment string in the plot file, or %PLOT_REGULAR if no special
 * comment is recognized.
 */

PlotType plot_type_from_string (const char *str)
{
    int i, ret = PLOT_REGULAR;

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (!strcmp(str + 2, ptinfo[i].pstr)) {
	    ret = ptinfo[i].ptype;
	    break;
	}
    }

    return ret;
}

static int write_plot_type_string (PlotType ptype, FILE *fp)
{
    int i, ret = 0;

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (ptype == ptinfo[i].ptype) {
	    fprintf(fp, "# %s\n", ptinfo[i].pstr);
	    ret = 1;
	    break;
	}
    }

    return ret;
}

/**
 * gnuplot_init:
 * @ptype: indication of the sort of plot to be made.
 * @fpp: pointer to stream to be opened.
 *
 * If we're in GUI mode: writes a unique temporary filename into
 * the internal variable #gretl_plotfile; opens plotfile for writing 
 * as @fpp; and writes initial lines into the output file to select 
 * the PNG terminal type, and direct gnuplot's ouput to a temporary
 * file in the gretl user directory.  
 *
 * If not in GUI mode, opens as @fpp the file %gpttmp.plt in the
 * gretl user directory.  
 *
 * This function is not used in batch mode.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gnuplot_init (PlotType ptype, FILE **fpp)
{
    int gui = gretl_in_gui_mode();
    char plotfile[FILENAME_MAX] = {0};

    if (gretl_looping()) {
	return E_OK;
    }

    /* 'gnuplot_path' is file-scope static var */
    if (*gnuplot_path == 0) {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    if (gui) {
	sprintf(plotfile, "%sgpttmp.XXXXXX", gretl_user_dir());
	if (mktemp(plotfile) == NULL) {
	    return E_FOPEN;
	}
    } else {
	sprintf(plotfile, "%sgpttmp.plt", gretl_user_dir());
    }

    set_gretl_plotfile(plotfile);

    *fpp = gretl_fopen(plotfile, "w");
    if (*fpp == NULL) {
	fprintf(stderr, "gnuplot_init: couldn't write to %s\n", plotfile);
	return E_FOPEN;
    } 

    if (gui) {
	fprintf(*fpp, "%s\n", get_gretl_png_term_line(ptype));
	fprintf(*fpp, "set output '%sgretltmp.png'\n", gretl_user_dir());
    }

    write_plot_type_string(ptype, *fpp);

#if GP_DEBUG
    fprintf(stderr, "gnuplot_init: set plotfile = '%s'\n", 
	    plotfile);
#endif

    return 0;
}

#ifdef ENABLE_NLS

#undef RECODE_DBG

static int recode_gnuplot_file (const char *fname)
{
    FILE *fp, *fq;
    const char *font;
    char oldline[512], newline[1024];
    char rname[FILENAME_MAX];
    int ttf = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	return 1;
    }

    strcpy(rname, fname);
    strcat(rname, "l2");

    fq = fopen(rname, "w");
    if (fq == NULL) {
	fclose(fp);
	return 1;
    }

    font = gretl_png_font();
    if (font != NULL && *font != '\0') {
	ttf = 1;
    }

    while (fgets(oldline, sizeof oldline, fp)) {
	if (isdigit((unsigned char) oldline[0])) {
	    fputs(oldline, fq);
	} else if (ttf) {
	    sprint_l2_to_html(newline, oldline, sizeof newline);
	    fputs(newline, fq);
#if RECODE_DBG
	    fprintf(stderr, "recode (sprint_l2_to_html):\n"
		    " original: '%s'\n modified: '%s'\n",
		    oldline, newline);
#endif
	} else {
	    fputs(oldline, fq);
	}
    }

    fclose(fp);
    fclose(fq);
	    
    remove(fname);
    rename(rname, fname);

    return 0;
}

#endif

/**
 * gnuplot_make_graph:
 *
 * Executes gnuplot, passing as an argument the gretl plotfile.
 *
 * Returns: the return value from the system command.
 */

int gnuplot_make_graph (void)
{
    int err = 0;
    char plotcmd[MAXLEN];

#ifdef ENABLE_NLS  
    if (use_latin_2() && gnuplot_has_ttf()) {
# if GP_DEBUG
	fprintf(stderr, "gnuplot_make_graph: calling recode_gnuplot_file()\n");
# endif
	recode_gnuplot_file(gretl_plotfile());
    } 
#endif

#ifdef WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", gretl_gnuplot_path(), gretl_plotfile());
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    sprintf(plotcmd, "%s%s \"%s\"", gretl_gnuplot_path(), 
	    (gretl_in_gui_mode())? "" : " -persist", gretl_plotfile());
    err = gretl_spawn(plotcmd);  
#endif

#if GP_DEBUG
    fprintf(stderr, "gnuplot_make_graph:\n"
	    " plotcmd='%s', err = %d\n", plotcmd, err);
#endif

    return err;
}

enum {
    GTITLE_VLS,
    GTITLE_RESID,
    GTITLE_AF,
    GTITLE_AFV
} graph_titles;

static int auto_plot_var (const char *s)
{
    if (strcmp(s, "t") == 0 ||
	strcmp(s, "annual") == 0 ||
	strcmp(s, "decdate") == 0 ||
	strcmp(s, "qtrs") == 0 ||
	strcmp(s, "months") == 0 ||
	strcmp(s, "hrs") == 0) {
	return 1;
    } else {
	return 0;
    }
}

static void make_gtitle (FILE *fp, int code, const char *s1, const char *s2)
{
    char title[128];
    char depvar[VNAMELEN];

    switch (code) {
    case GTITLE_VLS:
	sprintf(title, I_("%s versus %s (with least squares fit)"),
		s1, s2);
	break;
    case GTITLE_RESID:
	if (sscanf(s1, "residual for %15s", depvar) == 1) {
	    sprintf(title, I_("Regression residuals (= observed - fitted %s)"), 
		    depvar);
	}
	break;
    case GTITLE_AF:
	sprintf(title, I_("Actual and fitted %s"), s1);
	break;
    case GTITLE_AFV:
	if (s2 == NULL || auto_plot_var(s2)) {
	    sprintf(title, I_("Actual and fitted %s"), s1);
	} else {
	    sprintf(title, I_("Actual and fitted %s versus %s"), s1, s2);
	}
	break;
    default:
	*title = '\0';
	break;
    }

    if (*title != '\0') {
	fprintf(fp, "set title '%s'\n", title);
    }
}

static const char *front_strip (const char *s)
{
    while (*s) {
	if (isspace(*s) || *s == '{') {
	    s++;
	} else {
	    break;
	}
    }
	
    return s;
}

static void line_out (const char *s, int len, FILE *fp)
{
    char *p = malloc(len + 1);

    if (p != NULL) {
	*p = 0;
	strncat(p, s, len);
	fprintf(fp, "%s\n", front_strip(p));
	free(p);
    }
}

static void print_gnuplot_literal_lines (const char *s, FILE *fp)
{
    const char *p;

    p = s = front_strip(s);

    while (*s && *s != '}') {
	if (*s == ';') {
	    line_out(p, s - p, fp);
	    p = s + 1;
	}
	s++;
    }
}

static const char *series_name (const DATAINFO *pdinfo, int v)
{
    const char *ret = pdinfo->varname[v];

    if (pdinfo->varinfo != NULL && pdinfo->varinfo[v] != NULL) {
	ret = DISPLAYNAME(pdinfo, v);
	if (ret[0] == '\0') {
	    ret = pdinfo->varname[v];
	}
    } 

    return ret;
}

static int
get_gnuplot_output_file (FILE **fpp, GnuplotFlags flags, 
			 int *plot_count, int code)
{
    const char *plotfile = gretl_plotfile();
    int err = 0;

    *fpp = NULL;

    if ((flags & GP_FILE) && *plotfile != '\0') {
	*fpp = gretl_fopen(plotfile, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	}
    } else if ((flags & GP_BATCH) && plot_count != NULL) {  
	char fname[FILENAME_MAX];

	if (*plotfile == '\0' || strstr(plotfile, "gpttmp") != NULL) {
	    *plot_count += 1; 
	    sprintf(fname, "%sgpttmp%02d.plt", gretl_user_dir(), *plot_count);
	    set_gretl_plotfile(fname);
	} 
	plotfile = gretl_plotfile();
	*fpp = gretl_fopen(plotfile, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	}
    } else {
	/* note: gnuplot_init not used in batch mode */
	err = gnuplot_init(code, fpp);
    }

    return err;
}

static int 
get_ols_line (struct gnuplot_info *gpinfo, const int *list, 
	      double ***pZ, DATAINFO *pdinfo, char *ols_line)
{
    MODEL plotmod;
    char title[48];
    int olslist[4];
    double pval = 1.0;
    int err = 0;

#if GP_DEBUG
    fprintf(stderr, "gnuplot: doing get_ols_line\n");
#endif

    olslist[0] = 3;
    olslist[1] = list[1];
    olslist[2] = 0;
    olslist[3] = list[2];
	
    plotmod = lsq(olslist, pZ, pdinfo, OLS, OPT_A);
    err = plotmod.errcode;

    if (!err) {
	pval = t_pvalue_2(plotmod.coeff[1] / plotmod.sderr[1], plotmod.dfd);
	if (pval < .10) {
	    sprintf(title, "Y = %#.3g %c %#.3gX", plotmod.coeff[0],
		    (plotmod.coeff[1] > 0)? '+' : '-', 
		    fabs(plotmod.coeff[1]));
	}
    }

    gretl_push_c_numeric_locale();

    if (!err && pval < .10) {
	sprintf(ols_line, "%g + %g*x title '%s' w lines\n", 
		plotmod.coeff[0], plotmod.coeff[1], title);
	gpinfo->ols_ok = 1;
    }

    gretl_pop_c_numeric_locale();

    clear_model(&plotmod);

    return err;
}

static void 
print_x_range (struct gnuplot_info *gpinfo, const double *x, FILE *fp)
{
    if (gretl_isdummy(gpinfo->t1, gpinfo->t2, x)) {
	fputs("set xrange [-1:2]\n", fp);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", fp);
	gpinfo->xrange = 3;
    } else {
	double xmin, xmax;

	gretl_minmax(gpinfo->t1, gpinfo->t2, x, &xmin, &xmax);
	gpinfo->xrange = xmax - xmin;
	xmin -= gpinfo->xrange * .025;
	xmax += gpinfo->xrange * .025;
	fprintf(fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);
	gpinfo->xrange = xmax - xmin;
    }
}

/* two or more y vars plotted against some x: test to see if we want
   to use two y axes */

static void
check_for_yscale (struct gnuplot_info *gpinfo, const int *list, 
		  const double **Z)
{
    double ymin[6], ymax[6];
    double ratio;
    int i, j, oddcount;

#if GP_DEBUG
    fprintf(stderr, "gnuplot: doing check_for_yscale\n");
#endif

    /* find minima, maxima of the y-axis vars */
    for (i=1; i<list[0]; i++) {
	gretl_minmax(gpinfo->t1, gpinfo->t2, Z[list[i]], &ymin[i], &ymax[i]);
    }

    gpinfo->yscale = 0;

    for (i=1; i<list[0]; i++) {
	oddcount = 0;
	for (j=1; j<list[0]; j++) {
	    if (j == i) {
		continue;
	    }
	    ratio = ymax[i] / ymax[j];
	    if (ratio > 5.0 || ratio < 0.2) {
		gpinfo->yscale = 1;
		oddcount++;
	    }
	}
	if (oddcount == list[0] - 2) {
	    gpinfo->oddman = i;
	    break;
	}
    }
}

static int 
print_gp_dummy_data (struct gnuplot_info *gpinfo, int *list,
		     const double **Z, const DATAINFO *pdinfo,
		     FILE *fp)
{
    double xx, yy;
    int i, t;

    for (i=0; i<2; i++) {
	for (t=gpinfo->t1; t<=gpinfo->t2; t++) {
	    xx = Z[list[2]][t];
	    if (na(xx)) {
		continue;
	    }
	    yy = (i)? gpinfo->yvar2[t-gpinfo->t1] : gpinfo->yvar1[t-gpinfo->t1];
	    if (na(yy)) {
		fprintf(fp, "%.8g ?\n", xx);
	    } else {
		fprintf(fp, "%.8g %.8g", xx, yy);
		if (!gpinfo->ts_plot) {
		    if (pdinfo->markers) {
			fprintf(fp, " # %s", pdinfo->S[t]);
		    } else if (dataset_is_time_series(pdinfo)) {
			char obs[OBSLEN];

			ntodate(obs, t, pdinfo);
			fprintf(fp, " # %s", obs);
		    }
		}
		fputc('\n', fp);
	    }
	}
	fputs("e\n", fp);
    }

    return 0;
}

static void 
maybe_print_panel_jot (int t, const DATAINFO *pdinfo, FILE *fp)
{
    char obs[OBSLEN];
    int maj, min;

    ntodate(obs, t, pdinfo);
    sscanf(obs, "%d:%d", &maj, &min);
    if (maj > 1 && min == 1) {
	fprintf(fp, "%g ?\n", t + 0.5);
    }
}

static void
print_gp_data (struct gnuplot_info *gpinfo, const int *list,
	       const double **Z, const DATAINFO *pdinfo,
	       FILE *fp)
{
    double offset = 0.0;
    int datlist[3];
    int i, t;

    /* multi impulse plot? calculate offset for lines */
    if ((gpinfo->flags & GP_IMPULSES) && list[gpinfo->lo] > 2) {
	offset = 0.10 * gpinfo->xrange / gpinfo->npoints;
    }	

    datlist[0] = 2;
    datlist[1] = list[gpinfo->lo];

    for (i=1; i<gpinfo->lo; i++) {
	double xoff = offset * (i - 1);

	datlist[2] = list[i];
	for (t=gpinfo->t1; t<=gpinfo->t2; t++) {
	    const char *label = NULL;
	    char obs[OBSLEN];
	    int t_miss;

	    if (!gpinfo->ts_plot && i == 1) {
		if (pdinfo->markers) {
		    label = pdinfo->S[t];
		} else if (dataset_is_time_series(pdinfo)) {
		    ntodate(obs, t, pdinfo);
		    label = obs;
		}
	    }

	    if (gpinfo->ts_plot && pdinfo->structure == STACKED_TIME_SERIES) {
		maybe_print_panel_jot(t, pdinfo, fp);
	    }

	    t_miss = printvars(fp, t, datlist, Z, label, xoff);

	    if ((gpinfo->flags & GP_GUI) && gpinfo->miss == 0) {
		gpinfo->miss = t_miss;
	    }
	}
	fputs("e\n", fp);
    }
}

static void 
gp_info_init (struct gnuplot_info *gpinfo, GnuplotFlags flags,
	      int lo, const char *literal, int t1, int t2)
{
    gpinfo->flags = flags;
    gpinfo->ts_plot = 1;   /* plotting against time on x-axis? */
    gpinfo->yscale = 0;    /* two y axis scales needed? */
    gpinfo->impulses = 0;  /* plotting with impulses? */
    gpinfo->lo = lo;
    gpinfo->ols_ok = 0;    /* plot automatic OLS line? */
    gpinfo->t1 = t1;
    gpinfo->t2 = t2;
    gpinfo->npoints = 0;
    gpinfo->miss = 0;
    gpinfo->oddman = 0;
    gpinfo->xrange = 0.0;
    gpinfo->toomany = 0;
    gpinfo->yformula = NULL;

    if (lo > MAX_PLOT_LINES + 1) {
	gpinfo->toomany = 1;
    }

    if (lo > 2 && lo < 7 && !(flags & GP_RESIDS) && !(flags & GP_FA)
	&& !(flags & GP_DUMMY)) {
	/* allow probe for using two y axes */
	gpinfo->yscale = 1;
    } 

    if (flags & GP_IMPULSES) {
	gpinfo->impulses = 1;
    }    

    if (flags & GP_FA && literal != NULL && 
	!strncmp(literal, "yformula: ", 10)) {
	/* fitted vs actual plot with fitted given by formula */
	gpinfo->yformula = literal + 10;
    }

    if (literal != NULL && strstr(literal, "set style data")) {
	gpinfo->flags |= GP_DATA_STYLE;
    }

    gpinfo->yvar1 = gpinfo->yvar2 = NULL;
}

#if GP_DEBUG
static void print_gnuplot_flags (GnuplotFlags flags)
{
    fprintf(stderr, "*** gnuplot() called with flags:\n");

    if (flags & GP_IMPULSES) {
	fprintf(stderr, " GP_IMPULSES\n");
    }
    if (flags & GP_RESIDS) {
	fprintf(stderr, " GP_RESIDS\n");
    }	
    if (flags & GP_FA) {
	fprintf(stderr, " GP_FA\n");
    }
    if (flags & GP_DUMMY) {
	fprintf(stderr, " GP_DUMMY\n");
    }
    if (flags & GP_BATCH) {
	fprintf(stderr, " GP_BATCH\n");
    }
    if (flags & GP_GUI) {
	fprintf(stderr, " GP_GUI\n");
    }
    if (flags & GP_OLS_OMIT) {
	fprintf(stderr, " GP_OLS_OMIT\n");
    }
    if (flags & GP_DATA_STYLE) {
	fprintf(stderr, " GP_DATA_STYLE\n");
    }
    if (flags & GP_FILE) {
	fprintf(stderr, " GP_FILE\n");
    }
}
#endif

static void set_withstr (GnuplotFlags flags, const int *lines, 
			 int i, char *str)
{
    int ltest = 0;

    if (flags & GP_DATA_STYLE) {
	*str = 0;
	return;
    }

    if (lines != NULL) {
	ltest = lines[(flags & GP_GUI)? i - 1 : 0];
    }

    if (ltest) {
	strcpy(str, "w lines");
    } else {
	strcpy(str, "w points");
    }
}

static void graph_list_adjust_sample (int *list, 
				      struct gnuplot_info *ginfo,
				      const double **Z)
{
    int t1min = ginfo->t1;
    int t2max = ginfo->t2;
    int t_ok;
    int i, t;

    for (t=t1min; t<=t2max; t++) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    if (!na(Z[list[i]][t])) {
		t_ok = 1;
		break;
	    }
	}
	if (t_ok) {
	    break;
	} 
	t1min++;
    }

    for (t=t2max; t>t1min; t--) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    if (!na(Z[list[i]][t])) {
		t_ok = 1;
		break;
	    }
	}
	if (t_ok) {
	    break;
	}
	t2max--;
    }

    if (t2max > t1min) {
	for (i=1; i<=list[0]; i++) {
	    int all_missing = 1;

	    for (t=t1min; t<=t2max; t++) {
		if (!na(Z[list[i]][t])) {
		    all_missing = 0;
		    break;
		}
	    }
	    if (all_missing) {
		gretl_list_delete_at_pos(list, i);
		i--;
	    }
	}
    }

    ginfo->t1 = t1min;
    ginfo->t2 = t2max;
    ginfo->lo = list[0];
}

/**
 * gnuplot:
 * @list: list of variables to plot, by ID number.
 * @lines: vector of 1s and 0s to indicate whether variables should
 * be represented by lines or not (or NULL).
 * @literal: commands to be passed to gnuplot.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @plot_count: pointer to count of graphs drawn so far.
 * @flags: bitwise OR of zero or more options from #GnuplotFlags
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, -1 if the gnuplot system
 * command fails, or 1 if there are missing data values.
 */

int gnuplot (int *list, const int *lines, const char *literal,
	     double ***pZ, DATAINFO *pdinfo, 
	     int *plot_count, GnuplotFlags flags)
{
    FILE *fp = NULL;
    char s1[MAXDISP] = {0};
    char s2[MAXDISP] = {0};
    char xlabel[MAXDISP] = {0};
    char withstr[16] = {0};
    char keystr[48] = {0};
    char ols_line[128] = {0};
    int xvar, i;

    struct gnuplot_info gpinfo;

    *gretl_errmsg = '\0';

#if GP_DEBUG
    print_gnuplot_flags(flags);
#endif

    /* below: did have "height 1 width 1 box" for win32,
       "width 1 box" otherwise */
    strcpy(keystr, "set key left top\n");

    if (get_gnuplot_output_file(&fp, flags, plot_count, PLOT_REGULAR)) {
	return E_FOPEN;
    } 

    gp_info_init(&gpinfo, flags, list[0], literal, pdinfo->t1, pdinfo->t2);

    if (gpinfo.impulses) {
	strcpy(withstr, "w i");
    }

    /* check whether we're doing a time-series plot */
    if (!strcmp(pdinfo->varname[list[gpinfo.lo]], "time")) {
	int pv;

	pv = plotvar(pZ, pdinfo);
	if (pv > 0) {
	    list[gpinfo.lo] = pv;
	} else {
	    fclose(fp);
	    return E_ALLOC;
	}
    } else if (!auto_plot_var(pdinfo->varname[list[gpinfo.lo]])) {
	gpinfo.ts_plot = 0;
    } 

    /* set x-axis label for non-time series plots */
    if (!gpinfo.ts_plot) {
	if (flags & GP_DUMMY) {
	    strcpy(xlabel, series_name(pdinfo, list[2])); 
	} else {
	    strcpy(xlabel, series_name(pdinfo, list[gpinfo.lo]));
	}
    }	

    /* add a simple regression line if appropriate */
    if (!gpinfo.impulses && !(flags & GP_OLS_OMIT) && gpinfo.lo == 2 && 
	!gpinfo.ts_plot && !(flags & GP_RESIDS)) {
	get_ols_line(&gpinfo, list, pZ, pdinfo, ols_line);
	if (gpinfo.ols_ok) {
	    fprintf(fp, "# X = '%s' (%d)\n", pdinfo->varname[list[2]], list[2]);
	    fprintf(fp, "# Y = '%s' (%d)\n", pdinfo->varname[list[1]], list[1]);
	}
    }

    /* adjust sample range, and reject if it's empty */
    graph_list_adjust_sample(list, &gpinfo, (const double **) *pZ);
    if (gpinfo.t1 == gpinfo.t2 || gpinfo.lo == 1) {
	fclose(fp);
	return GRAPH_NO_DATA;
    }

    /* separation by dummy: create special vars */
    if (flags & GP_DUMMY) { 
	if (gpinfo.lo != 3 ||
	    factorized_vars(&gpinfo, (const double **) *pZ, list[1], list[3])) {
	    fclose(fp);
	    return -1;
	}
    } 

    gpinfo.npoints = gpinfo.t2 - gpinfo.t1 + 1;

    /* special tics for short time series plots */
    if (gpinfo.ts_plot) {
	if (gpinfo.toomany) {
	    fprintf(fp, "# multiple timeseries %d\n", pdinfo->pd);
	} else {
	    fprintf(fp, "# timeseries %d\n", pdinfo->pd);
	}
	if (pdinfo->pd == 4 && (gpinfo.t2 - gpinfo.t1) / 4 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 4\n", fp);
	}
	if (pdinfo->pd == 12 && (gpinfo.t2 - gpinfo.t1) / 12 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 12\n", fp);
	}
    } else if (gpinfo.toomany) {
	fputs("# multiple data series\n", fp);
    }

    fprintf(fp, "set xlabel '%s'\n", xlabel);
    fputs("set xzeroaxis\n", fp); 
    fputs("set missing \"?\"\n", fp);

    if (gpinfo.lo == 2) {
	/* only two variables */
	if (gpinfo.ols_ok) {
	    fputs(auto_ols_string, fp);
	    if (flags & GP_FA) {
		make_gtitle(fp, GTITLE_AFV, series_name(pdinfo, list[1]), 
			    series_name(pdinfo, list[2]));
	    } else {
		make_gtitle(fp, GTITLE_VLS, series_name(pdinfo, list[1]), 
			    xlabel);
	    }
	}
	if (flags & GP_RESIDS && !(flags & GP_DUMMY)) { 
	    make_gtitle(fp, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	    fprintf(fp, "set ylabel '%s'\n", I_("residual"));
	} else {
	    fprintf(fp, "set ylabel '%s'\n", series_name(pdinfo, list[1]));
	}
	if (!gpinfo.ols_ok) {
	    strcpy(keystr, "set nokey\n");
	}
    } else if ((flags & GP_RESIDS) && (flags & GP_DUMMY)) { 
	make_gtitle(fp, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	fprintf(fp, "set ylabel '%s'\n", I_("residual"));
    } else if (flags & GP_FA) {
	if (list[3] == pdinfo->v - 1) { 
	    /* x var is just time or index: is this always right? */
	    make_gtitle(fp, GTITLE_AF, series_name(pdinfo, list[2]), NULL);
	} else {
	    make_gtitle(fp, GTITLE_AFV, series_name(pdinfo, list[2]), 
			series_name(pdinfo, list[3]));
	}
	fprintf(fp, "set ylabel '%s'\n", series_name(pdinfo, list[2]));
    } 

    if (gpinfo.toomany) {
	strcpy(keystr, "set key outside\n");
    }

    fputs(keystr, fp);

    gretl_push_c_numeric_locale();

    xvar = (flags & GP_DUMMY)? list[gpinfo.lo - 1] : list[gpinfo.lo];
    print_x_range(&gpinfo, (const double *) (*pZ)[xvar], fp);

    if (gpinfo.yscale) { 
	check_for_yscale(&gpinfo, list, (const double **) *pZ);
	if (gpinfo.yscale) {
	    fputs("set ytics nomirror\n", fp);
	    fputs("set y2tics\n", fp);
	}
    }

#if GP_DEBUG    
    fprintf(stderr, "literal = '%s', yformula = '%s'\n", literal,
	    gpinfo.yformula);
#endif

    if (gpinfo.yformula != NULL) {
	/* cut out the "dummy" yvar that is in fact represented
	   by a formula rather than raw data */
	list[1] = list[2];
	list[2] = list[3];
	list[0] = gpinfo.lo = 2;
    } else if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    /* now print the 'plot' lines */
    fputs("plot \\\n", fp);
    if (gpinfo.yscale) {
	for (i=1; i<gpinfo.lo; i++) {
	    fprintf(fp, "'-' using 1:2 axes %s title '%s (%s)' %s%s",
		    (i == gpinfo.oddman)? "x1y2" : "x1y1",
		    series_name(pdinfo, list[i]), 
		    (i == gpinfo.oddman)? I_("right") : I_("left"),
		    (gpinfo.impulses)? "w impulses" : 
		    (gpinfo.ts_plot)? "w lines" : "w points",
		    (i == gpinfo.lo - 1)? "\n" : " , \\\n");
	}
    } else if (flags & GP_DUMMY) { 
	strcpy(s1, (flags & GP_RESIDS)? I_("residual") : 
	       series_name(pdinfo, list[1]));
	strcpy(s2, series_name(pdinfo, list[3]));
	fprintf(fp, " '-' using 1:2 title '%s (%s=1)', \\\n", s1, s2);
	fprintf(fp, " '-' using 1:2 title '%s (%s=0)'\n", s1, s2);
    } else if (gpinfo.yformula != NULL) {
	fprintf(fp, " '-' using 1:2 title '%s' w points , \\\n", I_("actual"));	
	fprintf(fp, "%s title '%s' w lines\n", gpinfo.yformula, I_("fitted"));
    } else if (flags & GP_FA) {
	set_withstr(flags, lines, 1, withstr);
	fprintf(fp, " '-' using 1:2 title '%s' %s lt 2, \\\n", I_("fitted"), withstr);
	set_withstr(flags, lines, 2, withstr);
	fprintf(fp, " '-' using 1:2 title '%s' %s lt 1\n", I_("actual"), withstr);	
    } else {
	for (i=1; i<gpinfo.lo; i++)  {
	    if (gpinfo.lo == 2) {
		*s1 = '\0';
	    } else {
		strcpy(s1, series_name(pdinfo, list[i]));
	    }
	    if (!gpinfo.impulses) { 
		set_withstr(gpinfo.flags, lines, i, withstr);
	    }
	    fprintf(fp, " '-' using 1:2 title '%s' %s", s1, withstr);
	    if (i < gpinfo.lo - 1 || gpinfo.ols_ok) {
	        fputs(" , \\\n", fp); 
	    } else {
	        fputc('\n', fp);
	    }
	}
    } 

    if (*ols_line != '\0') {
        fputs(ols_line, fp);
    }

    /* print the data to be graphed */
    if (flags & GP_DUMMY) {
	print_gp_dummy_data(&gpinfo, list, (const double **) *pZ, pdinfo,
			    fp);
	free(gpinfo.yvar1);
	free(gpinfo.yvar2);
    } else {
	print_gp_data(&gpinfo, list, (const double **) *pZ, pdinfo, fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (!(flags & GP_BATCH)) {
	if (gnuplot_make_graph()) {
	    gpinfo.miss = -1;
	}
    }

    return gpinfo.miss;
}

/**
 * multi_scatters:
 * @list: list of variables to plot, by ID number.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @plot_count: count of graphs shown to date.
 * @flags: option flags.
 *
 * Writes a gnuplot plot file to display up to 6 small graphs
 * based on the variables in @list, and calls gnuplot to make 
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_scatters (const int *list, double ***pZ, 
		    DATAINFO *pdinfo, int *plot_count, 
		    GnuplotFlags flags)
{
    int xvar = 0, yvar = 0, timevar = 0;
    int *plotlist = NULL;
    int pos, nplots = 0;
    FILE *fp = NULL;
    int i, t, err = 0;

    pos = gretl_list_separator_position(list);

    if (pos == 0) {
	/* plot against time or index */
	timevar = plotvar(pZ, pdinfo);
	if (timevar < 0) {
	    return E_ALLOC;
	}
	plotlist = gretl_list_copy(list);
	flags |= GP_LINES;
    } else if (pos > 2) { 
	/* plot several yvars against one xvar */
	plotlist = gretl_list_new(pos - 1);
	xvar = list[list[0]];
    } else {       
	/* plot one yvar against several xvars */
	plotlist = gretl_list_new(list[0] - pos);
	yvar = list[1];
    }

    if (plotlist == NULL) {
	return E_ALLOC;
    }

    if (yvar) {
	for (i=1; i<=plotlist[0]; i++) {
	   plotlist[i] = list[i + pos]; 
	}
    } else if (xvar) {
	for (i=1; i<pos; i++) {
	   plotlist[i] = list[i]; 
	}
    }

    /* max 6 plots */
    if (plotlist[0] > 6) {
	plotlist[0] = 6;
    }

    nplots = plotlist[0];
    gp_small_font_size = (nplots > 4)? 6 : 0;

    if (get_gnuplot_output_file(&fp, flags, plot_count, PLOT_MULTI_SCATTER)) {
	return E_FOPEN;
    }

    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (timevar) {
	double startdate = (*pZ)[timevar][pdinfo->t1];
	int jump = (pdinfo->t2 - pdinfo->t1 + 1) / (2 * pdinfo->pd);

	fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);
    } else {
	fputs("set noxtics\nset noytics\n", fp);
    }

    for (i=0; i<nplots; i++) {  
	int pv = plotlist[i+1];

	if (nplots <= 4) {
	    fputs("set size 0.45,0.5\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.5,0.5\n", fp);
	    else if (i == 2) fputs("0.0,0.0\n", fp);
	    else if (i == 3) fputs("0.5,0.0\n", fp);
	} else {
	    fputs("set size 0.31,0.45\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.32,0.5\n", fp);
	    else if (i == 2) fputs("0.64,0.5\n", fp);
	    else if (i == 3) fputs("0.0,0.0\n", fp);
	    else if (i == 4) fputs("0.32,0.0\n", fp);
	    else if (i == 5) fputs("0.64,0.0\n", fp);
	}

	if (timevar) {
	    fputs("set noxlabel\n", fp);
	    fputs("set noylabel\n", fp);
	    fprintf(fp, "set title '%s'\n", pdinfo->varname[pv]);
	} else {
	    fprintf(fp, "set xlabel '%s'\n",
		    (yvar)? pdinfo->varname[pv] :
		    pdinfo->varname[xvar]);
	    fprintf(fp, "set ylabel '%s'\n", 
		    (yvar)? pdinfo->varname[yvar] :
		    pdinfo->varname[pv]);
	}

	fputs("plot '-' using 1:2", fp);
	if (flags & GP_LINES) {
	    fputs(" with lines", fp);
	}
	fputc('\n', fp);

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    double xx;
	    int m;

	    m = (yvar)? plotlist[i+1] : (xvar)? xvar : timevar;
	    xx = (*pZ)[m][t];

	    if (na(xx)) {
		fputs("? ", fp);
	    } else {
		fprintf(fp, "%.8g ", xx);
	    }

	    m = (yvar)? yvar : pv;
	    xx = (*pZ)[m][t];

	    if (na(xx)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.8g\n", xx);
	    }
	}

	fputs("e\n", fp);
    } 

    gretl_pop_c_numeric_locale();

    fputs("set nomultiplot\n", fp);

    fclose(fp);

    if (!(flags & GP_BATCH)) {
	err = gnuplot_make_graph();
    }

    free(plotlist);

    return err;
}

static int get_3d_output_file (FILE **fpp)
{
    char fname[MAXLEN];
    int err = 0;

    sprintf(fname, "%sgpttmp.plt", gretl_user_dir());
    *fpp = gretl_fopen(fname, "w");

    if (*fpp == NULL) {
	err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
    }

    return err;
}

static void 
maybe_add_surface (const int *list, double ***pZ, DATAINFO *pdinfo, 
		   unsigned int flags, char *surface)
{
    MODEL smod;
    double umin, umax, vmin, vmax;
    int olslist[5];

    olslist[0] = 4;
    olslist[1] = list[3];
    olslist[2] = 0;
    olslist[3] = list[2];
    olslist[4] = list[1];

    gretl_minmax(pdinfo->t1, pdinfo->t2, (*pZ)[list[2]], &umin, &umax);
    gretl_minmax(pdinfo->t1, pdinfo->t2, (*pZ)[list[1]], &vmin, &vmax);

    smod = lsq(olslist, pZ, pdinfo, OLS, OPT_A);

    if (!smod.errcode && !na(smod.fstt) &&
	(fdist(smod.fstt, smod.dfn, smod.dfd) < .10 || flags & GP_FA)) {
	double uadj = (umax - umin) * 0.02;
	double vadj = (vmax - vmin) * 0.02;

	sprintf(surface, "[u=%g:%g] [v=%g:%g] "
		"%g+(%g)*u+(%g)*v title '', ", 
		umin - uadj, umax + uadj, 
		vmin - vadj, vmax + vadj,
		smod.coeff[0], smod.coeff[1],
		smod.coeff[2]);
    } 

    clear_model(&smod);
}

/**
 * gnuplot_3d:
 * @list: list of variables to plot, by ID number: Y, X, Z
 * @literal: literal command(s) to pass to gnuplot (or NULL)
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @plot_count: pointer to counter variable for plots (unused)
 * @flags: unused at present.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (int *list, const char *literal,
		double ***pZ, DATAINFO *pdinfo,  
		int *plot_count, GnuplotFlags flags)
{
    FILE *fq = NULL;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int orig_t1 = pdinfo->t1, orig_t2 = pdinfo->t2;
    int lo = list[0];
    int datlist[4];
    char surface[128] = {0};

    if (lo != 3) {
	fprintf(stderr, "gnuplot_3d needs three variables (only)\n");
	return -1;
    }

    if (get_3d_output_file(&fq)) {
	return E_FOPEN;
    }

    varlist_adjust_sample(list, &t1, &t2, (const double **) *pZ);

    /* if resulting sample range is empty, complain */
    if (t2 == t1) {
	fclose(fq);
	return GRAPH_NO_DATA;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    gretl_push_c_numeric_locale();

    maybe_add_surface(list, pZ, pdinfo, flags, surface);

    fprintf(fq, "set xlabel '%s'\n", series_name(pdinfo, list[2]));
    fprintf(fq, "set ylabel '%s'\n", series_name(pdinfo, list[1]));
    fprintf(fq, "set zlabel '%s'\n", series_name(pdinfo, list[3]));

    fputs("set missing \"?\"\n", fq);

    if (literal != NULL && *literal != 0) {
	print_gnuplot_literal_lines(literal, fq);
    }

#ifdef WIN32
    fprintf(fq, "splot %s'-' title '' w p lt 3\n", surface);
#else
    fprintf(fq, "splot %s'-' title ''\n", surface);
#endif

    datlist[0] = 3;
    datlist[1] = list[2];
    datlist[2] = list[1];
    datlist[3] = list[3];

    for (t=t1; t<=t2; t++) {
	const char *label = NULL;

	if (pdinfo->markers) {
	    label = pdinfo->S[t];
	}
	printvars(fq, t, datlist, (const double **) *pZ, label, 0.0);
    }	
    fputs("e\n", fq);

    gretl_pop_c_numeric_locale();

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    fclose(fq);

    return 0;
}

static void print_freq_test_label (char *s, const char *format, 
				   double v, double pv)
{
    gretl_pop_c_numeric_locale();
    sprintf(s, format, v, pv);
    gretl_push_c_numeric_locale();
}

static void print_freq_dist_label (char *s, int dist, double x, double y)
{
    int dcomma = 0;
#ifdef ENABLE_NLS
    char test[8];

    gretl_pop_c_numeric_locale();
    sprintf(test, "%g", 0.5);
    if (strchr(test, ',')) {
	dcomma = 1;
    }
#endif

    if (dist == DIST_NORMAL) {
	sprintf(s, "N(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    } else if (dist == DIST_GAMMA) {
	sprintf(s, "gamma(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    }

#ifdef ENABLE_NLS
    gretl_push_c_numeric_locale();
#endif
}

/* Below: a fix for the case where the y-range is by default
   degenerate, in which case gnuplot produces a graph OK, but
   issues a warning and returns non-zero.
*/

static void maybe_set_yrange (FreqDist *freq, double lambda, FILE *fp)
{
    double ymin = 1.0e+20;
    double ymax = -1.0e+20;
    int i;

    for (i=0; i<freq->numbins; i++) { 
	if (freq->f[i] > ymax) {
	    ymax = freq->f[i];
	}
	if (freq->f[i] < ymin) {
	    ymin = freq->f[i];
	}	
    }

    if (ymax == ymin) {
	fprintf(fp, "set yrange [%g:%g]\n", ymax * lambda * 0.99, 
		ymax * lambda * 1.01);
    } else {
	fprintf(fp, "set yrange [0.0:%g]\n", ymax * lambda * 1.1);
    }	
}

/**
 * plot_freq:
 * @freq: pointer to frequency distribution struct.
 * @dist: probability distribution code.
 *
 * Plot the actual frequency distribution for a variable versus a
 * theoretical distribution: Gaussian, gamma or none.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot_freq (FreqDist *freq, DistCode dist)
{
    double alpha = 0.0, beta = 0.0, lambda = 1.0;
    FILE *fp = NULL;
    int i, K = freq->numbins;
    char withstr[16] = {0};
    char label[80] = {0};
    double plotmin = 0.0, plotmax = 0.0;
    double barwidth, barskip;
    int plottype, use_boxes = 1;
    int err;

    if (K == 0) {
	return 1;
    }

    if (dist == DIST_NORMAL) {
	plottype = PLOT_FREQ_NORMAL;
    } else if (dist == DIST_GAMMA) {
	plottype = PLOT_FREQ_GAMMA;
    } else {
	plottype = PLOT_FREQ_SIMPLE;
    }

    if ((err = gnuplot_init(plottype, &fp))) {
	return err;
    }  

#if GP_DEBUG
    fprintf(stderr, "*** plot_freq called\n");
#endif  

    barwidth = freq->endpt[K-1] - freq->endpt[K-2];
    barskip = 0.005 * (freq->endpt[K] - freq->endpt[0]);

    if (K > 16) {
	barskip /= 2.0;
    }

    gretl_push_c_numeric_locale();

    if (dist) {
	lambda = 1.0 / (freq->n * barwidth - barskip);

	if (dist == DIST_NORMAL) {
	    fputs("# literal lines = 4\n", fp);
	    fprintf(fp, "sigma = %g\n", freq->sdx);
	    fprintf(fp, "mu = %g\n", freq->xbar);

	    plotmin = freq->endpt[0] - barwidth;
	    if (plotmin > freq->xbar - 3.3 * freq->sdx) {
		plotmin = freq->xbar - 3.3 * freq->sdx;
	    }
	    plotmax = freq->endpt[K-1] + barwidth;
	    if (plotmax < freq->xbar + 3.3 * freq->sdx) {
		plotmax = freq->xbar + 3.3 * freq->sdx;
	    }

	    if (!na(freq->test)) {
		fprintf(fp, "set label \"%s:\" at graph .03, graph .97%s\n",
			I_("Test statistic for normality"),
			label_front());
		print_freq_test_label(label, I_("Chi-squared(2) = %.3f pvalue = %.5f"), 
				      freq->test, chisq(freq->test, 2));
		fprintf(fp, "set label '%s' at graph .03, graph .93%s\n", 
			label, label_front());
	    }	
	} else if (dist == DIST_GAMMA) {
	    double var = freq->sdx * freq->sdx;

	    /* scale param = variance/mean */
	    beta = var / freq->xbar;
	    /* shape param = mean/scale */
	    alpha = freq->xbar / beta;

	    fputs("# literal lines = 4\n", fp);
	    fprintf(fp, "beta = %g\n", beta);
	    fprintf(fp, "alpha = %g\n", alpha);
	    plotmin = 0.0;
	    plotmax = freq->xbar + 4.0 * freq->sdx;

	    if (!na(freq->test)) {
		fprintf(fp, "set label '%s:' at graph .03, graph .97%s\n",
			I_("Test statistic for gamma"),
			label_front());
		print_freq_test_label(label, I_("z = %.3f pvalue = %.5f"), 
				      freq->test, normal_pvalue_2(freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93%s\n", 
			label, label_front());
	    }	
	}

	/* adjust min, max if needed */
	if (freq->midpt[0] < plotmin) {
	    plotmin = freq->midpt[0];
	}
	if (freq->midpt[K-1] > plotmax) {
	    plotmax = freq->midpt[K-1];
	}

	fprintf(fp, "set xrange [%.7g:%.7g]\n", plotmin, plotmax);
	fputs("set key right top\n", fp);
    } else { 
	/* plain frequency plot (no theoretical distribution shown) */
	lambda = 1.0 / freq->n;
	plotmin = freq->midpt[0] - barwidth;
	plotmax = freq->midpt[K-1] + barwidth;
	fprintf(fp, "set xrange [%.7g:%.7g]\n", plotmin, plotmax);
	maybe_set_yrange(freq, lambda, fp);
	fputs("set nokey\n", fp);
    }

    fprintf(fp, "set xlabel '%s'\n", freq->varname);
    if (dist) {
	fprintf(fp, "set ylabel '%s'\n", I_("Density"));
    } else {
	fprintf(fp, "set ylabel '%s'\n", I_("Relative frequency"));
    }

    if (isnan(lambda)) {
	if (fp != NULL) {
	    fclose(fp);
	}
	return 1;
    }

    if (freq->discrete) {
	fprintf(fp, "set boxwidth 0.75\n");
    }

    /* plot instructions */
    if (use_boxes) {
	if (gnuplot_has_style_fill()) {
	    fputs("set style fill solid 0.5\n", fp);
	}
	strcpy(withstr, "w boxes");
    } else {
	strcpy(withstr, "w impulses");
    }

    if (!dist) {
	fprintf(fp, "plot '-' using 1:2 %s\n", withstr);
    } else if (dist == DIST_NORMAL) {
	print_freq_dist_label(label, dist, freq->xbar, freq->sdx);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' %s , \\\n"
		"1.0/(sqrt(2.0*pi)*sigma)*exp(-.5*((x-mu)/sigma)**2) "
		"title '%s' w lines\n",
		freq->varname, withstr, label);
    } else if (dist == DIST_GAMMA) {
	print_freq_dist_label(label, dist, alpha, beta);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' %s ,\\\n"
		"x**(alpha-1.0)*exp(-x/beta)/(exp(lgamma(alpha))*(beta**alpha)) "
		"title '%s' w lines\n",
		freq->varname, withstr, label); 
    }

    for (i=0; i<K; i++) { 
	fprintf(fp, "%.8g %.8g\n", freq->midpt[i], lambda * freq->f[i]);
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    if (fp != NULL) {
	fclose(fp);
    }

    return gnuplot_make_graph();
}

int plot_fcast_errs (int t1, int t2, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     int time_series)
{
    FILE *fp = NULL;
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int do_errs = (maxerr != NULL);
    int t, n, err;

    /* don't graph empty portion of forecast */
    for (t=t2; t>=t1; t--) {
	if (na(depvar[t]) && na(yhat[t])) {
	    t2--;
	} else {
	    break;
	}
    }

    n = t2 - t1 + 1;

    if (n < 3) {
	/* won't draw a graph for 2 datapoints or less */
	return 1;
    }

    if ((err = gnuplot_init(PLOT_FORECAST, &fp))) {
	return err;
    }    

    /* check that we have any values for the actual var */
    for (t=t1; t<=t2; t++) {
	if (!na(depvar[t])) {
	    depvar_present = 1;
	    break;
	}
    }

    fputs("# forecasts with 95 pc conf. interval\n", fp);

    gretl_minmax(t1, t2, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);

    gretl_pop_c_numeric_locale();


    fputs("set missing \"?\"\n", fp);

    if (time_series) {
	fprintf(fp, "# timeseries %d\n", time_series);
    } else if (n < 33) {
	fputs("set xtics 1\n", fp);
    }

    fputs("set key left top\nplot \\\n", fp);
    if (depvar_present) {
	fprintf(fp, "'-' using 1:2 title '%s' w lines , \\\n",
		varname);
    }
    fprintf(fp, "'-' using 1:2 title '%s' w lines", I_("fitted"));
    if (do_errs) {
	fprintf(fp, " , \\\n'-' using 1:2:3 title '%s' w errorbars\n",
		I_("95 percent confidence interval"));
    } else {
	fputc('\n', fp);
    }

    gretl_push_c_numeric_locale();

    if (depvar_present) {
	for (t=t1; t<=t2; t++) {
	    if (na(depvar[t])) {
		fprintf(fp, "%.8g ?\n", obs[t]);
	    } else {
		fprintf(fp, "%.8g %.8g\n", obs[t], depvar[t]);
	    }
	}
	fputs("e\n", fp);
    }

    for (t=t1; t<=t2; t++) {
	if (na(yhat[t])) {
	    fprintf(fp, "%.8g ?\n", obs[t]);
	} else {
	    fprintf(fp, "%.8g %.8g\n", obs[t], yhat[t]);
	}
    }
    fputs("e\n", fp);

    if (do_errs) {
	for (t=t1; t<=t2; t++) {
	    if (na(yhat[t]) || na(maxerr[t])) {
		fprintf(fp, "%.8g ? ?\n", obs[t]);
	    } else {
		fprintf(fp, "%.8g %.8g %.8g\n", obs[t], yhat[t], maxerr[t]);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int garch_resid_plot (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    double *h, *obs;
    double sd2;
    int t, pv;
    int err;

    h = gretl_model_get_data(pmod, "garch_h");
    if (h == NULL) {
	return E_DATA;
    }

    if ((err = gnuplot_init(PLOT_GARCH, &fp))) {
	return err;
    }

    pv = plotvar(pZ, pdinfo);
    if (pv > 0) {
	obs = (*pZ)[pv];
    } else {
	fclose(fp);
	return E_ALLOC;
    }

    fprintf(fp, "set key left top\n"
	    "plot \\\n'-' using 1:2 title '%s' w lines , \\\n"
	    "'-' using 1:2 title '%s' w lines lt 2, \\\n" 
	    "'-' using 1:2 notitle w lines lt 2\n", 
	    I_("residual"), I_("+- sqrt(h(t))"));

    gretl_push_c_numeric_locale();

    for (t=pmod->t1; t<=pmod->t2; t++) {
	fprintf(fp, "%.8g %.8g\n", obs[t], pmod->uhat[t]);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = -sqrt(h[t]);
	fprintf(fp, "%.8g %.8g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = sqrt(h[t]);
	fprintf(fp, "%.8g %.8g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

void free_plotspec (GPT_SPEC *spec)
{
    int i;

    if (spec == NULL) return;

    if (spec->lines != NULL) free(spec->lines);
    if (spec->data != NULL) free(spec->data);
    if (spec->reglist != NULL) free(spec->reglist);

    for (i=0; i<4; i++) {
	if (spec->literal[i] != NULL) {
	    free(spec->literal[i]);
	}
    }

    if (spec->markers != NULL) {
	for (i=0; i<spec->n_markers; i++) {
	    free(spec->markers[i]);
	}
	free(spec->markers);
    }

    if (spec->labeled != NULL) {
	free(spec->labeled);
    }

    free(spec);
}

int get_termstr (const GPT_SPEC *spec, char *termstr)
{
    int cmds = 0;

    if (!strcmp(spec->termtype, "postscript color")) {
	strcpy(termstr, "postscript eps color"); 
    } else if (!strcmp(spec->termtype, "postscript")) {
	strcpy(termstr, "postscript eps"); 
    } else if (!strcmp(spec->termtype, "PDF")) {
	strcpy(termstr, "pdf");
    } else if (!strcmp(spec->termtype, "fig")) {
	strcpy(termstr, "fig");
    } else if (!strcmp(spec->termtype, "latex")) {
	strcpy(termstr, "latex");
    } else if (!strcmp(spec->termtype, "png")) { 
	const char *png_str = 
	    get_gretl_png_term_line(spec->code);

	strcpy(termstr, png_str + 9);
    } else if (!strcmp(spec->termtype, "emf color")) {
	const char *emf_str = 
	    get_gretl_emf_term_line(spec->code, 1);

	strcpy(termstr, emf_str + 9);
    } else if (!strcmp(spec->termtype, "plot commands")) { 
	cmds = 1;
    } else {
	strcpy(termstr, spec->termtype);
    }

    return cmds;
}

static char *escape_quotes (const char *s)
{
    if (strchr(s, '"') == NULL) {
	return NULL;
    } else {
	int qcount = 0;
	char *ret, *r;
	const char *p = s;

	while (*p) {
	    if (*p == '"') qcount++;
	    p++;
	}

	ret = malloc(strlen(s) + 1 + qcount);
	if (ret == NULL) {
	    return NULL;
	}

	r = ret;
	while (*s) {
	    if (*s == '"') {
		*r++ = '\\';
		*r++ = '"';
	    } else {
		*r++ = *s;
	    }
	    s++;
	}
	*r = 0;

	return ret;
    }
}

static void 
gp_string (FILE *fp, const char *fmt, const char *s, int png)
{
#ifdef ENABLE_NLS  
    if (png && use_latin_2()) {
	char htmlstr[128];

	sprint_l2_to_html(htmlstr, s, sizeof htmlstr);
	fprintf(fp, fmt, htmlstr);
    } else {
	fprintf(fp, fmt, s); 
    }
#else
    fprintf(fp, fmt, s);
#endif
}

const char *gp_justification_string (int j)
{
    if (j == GP_JUST_LEFT) {
	return "left";
    } else if (j == GP_JUST_CENTER) {
	return "center";
    } else if (j == GP_JUST_RIGHT) {
	return "right";
    } else {
	return "left";
    }
}

static void 
print_plot_labelspec (const GPT_LABEL *lbl, int png, FILE *fp)
{
    char *label = escape_quotes(lbl->text);

    gp_string(fp, "set label \"%s\" ", (label != NULL)? 
	      label : lbl->text, png);

    gretl_push_c_numeric_locale();

    fprintf(fp, "at %g,%g %s%s\n", 
	    lbl->pos[0], lbl->pos[1],
	    gp_justification_string(lbl->just),
	    label_front());

    gretl_pop_c_numeric_locale();

    if (label != NULL) {
	free(label);
    }
}

static int print_data_labels (const GPT_SPEC *spec, FILE *fp)
{
    const double *x, *y;
    double xrange, yrange;
    double yoff;
    int t;

    if (spec->n_lines > 2 || 
	spec->lines[0].ncols != 2 || 
	spec->markers == NULL) {
	return 1;
    }

    x = spec->data;
    y = x + spec->nobs;

    xrange = spec->range[0][1] - spec->range[0][0];
    yrange = spec->range[1][1] - spec->range[1][0];

    if (xrange == 0.0 || yrange == 0.0) {
	double ymin = 1.0e+16, ymax = -1.0e+16;
	double xmin = 1.0e+16, xmax = -1.0e+16;

	for (t=0; t<spec->nobs; t++) {
	    if (!na(x[t]) && !na(y[t])) {
		if (yrange == 0.0) {
		    if (y[t] < ymin) {
			ymin = y[t];
		    } else if (y[t] > ymax) {
			ymax = y[t];
		    }
		}
		if (xrange == 0.0) {
		    if (x[t] < xmin) {
			xmin = x[t];
		    } else if (x[t] > xmax) {
			xmax = x[t];
		    }
		}		
	    }
	}
	if (yrange == 0.0) {
	    yrange = ymax - ymin;
	}
	if (xrange == 0.0) {
	    xrange = xmax - xmin;
	}
    }   

    yoff = 0.03 * yrange;

    fputs("# printing data labels\n", fp);

    for (t=0; t<spec->nobs; t++) {
	double xoff = 0.0;

	if (spec->labeled != NULL && !spec->labeled[t]) {
	    /* printing only specified labels */
	    continue;
	}

	if (!na(x[t]) && !na(y[t])) {
	    if (x[t] > .90 * xrange) {
		xoff = -.02 * xrange;
	    }
	    fprintf(fp, "set label '%s' at %.8g,%.8g\n", spec->markers[t],
		    x[t] + xoff, y[t] + yoff);
	}
    }

    return 0;
}

int print_plotspec_details (const GPT_SPEC *spec, FILE *fp)
{
    int i, k, t;
    int png = get_png_output(spec);
    int started_data_lines = 0;
    int n_lines = spec->n_lines;
    double *x[4];
    int any_y2 = 0;
    int miss = 0;

    if (!string_is_blank(spec->titles[0])) {
	if ((spec->flags & GPTSPEC_OLS_HIDDEN) && 
	    is_auto_ols_string(spec->titles[0])) {
	    n_lines--;
	} else {
	    gp_string(fp, "set title '%s'\n", spec->titles[0], png);
	}
    }

    if (!string_is_blank(spec->titles[1])) {
	gp_string(fp, "set xlabel '%s'\n", spec->titles[1], png);
    }

    if (!string_is_blank(spec->titles[2])) {
	gp_string(fp, "set ylabel '%s'\n", spec->titles[2], png);
    }

    if ((spec->flags & GPTSPEC_Y2AXIS) && !string_is_blank(spec->titles[3])) {
	gp_string(fp, "set y2label '%s'\n", spec->titles[3], png);
    }

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	if (!string_is_blank(spec->labels[i].text)) {
	    print_plot_labelspec(&(spec->labels[i]), png, fp);
	}
    }

    if (spec->xzeroaxis) {
	fputs("set xzeroaxis\n", fp);
    }

    fputs("set missing \"?\"\n", fp);

    if (strcmp(spec->keyspec, "none") == 0) {
	fputs("set nokey\n", fp);
    } else {
	fprintf(fp, "set key %s\n", spec->keyspec);
    }

    k = (spec->flags & GPTSPEC_Y2AXIS)? 3 : 2;

    gretl_push_c_numeric_locale();

    for (i=0; i<k; i++) {
	if (na(spec->range[i][0]) || na(spec->range[i][1]) ||
	    spec->range[i][0] == spec->range[i][1]) {
	    continue;
	}
	fprintf(fp, "set %srange [%.7g:%.7g]\n",
		(i == 0)? "x" : (i == 1)? "y" : "y2",
		spec->range[i][0], spec->range[i][1]);
    }

    gretl_pop_c_numeric_locale();

    /* customized xtics? */
    if (!string_is_blank(spec->xtics)) {
	fprintf(fp, "set xtics %s\n", spec->xtics);
    }
    if (!string_is_blank(spec->mxtics)) {
	fprintf(fp, "set mxtics %s\n", spec->mxtics);
    }

    if (spec->flags & GPTSPEC_Y2AXIS) {
	/* using two y axes */
	fputs("set ytics nomirror\n", fp);
	fputs("set y2tics\n", fp);
    } else if (spec->flags & GPTSPEC_MINIMAL_BORDER) {
	/* suppressing part of border */
	fputs("set border 3\n", fp);
	if (string_is_blank(spec->xtics)) {
	    fputs("set xtics nomirror\n", fp);
	}
	fputs("set ytics nomirror\n", fp);
    }

    /* in case of plots that are editable (see gui client), it is
       important to write out the comment string that identifies the
       sort of graph, so that it will be recognized by type when
       it is redisplayed */

    write_plot_type_string(spec->code, fp);

    for (i=0; i<4; i++) {
	if (spec->literal[i] != NULL && *spec->literal[i] != '\0') {
	    fprintf(fp, "%s\n", spec->literal[i]);
	}
    }

    if (spec->flags & GPTSPEC_AUTO_OLS) {
	fputs(auto_ols_string, fp);
    }

    if ((spec->code == PLOT_FREQ_SIMPLE ||
	 spec->code == PLOT_FREQ_NORMAL ||
	 spec->code == PLOT_FREQ_GAMMA) && gnuplot_has_style_fill()) {
	fputs("set style fill solid 0.5\n", fp);
    }  

    if (spec->flags & GPTSPEC_ALL_MARKERS) {
	print_data_labels(spec, fp);
    }

    fputs("plot \\\n", fp);

    for (i=0; i<n_lines; i++) {
	if ((spec->flags & GPTSPEC_Y2AXIS) && spec->lines[i].yaxis != 1) {
	    any_y2 = 1;
	    break;
	}
    }

    for (i=0; i<n_lines; i++) {

	if (strcmp(spec->lines[i].scale, "NA")) {
	    if (!strcmp(spec->lines[i].scale, "1.0")) {
		fputs("'-' using 1", fp);
		for (k=2; k<=spec->lines[i].ncols; k++) {
		    fprintf(fp, ":%d", k);
		}
		fputc(' ', fp);
	    } else {
		fprintf(fp, "'-' using 1:($2*%s) ", 
			spec->lines[i].scale);
	    }
	} else {
	    fprintf(fp, "%s ", spec->lines[i].formula); 
	}

	if ((spec->flags & GPTSPEC_Y2AXIS) && spec->lines[i].yaxis != 1) {
	    fprintf(fp, "axes x1y%d ", spec->lines[i].yaxis);
	}

	gp_string(fp, "title '%s", spec->lines[i].title, png);

	if (any_y2) {
	    if (spec->lines[i].yaxis == 1) {
		fprintf(fp, " (%s)' ", I_("left"));
	    } else {
		fprintf(fp, " (%s)' ", I_("right"));
	    }
	} else {
	    fputs("' ", fp);
	}

	fprintf(fp, "w %s", spec->lines[i].style);
	if (spec->lines[i].type != 0) {
	    fprintf(fp, " lt %d", spec->lines[i].type);
	}

	if (i == n_lines - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp);
	}
    } 

    miss = 0;

    /* supply the data to gnuplot inline */

    gretl_push_c_numeric_locale();

    for (i=0; i<spec->n_lines; i++) { 
	int j;

	if (spec->lines[i].ncols == 0) {
	    continue;
	}

	if (!started_data_lines) {
	    x[0] = spec->data;
	    x[1] = x[0] + spec->nobs;
	    started_data_lines = 1;
	} 

	x[2] = x[1] + spec->nobs;
	x[3] = x[2] + spec->nobs;

	for (t=0; t<spec->nobs; t++) {
	    if (na(x[0][t])) {
		fputs("? ", fp);
		miss = 1;
	    } else {
		fprintf(fp, "%.8g ", x[0][t]);
	    }

	    for (j=1; j<spec->lines[i].ncols; j++) {
		if (na(x[j][t])) {
		    fputs("? ", fp);
		    miss = 1;
		} else {
		    fprintf(fp, "%.8g ",  x[j][t]);
		}
	    }

	    if (spec->markers != NULL && i == 0) {
		fprintf(fp, " # %s", spec->markers[t]);
	    }

	    fputc('\n', fp);
	}

	fputs("e\n", fp);

	x[1] += (spec->lines[i].ncols - 1) * spec->nobs;
    }

    gretl_pop_c_numeric_locale();

    return miss;
}

/* ship out a plot struct, to gnuplot or file */

int go_gnuplot (GPT_SPEC *spec, char *fname)
{
    FILE *fp = NULL;
    int dump = 0;
    int err = 0, miss;
    char termstr[72];

    dump = get_termstr(spec, termstr);

    if (dump) {  
	/* dump of gnuplot commands to named file */
	if (fname == NULL) {
	    return 1;  /* impossible */
	}
	fp = gretl_fopen(fname, "w");
	if (fp == NULL) {
	    return 1;
	}
    } else {     
	/* output to gnuplot, for screen or other "term" */
	if (spec->fp == NULL) {
	    fp = gretl_fopen(gretl_plotfile(), "w");
	}
	if (fp == NULL) {
	    return 1;
	}

	if (fname != NULL) { 
#ifdef ENABLE_NLS
	    fprint_gnuplot_encoding(termstr, fp);
#endif 
	    /* file, not screen display */
	    fprintf(fp, "set term %s\n", termstr);
	    fprintf(fp, "set output '%s'\n", fname);
	}
    }

    if (strstr(termstr, "png")) {
	set_png_output(spec);
    }

    miss = print_plotspec_details(spec, fp);
    fflush(fp);

    if (dump) {
	/* we're finished */
	fclose(fp);
    }
    
    if (!dump) {
	char plotcmd[MAXLEN];
#ifdef WIN32
	int winshow = 0;

	if (fname == NULL) { 
	    /* sending plot to screen */
	    fputs("pause -1\n", fp);
	    winshow = 1;
	} 
#endif
	fclose(fp);
	spec->fp = NULL;
	sprintf(plotcmd, "\"%s\" \"%s\"", gretl_gnuplot_path(), gretl_plotfile());
#ifdef WIN32
	if (winshow) {
	    err = (WinExec(plotcmd, SW_SHOWNORMAL) < 32);
	} else {
	    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
	}
#else
	if (gretl_spawn(plotcmd)) err = 1;
#endif 
    }

    if (miss) {
	err = 2;
    }

    return err;
}

int 
rmplot (const int *list, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    int (*range_mean_graph) (int, const double **, const DATAINFO *, PRN *);
    void *handle = NULL;
    int err;

    range_mean_graph = get_plugin_function("range_mean_graph", &handle);
    if (range_mean_graph == NULL) {
        return 1;
    }

    err = range_mean_graph(list[1], Z, pdinfo, prn);

    close_plugin(handle);

    if (!err && !gretl_in_batch_mode() && !gretl_looping()) {
        err = gnuplot_make_graph();
    }

    return err;
}

int 
hurstplot (const int *list, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    int (*hurst_exponent) (int, const double **, const DATAINFO *, PRN *);
    void *handle = NULL;
    int err;

    hurst_exponent = get_plugin_function("hurst_exponent", &handle);
    if (hurst_exponent == NULL) {
        return 1;
    }

    err = hurst_exponent(list[1], Z, pdinfo, prn);

    close_plugin(handle);

    if (!err && !gretl_in_batch_mode() && !gretl_looping()) {
        err = gnuplot_make_graph();
    } 

    return err;
}

static void get_x_and_y_sizes (int n, int *x, int *y)
{
    if (n == 2) {
	*x = 2;
	*y = 1;
    } else if (n == 3 || n == 4) {
	*x = 2;
	*y = 2;
    } else if (n == 5 || n == 6) {
	*x = 3;
	*y = 2;
    } else if (n > 6 && n < 10) {
	*x = 3;
	*y = 3;
    } else {
	*x = 0;
	*y = 0;
    }
}

/* panel: plot one variable as a time series, with separate plots for
   each cross-sectional unit 
*/

int 
gretl_panel_ts_plot (const int *list, const double **Z, DATAINFO *pdinfo) 
{
    FILE *fp = NULL;
    int i, j, k;
    int t, t0;
    int xnum, ynum;
    float xfrac, yfrac;
    float xorig = 0.0;
    float yorig;
    int err = 0;

    int T = pdinfo->pd;
    int nunits = pdinfo->n / T;

    get_x_and_y_sizes(nunits, &xnum, &ynum);
    if (xnum == 0 || ynum == 0) {
	return E_DATA;
    }

    err = gnuplot_init(PLOT_PANEL, &fp);
    if (err) {
	return err;
    }

    xfrac = 1.0 / xnum;
    yfrac = 1.0 / ynum;

    fputs("set key top left\n", fp);
    fputs("set multiplot\n", fp);
    fprintf(fp, "set xlabel '%s'\n", _("time"));
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    if (yfrac > 1.4 * xfrac) {
	yfrac = 1.4 * xfrac;
    }
    fprintf(fp, "set size %g,%g\n", xfrac, yfrac);

    k = 0;
    t0 = 0;

    for (i=0; i<xnum; i++) {

	yorig = 1.0 - yfrac;

	for (j=0; j<ynum; j++) {
	    int vj = list[1];

	    if (k == nunits) {
		break;
	    }

	    fprintf(fp, "set origin %g,%g\n", xorig, yorig);
	    fprintf(fp, "set title '%s (%d)'\n", pdinfo->varname[vj], k+1);
	    fputs("plot \\\n'-' using 1:2 notitle w lines\n", fp);

	    for (t=0; t<T; t++) {
		fprintf(fp, "%d %.8g\n", t+1, Z[vj][t+t0]);
	    }
	    fputs("e\n", fp);

	    k++;
	    if (k == nunits) {
		break;
	    }

	    t0 += T;
	    yorig -= yfrac;
	}

	if (k == nunits) {
	    break;
	}	

	xorig += xfrac;
    }

    fputs("unset multiplot\n", fp);
    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int 
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, int periods,
				 const double **Z,
				 const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    int confint = 0;
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int t, err;

    resp = gretl_VAR_get_impulse_response(var, targ, shock, periods,
					  Z, pdinfo);
    if (resp == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(resp) > 1) {
	confint = 1;
    }

    err = gnuplot_init((confint)? PLOT_IRFBOOT : PLOT_REGULAR, &fp);
    if (err) {
	gretl_matrix_free(resp);
	return err;
    }

    vtarg = gretl_VAR_get_variable_number(var, targ);
    vshock = gretl_VAR_get_variable_number(var, shock);

    if (!confint) {
	fputs("# impulse response plot\n", fp);
    }

    if (confint) {
	fputs("set key top left\n", fp);
	sprintf(title, I_("response of %s to a shock in %s, "
			  "with bootstrap confidence interval"),
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    } else {
	fputs("set nokey\n", fp);
	sprintf(title, I_("response of %s to a shock in %s"), 
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    }

    fprintf(fp, "set xlabel '%s'\n", _("periods"));
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n", title);

    if (confint) {
	fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines,\\\n", 
		I_("point estimate"));
	fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n",
		I_("0.025 and 0.975 quantiles"));
    } else {
	fputs("plot \\\n'-' using 1:2 w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    for (t=0; t<periods; t++) {
	fprintf(fp, "%d %.8g\n", t+1, gretl_matrix_get(resp, t, 0));
    }
    fputs("e\n", fp);

    if (confint) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.8g %.8g %.8g\n", t+1, 
		    gretl_matrix_get(resp, t, 0),
		    gretl_matrix_get(resp, t, 1),
		    gretl_matrix_get(resp, t, 2));
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);
    gretl_matrix_free(resp);

    return gnuplot_make_graph();
}

int 
gretl_VAR_plot_multiple_irf (GRETL_VAR *var, int periods,
			     const double **Z,
			     const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    int confint = 0;
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int t, err, i, j;

    int n = var->neqns;
    float plot_fraction = 1.0 / n;
    float xorig = 0.0;
    float yorig;

    gp_small_font_size = (n == 4)? 6 : 0;

    resp = gretl_VAR_get_impulse_response(var, 1, 1, periods, Z, pdinfo);
    if (resp == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(resp) > 1) {
	confint = 1;
    }

    err = gnuplot_init(PLOT_MULTI_IRF, &fp);
    if (err) {
	gretl_matrix_free(resp);
	return err;
    }

    if (!confint) {
	fputs("set nokey\n", fp);
    } else {
	fputs("set key top left\n", fp);
    }
    fputs("set multiplot\n", fp);
    fprintf(fp, "set xlabel '%s'\n", _("periods"));
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, "set size %g,%g\n", plot_fraction, plot_fraction);

    for (i=0; i<n; i++) {

	yorig = 1.0 - plot_fraction;
	vtarg = gretl_VAR_get_variable_number(var, i);

	for (j=0; j<n; j++) {

	    fprintf(fp, "set origin %g,%g\n", xorig, yorig);
	    resp = gretl_VAR_get_impulse_response(var, i, j, periods, Z, pdinfo);
	    if (resp == NULL) {
		return E_ALLOC;
	    }

	    vshock = gretl_VAR_get_variable_number(var, j);
	    sprintf(title, "%s -> %s", pdinfo->varname[vshock], pdinfo->varname[vtarg]);
	    fprintf(fp, "set title '%s'\n", title);

	    if (confint) {
		fputs("plot \\\n'-' using 1:2 notitle w lines,\\\n", fp); 
		fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
	    } else {
		fputs("plot \\\n'-' using 1:2 w lines\n", fp);
	    }

	    for (t=0; t<periods; t++) {
		fprintf(fp, "%d %.8g\n", t+1, gretl_matrix_get(resp, t, 0));
	    }
	    fputs("e\n", fp);

	    if (confint) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.8g %.8g %.8g\n", t+1, 
			    gretl_matrix_get(resp, t, 0),
			    gretl_matrix_get(resp, t, 1),
			    gretl_matrix_get(resp, t, 2));
		}
		fputs("e\n", fp);
	    }
	    
	    yorig -= plot_fraction;
	}

	xorig += plot_fraction;
    }

    fputs("unset multiplot\n", fp);
    gretl_pop_c_numeric_locale();

    fclose(fp);
    gretl_matrix_free(resp);

    return gnuplot_make_graph();
}

int gretl_VAR_residual_plot (const GRETL_VAR *var, 
			     double ***pZ, DATAINFO *pdinfo)
{
    const gretl_matrix *E;
    FILE *fp = NULL;
    int nvars, nobs, pv;
    int i, v, t, t1, err;

    E = gretl_VAR_get_residual_matrix(var);
    if (E == NULL) {
	return E_DATA;
    }

    t1 = gretl_VAR_get_t1(var);

    err = gnuplot_init(PLOT_REGULAR, &fp);
    if (err) {
	return err;
    }

    nvars = gretl_matrix_cols(E);
    nobs = gretl_matrix_rows(E);

    fputs("# VAR residual plot\n", fp);
    fputs("set key top left\n", fp);
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n", I_("VAR residuals"));

    fputs("plot \\\n", fp);
    for (i=0; i<nvars; i++) {
	v = gretl_VAR_get_variable_number(var, i);
	fprintf(fp, "'-' using 1:2 title '%s' w lines", pdinfo->varname[v]);
	if (i == nvars - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(",\\\n", fp); 
	}
    }

    pv = plotvar(pZ, pdinfo);
	
    gretl_push_c_numeric_locale();

    for (i=0; i<nvars; i++) {
	for (t=0; t<nobs; t++) {
	    double eti = gretl_matrix_get(E, t, i);

	    if (pv > 0) {
		fprintf(fp, "%g %.8g\n", (*pZ)[pv][t+t1], eti);
	    } else {
		fprintf(fp, "%d %.8g\n", t+1, eti);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int gretl_VAR_residual_mplot (const GRETL_VAR *var, 
			      double ***pZ, DATAINFO *pdinfo)
{
    const gretl_matrix *E;
    FILE *fp = NULL;
    double startdate;
    double xmin, xmax, xrange;
    int nvars, nobs;
    int i, v, t, t1;
    int jump, timevar;
    int err = 0;

    E = gretl_VAR_get_residual_matrix(var);
    if (E == NULL) {
	return E_DATA;
    }

    nvars = gretl_matrix_cols(E);
    if (nvars > 6) {
	return 1;
    }

    timevar = plotvar(pZ, pdinfo);
    if (timevar < 0) {
	return E_ALLOC;
    }

    nobs = gretl_matrix_rows(E);
    t1 = gretl_VAR_get_t1(var);

    err = gnuplot_init(PLOT_MULTI_SCATTER, &fp);
    if (err) {
	return err;
    }

    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    startdate = (*pZ)[timevar][t1];
    jump = nobs / (2 * pdinfo->pd);
    fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);

    gretl_minmax(t1, t1 + nobs - 1, (*pZ)[timevar], &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;
    fprintf(fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);	

    for (i=0; i<nvars; i++) { 

	if (nvars <= 4) {
	    fputs("set size 0.45,0.5\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.5,0.5\n", fp);
	    else if (i == 2) fputs("0.0,0.0\n", fp);
	    else if (i == 3) fputs("0.5,0.0\n", fp);
	} else {
	    fputs("set size 0.31,0.45\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.32,0.5\n", fp);
	    else if (i == 2) fputs("0.64,0.5\n", fp);
	    else if (i == 3) fputs("0.0,0.0\n", fp);
	    else if (i == 4) fputs("0.32,0.0\n", fp);
	    else if (i == 5) fputs("0.64,0.0\n", fp);
	}

	fputs("set noxlabel\n", fp);
	fputs("set noylabel\n", fp);
	v = gretl_VAR_get_variable_number(var, i);
	fprintf(fp, "set title '%s'\n", pdinfo->varname[v]);

	fputs("plot '-' using 1:2 with lines\n", fp);

	for (t=0; t<nobs; t++) {
	    double xx;

	    xx = (*pZ)[timevar][t+t1];
	    fprintf(fp, "%.8g\t", xx);
	    xx = gretl_matrix_get(E, t, i);
	    if (na(xx)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.8g\n", xx);
	    }
	}

	fputs("e\n", fp);
    } 

    gretl_pop_c_numeric_locale();
    fputs("set nomultiplot\n", fp);
    fclose(fp);

    return gnuplot_make_graph();
}

int gretl_VAR_roots_plot (GRETL_VAR *var)
{
    const gretl_matrix *lam;
    FILE *fp = NULL;
    double x, y;
    double px, py;
    int i, n, err;

    lam = gretl_VAR_get_roots(var);
    if (lam == NULL) {
	return E_ALLOC;
    }

    err = gnuplot_init(PLOT_VAR_ROOTS, &fp);
    if (err) {
	return err;
    }

    n = gretl_matrix_rows(lam);

    fprintf(fp, "set title '%s'\n", 
	    I_("VAR inverse roots in relation to the unit circle"));
    fputs("# literal lines = 8\n", fp);
    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1 w lines , \\\n"
	  "'-' w points pt 7\n", fp);

    gretl_push_c_numeric_locale();
    
    for (i=0; i<n; i++) {
        x = gretl_matrix_get(lam, i, 0);
        y = gretl_matrix_get(lam, i, 1);
	/* in polar form */
	px = atan2(y, x);
	py = sqrt(x * x + y * y);
	fprintf(fp, "%.8f %.8f # %.4f,%.4f\n", px, py, x, y);
    }

    gretl_pop_c_numeric_locale();

    fputs("e\n", fp);
    fclose(fp);

    return gnuplot_make_graph();
}

int confidence_ellipse_plot (gretl_matrix *V, double *b, double t, double c,
			     const char *iname, const char *jname)
{
    FILE *fp = NULL;
    double maxerr[2];
    double xcoeff[2];
    double ycoeff[2];
    double *e = NULL;
    int err = 0;

    maxerr[0] = t * sqrt(gretl_matrix_get(V, 0, 0));
    maxerr[1] = t * sqrt(gretl_matrix_get(V, 1, 1));

    err = gretl_invert_symmetric_matrix(V);
    if (err) {
	return err;
    }

    e = gretl_symmetric_matrix_eigenvals(V, 1, &err);
    if (err) {
	return err;
    }

    e[0] = sqrt(1.0 / e[0] * c);
    e[1] = sqrt(1.0 / e[1] * c);

    xcoeff[0] = e[0] * gretl_matrix_get(V, 0, 0);
    xcoeff[1] = e[1] * gretl_matrix_get(V, 0, 1);

    ycoeff[0] = e[0] * gretl_matrix_get(V, 1, 0);
    ycoeff[1] = e[1] * gretl_matrix_get(V, 1, 1);

    free(e);

    err = gnuplot_init(PLOT_ELLIPSE, &fp);
    if (err) {
	return err;
    }

    fprintf(fp, "set title '%s'\n",
	    /* xgettext:no-c-format */
	    I_("95% confidence ellipse and 95% marginal intervals"));
    fputs("# literal lines = 9\n", fp);
    fputs("set parametric\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);

    fprintf(fp, "set xlabel '%s'\n", iname);
    fprintf(fp, "set ylabel '%s'\n", jname);
    fprintf(fp, "set label '%.3g, %.3g' at ", b[0], b[1]);

    gretl_push_c_numeric_locale();

    fprintf(fp, "%g,%g point lt 2 pt 1 offset 3,3\n", b[0], b[1]);

    fprintf(fp, "x(t) = %g*cos(t)%+g*sin(t)%+g\n", xcoeff[0], xcoeff[1], b[0]);
    fprintf(fp, "y(t) = %g*cos(t)%+g*sin(t)%+g\n", ycoeff[0], ycoeff[1], b[1]);

    fputs("plot x(t), y(t) title '', \\\n", fp);
    fprintf(fp, "%g, y(t) title '' w lines lt 2, \\\n", b[0] - maxerr[0]);
    fprintf(fp, "%g, y(t) title '' w lines lt 2, \\\n", b[0] + maxerr[0]);
    fprintf(fp, "x(t), %g title '' w lines lt 2, \\\n", b[1] - maxerr[1]);
    fprintf(fp, "x(t), %g title '' w lines lt 2\n", b[1] + maxerr[1]);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int is_auto_ols_string (const char *s)
{
    if (strstr(s, "automatic OLS")) return 1;
    if (strstr(s, I_("with least squares fit"))) return 1;
    return 0;
}

static char gnuplot_pallette[4][8] = {
    "xff0000", 
    "x0000ff", 
    "x00cc00",  /* full-intensity green is not very legible */
    "xaabbcc"
};

const char *get_gnuplot_pallette (int i, PlotType ptype)
{
#if GP_DEBUG
    fprintf(stderr, "get_gnuplot_pallette: i=%d, ptype=%d\n",
	    i, ptype);
#endif
    if (i == 0 && (ptype == PLOT_FREQ_SIMPLE ||
		   ptype == PLOT_FREQ_NORMAL || 
		   ptype == PLOT_FREQ_GAMMA)) {
	return gnuplot_pallette[3];
    } else if (i >= 0 && i < 3) {
	return gnuplot_pallette[i];
    } else {
	return "";
    }
}

static int colstr_is_valid (const char *colstr)
{
    int i;
    const char *ok = "0123456789abcdef";

    if (*colstr != 'x') {
	return 0;
    }
    if (strlen(colstr) != 7) {
	return 0;
    }
    for (i=1; i<7; i++) {
	if (strchr(ok, colstr[i]) == NULL) {
	    return 0;
	}
    }

    return 1;
}

void set_gnuplot_pallette (int i, const char *colstr)
{
    if (i >= 0 && i <= 3 && colstr_is_valid(colstr)) {
	strcpy(gnuplot_pallette[i], colstr);
    } else {
	fprintf(stderr, "Invalid color spec, '%s'\n", colstr);
    }
}

#ifdef ENABLE_NLS

void pprint_gnuplot_encoding (const char *termstr, PRN *prn)
{
    if (strstr(termstr, "postscript")) {
	const char *enc = get_gnuplot_charset();

	if (enc != NULL) {
	    pprintf(prn, "set encoding %s\n", enc);
	}
    }
}

void fprint_gnuplot_encoding (const char *termstr, FILE *fp)
{
    if (strstr(termstr, "postscript")) {
	const char *enc = get_gnuplot_charset();

	if (enc != NULL) {
	    fprintf(fp, "set encoding %s\n", enc);
	}
    }
}

#endif 

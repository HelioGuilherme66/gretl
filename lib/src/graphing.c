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

/* graphing.c for gretl */

#include "libgretl.h"
#include "var.h"
#include "system.h"
#include "libset.h"
#include "matrix_extra.h"
#include "forecast.h"
#include "plotspec.h"
#include "usermat.h"

#include <unistd.h>
#include <glib.h>

#define GP_DEBUG 0

#ifdef WIN32
# include <windows.h>
#else
# include <signal.h>
# if HAVE_SYS_WAIT_H
#  include <sys/wait.h>
# endif
# ifndef WEXITSTATUS
#  define WEXITSTATUS(stat_val) ((unsigned)(stat_val) >> 8)
# endif
# ifndef WIFEXITED
#  define WIFEXITED(stat_val) (((stat_val) & 255) == 0)
# endif
#endif /* ! _WIN32 */

static char gnuplot_path[MAXLEN];
static int gp_small_font_size;

typedef struct gnuplot_info_ gnuplot_info;

struct gnuplot_info_ {
    GptFlags flags;
    FitType fit;
    int *list;
    int t1;
    int t2;
    double xrange;
    char xtics[64];
    char fmt[16];
    FILE *fp;
    const char *yformula;
    const double *x;
    double *yvar1;
    double *yvar2;
};

#define MAX_LETTERBOX_LINES 8

#define ts_plot(g)      ((g)->flags & GPT_TS)
#define use_impulses(g) ((g)->flags & GPT_IMPULSES)
#define use_lines(g)    ((g)->flags & GPT_LINES)

#if GP_DEBUG
static void print_gnuplot_flags (int flags, int revised);
#endif

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
    { PLOT_H_TEST,         "sampling distribution" },
    { PLOT_PROB_DIST,      "probability distribution" },
    { PLOT_TRI_GRAPH,      "TRAMO / X12A tri-graph" },
    { PLOT_VAR_ROOTS,      "VAR inverse roots plot" },
    { PLOT_ELLIPSE,        "confidence ellipse plot" },
    { PLOT_MULTI_IRF,      "multiple impulse responses" },
    { PLOT_PANEL,          "multiple panel plots" },
    { PLOT_BI_GRAPH,       "double time-series plot" },
    { PLOT_MANY_TS,        "multiple timeseries" },
    { PLOT_RQ_TAU,         "tau sequence plot" },
    { PLOT_BOXPLOTS,       "boxplots" },
    { PLOT_CURVE,          "curve" },
    { PLOT_QQ,             "QQ plot" },
    { PLOT_USER,           "user-defined plot" },
    { PLOT_TYPE_MAX,       NULL }
};

static void graph_list_adjust_sample (int *list, 
				      gnuplot_info *ginfo,
				      const double **Z);
static void clear_gpinfo (gnuplot_info *gi);
static void make_time_tics (gnuplot_info *gi,
			    const DATAINFO *pdinfo,
			    int many, char *xlabel,
			    PRN *prn);
    
#ifndef WIN32

#define SPAWN_DEBUG 0

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
    int child_pid = 0, sinp = 0, serr = 0;
    GError *error = NULL;
    gchar *argv[] = {
	NULL,
	NULL
    };

    if (*gnuplot_path == 0) {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    argv[0] = gnuplot_path;

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
	char errbuf[128];
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
# if SPAWN_DEBUG
		fprintf(stderr, "%s\n", errbuf);
#endif
		if (strstr(cmd, "font") != NULL) {
		    ret = 1;
		}
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

#endif /* !WIN32 */

#define gp_interactive(f) (!(f & GPT_BATCH))

static GptFlags get_gp_flags (gretlopt opt, int k, FitType *f)
{
    GptFlags flags = 0;

    if (gretl_in_batch_mode() || (opt & OPT_U)) {
	flags |= GPT_BATCH;
    }

    if (opt & OPT_G) {
	flags |= GPT_GUI;
    } 

    if (opt & OPT_R) {
	flags |= GPT_RESIDS;
    } else if (opt & OPT_F) {
	flags |= GPT_FA;
    }

    if (opt & OPT_M) {
	flags |= GPT_IMPULSES;
    } else if (opt & OPT_O) {
	flags |= GPT_LINES;
    }

    if (opt & OPT_Z) {
	flags |= GPT_DUMMY;
    } else if (opt & OPT_C) {
	flags |= GPT_XYZ;
    } else {
	if (opt & OPT_S) {
	    flags |= GPT_FIT_OMIT;
	}
	if (opt & OPT_T) {
	    flags |= GPT_IDX;
	}
    }

    if ((k == 2 || (k == 1 && (flags & GPT_IDX))) && !(opt & OPT_S)) {
	/* OPT_S suppresses auto-fit */
	if (opt & OPT_I) {
	    *f = PLOT_FIT_INVERSE;
	} else if (opt & OPT_Q) {
	    *f = PLOT_FIT_QUADRATIC;
	} else if (opt & OPT_L) {
	    *f = PLOT_FIT_LOESS;
	} else if (opt & OPT_N) {
	    *f = PLOT_FIT_OLS;
	}
    }

#if GP_DEBUG
    print_gnuplot_flags(flags, 0);
#endif

    return flags;
}

static void printvars (FILE *fp, int t, const int *list, const double **Z,
		       const double *x, const char *label, double offset)
{
    double xt;
    int i;

    if (x != NULL) {
	xt = x[t] + offset;
	fprintf(fp, "%.10g ", xt);
    }

    for (i=1; i<=list[0]; i++) {
	xt = Z[list[i]][t];
	if (na(xt)) {
	    fputs("? ", fp);
	} else {
	    if (x == NULL && i == 1) { 
		/* the x variable */
		xt += offset;
	    }
	    fprintf(fp, "%.10g ", xt);
	}
    }

    if (label != NULL) {
	fprintf(fp, "# %s", label);
    }

    fputc('\n', fp);
}

static int factorized_vars (gnuplot_info *gi, const double **Z)
{
    int ynum, dum;
    int fn, t, i;

    if (gi->list[0] != 3 || !gretl_isdummy(gi->t1, gi->t2, Z[gi->list[3]])) {
	gretl_errmsg_set(_("You must supply three variables, the last of "
			   "which is a dummy variable\n(with values 1 or 0)\n"));
	return E_DATA;
    }

    ynum = gi->list[1];
    dum = gi->list[3];

    fn = gi->t2 - gi->t1 + 1;

    gi->yvar1 = malloc(fn * sizeof *gi->yvar1);
    if (gi->yvar1 == NULL) {
	return E_ALLOC;
    }

    gi->yvar2 = malloc(fn * sizeof *gi->yvar2);
    if (gi->yvar2 == NULL) {
	free(gi->yvar1);
	return E_ALLOC;
    }

    i = 0;
    for (t=gi->t1; t<=gi->t2; t++) {
	if (na(Z[ynum][t])) {
	    gi->yvar1[i] = NADBL;
	    gi->yvar2[i] = NADBL;
	} else {
	    if (Z[dum][t] == 1.0) {
		gi->yvar1[i] = Z[ynum][t];
		gi->yvar2[i] = NADBL;
	    } else {
		gi->yvar1[i] = NADBL;
		gi->yvar2[i] = Z[ynum][t];
	    }
	}
	i++;
    }

    return 0;
}

#ifdef WIN32

int gnuplot_has_ttf (int reset)
{
    /* we know the gnuplot supplied with gretl for win32
       does TrueType fonts */
    return 1;
}

int gnuplot_pdf_terminal (void)
{
    return GP_PDF_CAIRO;
}

int gnuplot_png_terminal (void)
{
    return GP_PNG_CAIRO;
}
   
int gnuplot_has_latin5 (void)
{
    /* ... and that it supports ISO-8859-9 */
    return 1;
}

int gnuplot_has_cp1254 (void)
{
    /* ... and that it doesn't support CP1254 */
    return 0;
}

int gnuplot_has_bbox (void)
{
    /* ... and that it supports bounding box info */
    return 1;
}

int gnuplot_has_utf8 (void)
{
    /* ... and that it supports "set encoding utf8" */
    return 1;
}

#else /* !WIN32 */

int gnuplot_has_ttf (int reset)
{
    static int err = -1; 

    if (err == -1 || reset) {
	/* if we have cairo we know we (should be!) OK */
        err = gnuplot_test_command("set term pngcairo");
	if (err) {
	    /* otherwise (libgd) try some plausible ttf fonts */
	    err = gnuplot_test_command("set term png font Vera 8");
	    if (err) {
		err = gnuplot_test_command("set term png font luxisr 8");
	    }
	    if (err) {
		err = gnuplot_test_command("set term png font arial 8");
	    }
	}
    }

    return !err;
}

int gnuplot_has_latin5 (void)
{
    static int err = -1; 

    /* not OK in gnuplot 4.2.0 */

    if (err == -1) {
	err = gnuplot_test_command("set encoding iso_8859_9");
    }

    return !err;
}

int gnuplot_has_cp1254 (void)
{
    static int err = -1; 

    /* not OK in gnuplot 4.2.0 */

    if (err == -1) {
	err = gnuplot_test_command("set encoding cp1254");
    }

    return !err;
}

int gnuplot_pdf_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
	int err = gnuplot_test_command("set term pdfcairo");

	if (!err) {
	    ret = GP_PDF_CAIRO;
	} else {
	    err = gnuplot_test_command("set term pdf");
	    if (!err) {
		ret = GP_PDF_PDFLIB;
	    } else {
		ret = GP_PDF_NONE;
	    }
	}
    }

    return ret;
}

static int gnuplot_has_x11 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term x11");
    }

    return !err;
}

int gnuplot_has_wxt (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term wxt");
    }

    return !err;
}

int gnuplot_png_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
	int err = gnuplot_test_command("set term pngcairo");

	if (!err) {
	    fprintf(stderr, "gnuplot: using pngcairo driver\n");
	    ret = GP_PNG_CAIRO;
	} else {
	    /* try the old-style command: if it fails, we have 
	       the libgd driver, we hope! */
	    err = gnuplot_test_command("set term png color");
	    if (!err) {
		fprintf(stderr, "gnuplot: got old png driver\n");
		ret = GP_PNG_OLD;
	    } else {
		fprintf(stderr, "gnuplot: using libgd png driver\n");
		err = gnuplot_test_command("set term png truecolor");
		ret = (err)? GP_PNG_GD1 : GP_PNG_GD2;
	    }
	}
    }

    return ret;
}

int gnuplot_has_bbox (void)
{
    static int err = -1;

    /* not OK in gnuplot 4.2.0 */

    if (err == -1) {
	err = gnuplot_test_command("set term png ; "
				   "set output '/dev/null' ; "
				   "plot x ; print GPVAL_TERM_XMIN");
    }

    return !err;    
}

int gnuplot_has_utf8 (void)
{
    static int err = -1;

    /* not OK in gnuplot 4.2.0 */

    if (err == -1) {
	err = gnuplot_test_command("set encoding utf8");
    }

    return !err;    
}

#endif /* !WIN32 */

static int gnuplot_png_use_aa = 1;

void gnuplot_png_set_use_aa (int s)
{
    gnuplot_png_use_aa = s;
}

/* apparatus for handling plot colors */

static const gretlRGB default_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 }, /* full-intensity green is not very legible */
    { 0xbf, 0x25, 0xb2 },
    { 0x8f, 0xaa, 0xb3 },
    { 0xff, 0xa5, 0x00 },
    { 0x5f, 0x6b, 0x84 },  /* box fill */
    { 0xdd, 0xdd, 0xdd },  /* shade fill */    
};

static gretlRGB user_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 },
    { 0xbf, 0x25, 0xb2 },
    { 0x8f, 0xaa, 0xb3 },
    { 0xff, 0xa5, 0x00 },
    { 0x5f, 0x6b, 0x84 },
    { 0xdd, 0xdd, 0xdd }    
};

void print_rgb_hash (char *s, const gretlRGB *color)
{
    sprintf(s, "#%02x%02x%02x", color->r, color->g, color->b);
}

void gretl_rgb_get (gretlRGB *color, const char *s)
{
    int n, r, g, b;

    n = sscanf(s, "#%2x%2x%2x", &r, &g, &b);

    if (n == 3) {
	color->r = r;
	color->g = g;
	color->b = b;
    } else {
	color->r = color->g = color->b = 0;
    }
}

void print_palette_string (char *s)
{
    char colstr[8];
    int i;

    *s = '\0';

    for (i=0; i<N_GP_COLORS; i++) {
	sprintf(colstr, "x%02x%02x%02x", user_color[i].r, user_color[i].g, 
		user_color[i].b);
	strcat(s, colstr);
	if (i < N_GP_COLORS - 1) {
	    strcat(s, " ");
	}
    }
}

const gretlRGB *get_graph_color (int i)
{
    return (i >= 0 && i < N_GP_COLORS)? &user_color[i] : NULL;
}

void set_graph_palette (int i, gretlRGB color)
{
    if (i >= 0 && i < N_GP_COLORS) {
	user_color[i] = color;
    } else {
	fprintf(stderr, "Out of bounds color index %d\n", i);
    }
}

void set_graph_palette_from_string (int i, const char *s)
{
    int err = 0;

    if (i >= 0 && i < N_GP_COLORS) {
	unsigned int x[3];

	if (sscanf(s + 1, "%02x%02x%02x", &x[0], &x[1], &x[2]) == 3) {
	    user_color[i].r = x[0];
	    user_color[i].g = x[1];
	    user_color[i].b = x[2];
	} else {
	    err = 1;
	}
    } else {
	err = 1;
    }

    if (err) {
	fprintf(stderr, "Error in set_graph_palette_from_string(%d, '%s')\n", 
		i, s);
    }
}

void graph_palette_reset (int i)
{
    if (i == SHADECOLOR) {
	user_color[SHADECOLOR] = default_color[SHADECOLOR];
    } else if (i == BOXCOLOR) {
	user_color[BOXCOLOR] = default_color[BOXCOLOR];
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    user_color[i] = default_color[i];
	}
    }
}

/* Given a string @s such as "Sans 8" or "Bodoni MT 12", write
   the name part into @name and the point-size part into @psz.
   Return 2 if we got both a name and a size, 1 if we just
   got a name, 0 if we got nothing.
*/

int split_graph_fontspec (const char *s, char *name, int *psz)
{
    int i, k = 0, n = strlen(s);
    int nf = 0;

    for (i=n-1; i>0; i--) {
	if (isdigit(s[i])) k++;
	else break;
    }

    if (k > 0) {
	/* got a size */
	char ptstr[8];

	*ptstr = *name = '\0';
	strncat(ptstr, s + n - k, k);
	*psz = atoi(ptstr);
	strncat(name, s, n - k - 1);
	nf = 2;
    } else if (*s != '\0') {
	nf = 1;
	strcpy(name, s);
    }

    return nf;
}

#define USE_SMALL_FONT(t) (t == PLOT_MULTI_IRF || \
			   t == PLOT_MULTI_SCATTER || \
			   t == PLOT_PANEL)

static void 
write_gnuplot_font_string (char *fstr, PlotType ptype, int pngterm)
{
    const char *grfont = gretl_png_font();

    if (*grfont == '\0') {
	grfont = getenv("GRETL_PNG_GRAPH_FONT");
    }

    if (grfont == NULL || *grfont == '\0') {
	*fstr = '\0';
	return;
    }

    if (pngterm == GP_PNG_CAIRO) {
	char fname[128];
	int nf, fsize = 0;

	nf = split_graph_fontspec(grfont, fname, &fsize);
	if (nf == 2) {
	    if (USE_SMALL_FONT(ptype) && gp_small_font_size > 0) {
		fprintf(stderr, "Doing small font\n");
		sprintf(fstr, " font \"%s,%d\"", fname, gp_small_font_size);
	    } else {
		sprintf(fstr, " font \"%s,%d\"", fname, fsize);
	    }
	} else if (nf == 1) {
	    sprintf(fstr, " font \"%s\"", fname);
	}
    } else {
	int shrink = 0;

	if (USE_SMALL_FONT(ptype) && gp_small_font_size > 0) {
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
}

#ifndef WIN32

static void 
write_old_gnuplot_font_string (char *fstr, PlotType ptype)
{
    if (USE_SMALL_FONT(ptype)) {
	strcpy(fstr, " tiny");
    } else {
	strcpy(fstr, " small");
    }
}

#endif

/* requires gnuplot 4.2 or higher */

void write_plot_line_styles (int ptype, FILE *fp)
{
    char cstr[8];
    int i;
    
    if (frequency_plot_code(ptype)) {
	print_rgb_hash(cstr, &user_color[BOXCOLOR]);
	fprintf(fp, "set style line 1 lc rgb \"%s\"\n", cstr);
	fputs("set style line 2 lc rgb \"#000000\"\n", fp);
    } else if (ptype == PLOT_RQ_TAU) {
	fputs("set style line 1 lc rgb \"#000000\"\n", fp);
	for (i=1; i<BOXCOLOR; i++) {
	    print_rgb_hash(cstr, &user_color[i]);
	    fprintf(fp, "set style line %d lc rgb \"%s\"\n", i+1, cstr);
	}
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    print_rgb_hash(cstr, &user_color[i]);
	    fprintf(fp, "set style line %d lc rgb \"%s\"\n", i+1, cstr);
	}
	print_rgb_hash(cstr, &user_color[SHADECOLOR]);
	fprintf(fp, "set style line %d lc rgb \"%s\"\n", 
		SHADECOLOR + 1, cstr);
    }

    fputs("set style increment user\n", fp);
}

/* Get gnuplot to print the dimensions of a PNG plot, in terms
   of both pixels and data bounds, if gnuplot supports this.
*/

void print_plot_bounding_box_request (FILE *fp)
{
    fprintf(fp, "set print '%sgretltmp.png.bounds'\n", gretl_dotdir());
    fputs("print \"pixel_bounds: \", GPVAL_TERM_XMIN, GPVAL_TERM_XMAX, "
	  "GPVAL_TERM_YMIN, GPVAL_TERM_YMAX\n", fp);
    fputs("print \"data_bounds: \", GPVAL_X_MIN, GPVAL_X_MAX, "
	  "GPVAL_Y_MIN, GPVAL_Y_MAX\n", fp);
}

static void do_plot_bounding_box (void)
{
    FILE *fp = fopen(gretl_plotfile(), "a");

    if (fp != NULL) {
	print_plot_bounding_box_request(fp);
	fclose(fp);
    }
}

/**
 * get_gretl_png_term_line:
 * @ptype: indication of the sort of plot to be made, which
 * may made a difference to the color palette chosen.
 * @flags: plot option flags.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the PNG "terminal".  Checks the environment for setting of 
 * %GRETL_PNG_GRAPH_FONT.  Also appends a color-specification string 
 * if the gnuplot PNG driver supports this.
 *
 * Returns: the terminal string, "set term png ..."
 */

const char *get_gretl_png_term_line (PlotType ptype, GptFlags flags)
{
    static char png_term_line[256];
    char truecolor_string[12] = {0};
    char font_string[128];
    char size_string[16];
    int gpttf, pngterm = 0;

    *font_string = 0;
    *size_string = 0;

    pngterm = gnuplot_png_terminal();

#ifdef WIN32
    gpttf = 1;
#else
    gpttf = gnuplot_has_ttf(0);
#endif

    if (pngterm == GP_PNG_GD2 && gnuplot_png_use_aa) {
	strcpy(truecolor_string, " truecolor");
    }    

    /* plot font setup */
    if (gpttf) {
	write_gnuplot_font_string(font_string, ptype, pngterm);
    } 

#ifndef WIN32
    if (!gpttf) {
	write_old_gnuplot_font_string(font_string, ptype);
    }
#endif

    if (flags & GPT_LETTERBOX) {
	strcpy(size_string, " size 680,400");
    } else if (ptype == PLOT_VAR_ROOTS) {
	strcpy(size_string, " size 480,480");
    } else if (ptype == PLOT_QQ) {
	strcpy(size_string, " size 480,480");
    }

    if (pngterm == GP_PNG_CAIRO) {
	sprintf(png_term_line, "set term pngcairo%s%s",
		font_string, size_string);
	strcat(png_term_line, "\nset encoding utf8"); /* FIXME? */
    } else {
	sprintf(png_term_line, "set term png%s%s%s",
		truecolor_string, font_string, size_string); 
    }

#if GP_DEBUG
    fprintf(stderr, "png term line:\n'%s'\n", png_term_line);
#endif

    return png_term_line;
}

static void png_font_to_emf (const char *pngfont, char *emfline)
{
    char name[128];
    int pt;

    if (split_graph_fontspec(pngfont, name, &pt) == 2) {
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
    static char tline[256];
    const char *grfont = NULL;
    
    strcpy(tline, "set term emf ");

    if (color) {
	strcat(tline, "color ");
    } else {
	strcat(tline, "mono dash ");
    }

    /* font spec */
    grfont = gretl_png_font();
    if (grfont != NULL && *grfont != '\0') {
	png_font_to_emf(grfont, tline);
    }

    return tline;
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
    int i, len, ret = PLOT_REGULAR;

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	len = strlen(ptinfo[i].pstr);
	if (!strncmp(str + 2, ptinfo[i].pstr, len)) {
	    ret = ptinfo[i].ptype;
	    break;
	}
    }

    return ret;
}

int write_plot_type_string (PlotType ptype, FILE *fp)
{
    int i, ret = 0;

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (ptype == ptinfo[i].ptype) {
	    fprintf(fp, "# %s\n", ptinfo[i].pstr);
	    ret = 1;
	    break;
	}
    }

    if (get_local_decpoint() == ',') {
	fputs("set decimalsign ','\n", fp);
    }

    return ret;
}

#define PDF_CAIRO_STRING "set term pdfcairo font \"sans,5\""

static void print_term_string (int tt, FILE *fp)
{
    const char *tstr = NULL;

    if (tt == GP_TERM_EPS) {
	tstr = "set term postscript eps mono";
    } else if (tt == GP_TERM_PDF) {
	if (gnuplot_pdf_terminal() == GP_PDF_CAIRO) {
	    tstr = PDF_CAIRO_STRING;
	} else {
	    tstr = "set term pdf";
	}
    } else if (tt == GP_TERM_PNG) {
	tstr = get_gretl_png_term_line(PLOT_REGULAR, 0);
    } else if (tt == GP_TERM_EMF) {
	tstr = get_gretl_emf_term_line(PLOT_REGULAR, 1);
    } else if (tt == GP_TERM_FIG) {
	tstr = "set term fig";
    } else if (tt == GP_TERM_SVG) {
	tstr = "set term svg";
    }

    if (tstr != NULL) {
	fprintf(fp, "%s\n", tstr);
	if (tt != GP_TERM_EPS) {
	    write_plot_line_styles(PLOT_REGULAR, fp);
	}
    }
}

static int command_index_from_plot_type (PlotType p)
{
    if (p == PLOT_MULTI_SCATTER) {
	return SCATTERS;
    } else if (p == PLOT_BOXPLOTS) {
	return BXPLOT;
    } else if (p == PLOT_FORECAST) {
	return FCAST;
    } else {
	/* all other cases */
	return GNUPLOT;
    }
}

static int gretl_plot_count;
static int this_term_type;

/* recorder for filename given via --output=foo */
static char gnuplot_outname[FILENAME_MAX];

static int set_term_type_from_fname (const char *fname)
{
    if (has_suffix(fname, ".eps")) {
	this_term_type = GP_TERM_EPS;
    } else if (has_suffix(fname, ".ps")) {
	this_term_type = GP_TERM_EPS;
    } else if (has_suffix(fname, ".pdf")) {
	this_term_type = GP_TERM_PDF;
    } else if (has_suffix(fname, ".png")) {
	this_term_type = GP_TERM_PNG;
    } else if (has_suffix(fname, ".fig")) {
	this_term_type = GP_TERM_FIG;
    } else if (has_suffix(fname, ".emf")) {
	this_term_type = GP_TERM_EMF;
    } else if (has_suffix(fname, ".svg")) {
	this_term_type = GP_TERM_SVG;
    } 

    return this_term_type;
}

int specified_gp_output_format (void)
{
    return this_term_type;
}

void reset_plot_count (void)
{
    gretl_plot_count = 0;
}

static int make_temp_plot_name (char *fname)
{
    sprintf(fname, "%sgpttmp.XXXXXX", gretl_dotdir());
    return (mktemp(fname) == NULL)? E_FOPEN : 0;
}

static void print_set_output (const char *path, FILE *fp)
{
#ifdef WIN32
    char *s, buf[FILENAME_MAX];

    if (path == NULL) {
	strcpy(buf, gretl_dotdir());
    } else {
	strcpy(buf, path);
    }

    s = buf;
    while (*s) {
	if (*s == '\\') *s = '/';
	s++;
    }

    if (path == NULL) {
	fprintf(fp, "set output \"%sgretltmp.png\"\n", buf);
    } else {
	fprintf(fp, "set output \"%s\"\n", buf);
    }
#else
    if (path == NULL) {
	fprintf(fp, "set output '%sgretltmp.png'\n", gretl_dotdir());
    } else {
	fprintf(fp, "set output '%s'\n", path);
    }
#endif
}

static FILE *gp_set_up_batch (char *fname, PlotType ptype, 
			      GptFlags flags, int *err)
{
    int ci = command_index_from_plot_type(ptype);
    const char *optname = get_optval_string(ci, OPT_U);
    int fmt = GP_TERM_NONE;
    FILE *fp = NULL;

    if (optname != NULL && *optname != '\0') {
	/* user gave --output=<filename> */
	if (!strcmp(optname, "display")) {
	    /* switch to interactive mode */
	    return NULL;
	}
	fmt = set_term_type_from_fname(optname);
	if (fmt) {
	    /* input needs processing */
	    strcpy(gnuplot_outname, optname);
	    gretl_maybe_prepend_dir(gnuplot_outname);
	    make_temp_plot_name(fname);
	} else {
	    /* just passing commands through */
	    strcpy(fname, optname);
	    gretl_maybe_prepend_dir(fname);
	    this_term_type = GP_TERM_PLT;
	}
    } else {
	/* auto-constructed gnuplot commands filename */
	sprintf(fname, "%sgpttmp%02d.plt", gretl_workdir(), 
		++gretl_plot_count);
	this_term_type = GP_TERM_PLT;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
	if (*gnuplot_outname != '\0') {
	    /* write terminal/output lines */
	    print_term_string(fmt, fp);
	    print_set_output(gnuplot_outname, fp);
	}
    }

    return fp;
}

static FILE *gp_set_up_interactive (char *fname, PlotType ptype, 
				    GptFlags flags, int *err)
{
    int gui = gretl_in_gui_mode();
    FILE *fp = NULL;

    if (gui) {
	/* the filename should be unique */
	make_temp_plot_name(fname);
    } else {
	/* gretlcli: no need for uniqueness */
	sprintf(fname, "%sgpttmp.plt", gretl_dotdir());
    }

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
	if (gui) {
	    /* set up for PNG output */
	    fprintf(fp, "%s\n", get_gretl_png_term_line(ptype, flags));
	    print_set_output(NULL, fp);
	}
	write_plot_type_string(ptype, fp);
	write_plot_line_styles(ptype, fp);
    }

    return fp;
}

#ifndef WIN32

static int gnuplot_too_old (void)
{
    static int gperr = -1;

    if (gperr < 0) {
	gperr = gnuplot_test_command("set style line 2 lc rgb \"#0000ff\"");
    }

    if (gperr) {
	gretl_errmsg_set("Gnuplot is too old: must be >= version 4.2.0");
    }

    return gperr;
}

#endif

/* Open stream into which gnuplot commands will be written.

   When GPT_BATCH is set, we're either just dumping a 
   gnuplot command file for the user to process, or possibly
   generating output such as EPS, PDF, etc., in response to
   the --output=filename option.

   When GPT_BATCH is not set we're handling a graph that
   was set up interactively, either via the gretl GUI or
   in interactive mode in gretlcli.  We're going to display
   this graph, either as PNG in a gretl window or via
   gnuplot itself.

   Depending on the prospective use of the stream, we
   may write some initializations into it, the primary
   case being when we're going to produce PNG output
   for display in the GUI.
*/

static FILE *open_gp_stream (PlotType ptype, GptFlags flags, int *err)
{
    char fname[FILENAME_MAX] = {0};
    int batch = (flags & GPT_BATCH);
    FILE *fp = NULL;

    /* ensure we have 'gnuplot_path' in place (file-scope static var) */
    if (*gnuplot_path == '\0') {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

#ifndef WIN32
    if (gnuplot_too_old()) {
	*err = E_EXTERNAL;
	return NULL;
    }
#endif

    this_term_type = GP_TERM_NONE;
    *gnuplot_outname = '\0';

    if (batch) {
	fp = gp_set_up_batch(fname, ptype, flags, err);
    }

    if (!batch || (fp == NULL && !*err)) {
	fp = gp_set_up_interactive(fname, ptype, flags, err);
    }

    if (fp == NULL && *fname) {
	fprintf(stderr, "open_gp_stream: couldn't write to %s\n", fname);
    }

#if GPDEBUG
    fprintf(stderr, "open_gp_stream: '%s'\n", gretl_plotfile());
#endif

    return fp;
}

/**
 * get_gnuplot_batch_stream:
 * @ptype: indication of the sort of plot to be made.
 * @err: location to receive error code.
 *
 * Returns: writable stream on success, %NULL on failure.
 */

FILE *get_gnuplot_batch_stream (PlotType ptype, int *err)
{
    return open_gp_stream(ptype, GPT_BATCH, err);
}

/**
 * get_plot_input_stream:
 * @ptype: indication of the sort of plot to be made.
 * @err: location to receive error code.
 *
 * If we're in GUI mode: writes a unique temporary filename into
 * the internal variable #gretl_plotfile; opens plotfile for writing;
 * and writes initial lines into the output file to select 
 * the PNG terminal type and direct gnuplot's output to a temporary
 * file in the gretl user directory.  
 *
 * If not in GUI mode, opens the file %gpttmp.plt in the gretl
 * user directory.  
 *
 * This function is not for use in batch mode.
 *
 * Returns: writable stream on success, %NULL on failure.
 */

FILE *get_plot_input_stream (PlotType ptype, int *err)
{
    return open_gp_stream(ptype, 0, err);
}

/**
 * gnuplot_cleanup:
 *
 * Removes any temporary gnuplot input file written in
 * the user's dot directory.
 */

void gnuplot_cleanup (void)
{
    const char *p, *fname = gretl_plotfile();

    p = strstr(fname, "gpttmp");

    if (p != NULL) {
	int pnum;

	if (sscanf(p, "gpttmp%d.plt", &pnum) == 0) {
	    gretl_remove(fname);
	}
    }
}

static int graph_file_written;

int graph_written_to_file (void)
{
    return graph_file_written;
}

static void remove_old_png (char *buf)
{
    sprintf(buf, "%sgretltmp.png", gretl_dotdir());
    remove(buf);
}

/**
 * gnuplot_make_graph:
 *
 * Executes gnuplot to produce a graph: in the gretl GUI
 * in interactive mode this will be a PNG file. In the
 * CLI program in interactive mode there will be a direct
 * call to gnuplot to display the graph. In batch mode
 * the type of file written depends on the options selected
 * by the user.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gnuplot_make_graph (void)
{
    char buf[MAXLEN];
    const char *fname = gretl_plotfile();
    int gui = gretl_in_gui_mode();
    int fmt, err = 0;

    graph_file_written = 0;

    fmt = specified_gp_output_format();

    if (fmt == GP_TERM_PLT) {
	/* no-op: just the plot commands are wanted */
	graph_file_written = 1;
	return 0;
    } else if (fmt == GP_TERM_PDF) {
	/* can we do this? */
	if (gnuplot_pdf_terminal() == GP_PDF_NONE) {
	    gretl_errmsg_set(_("Gnuplot does not support PDF output "
			       "on this system"));
	    return E_EXTERNAL;
	}
    } else if (fmt == GP_TERM_NONE && gui) {
	if (gnuplot_has_bbox()) {
	    do_plot_bounding_box();
	}
	/* ensure we don't get stale output */
	remove_old_png(buf);
    }

#ifdef WIN32
    sprintf(buf, "\"%s\" \"%s\"", gretl_gnuplot_path(), fname);
    err = winfork(buf, NULL, SW_SHOWMINIMIZED, 0);
#else
    if (gui || fmt) {
	sprintf(buf, "%s \"%s\"", gretl_gnuplot_path(), fname);
    } else {
	/* gretlcli, interactive */
	sprintf(buf, "%s -persist \"%s\"", gretl_gnuplot_path(), fname);
    }
    err = gretl_spawn(buf);  
#endif

#if GP_DEBUG
    fprintf(stderr, "gnuplot_make_graph:\n"
	    " command='%s', err = %d\n", buf, err);
#endif

    if (fmt) {
	if (err) {
	    /* leave the bad file for diagnostic purposes */
	    fprintf(stderr, "err = %d: bad file is '%s'\n", err, fname);
	} else {
	    /* remove the temporary input file */
	    remove(fname);
	    set_gretl_plotfile(gnuplot_outname);
	    graph_file_written = 1;
	}	
    }

    return err;
}

enum {
    GTITLE_VLS,
    GTITLE_RESID,
    GTITLE_AF,
    GTITLE_AFV
} graph_titles;

static void make_gtitle (gnuplot_info *gi, int code, 
			 const char *s1, const char *s2)
{
    char depvar[VNAMELEN];
    char title[128];

    switch (code) {
    case GTITLE_VLS:
	if (gi->fit == PLOT_FIT_OLS) {
	    sprintf(title, _("%s versus %s (with least squares fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_INVERSE) {
	    sprintf(title, _("%s versus %s (with inverse fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    sprintf(title, _("%s versus %s (with quadratic fit)"),
		    s1, s2);
	}	    
	break;
    case GTITLE_RESID:
	if (sscanf(s1, "residual for %15s", depvar) == 1) {
	    sprintf(title, _("Regression residuals (= observed - fitted %s)"), 
		    depvar);
	}
	break;
    case GTITLE_AF:
	sprintf(title, _("Actual and fitted %s"), s1);
	break;
    case GTITLE_AFV:
	if (s2 == NULL || (gi->flags & GPT_TS)) {
	    sprintf(title, _("Actual and fitted %s"), s1);
	} else {
	    sprintf(title, _("Actual and fitted %s versus %s"), s1, s2);
	}
	break;
    default:
	*title = '\0';
	break;
    }

    if (*title != '\0') {
	fprintf(gi->fp, "set title \"%s\"\n", title);
    }
}

static void print_axis_label (char axis, const char *s, FILE *fp)
{
    if (strchr(s, '\'')) {
	fprintf(fp, "set %clabel \"%s\"\n", axis, s);
    } else {
	fprintf(fp, "set %clabel '%s'\n", axis, s);
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

static int loess_plot (gnuplot_info *gi, const char *literal,
		       const double **Z, const DATAINFO *pdinfo)
{
    gretl_matrix *y = NULL;
    gretl_matrix *x = NULL;
    gretl_matrix *yh = NULL;
    int xno, yno = gi->list[1];
    const double *yvar = Z[yno];
    const double *xvar;
    FILE *fp = NULL;
    char title[96];
    int t, T, d = 1;
    double q = 0.5;
    int err = 0;

    if (gi->x != NULL) {
	xno = 0;
	xvar = gi->x;
    } else {
	xno = gi->list[2];
	xvar = Z[xno];
    }

    graph_list_adjust_sample(gi->list, gi, Z);
    if (gi->t1 == gi->t2 || gi->list[0] != 2) {
	return GRAPH_NO_DATA;
    }

    fp = open_gp_stream(PLOT_REGULAR, gi->flags, &err);
    if (err) {
	return E_FOPEN;
    } 

    err = gretl_plotfit_matrices(yvar, xvar, PLOT_FIT_LOESS,
				 gi->t1, gi->t2, &y, &x);

    if (!err) {
	err = sort_pairs_by_x(x, y, NULL, NULL); /* markers! */
    }

    if (!err) {
	yh = loess_fit(x, y, d, q, OPT_R, &err);
    }

    if (err) {
	fclose(fp);
	goto bailout;
    }

    if (xno > 0) {
	const char *s1 = var_get_graph_name(pdinfo, yno);
	const char *s2 = var_get_graph_name(pdinfo, xno);

	sprintf(title, _("%s versus %s (with loess fit)"), s1, s2);
	print_keypos_string(GP_KEY_LEFT_TOP, fp);
	fprintf(fp, "set title \"%s\"\n", title);
	print_axis_label('y', s1, fp);
	print_axis_label('x', s2, fp);
    } else {
	print_keypos_string(GP_KEY_LEFT_TOP, fp);
	print_axis_label('y', var_get_graph_name(pdinfo, yno), fp);
    }

    print_auto_fit_string(PLOT_FIT_LOESS, fp);

    if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 title '' w points, \\\n", fp);
    sprintf(title, _("loess fit, d = %d, q = %g"), d, q);
    fprintf(fp, " '-' using 1:2 title \"%s\" w lines\n", title);

    T = gretl_vector_get_length(yh);

    gretl_push_c_numeric_locale();

    for (t=0; t<T; t++) {
	fprintf(fp, "%.10g %.10g\n", x->val[t], y->val[t]);
    }
    fputs("e\n", fp);

    for (t=0; t<T; t++) {
	fprintf(fp, "%.10g %.10g\n", x->val[t], yh->val[t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    err = gnuplot_make_graph();

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(x);
    gretl_matrix_free(yh);
    clear_gpinfo(gi);

    return err;
} 

static int get_fitted_line (gnuplot_info *gi, 
			    const double **Z, const DATAINFO *pdinfo, 
			    char *targ)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *V = NULL;
    const double *yvar, *xvar = NULL;
    double x0, s2, *ps2 = NULL;
    int err = 0;

    if (gi->x != NULL && (pdinfo->pd == 1 || pdinfo->pd == 4 || pdinfo->pd == 12)) {
	x0 = gi->x[gi->t1];
    } else {
	xvar = Z[gi->list[2]];
	x0 = NADBL;
    }

    yvar = Z[gi->list[1]];

    if (gi->fit == PLOT_FIT_NONE) {
	/* Doing first-time automatic OLS: we want to check for
	   statistical significance of the slope coefficient
	   to see if it's worth drawing the fitted line, so
	   we have to allocate the variance matrix; otherwise
	   we only need the coefficients.
	*/
	V = gretl_matrix_alloc(2, 2);
	if (V == NULL) {
	    return E_ALLOC;
	}
	ps2 = &s2;
    }

    err = gretl_plotfit_matrices(yvar, xvar, gi->fit,
				 pdinfo->t1, pdinfo->t2,
				 &y, &X);

    if (!err) {
	int k = (gi->fit == PLOT_FIT_QUADRATIC)? 3 : 2;

	b = gretl_column_vector_alloc(k);
	if (b == NULL) {
	    err = E_ALLOC;
	}
    }
    
    if (!err) {
	err = gretl_matrix_ols(y, X, b, V, NULL, ps2);
    }

    if (!err && gi->fit == PLOT_FIT_NONE) {
	/* the "automatic" case */
	double pv, v = gretl_matrix_get(V, 1, 1);
	int T = gretl_vector_get_length(y);

	pv = student_pvalue_2(T - 2, b->val[1] / sqrt(v));
	/* show the line if the two-tailed p-value for the slope coeff
	   is less than 0.1, otherwise discard it */
	if (pv < 0.10) {
	    gi->fit = PLOT_FIT_OLS;
	}
    }
	    
    if (!err && gi->fit != PLOT_FIT_NONE) {
	char title[MAXTITLE], formula[GP_MAXFORMULA];
	double pd = pdinfo->pd;

	set_plotfit_line(title, formula, gi->fit, b->val, x0, pd);
	sprintf(targ, "%s title \"%s\" w lines\n", formula, title);
	gi->flags |= GPT_AUTO_FIT;
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

/* support the "fit" options for a single time-series plot */

static int time_fit_plot (gnuplot_info *gi, const char *literal,
			  const double **Z, const DATAINFO *pdinfo)
{
    int yno = gi->list[1];
    const double *yvar = Z[yno];
    FILE *fp = NULL;
    char fitline[128] = {0};
    PRN *prn;
    int t, err = 0;

    if (gi->x == NULL) {
	return E_DATA;
    }

    graph_list_adjust_sample(gi->list, gi, Z);
    if (gi->t1 == gi->t2) {
	return GRAPH_NO_DATA;
    }

    err = get_fitted_line(gi, Z, pdinfo, fitline);
    if (err) {
	return err;
    }

    gi->flags |= GPT_LETTERBOX;

    fp = open_gp_stream(PLOT_REGULAR, gi->flags, &err);
    if (err) {
	return err;
    } 

    prn = gretl_print_new_with_stream(fp);

    if (prn != NULL) {
	make_time_tics(gi, pdinfo, 0, NULL, prn);
	gretl_print_detach_stream(prn);
	gretl_print_destroy(prn);
    }

    print_keypos_string(GP_KEY_LEFT_TOP, fp);
    print_axis_label('y', var_get_graph_name(pdinfo, yno), fp);

    print_auto_fit_string(gi->fit, fp);

    if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 title '' w lines, \\\n", fp);
    
    gretl_push_c_numeric_locale();

    fprintf(fp, " %s", fitline);

    for (t=gi->t1; t<=gi->t2; t++) {
	fprintf(fp, "%.10g %.10g\n", gi->x[t], yvar[t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    err = gnuplot_make_graph();

    clear_gpinfo(gi);

    return err;
} 

static int check_tic_labels (double vmin, double vmax,
			     gnuplot_info *gi)
{
    char s1[32], s2[32];
    int d, err = 0;

    for (d=6; d<12; d++) {
	sprintf(s1, "%.*g", d, vmin);
	sprintf(s2, "%.*g", d, vmax);
	if (strcmp(s1, s2)) {
	    break;
	}
    }

    if (d > 6) {
	sprintf(gi->fmt, "%% .%dg", d+1);
	sprintf(gi->xtics, "%.*g %#.6g", d+1, vmin, 
		(vmax - vmin)/ 4.0);
    }

    return err;
}

static void check_y_tics (gnuplot_info *gi, const double **Z,
			  FILE *fp)
{
    double ymin, ymax;

    *gi->fmt = '\0';

    gretl_minmax(gi->t1, gi->t2, Z[gi->list[1]], &ymin, &ymax);
    check_tic_labels(ymin, ymax, gi);

    if (*gi->fmt != '\0') {
	fprintf(fp, "set format y \"%s\"\n", gi->fmt);
    }
}

/* Find the minimum and maximum x-axis values and construct the gnuplot
   x-range.  We have to be a bit careful here to include only values
   that will actually display on the plot, i.e. x-values that are
   accompanied by at least one non-missing y-axis value.  We also have
   to avoid creating an "empty x range" that will choke gnuplot.
*/

static void 
print_x_range_from_list (gnuplot_info *gi, const double **Z, const int *list)
{
    int vx, k, dk = -1;
    const double *x;

    if (gi->flags & GPT_DUMMY) {
	dk = list[0];
	k = dk - 1;
    } else {
	k = list[0];
    }

    vx = list[k];
    x = Z[vx];
    
    if (gretl_isdummy(gi->t1, gi->t2, x)) {
	fputs("set xrange [-1:2]\n", gi->fp);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", gi->fp);
	gi->xrange = 3;
    } else {
	double xmin, xmin0 = NADBL;
	double xmax, xmax0 = NADBL;
	int t, i, vy, obs_ok;

	for (t=gi->t1; t<=gi->t2; t++) {
	    obs_ok = 0;
	    if (!na(x[t]) && (dk < 0 || !na(Z[dk][t]))) {
		for (i=1; i<k; i++) {
		    vy = list[i];
		    if (!na(Z[vy][t])) {
			/* got x obs and at least one y obs */
			obs_ok = 1;
			break;
		    }
		}
	    }
	    if (obs_ok) {
		if (na(xmin0) || x[t] < xmin0) {
		    xmin0 = x[t];
		}
		if (na(xmax0) || x[t] > xmax0) {
		    xmax0 = x[t];
		}
	    }
	}
		    
	gi->xrange = xmax0 - xmin0;

	if (gi->xrange == 0.0) {
	    /* construct a non-empty range */
	    xmin = xmin0 - 0.5;
	    xmax = xmin0 + 0.5;
	} else {
	    xmin = xmin0 - gi->xrange * .025;
	    if (xmin0 >= 0.0 && xmin < 0.0) {
		xmin = 0.0;
	    }
	    xmax = xmax0 + gi->xrange * .025;
	}

	fprintf(gi->fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
	gi->xrange = xmax - xmin;
	check_tic_labels(xmin0, xmax0, gi);
    }
}

static void 
print_x_range (gnuplot_info *gi, const double *x)
{
    if (gretl_isdummy(gi->t1, gi->t2, x)) {
	fputs("set xrange [-1:2]\n", gi->fp);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", gi->fp);
	gi->xrange = 3;
    } else {
	double xmin0, xmin, xmax0, xmax;

	gretl_minmax(gi->t1, gi->t2, x, &xmin0, &xmax0);
	gi->xrange = xmax0 - xmin0;
	xmin = xmin0 - gi->xrange * .025;
	if (xmin0 >= 0.0 && xmin < 0.0) {
	    xmin = 0.0;
	}
	xmax = xmax0 + gi->xrange * .025;
	fprintf(gi->fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
	gi->xrange = xmax - xmin;
    }
}

/* two or more y vars plotted against some x: test to see if we want
   to use two y axes */

static void
check_for_yscale (gnuplot_info *gi, const double **Z, int *oddman)
{
    double ymin[6], ymax[6];
    double ratio;
    int i, j, oddcount;

#if GP_DEBUG
    fprintf(stderr, "gnuplot: doing check_for_yscale\n");
#endif

    /* find minima, maxima of the y-axis vars */
    for (i=1; i<gi->list[0]; i++) {
	gretl_minmax(gi->t1, gi->t2, Z[gi->list[i]], 
		     &ymin[i-1], &ymax[i-1]);
    }

    gi->flags &= ~GPT_Y2AXIS;

    for (i=1; i<gi->list[0]; i++) {
	oddcount = 0;
	for (j=1; j<gi->list[0]; j++) {
	    if (j == i) {
		continue;
	    }
	    ratio = ymax[i-1] / ymax[j-1];
	    if (ratio > 5.0 || ratio < 0.2) {
		gi->flags |= GPT_Y2AXIS;
		oddcount++;
	    }
	}
	if (oddcount == gi->list[0] - 2) {
	    /* series at list position i differs considerably in scale
	       from all the others in the list */
	    *oddman = i;
	    break;
	}
    }

    if (*oddman == 0) {
	gi->flags &= ~GPT_Y2AXIS;
    }
}

static int print_gp_dummy_data (gnuplot_info *gi, 
				const double **Z, 
				const DATAINFO *pdinfo)
{
    double xx, yy;
    int i, s, t;

    for (i=0; i<2; i++) {
	for (t=gi->t1; t<=gi->t2; t++) {
	    s = t - gi->t1;
	    if (gi->x != NULL) {
		xx = gi->x[t];
	    } else {
		xx = Z[gi->list[2]][t];
		if (na(xx)) {
		    continue;
		}
	    }
	    yy = (i > 0)? gi->yvar2[s] : gi->yvar1[s];
	    if (na(yy)) {
		fprintf(gi->fp, "%.10g ?\n", xx);
	    } else {
		fprintf(gi->fp, "%.10g %.10g", xx, yy);
		if (!(gi->flags & GPT_TS)) {
		    if (pdinfo->markers) {
			fprintf(gi->fp, " # %s", pdinfo->S[t]);
		    } else if (dataset_is_time_series(pdinfo)) {
			char obs[OBSLEN];

			ntodate(obs, t, pdinfo);
			fprintf(gi->fp, " # %s", obs);
		    }
		}
		fputc('\n', gi->fp);
	    }
	}
	fputs("e\n", gi->fp);
    }

    return 0;
}

/* for printing panel time-series graph: insert a discontinuity
   between the panel units */

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

/* sanity check for totally empty graph */

static int 
all_graph_data_missing (const int *list, int t, const double **Z)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (!na(Z[list[i]][t])) {
	    return 0;
	}
    }

    return 1;
}

static void print_gp_data (gnuplot_info *gi, const double **Z, 
			   const DATAINFO *pdinfo)
{
    int n = gi->t2 - gi->t1 + 1;
    double offset = 0.0;
    int datlist[3];
    int ynum = 2;
    int nomarkers = 0;
    int i, t;

    /* multi impulse plot? calculate offset for lines */
    if (use_impulses(gi) && gi->list[0] > 2) { 
	offset = 0.10 * gi->xrange / n;
    }

    if (gi->x != NULL) {
	datlist[0] = 1;
	ynum = 1;
    } else {
	datlist[0] = 2;
	datlist[1] = gi->list[gi->list[0]];
    }

    if (use_impulses(gi) || use_lines(gi)) {
	nomarkers = 1;
    }

    /* loop across the variables, printing x then y[i] for each i */

    for (i=1; i<gi->list[0]; i++) {
	double xoff = offset * (i - 1);

	datlist[ynum] = gi->list[i];

	for (t=gi->t1; t<=gi->t2; t++) {
	    const char *label = NULL;
	    char obs[OBSLEN];

	    if (gi->x == NULL && all_graph_data_missing(gi->list, t, Z)) {
		continue;
	    }

	    if (!(gi->flags & GPT_TS) && i == 1) {
		if (pdinfo->markers) {
		    label = pdinfo->S[t];
		} else if (!nomarkers && dataset_is_time_series(pdinfo)) {
		    ntodate(obs, t, pdinfo);
		    label = obs;
		}
	    }

	    if ((gi->flags & GPT_TS) && pdinfo->structure == STACKED_TIME_SERIES) {
		maybe_print_panel_jot(t, pdinfo, gi->fp);
	    }

	    printvars(gi->fp, t, datlist, Z, gi->x, label, xoff);
	}

	fputs("e\n", gi->fp);
    }
}

static int
gpinfo_init (gnuplot_info *gi, gretlopt opt, const int *list, 
	     const char *literal, int t1, int t2)
{
    int l0 = list[0];

    gi->fit = PLOT_FIT_NONE;
    gi->flags = get_gp_flags(opt, l0, &gi->fit);

    if (gi->fit == PLOT_FIT_NONE) {
	gi->flags |= GPT_TS; /* may be renounced later */
    }

    gi->t1 = t1;
    gi->t2 = t2;
    gi->xrange = 0.0;
    gi->xtics[0] = '\0';
    gi->fmt[0] = '\0';
    gi->yformula = NULL;
    gi->fp = NULL;

    gi->x = NULL;
    gi->yvar1 = NULL;
    gi->yvar2 = NULL;
    gi->list = NULL;

    if (l0 < 2 && !(gi->flags & GPT_IDX)) {
	return E_ARGS;
    }

    if ((gi->flags & GPT_DUMMY) && (gi->flags & GPT_IDX)) {
	return E_BADOPT;
    }

    gi->list = gretl_list_copy(list);
    if (gi->list == NULL) {
	return E_ALLOC;
    }

    if ((l0 > 2 || (l0 > 1 && (gi->flags & GPT_IDX))) && 
	 l0 < 7 && !(gi->flags & GPT_RESIDS) && !(gi->flags & GPT_FA)
	&& !(gi->flags & GPT_DUMMY)  & !(opt & OPT_Y)) { /* FIXME GPT_XYZ */
	/* allow probe for using two y axes */
#if GP_DEBUG
	fprintf(stderr, "l0 = %d, setting y2axis probe\n", l0);
#endif
	gi->flags |= GPT_Y2AXIS;
    } 

    if ((gi->flags & GPT_FA) && literal != NULL && 
	!strncmp(literal, "yformula: ", 10)) {
	/* fitted vs actual plot with fitted given by formula */
	gi->yformula = literal + 10;
    }

    if (literal != NULL && strstr(literal, "set style data")) {
	gi->flags |= GPT_DATA_STYLE;
    }

#if GP_DEBUG
    print_gnuplot_flags(gi->flags, 1);
#endif

    return 0;
}

static void clear_gpinfo (gnuplot_info *gi)
{
    free(gi->yvar1);
    free(gi->yvar2);
    free(gi->list);

    if (gi->fp != NULL) {
	fclose(gi->fp);
    }
}

#if GP_DEBUG
static void print_gnuplot_flags (int flags, int revised)
{
    if (revised) {
	fprintf(stderr, "*** gnuplot flags after initial revision:\n");
    } else {
	fprintf(stderr, "*** gnuplot() called with flags:\n");
    }

    if (flags & GPT_IMPULSES) {
	fprintf(stderr, " GPT_IMPULSES\n");
    }
    if (flags & GPT_LINES) {
	fprintf(stderr, " GPT_LINES\n");
    }
    if (flags & GPT_RESIDS) {
	fprintf(stderr, " GPT_RESIDS\n");
    }	
    if (flags & GPT_FA) {
	fprintf(stderr, " GPT_FA\n");
    }
    if (flags & GPT_DUMMY) {
	fprintf(stderr, " GPT_DUMMY\n");
    }
    if (flags & GPT_XYZ) {
	fprintf(stderr, " GPT_XYZ\n");
    }
    if (flags & GPT_BATCH) {
	fprintf(stderr, " GPT_BATCH\n");
    }
    if (flags & GPT_GUI) {
	fprintf(stderr, " GPT_GUI\n");
    }
    if (flags & GPT_FIT_OMIT) {
	fprintf(stderr, " GPT_FIT_OMIT\n");
    }
    if (flags & GPT_DATA_STYLE) {
	fprintf(stderr, " GPT_DATA_STYLE\n");
    }
    if (flags & GPT_IDX) {
	fprintf(stderr, " GPT_IDX\n");
    }
    if (flags & GPT_TS) {
	fprintf(stderr, " GPT_TS\n");
    }
    if (flags & GPT_Y2AXIS) {
	fprintf(stderr, " GPT_Y2AXIS\n");
    }
    if (flags & GPT_AUTO_FIT) {
	fprintf(stderr, " GPT_AUTO_FIT\n");
    }
    if (flags & GPT_FIT_HIDDEN) {
	fprintf(stderr, " GPT_FIT_HIDDEN\n");
    }
}
#endif

static void set_lwstr (const DATAINFO *pdinfo, int v, char *s)
{
    int w = var_get_linewidth(pdinfo, v);

    if (w > 1) {
	sprintf(s, " lw %d", w);
    } else {
	*s = '\0';
    }
}

static void set_withstr (GptFlags flags, char *str)
{
    if (flags & GPT_DATA_STYLE) {
	*str = 0;
    } else if (flags & GPT_LINES) {
	strcpy(str, "w lines");
    } else {
	strcpy(str, "w points");
    }
}

static void graph_list_adjust_sample (int *list, 
				      gnuplot_info *ginfo,
				      const double **Z)
{
    int t1min = ginfo->t1;
    int t2max = ginfo->t2;
    int t_ok;
    int i, t, vi;

    for (t=t1min; t<=t2max; t++) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && !na(Z[vi][t])) {
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
	    vi = list[i];
	    if (vi > 0 && !na(Z[vi][t])) {
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

	    vi = list[i];
	    for (t=t1min; t<=t2max; t++) {
		if (!na(Z[vi][t])) {
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
}

static int maybe_add_plotx (gnuplot_info *gi, 
			    const DATAINFO *pdinfo)
{
    int k = gi->list[0];
    int add0 = 0;

    /* are we really doing a time-series plot? */
    if (k > 1 && !strcmp(pdinfo->varname[gi->list[k]], "time")) {
	; /* yes */
    } else if (gi->flags & GPT_IDX) {
	add0 = 1; /* yes */
    } else {
	/* no: get out */
	gi->flags &= ~GPT_TS;
	return 0;
    }

    gi->x = gretl_plotx(pdinfo);
    if (gi->x == NULL) {
	return E_ALLOC;
    }

    /* a bit ugly, but add a dummy list entry for
       the 'virtual' plot variable */
    if (add0) {
	gretl_list_append_term(&gi->list, 0);
	if (gi->list == NULL) {
	    return E_ALLOC;
	} 
    }

#if GP_DEBUG
    fprintf(stderr, "maybe_add_plotx: gi->x at %p\n", 
	    (void *) gi->x);
    printlist(gi->list, "gi->list");
#endif

    return 0;
}

void gnuplot_missval_string (FILE *fp)
{
    fputs("set datafile missing \"?\"\n", fp);
}

static void graph_month_name (char *mname, int m)
{
    struct tm mt;

    mt.tm_sec = 0;
    mt.tm_min = 0;
    mt.tm_hour = 0;
    mt.tm_mday = 1;
    mt.tm_mon = m - 1;
    mt.tm_year = 100;

    strftime(mname, 7, "%b", &mt);
    mname[4] = '\0';
}

/* for short daily time-series plots: write month names
   into the xtics */

static void make_named_month_tics (const gnuplot_info *gi, double yrs, 
				   PRN *prn)
{
    double t0 = gi->x[gi->t1];
    double t1 = gi->x[gi->t2];
    double x, tw = 1.0/12;
    int i, m, n = 0;
    char mname[8];
    int notfirst = 0;
    int scale = (int) (yrs * 1.5);

    if (scale == 0) {
	scale = 1;
    }

    t0 += (1.0 - (t0 - floor(t0)) * 12.0) / 12.0;
    for (x=t0; x<t1; x+=tw) n++;

    x = (t0 - floor(t0)) * 12;
    m = 1 + ((x - floor(x) > .8)? ceil(x) : floor(x));
    if (m > 12) m -= 12;

    pputs(prn, "# literal lines = 1\n"); 
    pputs(prn, "set xtics ("); 
    x = t0;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	if (m == 1) {
	    if (notfirst) {
		pputs(prn, ", ");
	    }
	    pprintf(prn, "\"%4.0f\" %.10g", x, x);
	    notfirst = 1;
	} else if ((scale == 1) || (m % scale == 1)) {
	    graph_month_name(mname, m);
	    if (notfirst) {
		pputs(prn, ", ");
	    }
	    pprintf(prn, "\"%s\" %.10g", mname, x);
	    notfirst = 1;
	}
	m++;
	x += tw;
	if (m > 12) m -= 12;
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

/* Below: we're making a combined time series plot for panel data.
   That is, time series for unit 1, followed by time series for unit
   2, etc.  We'd like to show tic marks to represent the start of each
   unit's time series, but we have to watch out for the case where
   there are "too many" units -- we don't want a dense fudge of marks
   on the x-axis.  In that case we put a tic mark only for every k'th
   unit.
*/

static void make_panel_unit_tics (const DATAINFO *pdinfo, 
				  gnuplot_info *gi, 
				  PRN *prn)
{
    int maxtics, ticskip;
    double ntics;
    int printed;
    int u, t, n;

    pputs(prn, "# literal lines = 1\n"); 
    pputs(prn, "set xtics ("); 

    gretl_push_c_numeric_locale();

    /* how many panel units are included in the plot? */
    maxtics = pdinfo->paninfo->unit[gi->t2] - 
	pdinfo->paninfo->unit[gi->t1] + 1;

    ntics = maxtics;
    while (ntics > 20) {
	ntics /= 1.5;
    }

    ticskip = maxtics / ceil(ntics);

    if (ticskip == 1 && ntics < maxtics) {
	/* otherwise we'll get an incomplete scale */
	ntics = maxtics;
    }

    n = printed = 0;
    for (t=gi->t1; t<=gi->t2 && printed<ntics; t++) {
	u = pdinfo->paninfo->unit[t];
	if (t == gi->t1 || u != pdinfo->paninfo->unit[t-1]) {
	    if (n % ticskip == 0) {
		pprintf(prn, "\"%d\" %.10g", u + 1, gi->x[t]);
		if (++printed < ntics) {
		    pputs(prn, ", ");
		}
	    }
	    n++;
	}
    } 

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

/* panel plot includes two or more time series, end-to-end */

static int panel_plot (const DATAINFO *pdinfo, int t1, int t2)
{
    if (dataset_is_panel(pdinfo)) {
	if (pdinfo->paninfo->unit[t2] != pdinfo->paninfo->unit[t1]) {
	    return 1;
	}
    }

    return 0;
}

static void make_calendar_tics (const DATAINFO *pdinfo,
				const gnuplot_info *gi,
				PRN *prn)
{
    int T = gi->t2 - gi->t1 + 1;
    double yrs;

    if (pdinfo->pd == 52) {
	yrs = T / 52.0;
    } else {
	yrs = T / (pdinfo->pd * 52.0);
    }

    if (yrs <= 3) {
	make_named_month_tics(gi, yrs, prn);
    } else if (yrs < 6) {
	/* don't show ugly "fractions of years" */
	pputs(prn, "set xtics 1\n");
	if (yrs < 3) {
	    /* put monthly minor tics */
	    pputs(prn, "set mxtics 12\n");
	} else if (yrs < 5) {
	    /* quarterly minor tics */
	    pputs(prn, "set mxtics 4\n");
	}
    }
}

/* special tics for time series plots */

static void make_time_tics (gnuplot_info *gi,
			    const DATAINFO *pdinfo,
			    int many, char *xlabel,
			    PRN *prn)
{
    if (many) {
	pprintf(prn, "# multiple timeseries %d\n", pdinfo->pd);
    } else {
	pprintf(prn, "# timeseries %d", pdinfo->pd);
	gi->flags |= GPT_LETTERBOX;
	pputs(prn, " (letterbox)\n");
    } 

    if (pdinfo->pd == 4 && (gi->t2 - gi->t1) / 4 < 8) {
	pputs(prn, "set xtics nomirror 0,1\n"); 
	pputs(prn, "set mxtics 4\n");
    } else if (pdinfo->pd == 12 && (gi->t2 - gi->t1) / 12 < 8) {
	pputs(prn, "set xtics nomirror 0,1\n"); 
	pputs(prn, "set mxtics 12\n");
    } else if (dated_daily_data(pdinfo) || dated_weekly_data(pdinfo)) {
	make_calendar_tics(pdinfo, gi, prn);
    } else if (panel_plot(pdinfo, gi->t1, gi->t2)) {
	make_panel_unit_tics(pdinfo, gi, prn);
	if (xlabel != NULL) {
	    strcpy(xlabel, _("time series by group"));
	}
    }
}

/* Respond to use of the option --matrix=<matname> in the gnuplot
   command, or create a plot directly from a matrix and a plot list.
*/

int matrix_plot (gretl_matrix *m, const int *list, const char *literal, 
		 gretlopt opt)
{
    double **Z = NULL;
    DATAINFO *dinfo = NULL;
    int *plotlist = NULL;
    int err = 0;

    if (m == NULL) {
	/* we need to find the matrix by name */
	const char *mname = get_optval_string(GNUPLOT, OPT_X);
    
	if (mname == NULL) {
	    err = E_DATA;
	} else {
	    m = get_matrix_by_name(mname);
	    if (m == NULL) {
		err = E_DATA;
	    }
	}
    }

    if (!err) {
	dinfo = gretl_dataset_from_matrix(m, list, &Z, &err);
    }

    if (err) {
	return err;
    }

    plotlist = gretl_consecutive_list_new(1, dinfo->v - 1);
    if (plotlist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	opt &= ~OPT_X;
	err = gnuplot(plotlist, literal, (const double **) Z, dinfo, opt);
    }

    destroy_dataset(Z, dinfo);   
    free(plotlist);

    return err;
}

#define fit_opts (OPT_I | OPT_L | OPT_Q | OPT_N)

/**
 * gnuplot:
 * @plotlist: list of variables to plot, by ID number.
 * @literal: commands to be passed to gnuplot.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: option flags.
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gnuplot (const int *plotlist, const char *literal,
	     const double **Z, const DATAINFO *pdinfo, 
	     gretlopt opt)
{
    PRN *prn = NULL;
    FILE *fp = NULL;
    int *list = NULL;
    char s1[MAXDISP] = {0};
    char s2[MAXDISP] = {0};
    char xlabel[MAXDISP] = {0};
    char withstr[16] = {0};
    char lwstr[8] = {0};
    char keystr[48] = {0};
    char fit_line[128] = {0};
    int oddman = 0;
    int many = 0;
    gnuplot_info gi;
    int i, err = 0;

    gretl_error_clear();

    err = incompatible_options(opt, fit_opts);
    if (err) {
	return err;
    }

    if ((opt & OPT_T) && (opt & fit_opts)) {
	if (plotlist[0] > 1 || !dataset_is_time_series(pdinfo)) {
	    return E_BADOPT;
	}
    }

    if (opt & OPT_X) {
	return matrix_plot(NULL, plotlist, literal, opt);
    }

#if GP_DEBUG
    printlist(plotlist, "gnuplot: plotlist");
#endif

    err = gpinfo_init(&gi, opt, plotlist, literal, 
		      pdinfo->t1, pdinfo->t2);
    if (err) {
	goto bailout;
    }

    err = maybe_add_plotx(&gi, pdinfo);
    if (err) {
	goto bailout;
    }

    if (gi.fit == PLOT_FIT_LOESS) {
	return loess_plot(&gi, literal, Z, pdinfo);
    }

    if (opt & OPT_T && (opt & fit_opts)) {
	return time_fit_plot(&gi, literal, Z, pdinfo);
    }

    if (gi.list[0] > MAX_LETTERBOX_LINES + 1) {
	many = 1;
    }

    /* convenience pointer */
    list = gi.list;

    /* below: did have "height 1 width 1 box" for win32,
       "width 1 box" otherwise */
    strcpy(keystr, "set key left top\n");

    if (gi.flags & GPT_IMPULSES) {
	strcpy(withstr, "w i");
    }

    /* set x-axis label for non-time series plots */
    if (!(gi.flags & GPT_TS)) {
	int v = (gi.flags & GPT_DUMMY)? list[2] : list[list[0]];

	strcpy(xlabel, var_get_graph_name(pdinfo, v));
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	goto bailout;
    }

    /* adjust sample range, and reject if it's empty */
    graph_list_adjust_sample(list, &gi, Z);
    if (gi.t1 == gi.t2 || list[0] < 2) {
	err = GRAPH_NO_DATA;
	goto bailout;
    }

    /* add a simple regression line if appropriate */
    if (!use_impulses(&gi) && !(gi.flags & GPT_FIT_OMIT) && list[0] == 2 && 
	!(gi.flags & GPT_TS) && !(gi.flags & GPT_RESIDS)) {
	get_fitted_line(&gi, Z, pdinfo, fit_line);
	pprintf(prn, "# X = '%s' (%d)\n", pdinfo->varname[list[2]], list[2]);
	pprintf(prn, "# Y = '%s' (%d)\n", pdinfo->varname[list[1]], list[1]);
    } 

    /* separation by dummy: create special vars */
    if (gi.flags & GPT_DUMMY) { 
	err = factorized_vars(&gi, Z);
	if (err) {
	    goto bailout;
	}
    } 

    /* special tics for time series plots */
    if (gi.flags & GPT_TS) {
	make_time_tics(&gi, pdinfo, many, xlabel, prn);
    }

    /* open file and dump the prn into it: we delay writing
       the file header till we know a bit more about the plot
    */
    fp = open_gp_stream(PLOT_REGULAR, gi.flags, &err);
    if (err) {
	gretl_print_destroy(prn);
	goto bailout;
    } 

    gi.fp = fp;
    fputs(gretl_print_get_buffer(prn), fp);
    gretl_print_destroy(prn);

    print_axis_label('x', xlabel, fp);
    fputs("set xzeroaxis\n", fp); 
    gnuplot_missval_string(fp);

    if (list[0] == 2) {
	/* only two variables */
	if (gi.flags & GPT_AUTO_FIT) {
	    print_auto_fit_string(gi.fit, fp);
	    if (gi.flags & GPT_FA) {
		make_gtitle(&gi, GTITLE_AFV, var_get_graph_name(pdinfo, list[1]), 
			    var_get_graph_name(pdinfo, list[2]));
	    } else {
		make_gtitle(&gi, GTITLE_VLS, var_get_graph_name(pdinfo, list[1]), 
			    xlabel);
	    }
	}
	if (gi.flags & GPT_RESIDS && !(gi.flags & GPT_DUMMY)) { 
	    make_gtitle(&gi, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	    fprintf(fp, "set ylabel '%s'\n", _("residual"));
	} else {
	    print_axis_label('y', var_get_graph_name(pdinfo, list[1]), fp);
	}
	if (!(gi.flags & GPT_AUTO_FIT)) {
	    strcpy(keystr, "set nokey\n");
	}
    } else if ((gi.flags & GPT_RESIDS) && (gi.flags & GPT_DUMMY)) { 
	make_gtitle(&gi, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	fprintf(fp, "set ylabel '%s'\n", _("residual"));
    } else if (gi.flags & GPT_FA) {
	if (list[3] == pdinfo->v - 1) { 
	    /* x var is just time or index: is this always right? */
	    make_gtitle(&gi, GTITLE_AF, var_get_graph_name(pdinfo, list[2]), NULL);
	} else {
	    make_gtitle(&gi, GTITLE_AFV, var_get_graph_name(pdinfo, list[2]), 
			var_get_graph_name(pdinfo, list[3]));
	}
	print_axis_label('y', var_get_graph_name(pdinfo, list[2]), fp);
    } 

    if (many) {
	strcpy(keystr, "set key outside\n");
    }

    fputs(keystr, fp);

    gretl_push_c_numeric_locale();

    if (gi.x != NULL) {
	print_x_range(&gi, gi.x);
    } else {
	print_x_range_from_list(&gi, Z, list);
    }

    if (*gi.fmt != '\0' && *gi.xtics != '\0') {
	/* remedial handling of broken tics */
	fprintf(fp, "set format x \"%s\"\n", gi.fmt);
	fprintf(fp, "set xtics %s\n", gi.xtics); 
    }

    if (gi.flags & GPT_Y2AXIS) { 
	check_for_yscale(&gi, Z, &oddman);
	if (gi.flags & GPT_Y2AXIS) {
	    fputs("set ytics nomirror\n", fp);
	    fputs("set y2tics\n", fp);
	}
    } else if (gi.yformula == NULL && list[0] == 2) {
	check_y_tics(&gi, Z, fp);
    }

#if GP_DEBUG    
    fprintf(stderr, "literal = '%s', yformula = '%s'\n", literal,
	    gi.yformula);
#endif

    if (gi.yformula != NULL) {
	/* cut out the "dummy" yvar that is in fact represented
	   by a formula rather than raw data */
	list[1] = list[2];
	list[2] = list[3];
	list[0] = 2;
    } else if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    /* now print the 'plot' lines */
    fputs("plot \\\n", fp);
    if (gi.flags & GPT_Y2AXIS) {
	for (i=1; i<list[0]; i++) {
	    set_lwstr(pdinfo, list[i], lwstr);
	    fprintf(fp, "'-' using 1:($2) axes %s title \"%s (%s)\" %s%s%s",
		    (i == oddman)? "x1y2" : "x1y1",
		    var_get_graph_name(pdinfo, list[i]), 
		    (i == oddman)? _("right") : _("left"),
		    (use_impulses(&gi))? "w impulses" : 
		    (gi.flags & GPT_TS)? "w lines" : "w points",
		    lwstr,
		    (i == list[0] - 1)? "\n" : ", \\\n");
	}
    } else if (gi.flags & GPT_DUMMY) { 
	strcpy(s1, (gi.flags & GPT_RESIDS)? _("residual") : 
	       var_get_graph_name(pdinfo, list[1]));
	strcpy(s2, var_get_graph_name(pdinfo, list[3]));
	fprintf(fp, " '-' using 1:($2) title \"%s (%s=1)\", \\\n", s1, s2);
	fprintf(fp, " '-' using 1:($2) title \"%s (%s=0)\"\n", s1, s2);
    } else if (gi.yformula != NULL) {
	fprintf(fp, " '-' using 1:($2) title \"%s\" w points , \\\n", _("actual"));	
	fprintf(fp, "%s title '%s' w lines\n", gi.yformula, _("fitted"));
    } else if (gi.flags & GPT_FA) {
	set_withstr(gi.flags, withstr);
	fprintf(fp, " '-' using 1:($2) title \"%s\" %s lt 2, \\\n", _("fitted"), withstr);
	fprintf(fp, " '-' using 1:($2) title \"%s\" %s lt 1\n", _("actual"), withstr);	
    } else {
	for (i=1; i<list[0]; i++)  {
	    set_lwstr(pdinfo, list[i], lwstr);
	    if (list[0] == 2) {
		*s1 = '\0';
	    } else {
		strcpy(s1, var_get_graph_name(pdinfo, list[i]));
	    }
	    if (!use_impulses(&gi)) { 
		set_withstr(gi.flags, withstr);
	    }
	    fprintf(fp, " '-' using 1:($2) title \"%s\" %s%s", s1, withstr, lwstr);
	    if (i < list[0] - 1 || (gi.flags & GPT_AUTO_FIT)) {
	        fputs(" , \\\n", fp); 
	    } else {
	        fputc('\n', fp);
	    }
	}
    } 

    if (*fit_line != '\0') {
        fputs(fit_line, fp);
    }

    /* print the data to be graphed */
    if (gi.flags & GPT_DUMMY) {
	print_gp_dummy_data(&gi, Z, pdinfo);
    } else {
	print_gp_data(&gi, Z, pdinfo);
    }

    /* flush stream */
    fclose(gi.fp);
    gi.fp = NULL;

    gretl_pop_c_numeric_locale();

    err = gnuplot_make_graph();

 bailout:

    clear_gpinfo(&gi);

    return err;
}

int theil_forecast_plot (const int *plotlist, const double **Z, 
			 const DATAINFO *pdinfo, gretlopt opt)
{
    FILE *fp = NULL;
    gnuplot_info gi;
    int vx, vy;
    int err = 0;

    gretl_error_clear();

    if (plotlist[0] != 2) {
	return E_DATA;
    }

    err = gpinfo_init(&gi, opt | OPT_S, plotlist, NULL, 
		      pdinfo->t1, pdinfo->t2);
    if (err) {
	goto bailout;
    }

    /* ensure the time-series flag is unset */
    gi.flags &= ~GPT_TS;

    graph_list_adjust_sample(gi.list, &gi, Z);
    if (gi.t1 == gi.t2) {
	err = GRAPH_NO_DATA;
	goto bailout;
    }

    fp = open_gp_stream(PLOT_REGULAR, gi.flags, &err);
    if (err) {
	goto bailout;
    } 

    gi.fp = fp;

    vx = gi.list[2];
    vy = gi.list[1];

    print_axis_label('x', var_get_graph_name(pdinfo, vx), fp);
    print_axis_label('y', var_get_graph_name(pdinfo, vy), fp);
	   
    fputs("set xzeroaxis\n", fp); 
    gnuplot_missval_string(fp);
    fputs("set key left top\n", fp);

    gretl_push_c_numeric_locale();

    print_x_range_from_list(&gi, Z, gi.list);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:($2) notitle w points , \\\n", fp);
    fprintf(fp, " x title \"%s\" w lines\n", _("actual = predicted"));

    print_gp_data(&gi, Z, pdinfo);

    fclose(gi.fp);
    gi.fp = NULL;

    gretl_pop_c_numeric_locale();

    err = gnuplot_make_graph();

 bailout:

    clear_gpinfo(&gi);

    return err;
}

/**
 * multi_scatters:
 * @list: list of variables to plot, by ID number.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: can include %OPT_L to use lines, %OPT_U to
 * direct output to a named file.
 *
 * Writes a gnuplot plot file to display up to 6 small graphs
 * based on the variables in @list, and calls gnuplot to make 
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_scatters (const int *list, const double **Z, 
		    const DATAINFO *pdinfo, gretlopt opt)
{
    GptFlags flags = 0;
    int xvar = 0, yvar = 0;
    const double *obs = NULL;
    int *plotlist = NULL;
    int pos, nplots = 0;
    FILE *fp = NULL;
    int i, t, err = 0;

    if (gretl_in_batch_mode()) {
	flags |= GPT_BATCH;
    }

    if (opt & OPT_L) {
	flags |= GPT_LINES;
    }

    pos = gretl_list_separator_position(list);

    if (pos == 0) {
	/* plot against time or index */
	obs = gretl_plotx(pdinfo);
	if (obs == NULL) {
	    return E_ALLOC;
	}
	plotlist = gretl_list_copy(list);
	flags |= GPT_LINES;
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

    fp = open_gp_stream(PLOT_MULTI_SCATTER, flags, &err);
    if (err) {
	return err;
    }

    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (obs != NULL) {
	double startdate = obs[pdinfo->t1];
	double enddate = obs[pdinfo->t2];
	int jump;

	fprintf(fp, "set xrange [%g:%g]\n", floor(startdate), ceil(enddate));

	if (pdinfo->pd == 1) {
	    jump = (pdinfo->t2 - pdinfo->t1 + 1) / 6;
	} else {
	    jump = (pdinfo->t2 - pdinfo->t1 + 1) / (4 * pdinfo->pd);
	}

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

	if (obs != NULL) {
	    fputs("set noxlabel\n", fp);
	    fputs("set noylabel\n", fp);
	    fprintf(fp, "set title '%s'\n", 
		    var_get_graph_name(pdinfo, pv));
	} else {
	    fprintf(fp, "set xlabel '%s'\n",
		    (yvar)? pdinfo->varname[pv] :
		    pdinfo->varname[xvar]);
	    fprintf(fp, "set ylabel '%s'\n", 
		    (yvar)? pdinfo->varname[yvar] :
		    pdinfo->varname[pv]);
	}

	fputs("plot '-' using 1:2", fp);
	if (flags & GPT_LINES) {
	    fputs(" with lines", fp);
	}
	fputc('\n', fp);

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    double xx;

	    xx = (yvar)? Z[pv][t] : (xvar)? Z[xvar][t] : obs[t];

	    if (na(xx)) {
		fputs("? ", fp);
	    } else {
		fprintf(fp, "%.10g ", xx);
	    }

	    xx = (yvar)? Z[yvar][t] : Z[pv][t];

	    if (na(xx)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.10g\n", xx);
	    }
	}

	fputs("e\n", fp);
    } 

    gretl_pop_c_numeric_locale();

    fputs("set nomultiplot\n", fp);

    fclose(fp);

    err = gnuplot_make_graph();

    free(plotlist);

    return err;
}

static int get_3d_output_file (FILE **fpp)
{
    char fname[MAXLEN];
    int err = 0;

    sprintf(fname, "%sgpttmp.plt", gretl_dotdir());
    *fpp = gretl_fopen(fname, "w");

    if (*fpp == NULL) {
	err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
    }

    return err;
}

static gchar *maybe_get_surface (const int *list, 
				 double ***pZ, DATAINFO *pdinfo, 
				 gretlopt opt)
{
    MODEL smod;
    double umin, umax, vmin, vmax;
    int olslist[5];
    gchar *ret = NULL;

    olslist[0] = 4;
    olslist[1] = list[3];
    olslist[2] = 0;
    olslist[3] = list[2];
    olslist[4] = list[1];

    gretl_minmax(pdinfo->t1, pdinfo->t2, (*pZ)[list[2]], &umin, &umax);
    gretl_minmax(pdinfo->t1, pdinfo->t2, (*pZ)[list[1]], &vmin, &vmax);

    smod = lsq(olslist, pZ, pdinfo, OLS, OPT_A);

    if (!smod.errcode && !na(smod.fstt) &&
	(snedecor_cdf_comp(smod.dfn, smod.dfd, smod.fstt) < .10 || (opt & OPT_F))) {
	double uadj = (umax - umin) * 0.02;
	double vadj = (vmax - vmin) * 0.02;

	ret = g_strdup_printf("[u=%g:%g] [v=%g:%g] "
			      "%g+(%g)*u+(%g)*v title ''", 
			      umin - uadj, umax + uadj, 
			      vmin - vadj, vmax + vadj,
			      smod.coeff[0], smod.coeff[1],
			      smod.coeff[2]);
    } 

    clear_model(&smod);

    return ret;
}

/**
 * gnuplot_3d:
 * @list: list of variables to plot, by ID number: Y, X, Z
 * @literal: literal command(s) to pass to gnuplot (or NULL)
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: unused at present.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (int *list, const char *literal,
		double ***pZ, DATAINFO *pdinfo,  
		gretlopt opt)
{
    FILE *fq = NULL;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int orig_t1 = pdinfo->t1, orig_t2 = pdinfo->t2;
    int lo = list[0];
    int datlist[4];
    int addstyle = 0;
    gchar *surface = NULL;

    if (lo != 3) {
	fprintf(stderr, "gnuplot_3d needs three variables (only)\n");
	return E_DATA;
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

#ifndef WIN32
    if (gnuplot_has_wxt()) {
	fputs("set term wxt\n", fq);
    } else if (gnuplot_has_x11()) {
	fputs("set term x11\n", fq);
    } else {
	fclose(fq);
	return E_EXTERNAL;
    }
#endif

    gretl_push_c_numeric_locale();

    /* try to ensure we don't get "invisible" green datapoints */
    fprintf(fq, "set style line 2 lc rgb \"#0000ff\"\n");
    addstyle = 1;
    
    print_axis_label('x', var_get_graph_name(pdinfo, list[2]), fq);
    print_axis_label('y', var_get_graph_name(pdinfo, list[1]), fq);
    print_axis_label('z', var_get_graph_name(pdinfo, list[3]), fq);

    gnuplot_missval_string(fq);

    if (literal != NULL && *literal != 0) {
	print_gnuplot_literal_lines(literal, fq);
    }

    surface = maybe_get_surface(list, pZ, pdinfo, opt);

    if (surface != NULL) {
	if (addstyle) {
	    fprintf(fq, "splot %s, \\\n'-' title '' w p ls 2\n", surface);
	} else {
	    fprintf(fq, "splot %s, \\\n'-' title '' w p lt 3\n", surface);
	}
	g_free(surface);
    } else {
	if (addstyle) {
	    fputs("splot '-' title '' w p ls 2\n", fq);
	} else {
	    fputs("splot '-' title '' w p lt 3\n", fq);
	}
    }

    datlist[0] = 3;
    datlist[1] = list[2];
    datlist[2] = list[1];
    datlist[3] = list[3];

    for (t=t1; t<=t2; t++) {
	const char *label = NULL;

	if (pdinfo->markers) {
	    label = pdinfo->S[t];
	}
	printvars(fq, t, datlist, (const double **) *pZ, NULL, label, 0.0);
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
    char test[8];

    gretl_pop_c_numeric_locale();

    sprintf(test, "%g", 0.5);
    if (strchr(test, ',')) {
	dcomma = 1;
    }

    if (dist == D_NORMAL) {
	sprintf(s, "N(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    } else if (dist == D_GAMMA) {
	sprintf(s, "gamma(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    }

    gretl_push_c_numeric_locale();
}

/* Below: a fix for the case where the y-range is by default
   degenerate, in which case gnuplot produces a graph OK, but
   issues a warning and returns non-zero.
*/

static void maybe_set_yrange (FreqDist *freq, double lambda, FILE *fp)
{
    double ymin = 1.0e+200;
    double ymax = -1.0e+200;
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
	fprintf(fp, "set yrange [%.10g:%.10g]\n", ymax * lambda * 0.99, 
		ymax * lambda * 1.01);
    } else {
	fprintf(fp, "set yrange [0.0:%.10g]\n", ymax * lambda * 1.1);
    }	
}

static double minskip (FreqDist *freq)
{
    double s, ms = freq->midpt[1] - freq->midpt[0];
    int i;

    for (i=2; i<freq->numbins; i++) {
	s = freq->midpt[i] - freq->midpt[i-1];
	if (s < ms) {
	    ms = s;
	}
    }

    return ms;
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
    char withstr[32] = {0};
    char label[80] = {0};
    double plotmin = 0.0, plotmax = 0.0;
    double barwidth;
    const double *endpt;
    int plottype, use_boxes = 1;
    int err = 0;

    if (K == 0) {
	return E_DATA;
    }

    if (K == 1) {
	gretl_errmsg_sprintf(_("'%s' is a constant"), freq->varname);
	return E_DATA;
    }

    if (dist == D_NORMAL) {
	plottype = PLOT_FREQ_NORMAL;
    } else if (dist == D_GAMMA) {
	plottype = PLOT_FREQ_GAMMA;
    } else {
	plottype = PLOT_FREQ_SIMPLE;
    }

    fp = get_plot_input_stream(plottype, &err);
    if (err) {
	return err;
    }  

#if GP_DEBUG
    fprintf(stderr, "*** plot_freq called\n");
#endif  

    if (freq->discrete) {
	endpt = freq->midpt;
	barwidth = minskip(freq); 
	use_boxes = 0;
    } else {
	/* equally sized bins, width to be determined */
	endpt = freq->endpt;
	barwidth = freq->endpt[K-1] - freq->endpt[K-2];
    }

    gretl_push_c_numeric_locale();

    if (dist) {
	lambda = 1.0 / (freq->n * barwidth);

	if (dist == D_NORMAL) {
	    fputs("# literal lines = 4\n", fp);
	    fprintf(fp, "sigma = %g\n", freq->sdx);
	    fprintf(fp, "mu = %g\n", freq->xbar);

	    plotmin = endpt[0] - barwidth;
	    if (plotmin > freq->xbar - 3.3 * freq->sdx) {
		plotmin = freq->xbar - 3.3 * freq->sdx;
	    }

	    plotmax = endpt[K-1] + barwidth;
	    if (plotmax < freq->xbar + 3.3 * freq->sdx) {
		plotmax = freq->xbar + 3.3 * freq->sdx;
	    }

	    if (!na(freq->test)) {
		fprintf(fp, "set label \"%s:\" at graph .03, graph .97 front\n",
			_("Test statistic for normality"));
		print_freq_test_label(label, _("Chi-squared(2) = %.3f pvalue = %.5f"), 
				      freq->test, chisq_cdf_comp(2, freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93 front\n", 
			label);
	    }	
	} else if (dist == D_GAMMA) {
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
		fprintf(fp, "set label '%s:' at graph .03, graph .97 front\n",
			_("Test statistic for gamma"));
		print_freq_test_label(label, _("z = %.3f pvalue = %.5f"), 
				      freq->test, normal_pvalue_2(freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93 front\n", 
			label);
	    }	
	}

	/* adjust min, max if needed */
	if (freq->midpt[0] < plotmin) {
	    plotmin = freq->midpt[0];
	}
	if (freq->midpt[K-1] > plotmax) {
	    plotmax = freq->midpt[K-1];
	}

	fprintf(fp, "set xrange [%.10g:%.10g]\n", plotmin, plotmax);
	fputs("set key right top\n", fp);
    } else { 
	/* plain frequency plot (no theoretical distribution shown) */
	lambda = 1.0 / freq->n;
	plotmin = freq->midpt[0] - barwidth;
	plotmax = freq->midpt[K-1] + barwidth;
	fprintf(fp, "set xrange [%.10g:%.10g]\n", plotmin, plotmax);
	maybe_set_yrange(freq, lambda, fp);
	fputs("set nokey\n", fp);
    }

    fprintf(fp, "set xlabel '%s'\n", freq->varname);
    if (dist) {
	fprintf(fp, "set ylabel '%s'\n", _("Density"));
    } else {
	fprintf(fp, "set ylabel '%s'\n", _("Relative frequency"));
    }

    if (isnan(lambda)) {
	if (fp != NULL) {
	    fclose(fp);
	}
	return 1;
    }

    /* plot instructions */
    if (use_boxes) {
	fputs("set style fill solid 0.6\n", fp);
	strcpy(withstr, "w boxes");
    } else {
	strcpy(withstr, "w impulses linewidth 3");
    }

    if (!dist) {
	fprintf(fp, "plot '-' using 1:2 %s\n", withstr);
    } else if (dist == D_NORMAL) {
	print_freq_dist_label(label, dist, freq->xbar, freq->sdx);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title \"%s\" %s, \\\n"
		"1.0/(sqrt(2.0*pi)*sigma)*exp(-.5*((x-mu)/sigma)**2) "
		"title \"%s\" w lines\n",
		freq->varname, withstr, label);
    } else if (dist == D_GAMMA) {
	print_freq_dist_label(label, dist, alpha, beta);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' %s, \\\n"
		"x**(alpha-1.0)*exp(-x/beta)/(exp(lgamma(alpha))*(beta**alpha)) "
		"title \"%s\" w lines\n",
		freq->varname, withstr, label); 
    }

    for (i=0; i<K; i++) { 
	fprintf(fp, "%.10g %.10g\n", freq->midpt[i], lambda * freq->f[i]);
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    if (fp != NULL) {
	fclose(fp);
    }

    return gnuplot_make_graph();
}

static void print_y_data (const double *x, const double *y,
			  int t1, int t2, FILE *fp)
{
    int t;

    for (t=t1; t<=t2; t++) {
	if (na(y[t])) {
	    fprintf(fp, "%.10g ?\n", x[t]);
	} else {
	    fprintf(fp, "%.10g %.10g\n", x[t], y[t]);
	}
    }
    fputs("e\n", fp);
}

enum {
    CONF_BARS,
    CONF_FILL,
    CONF_LOW,
    CONF_HIGH
};

static void print_confband_data (const double *x, const double *y,
				 const double *e, int t1, int t2, 
				 int mode, FILE *fp)
{
    int t;

    for (t=t1; t<=t2; t++) {
	if (na(y[t]) || na(e[t])) {
	    fprintf(fp, "%.10g ? ?\n", x[t]);
	} else if (mode == CONF_FILL) {
	    fprintf(fp, "%.10g %.10g %.10g\n", x[t], y[t] - e[t], y[t] + e[t]);
	} else if (mode == CONF_LOW) {
	    fprintf(fp, "%.10g %.10g\n", x[t], y[t] - e[t]);
	} else if (mode == CONF_HIGH) {
	    fprintf(fp, "%.10g %.10g\n", x[t], y[t] + e[t]);
	} else {
	    fprintf(fp, "%.10g %.10g %.10g\n", x[t], y[t], e[t]);
	} 
    }
    fputs("e\n", fp);
}

int plot_fcast_errs (const FITRESID *fr, const double *maxerr,
		     const DATAINFO *pdinfo, gretlopt opt)
{
    FILE *fp = NULL;
    const double *obs = NULL;
    GptFlags flags = 0;
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int use_fill = 0, use_lines = 0;
    int do_errs = (maxerr != NULL);
    char cistr[64];
    int t2 = fr->t2;
    int t1, yhmin;
    int t, n, err = 0;

    /* note: yhmin is the first obs at which to start plotting y-hat */
    if (do_errs) {
	t1 = fr->t0;
	yhmin = (opt & OPT_H)? fr->t0 : fr->t1;
    } else {
	t1 = (fr->t0 >= 0)? fr->t0 : 0;
	yhmin = t1;
    }

    /* don't graph empty trailing portion of forecast */
    for (t=fr->t2; t>=t1; t--) {
	if (na(fr->actual[t]) && na(fr->fitted[t])) {
	    t2--;
	} else {
	    break;
	}
    }

    n = t2 - t1 + 1;

    if (n < 3) {
	/* won't draw a graph for 2 data points or less */
	return 1;
    }

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if (opt & OPT_U) {
	/* respond to command-line --plot=fname option */
	flags = GPT_BATCH;
    }

    fp = open_gp_stream(PLOT_FORECAST, flags, &err);
    if (err) {
	return err;
    }    

    /* check that we have any values for the actual var */
    for (t=t1; t<=t2; t++) {
	if (!na(fr->actual[t])) {
	    depvar_present = 1;
	    break;
	}
    }

    if (do_errs) {
	if (opt & OPT_F) {
	    use_fill = 1;
	} else if (opt & OPT_L) {
	    use_lines = 1;
	}
    }

    gretl_minmax(t1, t2, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();
    fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
    gretl_pop_c_numeric_locale();

    gnuplot_missval_string(fp);

    if (dataset_is_time_series(pdinfo)) {
	fprintf(fp, "# timeseries %d\n", pdinfo->pd);
    } else if (n < 33) {
	fputs("set xtics 1\n", fp);
    }

    if (do_errs && !use_fill && !use_lines && n > 150) {
	use_fill = 1;
    }

    if (use_fill) {
	fprintf(fp, "set style fill solid 0.4\n");
    }

    fputs("set key left top\n", fp);
    fputs("plot \\\n", fp);

    if (do_errs) {
	sprintf(cistr, _("%g percent interval"), 100 * (1 - fr->alpha));
    }

    if (use_fill) {
	/* plot the confidence bands first so the other lines
	   come out on top */
	if (do_errs) {
	    fprintf(fp, "'-' using 1:2:3 title '%s' w filledcurve lt 3 , \\\n",
		    cistr);
	} 
	if (depvar_present) {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1 , \\\n",
		    fr->depvar);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines lt 2\n", _("forecast"));
    } else {
	/* plot confidence bands last */
	if (depvar_present) {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines , \\\n",
		    fr->depvar);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines", _("forecast"));
	if (do_errs) {
	    if (use_lines) {
		fprintf(fp, " , \\\n'-' using 1:2 title '%s' w lines , \\\n",
			cistr);
		fputs("'-' using 1:2 notitle '%s' w lines lt 3\n", fp);
	    } else {
		fprintf(fp, " , \\\n'-' using 1:2:3 title '%s' w errorbars\n",
			cistr);
	    }
	} else {
	    fputc('\n', fp);
	}
    }

    gretl_push_c_numeric_locale();

    /* write out the inline data, the order depending on whether
       or not we're using fill style for the confidence
       bands
    */

    if (use_fill) {
	if (do_errs) {
	    print_confband_data(obs, fr->fitted, maxerr, yhmin, t2, CONF_FILL, fp);
	}
	if (depvar_present) {
	    print_y_data(obs, fr->actual, t1, t2, fp);
	}
	print_y_data(obs, fr->fitted, yhmin, t2, fp);
    } else {
	if (depvar_present) {
	    print_y_data(obs, fr->actual, t1, t2, fp);
	}
	print_y_data(obs, fr->fitted, yhmin, t2, fp);
	if (do_errs) {
	    if (use_lines) {
		print_confband_data(obs, fr->fitted, maxerr, yhmin, t2, CONF_LOW, fp);
		print_confband_data(obs, fr->fitted, maxerr, yhmin, t2, CONF_HIGH, fp);
	    } else {
		print_confband_data(obs, fr->fitted, maxerr, yhmin, t2, CONF_BARS, fp);
	    }
	}	
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

#ifndef min
# define min(x,y) (((x)<(y))? (x):(y))
#endif

#ifndef max
# define max(x,y) (((x)>(y))? (x):(y))
#endif

int plot_tau_sequence (const MODEL *pmod, const DATAINFO *pdinfo,
		       int k)
{
    FILE *fp;
    gretl_matrix *tau = gretl_model_get_data(pmod, "rq_tauvec");
    gretl_matrix *B = gretl_model_get_data(pmod, "rq_sequence");
    double tau_i, bi, se, blo, bhi;
    double alpha, cval, tcrit, olsband;
    double ymin[2], ymax[2];
    gchar *tmp;
    int ntau, bcols;
    int i, j, err = 0;

    if (tau == NULL || B == NULL) {
	return E_DATA;
    }

    ntau = gretl_vector_get_length(tau);
    if (ntau == 0) {
	return E_DATA;
    }

    fp = get_plot_input_stream(PLOT_RQ_TAU, &err);
    if (err) {
	return err;
    } 

    bcols = gretl_matrix_cols(B);

    alpha = gretl_model_get_double(pmod, "rq_alpha");
    if (na(alpha)) {
	alpha = .05;
    }

    cval = 100 * (1 - alpha);
    tcrit = student_cdf_inverse(pmod->dfd, 1 - alpha/2);
    olsband = tcrit * pmod->sderr[k];

    /* Try to figure best placement of key */

    j = k * ntau;
    if (bcols == 3) {
	blo = gretl_matrix_get(B, j, 1);
	bhi = gretl_matrix_get(B, j, 2);
    } else {
	bi = gretl_matrix_get(B, j, 0);
	se = gretl_matrix_get(B, j, 1);
	blo = bi - tcrit * se;
	bhi = bi + tcrit * se;
    }
    ymin[0] = min(blo, pmod->coeff[k] - olsband);
    ymax[0] = max(bhi, pmod->coeff[k] + olsband);

    j += ntau - 1;
    if (bcols == 3) {
	blo = gretl_matrix_get(B, j, 1);
	bhi = gretl_matrix_get(B, j, 2);
    } else {
	bi = gretl_matrix_get(B, j, 0);
	se = gretl_matrix_get(B, j, 1);
	blo = bi - tcrit * se;
	bhi = bi + tcrit * se;
    }
    ymin[1] = min(blo, pmod->coeff[k] - olsband);
    ymax[1] = max(bhi, pmod->coeff[k] + olsband);

    fputs("set xrange [0.0:1.0]\n", fp);
    fputs("set xlabel 'tau'\n", fp);

    tmp = g_strdup_printf(_("Coefficient on %s"), 
			  var_get_graph_name(pdinfo, pmod->list[k+2]));
    fprintf(fp, "set title \"%s\"\n", tmp);
    g_free(tmp);

    fputs("set style fill solid 0.4\n", fp);

    if (ymax[0] < .88 * ymax[1]) {
	fputs("set key left top\n", fp);
    } else if (ymax[1] < .88 * ymax[0]) {
	fputs("set key right top\n", fp);
    } else if (ymin[0] < .88 * ymin[1]) {
	fputs("set key right bottom\n", fp);
    } else {
	fputs("set key left bottom\n", fp);
    }

    fputs("plot \\\n", fp);

    /* plot the rq confidence band first so the other lines
       come out on top */
    fputs("'-' using 1:2:3 notitle w filledcurve lt 3 , \\\n", fp);

    /* rq estimates */
    tmp = g_strdup_printf(_("Quantile estimates with %g%% band"), cval);
    fprintf(fp, "'-' using 1:2 title '%s' w lp lt 1 , \\\n", tmp);
    g_free(tmp);

    /* numeric output coming up! */
    gretl_push_c_numeric_locale();

    /* ols estimate plus (1 - alpha) band */
    tmp = g_strdup_printf(_("OLS estimate with %g%% band"), cval);
    fprintf(fp, "%g title '%s' w lines lt 2 , \\\n", pmod->coeff[k], tmp);
    g_free(tmp);
    fprintf(fp, "%g notitle w dots lt 2 , \\\n", pmod->coeff[k] + olsband);
    fprintf(fp, "%g notitle w dots lt 2\n", pmod->coeff[k] - olsband);

    /* write out the interval values */

    for (i=0, j=k*ntau; i<ntau; i++, j++) {
	tau_i = gretl_vector_get(tau, i);
	if (bcols == 3) {
	    blo = gretl_matrix_get(B, j, 1);
	    bhi = gretl_matrix_get(B, j, 2);
	} else {
	    bi = gretl_matrix_get(B, j, 0);
	    se = gretl_matrix_get(B, j, 1);
	    blo = bi - tcrit * se;
	    bhi = bi + tcrit * se;
	}
	fprintf(fp, "%.10g %.10g %.10g\n", tau_i, blo, bhi);
    }
    fputs("e\n", fp);

    for (i=0, j=k*ntau; i<ntau; i++, j++) {
	tau_i = gretl_vector_get(tau, i);
	bi = gretl_matrix_get(B, j, 0);
	fprintf(fp, "%.10g %.10g\n", tau_i, bi);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();    
}

int garch_resid_plot (const MODEL *pmod, const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    const double *obs;
    const double *h;
    double sd2;
    int t, err = 0;

    h = gretl_model_get_data(pmod, "garch_h");
    if (h == NULL) {
	return E_DATA;
    }

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    fp = get_plot_input_stream(PLOT_GARCH, &err);
    if (err) {
	return err;
    }

    fputs("set key left top\n", fp);

    fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines, \\\n"
	    "'-' using 1:2 title '%s' w lines lt 2, \\\n" 
	    "'-' using 1:2 notitle w lines lt 2\n", 
	    _("residual"), _("+- sqrt(h(t))"));

    gretl_push_c_numeric_locale();

    for (t=pmod->t1; t<=pmod->t2; t++) {
	fprintf(fp, "%.10g %.10g\n", obs[t], pmod->uhat[t]);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = -sqrt(h[t]);
	fprintf(fp, "%.10g %.10g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = sqrt(h[t]);
	fprintf(fp, "%.10g %.10g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int rmplot (const int *list, const double **Z, DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn)
{
    int (*range_mean_graph) (int, const double **, const DATAINFO *, 
			     gretlopt, PRN *);
    void *handle = NULL;
    int err;

    range_mean_graph = get_plugin_function("range_mean_graph", &handle);
    if (range_mean_graph == NULL) {
        return 1;
    }

    err = range_mean_graph(list[1], Z, pdinfo, opt, prn);

    close_plugin(handle);

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

static int panel_ytic_width (double ymin, double ymax)
{
    char s1[16], s2[16];
    int n1, n2;

    if (ymin < 0 && ymax > 0) {
	sprintf(s1, "% g", ymin);
	sprintf(s2, "% g", ymax);
    } else {
	sprintf(s1, "%g", ymin);
	sprintf(s2, "%g", ymax);
    }

    n1 = strlen(s1);
    n2 = strlen(s2);

    return (n1 > n2)? n1 : n2;
}

/* Panel: plot one variable as a time series, with separate plots for
   each cross-sectional unit.  By default we arrange the plots in a
   grid, but if OPT_V is given we make each plot full width and
   stack the plots vertically on the "page".
*/

int 
gretl_panel_ts_plot (const int *list, const double **Z, DATAINFO *pdinfo,
		     gretlopt opt) 
{
    FILE *fp = NULL;
    int i, j, k, v, t, t0;
    int w, xnum, ynum;
    float xfrac, yfrac;
    float yorig, xorig = 0.0;
    const double *y;
    double yt, ymin, ymax, incr;
    int nunits, T = pdinfo->pd;
    int err = 0;

    nunits = pdinfo->paninfo->unit[pdinfo->t2] -
	pdinfo->paninfo->unit[pdinfo->t1] + 1;
    
    if (opt & OPT_V) {
	xnum = 1;
	ynum = nunits;
    } else {
	get_x_and_y_sizes(nunits, &xnum, &ynum);
    }

    if (xnum == 0 || ynum == 0) {
	return E_DATA;
    }

    gp_small_font_size = (nunits > 4)? 7 : 0;

    fp = get_plot_input_stream(PLOT_PANEL, &err);
    if (err) {
	return err;
    }

    v = list[1];
    y = Z[v];
    gretl_minmax(pdinfo->t1, pdinfo->t2, y, &ymin, &ymax);
    w = panel_ytic_width(ymin, ymax);

    xfrac = 1.0 / xnum;
    yfrac = 1.0 / ynum;

    fputs("set key left top\n", fp);
    fputs("set datafile missing \"?\"\n", fp);
    fputs("set xtics nomirror\n", fp);
    fputs("set ytics nomirror\n", fp);
    fprintf(fp, "set format y \"%%%dg\"\n", w);
    fputs("set multiplot\n", fp);

    if (opt & OPT_V) {
	fputs("set noxlabel\n", fp);
    } else {
	fprintf(fp, "set xlabel '%s'\n", _("time"));
    }

    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    if (yfrac > 1.4 * xfrac) {
	yfrac = 1.4 * xfrac;
    }
    fprintf(fp, "set size %g,%g\n", xfrac, yfrac);

    k = 0;
    t0 = pdinfo->t1;

    for (i=0; i<xnum && k<nunits; i++) {

	yorig = 1.0 - yfrac;

	for (j=0; j<ynum && k<nunits; j++, k++) {

	    fprintf(fp, "set origin %g,%g\n", xorig, yorig);

	    if (opt & OPT_V) {
		gretl_minmax(t0, t0 + T - 1, y, &ymin, &ymax);
		incr = (ymax - ymin) / 2.0;
		fprintf(fp, "set ytics %g\n", incr);
		fprintf(fp, "set ylabel '%s (%d)'\n", pdinfo->varname[v], k+1);
	    } else {
		fprintf(fp, "set title '%s (%d)'\n", pdinfo->varname[v], k+1);
	    }

	    fputs("plot \\\n'-' using 1:($2) notitle w lines\n", fp);

	    for (t=0; t<T; t++) {
		yt = y[t+t0];
		if (na(yt)) {
		    fprintf(fp, "%d ?\n", t+1);
		} else {
		    fprintf(fp, "%d %.10g\n", t+1, yt);
		}
	    }
	    fputs("e\n", fp);

	    if (k < nunits) {
		t0 += T;
		yorig -= yfrac;
	    }
	}

	if (k < nunits) {
	    xorig += xfrac;
	}
    }

    fputs("unset multiplot\n", fp);
    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int 
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, int periods,
				 double alpha,
				 const double **Z,
				 const DATAINFO *pdinfo,
				 gretlopt opt)
{
    FILE *fp = NULL;
    int confint = 0;
    int use_fill = !(opt & OPT_E);
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int t, err = 0;

    if (alpha < 0.01 || alpha > 0.5) {
	return E_DATA;
    }

    resp = gretl_VAR_get_impulse_response(var, targ, shock, periods,
					  alpha, Z, pdinfo);
    if (resp == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(resp) > 1) {
	confint = 1;
    }

    fp = get_plot_input_stream((confint)? PLOT_IRFBOOT : PLOT_REGULAR, &err);
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
	fputs("set key left top\n", fp);
	sprintf(title, _("response of %s to a shock in %s, "
			 "with bootstrap confidence interval"),
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    } else {
	fputs("set nokey\n", fp);
	sprintf(title, _("response of %s to a shock in %s"), 
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    }

    fprintf(fp, "set xlabel '%s'\n", _("periods"));
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n", title);

    if (confint) {
	double ql = alpha / 2;
	double qh = 1.0 - ql;

	fputs("plot \\\n", fp);
	if (use_fill) {
	    sprintf(title, _("%g percent confidence band"), 100 * (1 - alpha));
	    fprintf(fp, "'-' using 1:2:3 title '%s' w filledcurve lt %d, \\\n", 
		    title, SHADECOLOR + 1);
	    fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1\n", _("point estimate"));
	} else {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n", 
		_("point estimate"));
	    sprintf(title, _("%g and %g quantiles"), ql, qh);
	    fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n", title);
	}
    } else {
	fputs("plot \\\n'-' using 1:2 w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    if (confint && use_fill) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.10g %.10g\n", t+1, 
		    gretl_matrix_get(resp, t, 1),
		    gretl_matrix_get(resp, t, 2));
	}
	fputs("e\n", fp);
    }

    for (t=0; t<periods; t++) {
	fprintf(fp, "%d %.10g\n", t+1, gretl_matrix_get(resp, t, 0));
    }
    fputs("e\n", fp);

    if (confint && !use_fill) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.10g %.10g %.10g\n", t+1, 
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
			     double alpha,
			     const double **Z,
			     const DATAINFO *pdinfo,
			     gretlopt opt)
{
    FILE *fp = NULL;
    int confint = 0;
    int use_fill = !(opt & OPT_E);
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int n = var->neqns;
    float plot_fraction = 1.0 / n;
    float xorig = 0.0;
    float yorig;
    int t, i, j;
    int err = 0;

    gp_small_font_size = (n == 4)? 6 : 0;

    resp = gretl_VAR_get_impulse_response(var, 1, 1, periods, alpha, Z, pdinfo);
    if (resp == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(resp) > 1) {
	confint = 1;
    }

    fp = get_plot_input_stream(PLOT_MULTI_IRF, &err);
    if (err) {
	gretl_matrix_free(resp);
	return err;
    }

    if (!confint) {
	fputs("set nokey\n", fp);
    } else {
	fputs("set key left top\n", fp);
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
	    resp = gretl_VAR_get_impulse_response(var, i, j, periods, alpha, Z, pdinfo);
	    if (resp == NULL) {
		return E_ALLOC;
	    }

	    vshock = gretl_VAR_get_variable_number(var, j);
	    sprintf(title, "%s -> %s", pdinfo->varname[vshock], pdinfo->varname[vtarg]);
	    fprintf(fp, "set title '%s'\n", title);

	    fputs("plot \\\n", fp);

	    if (confint && use_fill) {
		fprintf(fp, "'-' using 1:2:3 notitle w filledcurve lt %d, \\\n", 
			SHADECOLOR + 1);
		fputs("'-' using 1:2 notitle w lines lt 1\n", fp);
	    } else if (confint) {
		fputs("'-' using 1:2 notitle w lines, \\\n", fp); 
		fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
	    } else {
		fputs("'-' using 1:2 notitle w lines\n", fp);
	    }

	    if (confint && use_fill) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.10g %.10g\n", t+1, 
			    gretl_matrix_get(resp, t, 1),
			    gretl_matrix_get(resp, t, 2));
		}
		fputs("e\n", fp);
	    }		

	    for (t=0; t<periods; t++) {
		fprintf(fp, "%d %.10g\n", t+1, gretl_matrix_get(resp, t, 0));
	    }
	    fputs("e\n", fp);

	    if (confint && !use_fill) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.10g %.10g %.10g\n", t+1, 
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

int gretl_system_residual_plot (void *p, int ci, const DATAINFO *pdinfo)
{
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *E = NULL;
    FILE *fp = NULL;
    const double *obs;
    int nvars, nobs;
    int i, v, t, t1;
    int err = 0;

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) p;
	E = gretl_VAR_get_residual_matrix(var);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) p;
	E = sys->E;
    }

    if (E == NULL) {
	return E_DATA;
    }

    t1 = E->t1;

    fp = get_plot_input_stream(PLOT_REGULAR, &err);
    if (err) {
	return err;
    }

    obs = gretl_plotx(pdinfo);

    nvars = gretl_matrix_cols(E);
    nobs = gretl_matrix_rows(E);

    fputs("# system residual plot\n", fp);
    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (ci == VAR) {
	fprintf(fp, "set title '%s'\n", _("VAR residuals"));
    } else {
	fprintf(fp, "set title '%s'\n", _("System residuals"));
    }

    fputs("plot \\\n", fp);
    for (i=0; i<nvars; i++) {
	if (var != NULL) {
	    v = gretl_VAR_get_variable_number(var, i);
	} else {
	    v = system_get_depvar(sys, i);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines", pdinfo->varname[v]);
	if (i == nvars - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp); 
	}
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<nvars; i++) {
	for (t=0; t<nobs; t++) {
	    double eti = gretl_matrix_get(E, t, i);

	    if (obs != NULL) {
		fprintf(fp, "%g %.10g\n", obs[t+t1], eti);
	    } else {
		fprintf(fp, "%d %.10g\n", t+1, eti);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int gretl_system_residual_mplot (void *p, int ci, const DATAINFO *pdinfo) 
{
    const gretl_matrix *E = NULL;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    FILE *fp = NULL;
    const double *obs;
    double startdate;
    double xmin, xmax, xrange;
    int nvars, nobs, jump;
    int i, v, t, t1;
    int err = 0;

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) p;
	E = gretl_VAR_get_residual_matrix(var);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) p;
	E = sys->E;
    }

    if (E == NULL) {
	return E_DATA;
    }

    nvars = gretl_matrix_cols(E);
    if (nvars > 6) {
	return 1;
    }

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    nobs = gretl_matrix_rows(E);
    t1 = E->t1;

    fp = get_plot_input_stream(PLOT_MULTI_SCATTER, &err);
    if (err) {
	return err;
    }

    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    startdate = obs[t1];
    jump = nobs / (2 * pdinfo->pd);
    fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);

    gretl_minmax(t1, t1 + nobs - 1, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;
    fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);	

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
	if (var != NULL) {
	    v = gretl_VAR_get_variable_number(var, i);
	} else {
	    v = system_get_depvar(sys, i);
	}
	fprintf(fp, "set title '%s'\n", pdinfo->varname[v]);

	fputs("plot '-' using 1:2 with lines\n", fp);

	for (t=0; t<nobs; t++) {
	    double xx;

	    fprintf(fp, "%.10g\t", obs[t+t1]);
	    xx = gretl_matrix_get(E, t, i);
	    if (na(xx)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.10g\n", xx);
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
    int i, n, err = 0;

    lam = gretl_VAR_get_roots(var, &err);
    if (err) {
	return err;
    }

    fp = get_plot_input_stream(PLOT_VAR_ROOTS, &err);
    if (err) {
	return err;
    }

    n = gretl_matrix_rows(lam);

    fprintf(fp, "set title '%s'\n", 
	    _("VAR inverse roots in relation to the unit circle"));
    fputs("# literal lines = 8\n", fp);
    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1 w lines, \\\n'-' w points pt 7\n", fp);

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

/**
 * confidence_ellipse_plot:
 * @V: 2x2 covariance matrix.
 * @b: 2-vector containing point estimates
 * @tcrit: critical t-value for 1 - alpha confidence.
 * @Fcrit: critical F-value for 1 - alpha confidence.
 * @alpha: nominal non-coverage, as decimal.
 * @iname: name of first parameter.
 * @jname: name of second parameter.
 *
 * Plots a 95% confidence ellipse for the parameter estimates
 * in @b with covariance @V.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int confidence_ellipse_plot (gretl_matrix *V, double *b, 
			     double tcrit, double Fcrit, double alpha, 
			     const char *iname, const char *jname)
{
    FILE *fp = NULL;
    double maxerr[2];
    double xcoeff[2];
    double ycoeff[2];
    double cval = 100 * (1 - alpha);
    gretl_matrix *e = NULL;
    gchar *title;
    int err = 0;

    maxerr[0] = tcrit * sqrt(gretl_matrix_get(V, 0, 0));
    maxerr[1] = tcrit * sqrt(gretl_matrix_get(V, 1, 1));

    err = gretl_invert_symmetric_matrix(V);
    if (err) {
	return err;
    }

    e = gretl_symmetric_matrix_eigenvals(V, 1, &err);
    if (err) {
	return err;
    }

    e->val[0] = sqrt(1.0 / e->val[0] * Fcrit);
    e->val[1] = sqrt(1.0 / e->val[1] * Fcrit);

    xcoeff[0] = e->val[0] * gretl_matrix_get(V, 0, 0);
    xcoeff[1] = e->val[1] * gretl_matrix_get(V, 0, 1);

    ycoeff[0] = e->val[0] * gretl_matrix_get(V, 1, 0);
    ycoeff[1] = e->val[1] * gretl_matrix_get(V, 1, 1);

    gretl_matrix_free(e);

    fp = get_plot_input_stream(PLOT_ELLIPSE, &err);
    if (err) {
	return err;
    }

    title = g_strdup_printf(_("%g%% confidence ellipse and %g%% marginal intervals"),
			    cval, cval);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);

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

#define MAKKONEN_POS 0

/* Probability of non-exceedance of the kth value in a set of n
   rank-ordered values.  See L. Makkonen, 'Bringing Closure to the
   Plotting Position Controversy', Communications in Statistics -
   Theory and Methods, vol 37, January 2008, for an argument in favor
   of using k / (n + 1); but also see many uses of (k - 1/2) / n in
   the literature.  
*/

static double plotpos (int k, int n)
{
#if MAKKONEN_POS
    return k / (n + 1.0);
#else
    return (k - 0.5) / n;
#endif
}

static double quantile_interp (const double *y, int n,
			       double ftarg)
{
    double f, ret = NADBL;
    int i;

    for (i=0; i<n; i++) {
	f = plotpos(i+1, n);
	if (f >= ftarg) {
	    if (f > ftarg && i > 0) {
		double f0 = plotpos(i, n);
		double d = (ftarg - f0) / (f - f0);

		ret = (1-d) * y[i-1] + d * y[i];
	    } else {
		ret = y[i];
	    }
	    break;
	}
    }

    return ret;
}

static int qq_plot_two_series (const int *list, const double **Z,
			       const DATAINFO *pdinfo)
{
    double *x = NULL;
    double *y = NULL;
    double f, qx, qy;
    FILE *fp = NULL;
    int vx = list[1];
    int vy = list[2];
    int nx = 10, ny = 10;
    int i, n, err = 0;

    x = gretl_sorted_series(vx, Z, pdinfo, OPT_NONE, &nx, &err);

    if (!err) {
	y = gretl_sorted_series(vy, Z, pdinfo, OPT_NONE, &ny, &err);
	if (err) {
	    free(x);
	    x = NULL;
	} 
    }

    if (!err) {
	/* take the smaller sample as basis */
	n = (nx > ny)? ny : nx;
    }

    if (!err) {
	fp = get_plot_input_stream(PLOT_QQ, &err);
    }

    if (err) {
	free(x);
	free(y);
	return err;
    }   

    fputs("set title \"Q-Q plot\"\n", fp);
    fputs("set datafile missing '?'\n", fp);
    fputs("set key top left\n", fp);
    fprintf(fp, "set xlabel \"%s\"\n", var_get_graph_name(pdinfo, vx));
    fprintf(fp, "set ylabel \"%s\"\n", var_get_graph_name(pdinfo, vy));
    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    fputs(" x notitle w lines\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	f = plotpos(i+1, n);

	if (nx == ny) {
	    qx = x[i];
	    qy = y[i];
	} else if (nx == n) {
	    qx = x[i];
	    qy = quantile_interp(y, ny, f);
	} else {
	    qx = quantile_interp(x, nx, f);
	    qy = y[i];
	}

	if (!na(qx) && !na(qy)) {
	    fprintf(fp, "%.12g %.12g\n", qx, qy); 
	}
    } 

    fputs("e\n", fp);
    
    gretl_pop_c_numeric_locale();

    free(x);
    free(y);

    fclose(fp);

    return gnuplot_make_graph();    
}

static int normal_qq_plot (const int *list, const double **Z, 
			   const DATAINFO *pdinfo, gretlopt opt)
{
    int zscores = 0;
    double ym = 0, ys = 1;
    double p, qx, qy;
    double *y = NULL;
    FILE *fp = NULL;
    int v = list[1];
    int i, n = 20;
    int err = 0;

    y = gretl_sorted_series(v, Z, pdinfo, OPT_NONE, &n, &err);

    if (!err && y[0] == y[n-1]) {
	gretl_errmsg_sprintf(_("%s is a constant"), pdinfo->varname[v]);
	err = E_DATA;
    }

    if (err) {
	return err;
    } 

    if (opt & OPT_Z) {
	/* standardize the data */
	zscores = 1;
    }

    if (!(opt & OPT_R)) {
	ym = gretl_mean(0, n-1, y);
	ys = gretl_stddev(0, n-1, y);

	if (zscores) {
	    /* standardize y */
	    for (i=0; i<n; i++) {
		y[i] = (y[i] - ym) / ys;
	    }
	}
    }

    fp = get_plot_input_stream(PLOT_QQ, &err);
    if (err) {
	free(y);
	return err;
    }

    fprintf(fp, "set title \"Q-Q plot for %s\"\n", 
	    var_get_graph_name(pdinfo, v));
    fputs("set datafile missing '?'\n", fp);
    fputs("set xlabel \"Normal quantiles\"\n", fp);

    if (opt & OPT_R) {
	fputs("set nokey\n", fp);
	fputs("plot \\\n", fp);
	fputs(" '-' using 1:2 notitle w points\n", fp);
    } else {
	fputs("set key top left\n", fp);
	fputs("plot \\\n", fp);
	fputs(" '-' using 1:2 notitle w points, \\\n", fp);
	fputs(" x title \"y = x\" w lines\n", fp);
    }
    
    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	p = plotpos(i+1, n);
	/* empirical quantile */
	qy = y[i];
	/* normal quantile */
	qx = normal_critval(1 - p);
	if (!na(qx) && !zscores && !(opt & OPT_R)) {
	    qx = ys * qx + ym;
	}
	if (!na(qx) && !na(qy)) {
	    fprintf(fp, "%.12g %.12g\n", qx, qy); 
	}
    } 

    fputs("e\n", fp);
    
    gretl_pop_c_numeric_locale();

    free(y);
    fclose(fp);

    return gnuplot_make_graph();
}

int qq_plot (const int *list, const double **Z, 
	     const DATAINFO *pdinfo, gretlopt opt)
{
    int err;

    if (list[0] == 1) {
	/* one series against normal */
	err = normal_qq_plot(list, Z, pdinfo, opt);
    } else if (list[0] == 2) {
	/* two empirical series */
	err = qq_plot_two_series(list, Z, pdinfo);
    } else {
	err = E_DATA;
    }

    return err;
}

/**
 * xy_plot_with_control:
 * @list: list of variables by ID number: Y, X, control.
 * @literal: extra gnuplot commands or %NULL.
 * @Z: data array.
 * @pdinfo: pointer to dataset information.
 * @opt: use %OPT_G for GUI graph.
 *
 * Constructs a scatterplot of modified Y against modified X,
 * where the modification consists in taking the residuals from
 * OLS regression of the variable in question on the control variable,
 * a la Frisch-Waugh-Lovell.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int xy_plot_with_control (const int *list, const char *literal,
			  const double **Z, const DATAINFO *pdinfo,
			  gretlopt opt)
{
    int vy = list[1], vx = list[2], vz = list[3];
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int mlist[3] = {2, 0, 0};
    MODEL mod;
    DATAINFO *ginfo = NULL;
    double **gZ = NULL;
    int s, t, T;
    int err = 0;

    if (list == NULL || list[0] != 3) {
	return E_DATA;
    }

    /* make sure we have a uniform sample for the two regressions */

    for (t=t1; t<=t2; t++) {
	if (na(Z[vy][t]) || na(Z[vx][t]) || na(Z[vz][t])) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=t2; t>=t1; t--) {
	if (na(Z[vy][t]) || na(Z[vx][t]) || na(Z[vz][t])) {
	    t2--;
	} else {
	    break;
	}
    }  

    /* count the maximum usable observations */
    T = t2 - t1 + 1;
    if (T < 3) {
	return E_DATA;
    }

    /* recount the usable observations, allowing for NAs in reduced
       sample */
    for (t=t1; t<=t2; t++) {
	if (na(Z[vy][t]) || na(Z[vx][t]) || na(Z[vz][t])) {
	    T--;
	}
    }
    if (T < 3) {
	return E_DATA;
    }    

    /* create temporary dataset */

    ginfo = create_auxiliary_dataset(&gZ, 4, T);
    if (ginfo == NULL) {
	return E_ALLOC;
    }

    strcpy(ginfo->varname[1], pdinfo->varname[vy]);
    strcpy(ginfo->varname[2], pdinfo->varname[vx]);

    sprintf(ginfo->varinfo[1]->display_name, _("adjusted %s"), pdinfo->varname[vy]);
    sprintf(ginfo->varinfo[2]->display_name, _("adjusted %s"), pdinfo->varname[vx]);
    
    s = 0;
    for (t=t1; t<=t2; t++) {
	if (!na(Z[vy][t]) && !na(Z[vx][t]) && !na(Z[vz][t])) {
	    gZ[1][s] = Z[vy][t];
	    gZ[2][s] = Z[vx][t];
	    gZ[3][s] = Z[vz][t];
	    s++;
	}
    }

    /* regress Y on Z and save the residuals */

    mlist[1] = 1;
    mlist[2] = 3;
    mod = lsq(mlist, &gZ, ginfo, OLS, OPT_A);
    err = mod.errcode;
    if (err) {
	clear_model(&mod);
	goto bailout;
    } else {
	for (t=0; t<mod.nobs; t++) {
	    gZ[1][t] = mod.uhat[t];
	}
	clear_model(&mod);
    }

    /* regress X on Z and save the residuals */

    mlist[1] = 2;    
    mod = lsq(mlist, &gZ, ginfo, OLS, OPT_A);
    err = mod.errcode;
    if (err) {
	clear_model(&mod);
	goto bailout;
    } else {
	for (t=0; t<mod.nobs; t++) {
	    gZ[2][t] = mod.uhat[t];
	}
	clear_model(&mod);
    }

    /* call for scatter of purged Y against purged X */

    mlist[1] = 1;
    mlist[2] = 2;
    err = gnuplot(mlist, literal, (const double **) gZ, 
		  ginfo, opt | OPT_C);

 bailout:

    /* trash the temporary dataset */
    destroy_dataset(gZ, ginfo);

    return err;
}

int is_auto_fit_string (const char *s)
{
    /* FIXME? */
    if (strstr(s, "automatic fit")) return 1;
    if (strstr(s, _("with least squares fit"))) return 1;
    return 0;
}

/** 
 * gnuplot_process_file:
 * @opt: may include %OPT_U for output to specified file.
 * @prn: gretl printing struct.
 *
 * Respond to the "gnuplot" command with %OPT_D, to specify
 * that input should be taken from a user-created gnuplot
 * command file.
 *
 * Returns: 0 on success, or if ignored; otherwise error code.
 */

int gnuplot_process_file (gretlopt opt, PRN *prn)
{
    const char *inname = get_optval_string(GNUPLOT, OPT_D);
    FILE *fp, *fq;
    int err = 0;

    if (inname == NULL && *inname == '\0') {
	return E_DATA;
    }

    fp = gretl_fopen(inname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fq = open_gp_stream(PLOT_USER, GPT_BATCH, &err);

    if (err) {
	fclose(fp);
    } else {
	char line[1024];

	while (fgets(line, sizeof line, fp)) {
	    fputs(line, fq);
	}

	fclose(fp);
	fclose(fq);

	err = gnuplot_make_graph();
    }

    return err;
}


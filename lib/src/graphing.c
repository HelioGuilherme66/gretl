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

#define ts_plot(g) ((g)->flags & GPT_TS)
#define use_impulses(g) ((g)->flags & GPT_IMPULSES)

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
    { PLOT_TYPE_MAX,       NULL }
};

static void graph_list_adjust_sample (int *list, 
				      gnuplot_info *ginfo,
				      const double **Z);
static void clear_gpinfo (gnuplot_info *gi);
    
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

static GptFlags get_gp_flags (gretlopt opt, int k, FitType *f)
{
    GptFlags flags = 0;

    if (opt & OPT_B) {
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

    if (k == 2 && !(opt & OPT_S)) {
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
    int ynum = gi->list[1];
    int dum = gi->list[3];
    int fn, t, i;

    fn = gi->t2 - gi->t1 + 1;

    gi->yvar1 = malloc(fn * sizeof *gi->yvar1);
    if (gi->yvar1 == NULL) {
	return 1;
    }

    gi->yvar2 = malloc(fn * sizeof *gi->yvar2);
    if (gi->yvar2 == NULL) {
	free(gi->yvar1);
	return 1;
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
   
int gnuplot_has_style_fill (void)
{
    /* ... and that it does style fill */
    return 1;
}

static int gnuplot_uses_datafile_missing (void)
{
    /* yup */
    return 1;
}

const char *gnuplot_label_front_string (void)
{
    /* ... and that it handles "front" for labels */
    return " front";
}

int gnuplot_has_latin5 (void)
{
    /* ... and that it supports ISO-8859-9 */
    return 1;
}

int gnuplot_has_cp1250 (void)
{
    /* ... and that it supports CP1250 */
    return 1;
}

int gnuplot_has_cp1254 (void)
{
    /* ... and that it doesn't support CP1254 */
    return 0;
}

int gnuplot_has_rgb (void)
{
    /* ... and that it supports rgb line-color specs */
    return 1;
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

#else

int gnuplot_has_ttf (int reset)
{
    static int err = -1; 

    /* try a range of ttf fonts that might plausibly be installed
       with X11 */

    if (err == -1 || reset) {
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

static int gnuplot_has_size (void)
{
    static int err = -1; 
    
    if (err == -1) {
	err = gnuplot_test_command("set term png size 640,480");
    }

    return !err;
}

int gnuplot_has_latin5 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set encoding iso_8859_9");
    }

    return !err;
}

int gnuplot_has_cp1250 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set encoding cp1250");
    }

    return !err;
}

int gnuplot_has_cp1254 (void)
{
    static int err = -1; 

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

    if (err == -1) {
	err = gnuplot_test_command("set encoding utf8");
    }

    return !err;    
}

int gnuplot_has_style_fill (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set style fill solid");
    }

    return !err;
}

static int gnuplot_uses_datafile_missing (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set datafile missing \"?\"");
    }

    return !err;
}

int gnuplot_has_rgb (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set style line 2 lc rgb \"#0000ff\"");
    }

    return !err;
}

const char *gnuplot_label_front_string (void)
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

#endif /* !WIN32 */

static int gnuplot_png_use_aa = 1;

void gnuplot_png_set_use_aa (int s)
{
    gnuplot_png_use_aa = s;
}

/* apparatus for handling plot colors */

enum {
    OLD_PNG_COLOR,  /* plain old "color" option */
    GD_PNG_COLOR,   /* pre-4.2.0 libgd-based PNG color spec */
    RGB_LINE_COLOR  /* per-line rgb settings */
};

static const gretlRGB default_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 }, /* full-intensity green is not very legible */
    { 0xbf, 0x25, 0xb2 },
    { 0x8f, 0xaa, 0xb3 },
    { 0xff, 0xa5, 0x00 },
    { 0x5f, 0x6b, 0x84 }  /* color for box fill */
};

static gretlRGB user_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 },
    { 0xbf, 0x25, 0xb2 },
    { 0x8f, 0xaa, 0xb3 },
    { 0xff, 0xa5, 0x00 },
    { 0x5f, 0x6b, 0x84 }
};

static void print_rgb_x (char *s, gretlRGB color)
{
    sprintf(s, "x%02x%02x%02x", color.r, color.g, color.b);
}

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
    if (i == BOXCOLOR) {
	user_color[BOXCOLOR] = default_color[BOXCOLOR];
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    user_color[i] = default_color[i];
	}
    }
}

static int split_fontname (const char *s, char *name, int *psz)
{
    int i, k = 0, n = strlen(s);
    int nf = 0;

    for (i=n-1; i>0; i--) {
	if (isdigit(s[i])) k++;
	else break;
    }

    if (k > 0) {
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
write_gnuplot_font_string (char *fstr, const char *grfont, PlotType ptype,
			   int pngterm)
{
    if (pngterm == GP_PNG_CAIRO) {
	char fname[128];
	int nf, fsize = 0;

	nf = split_fontname(grfont, fname, &fsize);
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

/* we need this only if we don't have per-line rgb
   settings, which are in gnuplot 4.2 and higher */

static char *make_png_colorspec (char *targ, int ptype)
{
    char cstr[8];
    int i;

    /* background etc. */
    strcpy(targ, " xffffff x000000 x202020");

    if (frequency_plot_code(ptype)) {
	strcat(targ, " ");
	print_rgb_x(cstr, user_color[BOXCOLOR]);
	strcat(targ, cstr);
	strcat(targ, " x000000");
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    strcat(targ, " ");
	    print_rgb_x(cstr, user_color[i]);
	    strcat(targ, cstr);
	}
    }

    return targ;
}

/* we use this mechanism with gnuplot 4.2 and higher */

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
    }

    fputs("set style increment user\n", fp);
}

/* end colors apparatus */

/* Get gnuplot to print the dimensions of a PNG plot, in terms
   of both pixels and data bounds, if gnuplot supports this.
*/

void print_plot_bounding_box_request (FILE *fp)
{
    fprintf(fp, "set print '%sgretltmp.png.bounds'\n", gretl_dot_dir());
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
    char color_string[64];
    int gpcolors, gpttf = 1, gpsize = 1;
    int pngterm = 0;
    const char *grfont = NULL;

    *font_string = 0;
    *size_string = 0;
    *color_string = 0;

    pngterm = gnuplot_png_terminal();

#ifdef WIN32
    gpcolors = RGB_LINE_COLOR;
#else
    if (gnuplot_has_rgb()) {
	gpcolors = RGB_LINE_COLOR;
    } else if (pngterm == GP_PNG_OLD) {
	gpcolors = OLD_PNG_COLOR;
    } else {
	gpcolors = GD_PNG_COLOR;
    }

    gpttf = gnuplot_has_ttf(0);
    gpsize = gnuplot_has_size();
#endif

    if (pngterm == GP_PNG_GD2 && gnuplot_png_use_aa) {
	strcpy(truecolor_string, " truecolor");
    }    

    /* plot font setup */
    if (gpttf) {
	grfont = gretl_png_font();
	if (*grfont == 0) {
	    grfont = getenv("GRETL_PNG_GRAPH_FONT");
	}
	if (grfont != NULL && *grfont != 0) {
	    write_gnuplot_font_string(font_string, grfont, ptype, pngterm);
	}
    } 

#ifndef WIN32
    if (!gpttf) {
	write_old_gnuplot_font_string(font_string, ptype);
    }
#endif

    /* plot color setup */
    if (gpcolors == GD_PNG_COLOR) {
	make_png_colorspec(color_string, ptype);
    } else if (gpcolors == OLD_PNG_COLOR) {
	strcpy(color_string, " color");
    } else {
	/* handled via styles */
	*color_string = '\0';
    }

    if (gpsize) {
	if (flags & GPT_LETTERBOX) {
	    strcpy(size_string, " size 680,400");
	} else if (ptype == PLOT_VAR_ROOTS) {
	    strcpy(size_string, " size 480,480");
	}
    }

    if (pngterm == GP_PNG_CAIRO) {
	sprintf(png_term_line, "set term pngcairo%s%s",
		font_string, size_string);
	strcat(png_term_line, "\nset encoding utf8"); /* FIXME? */
    } else {
	sprintf(png_term_line, "set term png%s%s%s%s",
		truecolor_string, font_string, size_string, 
		color_string);
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

    if (split_fontname(pngfont, name, &pt) == 2) {
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

    return ret;
}

static int real_gnuplot_init (PlotType ptype, int flags, FILE **fpp)
{
    int gui = gretl_in_gui_mode();
    char plotfile[FILENAME_MAX] = {0};

    /* 'gnuplot_path' is file-scope static var */
    if (*gnuplot_path == 0) {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    if (gui) {
	sprintf(plotfile, "%sgpttmp.XXXXXX", gretl_dot_dir());
	if (mktemp(plotfile) == NULL) {
	    return E_FOPEN;
	}
    } else {
	sprintf(plotfile, "%sgpttmp.plt", gretl_dot_dir());
    }

    set_gretl_plotfile(plotfile);

    *fpp = gretl_fopen(plotfile, "w");
    if (*fpp == NULL) {
	fprintf(stderr, "gnuplot_init: couldn't write to %s\n", plotfile);
	return E_FOPEN;
    } 

    if (gui) {
	fprintf(*fpp, "%s\n", get_gretl_png_term_line(ptype, flags));
	fprintf(*fpp, "set output '%sgretltmp.png'\n", gretl_dot_dir());
    }

    write_plot_type_string(ptype, *fpp);

    if (gnuplot_has_rgb()) {
	write_plot_line_styles(ptype, *fpp);
    }

#if GP_DEBUG
    fprintf(stderr, "gnuplot_init: set plotfile = '%s'\n", 
	    plotfile);
#endif

    return 0;
}

/**
 * gnuplot_init:
 * @ptype: indication of the sort of plot to be made.
 * @fpp: pointer to stream to be opened.
 *
 * If we're in GUI mode: writes a unique temporary filename into
 * the internal variable #gretl_plotfile; opens plotfile for writing 
 * as @fpp; and writes initial lines into the output file to select 
 * the PNG terminal type, and direct gnuplot's output to a temporary
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
    return real_gnuplot_init(ptype, 0, fpp);
}

static int gretl_plot_count;

void reset_plot_count (void)
{
    gretl_plot_count = 0;
}

/* initialization for gnuplot output file in batch mode: this
   is wanted when drawing batch boxplots */

FILE *gnuplot_batch_init (int *err)
{
    const char *optname = get_optval_string(BXPLOT, OPT_U);
    char fname[FILENAME_MAX];
    FILE *fp = NULL;

    if (optname != NULL && *optname != '\0') {
	/* user gave --output=<filename> */
	strcpy(fname, optname);
	gretl_maybe_prepend_dir(fname);
	fp = gretl_fopen(fname, "w");
    } else {
	sprintf(fname, "%sgpttmp%02d.plt", gretl_work_dir(),
		++gretl_plot_count);
	fp = gretl_fopen(fname, "w");
    }

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
    }

    return fp;
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
    } else if (tt == GP_TERM_FIG) {
	tstr = "set term fig";
    }

    if (tstr != NULL) {
	fprintf(fp, "%s\n", tstr);
    }
}

/* Produce formatted graph output in batch mode: at present
   we only recognize EPS, PDF, PNG and XFig */

static int make_graph_special (const char *fname, int fmt)
{
    char plotcmd[MAXLEN];
    char tmp[FILENAME_MAX];
    char line[1024];
    FILE *fp, *fq;
    int err;

    if (fmt == GP_TERM_PDF && gnuplot_pdf_terminal() == GP_PDF_NONE) {
	strcpy(gretl_errmsg, _("Gnuplot does not support PDF output "
			       "on this system"));
	return E_EXTERNAL;
    }

    strcpy(tmp, fname);
    strcpy(strrchr(tmp, '.'), ".gp");

    /* file to contain augmented gnuplot input */
    fp = gretl_fopen(tmp, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    /* existing plot file */
    fq = gretl_fopen(fname, "r");
    if (fq == NULL) {
	fclose(fp);
	return E_FOPEN;
    }    

    /* write terminal/output lines */
    print_term_string(fmt, fp);
    fprintf(fp, "set output '%s'\n", fname);

    /* transcribe original plot content */
    while (fgets(line, sizeof line, fq)) {
	fputs(line, fp);
    }

    fclose(fq);
    fclose(fp);

#ifdef WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", gretl_gnuplot_path(), tmp);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    sprintf(plotcmd, "%s \"%s\"", gretl_gnuplot_path(), tmp);
    err = gretl_spawn(plotcmd);
#endif 

    /* remove the temporary input file */
    if (err) {
	fprintf(stderr, "err = %d: bad file is '%s'\n", err, tmp);
    } else {
	remove(tmp);
    }

    return err;
}

int specified_gp_output_format (void)
{
    const char *fname = gretl_plotfile();

    if (has_suffix(fname, ".eps")) {
	return GP_TERM_EPS;
    } else if (has_suffix(fname, ".pdf")) {
	return GP_TERM_PDF;
    } else if (has_suffix(fname, ".png")) {
	return GP_TERM_PNG;
    } else if (has_suffix(fname, ".fig")) {
	return GP_TERM_FIG;
    } else {
	return GP_TERM_NONE;
    }
}

/**
 * gnuplot_make_graph:
 *
 * Executes gnuplot to produce a grapg file in PNG format.
 *
 * Returns: the return value from the system command.
 */

int gnuplot_make_graph (void)
{
    const char *plotfile = gretl_plotfile();
    char plotcmd[MAXLEN];
    int fmt, err = 0;

    fmt = specified_gp_output_format();

    if (fmt != GP_TERM_NONE) {
	return make_graph_special(plotfile, fmt);
    }

    if (gretl_in_gui_mode() && gnuplot_has_bbox()) {
	do_plot_bounding_box();
    }

#ifdef WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", gretl_gnuplot_path(), plotfile);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    sprintf(plotcmd, "%s%s \"%s\"", gretl_gnuplot_path(), 
	    (gretl_in_gui_mode())? "" : " -persist", plotfile);
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

static int
get_gnuplot_output_file (FILE **fpp, GptFlags flags, int code)
{
    const char *plotfile = gretl_plotfile();
    int err = 0;

    *fpp = NULL;

    if ((flags & GPT_FILE) && *plotfile != '\0') {
	*fpp = gretl_fopen(plotfile, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	}
    } else if (flags & GPT_BATCH) {
	const char *optname = get_optval_string(GNUPLOT, OPT_U);
	char fname[FILENAME_MAX];

	if (optname != NULL && *optname != '\0') {
	    /* user gave --output=<filename> */
	    strcpy(fname, optname);
	    gretl_maybe_prepend_dir(fname);
	} else {
	    sprintf(fname, "%sgpttmp%02d.plt", gretl_work_dir(), 
		    ++gretl_plot_count);
	}
	*fpp = gretl_fopen(fname, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	} else {
	    set_gretl_plotfile(fname);
	}
    } else {
	/* note: gnuplot_init is not used in batch mode */
	err = real_gnuplot_init(code, flags, fpp);
    }

#if GPDEBUG
    fprintf(stderr, "get_gnuplot_output_file: '%s'\n", gretl_plotfile());
#endif

    return err;
}

static int 
loess_plot (gnuplot_info *gi, const double **Z, const DATAINFO *pdinfo)
{
    gretl_matrix *y = NULL;
    gretl_matrix *x = NULL;
    gretl_matrix *yh = NULL;
    int yno = gi->list[1];
    int xno = gi->list[2];
    const char *s1, *s2;
    FILE *fp = NULL;
    char title[96];
    int t, T, d = 1;
    double q = 0.5;
    int err;

    graph_list_adjust_sample(gi->list, gi, Z);
    if (gi->t1 == gi->t2 || gi->list[0] != 2) {
	return GRAPH_NO_DATA;
    }

    if (get_gnuplot_output_file(&fp, gi->flags, PLOT_REGULAR)) {
	return E_FOPEN;
    } 

    err = gretl_plotfit_matrices(yno, xno, PLOT_FIT_LOESS, Z,
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

    s1 = var_get_graph_name(pdinfo, yno);
    s2 = var_get_graph_name(pdinfo, xno);

    sprintf(title, _("%s versus %s (with loess fit)"), s1, s2);
    print_keypos_string(GP_KEY_LEFT_TOP, fp);
    fprintf(fp, "set title \"%s\"\n", title);
    print_axis_label('y', s1, fp);
    print_axis_label('x', s2, fp);
    print_auto_fit_string(PLOT_FIT_LOESS, fp);

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

    if (!(gi->flags & GPT_BATCH)) {
	err = gnuplot_make_graph();
    }    

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
    int yno = gi->list[1];
    int xno = gi->list[2];
    double s2, *ps2 = NULL;
    FitType f = gi->fit;
    char title[72];
    int k, err;

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
	f = PLOT_FIT_OLS;
    }

    /* columns needed in X matrix */
    k = (f == PLOT_FIT_QUADRATIC)? 3 : 2;

    err = gretl_plotfit_matrices(yno, xno, f, Z,
				 pdinfo->t1, pdinfo->t2,
				 &y, &X);

    if (!err) {
	b = gretl_column_vector_alloc(k);
	if (b == NULL) {
	    err = E_ALLOC;
	}
    }
    
    if (!err) {
	err = gretl_matrix_ols(y, X, b, V, NULL, ps2);
    }

    if (!err) {
	double *c = b->val;

	if (gi->fit == PLOT_FIT_NONE) {
	    /* the "automatic" case */
	    double pv, v = gretl_matrix_get(V, 1, 1);
	    int T = gretl_vector_get_length(y);

	    pv = student_pvalue_2(T - k, c[1] / sqrt(v));
	    /* show the line if the p-value for the slope coeff is
	       less than 0.1, otherwise discard it */
	    if (pv < .10) {
		sprintf(title, "Y = %#.3g %c %#.3gX", b->val[0],
			(c[1] > 0)? '+' : '-', fabs(c[1]));
		gretl_push_c_numeric_locale();
		sprintf(targ, "%.10g + %.10g*x title '%s' w lines\n", 
			c[0], c[1], title);
		gretl_pop_c_numeric_locale();
		gi->fit = PLOT_FIT_OLS;
	    }
	} else if (gi->fit == PLOT_FIT_OLS) {
	    sprintf(title, "Y = %#.3g %c %#.3gX", c[0],
		    (c[1] > 0)? '+' : '-', fabs(c[1]));
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%.10g + %.10g*x title '%s' w lines\n", 
		    c[0], c[1], title);
	    gretl_pop_c_numeric_locale();
	} else if (gi->fit == PLOT_FIT_INVERSE) {
	    sprintf(title, "Y = %#.3g %c %#.3g/X", c[0],
		    (c[1] > 0)? '+' : '-', fabs(c[1]));
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%.10g + %.10g/x title '%s' w lines\n", 
		    c[0], c[1], title);
	    gretl_pop_c_numeric_locale();
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    sprintf(title, "Y = %#.3g %c %#.3gX %c %#.3gX^2", c[0],
		    (c[1] > 0)? '+' : '-', fabs(c[1]),
		    (c[2] > 0)? '+' : '-', fabs(c[2])),
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%.10g + %.10g*x + %.10g*x**2 title '%s' w lines\n", 
		    c[0], c[1], c[2], title);
	    gretl_pop_c_numeric_locale();
	}

	if (gi->fit != PLOT_FIT_NONE) {
	    gi->flags |= GPT_AUTO_FIT;
	}
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

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

static void
print_gp_data (gnuplot_info *gi, const double **Z, 
	       const DATAINFO *pdinfo)
{
    int n = gi->t2 - gi->t1 + 1;
    double offset = 0.0;
    int datlist[3];
    int ynum = 2;
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
		} else if (dataset_is_time_series(pdinfo)) {
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
    if (flags & GPT_FILE) {
	fprintf(stderr, " GPT_FILE\n");
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
    if (gnuplot_uses_datafile_missing()) {
	fputs("set datafile missing \"?\"\n", fp);
    } else {
	fputs("set missing \"?\"\n", fp);
    }
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

    maxtics = pdinfo->paninfo->unit[gi->t2] - 
	pdinfo->paninfo->unit[gi->t1] + 1;

    ntics = maxtics;
    while (ntics > 14) {
	ntics /= 1.5;
    }

    ticskip = maxtics / ceil(ntics);

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
    double yrs = (gi->t2 - gi->t1 + 1.0) / (pdinfo->pd * 52.0);

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
    int i, err;

    gretl_error_clear();

    err = incompatible_options(opt, OPT_T | OPT_I | OPT_L | OPT_Q | OPT_N);
    if (err) {
	return err;
    }

#if GP_DEBUG
    printlist(plotlist, "gnuplot: plotlist");
#endif

    err = gpinfo_init(&gi, opt, plotlist, literal, 
		      pdinfo->t1, pdinfo->t2);
    if (err) {
	goto bailout;
    }

    if (gi.fit == PLOT_FIT_LOESS) {
	return loess_plot(&gi, Z, pdinfo);
    }

    if (gi.list[0] > MAX_LETTERBOX_LINES + 1) {
	many = 1;
    }

    err = maybe_add_plotx(&gi, pdinfo);
    if (err) {
	goto bailout;
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
	if (list[0] != 3 || factorized_vars(&gi, Z)) {
	    err = E_DATA;
	    goto bailout;
	}
    } 

    /* special tics for time series plots */
    if (gi.flags & GPT_TS) {
	if (many) {
	    pprintf(prn, "# multiple timeseries %d\n", pdinfo->pd);
	} else {
	    int gpsize = 1;

#ifndef WIN32
	    gpsize = gnuplot_has_size();
#endif
	    pprintf(prn, "# timeseries %d", pdinfo->pd);
	    if (gpsize) {
		gi.flags |= GPT_LETTERBOX;
		pputs(prn, " (letterbox)\n");
	    } else {
		pputc(prn, '\n');
	    }
	} 
	if (pdinfo->pd == 4 && (gi.t2 - gi.t1) / 4 < 8) {
	    pputs(prn, "set xtics nomirror 0,1\n"); 
	    pputs(prn, "set mxtics 4\n");
	} else if (pdinfo->pd == 12 && (gi.t2 - gi.t1) / 12 < 8) {
	    pputs(prn, "set xtics nomirror 0,1\n"); 
	    pputs(prn, "set mxtics 12\n");
	} else if (dated_daily_data(pdinfo)) {
	    make_calendar_tics(pdinfo, &gi, prn);
	} else if (panel_plot(pdinfo, gi.t1, gi.t2)) {
	    make_panel_unit_tics(pdinfo, &gi, prn);
	    strcpy(xlabel, _("time series by group"));
	}
    }

    /* open file and dump the prn into it: we delaying writing
       the file header till we know a bit more about the plot
    */
    if (get_gnuplot_output_file(&fp, gi.flags, PLOT_REGULAR)) {
	err = E_FOPEN;
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

    if (!(gi.flags & GPT_BATCH) || specified_gp_output_format()) {
	err = gnuplot_make_graph();
    }

 bailout:

    clear_gpinfo(&gi);

    return err;
}

/**
 * multi_scatters:
 * @list: list of variables to plot, by ID number.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: can include %OPT_L to use lines.
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

    if (opt & OPT_B) {
	flags = GPT_BATCH;
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

    if (get_gnuplot_output_file(&fp, flags, PLOT_MULTI_SCATTER)) {
	return E_FOPEN;
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

    if (!(flags & GPT_BATCH)) {
	err = gnuplot_make_graph();
    }

    free(plotlist);

    return err;
}

static int get_3d_output_file (FILE **fpp)
{
    char fname[MAXLEN];
    int err = 0;

    sprintf(fname, "%sgpttmp.plt", gretl_dot_dir());
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
		   gretlopt opt, char *surface)
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
	(snedecor_cdf_comp(smod.dfn, smod.dfd, smod.fstt) < .10 || (opt & OPT_F))) {
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

#ifndef WIN32
    if (gnuplot_has_x11()) {
	/* wxt is too slow/jerky? */
	fputs("set term x11\n", fq);
    }
#endif

    gretl_push_c_numeric_locale();

    maybe_add_surface(list, pZ, pdinfo, opt, surface);
    
    print_axis_label('x', var_get_graph_name(pdinfo, list[2]), fq);
    print_axis_label('y', var_get_graph_name(pdinfo, list[1]), fq);
    print_axis_label('z', var_get_graph_name(pdinfo, list[3]), fq);

    gnuplot_missval_string(fq);

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
#ifdef ENABLE_NLS
    char test[8];

    gretl_pop_c_numeric_locale();
    sprintf(test, "%g", 0.5);
    if (strchr(test, ',')) {
	dcomma = 1;
    }
#endif

    if (dist == D_NORMAL) {
	sprintf(s, "N(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    } else if (dist == D_GAMMA) {
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
    int err;

    if (K == 0) {
	return E_DATA;
    }

    if (K == 1) {
	sprintf(gretl_errmsg, _("'%s' is a constant"), freq->varname);
	return E_DATA;
    }

    if (dist == D_NORMAL) {
	plottype = PLOT_FREQ_NORMAL;
    } else if (dist == D_GAMMA) {
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
		fprintf(fp, "set label \"%s:\" at graph .03, graph .97%s\n",
			_("Test statistic for normality"),
			gnuplot_label_front_string());
		print_freq_test_label(label, _("Chi-squared(2) = %.3f pvalue = %.5f"), 
				      freq->test, chisq_cdf_comp(2, freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93%s\n", 
			label, gnuplot_label_front_string());
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
		fprintf(fp, "set label '%s:' at graph .03, graph .97%s\n",
			_("Test statistic for gamma"),
			gnuplot_label_front_string());
		print_freq_test_label(label, _("z = %.3f pvalue = %.5f"), 
				      freq->test, normal_pvalue_2(freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93%s\n", 
			label, gnuplot_label_front_string());
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
	if (gnuplot_has_style_fill()) {
	    fputs("set style fill solid 0.6\n", fp);
	}
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
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int use_fill = 0, use_lines = 0;
    int do_errs = (maxerr != NULL);
    char cistr[64];
    int t2 = fr->t2;
    int t1, yhmin;
    int t, n, err;

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

    if ((err = gnuplot_init(PLOT_FORECAST, &fp))) {
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
    int i, j, err;

    if (tau == NULL || B == NULL) {
	return E_DATA;
    }

    ntau = gretl_vector_get_length(tau);
    if (ntau == 0) {
	return E_DATA;
    }

    if ((err = gnuplot_init(PLOT_RQ_TAU, &fp))) {
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
    int t, err;

    h = gretl_model_get_data(pmod, "garch_h");
    if (h == NULL) {
	return E_DATA;
    }

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if ((err = gnuplot_init(PLOT_GARCH, &fp))) {
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

int rmplot (const int *list, const double **Z, DATAINFO *pdinfo, PRN *prn)
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

    err = gnuplot_init(PLOT_PANEL, &fp);
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
				 const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    int confint = 0;
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int t, err;

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

	fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines, \\\n", 
		_("point estimate"));
	sprintf(title, _("%g and %g quantiles"), ql, qh);
	fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n", title);
    } else {
	fputs("plot \\\n'-' using 1:2 w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    for (t=0; t<periods; t++) {
	fprintf(fp, "%d %.10g\n", t+1, gretl_matrix_get(resp, t, 0));
    }
    fputs("e\n", fp);

    if (confint) {
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

    resp = gretl_VAR_get_impulse_response(var, 1, 1, periods, alpha, Z, pdinfo);
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

	    if (confint) {
		fputs("plot \\\n'-' using 1:2 notitle w lines, \\\n", fp); 
		fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
	    } else {
		fputs("plot \\\n'-' using 1:2 w lines\n", fp);
	    }

	    for (t=0; t<periods; t++) {
		fprintf(fp, "%d %.10g\n", t+1, gretl_matrix_get(resp, t, 0));
	    }
	    fputs("e\n", fp);

	    if (confint) {
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
    int i, v, t, t1, err;

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

    err = gnuplot_init(PLOT_REGULAR, &fp);
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

    err = gnuplot_init(PLOT_MULTI_SCATTER, &fp);
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

    err = gnuplot_init(PLOT_VAR_ROOTS, &fp);
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

    err = gnuplot_init(PLOT_ELLIPSE, &fp);
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


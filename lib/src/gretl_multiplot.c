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

/* Support for the commands "gpbuild" and "gridplot", which produce
   plots using gnuplot's "multiplot" apparatus. This facility was
   added in May/June 2023.
*/

#include "libgretl.h"
#include "uservar.h"
#include "usermat.h"
#include "gretl_multiplot.h"

#define GRID_DEBUG 0

typedef struct {
    gretlopt flag;
    int *target;
    int min;
    int max;
    int def;
} mplot_option;

static gretl_array *mp_array;
static gretl_matrix *mp_layout;
static const char *array_name;
static int prior_array;
static int mp_fontsize = 10;
static int mp_width = 800;
static int mp_height = 600;
static int mp_rows;
static int mp_cols;
static int mp_collecting;

int gretl_multiplot_collecting (void)
{
    return mp_array != NULL && mp_collecting;
}

static const mplot_option mp_options[] = {
    { OPT_F, &mp_fontsize, 4, 24, 10 },
    { OPT_W, &mp_width,    200, 2048, 800 },
    { OPT_H, &mp_height,   200, 2048, 600 },
    { OPT_R, &mp_rows,     1, 12, 0 },
    { OPT_C, &mp_cols,     1, 12, 0 }
};

static int n_mp_options = G_N_ELEMENTS(mp_options);

static int set_multiplot_sizes (gretlopt opt)
{
    const mplot_option *mpo;
    int i, k, err = 0;

    for (i=0; i<n_mp_options && !err; i++) {
	mpo = &mp_options[i];
	if (opt & mpo->flag) {
	    k = get_optval_int(GRIDPLOT, mpo->flag, &err);
	    if (!err && (k < mpo->min || k > mpo->max)) {
		gretl_errmsg_set("gridplot: out-of-bounds option value");
		err = E_INVARG;
	    }
	    if (!err) {
		*mpo->target = k;
	    }
	}
    }

    return err;
}

static int initialize_mp_array (const char *param, DATASET *dset)
{
    int err = 0;

    mp_array = get_strings_array_by_name(param);
    if (mp_array != NULL) {
	prior_array = 1;
	gretl_array_void_content(mp_array);
    } else {
	gchar *s = g_strdup_printf("%s = array(0)", param);

	err = generate(s, dset, GRETL_TYPE_STRINGS, OPT_NONE, NULL);
	if (!err) {
	    mp_array = get_strings_array_by_name(param);
	}
	g_free(s);
    }

    if (mp_array != NULL) {
	array_name = param;
	mp_collecting = 1;
    }

    return err;
}

/* respond to a validated case of the --layout=matrix option */

static void set_mp_layout (const gretl_matrix *m)
{
    if (mp_layout != NULL) {
	gretl_matrix_free(mp_layout);
    }
    mp_layout = gretl_matrix_copy(m);
    if (mp_layout != NULL) {
	if (mp_rows != m->rows || mp_cols != m->cols) {
	    mp_rows = m->rows;
	    mp_cols = m->cols;
	}
    }
}

static void set_multiplot_defaults (void)
{
    int i;

    if (mp_layout != NULL) {
	gretl_matrix_free(mp_layout);
	mp_layout = NULL;
    }
    for (i=0; i<n_mp_options; i++) {
	*(mp_options[i].target) = mp_options[i].def;
    }
}

void gretl_multiplot_clear (int err)
{
    if (mp_array != NULL) {
	if (err) {
	    if (prior_array) {
		/* got a pre-existing uservar array */
		gretl_array_void_content(mp_array);
	    } else {
		/* scrub temporary uservar */
		user_var_delete_by_name(array_name, NULL);
	    }
	}
        mp_array = NULL;
    }

    array_name = NULL;
    prior_array = 0;
    mp_collecting = 0;
    set_multiplot_defaults();
}

/* called from interact.c: process_command_error() */

void gretl_multiplot_destroy (void)
{
    gretl_multiplot_clear(1);
}

/* This responds to the starting command for a "gpbuild" block */

int gretl_multiplot_start (const char *param, gretlopt opt,
			   DATASET *dset)
{
    int err = 0;

    if (mp_array == NULL) {
	err = initialize_mp_array(param, dset);
    } else {
        gretl_errmsg_set(_("gridplot: cannot be nested"));
        err = E_DATA;
    }

    return err;
}

static int invalid_mp_error (int ci)
{
    gretl_errmsg_sprintf(_("%s: invalid (multiplot) plot specification"),
			 gretl_command_word(ci));
    return E_INVARG;
}

/* Append a plot specification to the @multiplot array */

int gretl_multiplot_add_plot (gchar *buf)
{
    int err = 0;

    if (mp_array != NULL && buf != NULL) {
        if (strstr(buf, "set multiplot")) {
            err = invalid_mp_error(GPBUILD);
        } else {
	    gretl_array_append_string(mp_array, buf, 1);
        }
    } else {
	gretl_errmsg_set("gretl_multiplot_add_plot: failed");
        err = E_DATA;
    }

#if GRID_DEBUG
    fprintf(stderr, "gretl_multiplot_add_plot, err = %d\n", err);
#endif

    return err;
}

static int multiplot_set_grid (int n)
{
    int err = 0;

#if GRID_DEBUG
    fprintf(stderr, "multiplot_set_grid: n=%d, prior size %d x %d\n",
	    n, mp_rows, mp_cols);
#endif

    if (mp_rows == 0 && mp_cols == 0) {
	/* fully automatic grid */
	mp_rows = ceil(sqrt((double) n));
	mp_cols = ceil((double) n / mp_rows);
    } else if (mp_rows == 0) {
	/* automatic rows */
        if (mp_cols > n) {
            mp_cols = n;
            mp_rows = 1;
        } else {
            mp_rows = ceil((double) n / mp_cols);
        }
    } else if (mp_cols == 0) {
	/* automatic cols */
        if (mp_rows > n) {
            mp_rows = n;
            mp_cols = 1;
        } else {
            mp_cols = ceil((double) n / mp_rows);
        }
    } else if (mp_rows * mp_cols < n) {
	gretl_errmsg_sprintf("Specified grid (%d by %d) is too small "
			     "for %d sub-plots", mp_rows, mp_cols, n);
	err = E_INVARG;
    } else if (mp_rows * mp_cols > n) {
	int ar = ceil(sqrt((double) n));
	int ac = ceil((double) n / ar);

	if (mp_rows * mp_cols > ar * ac) {
	    gretl_errmsg_sprintf("Specified grid (%d by %d) is too big "
				 "for %d sub-plots", mp_rows, mp_cols, n);
	    err = E_INVARG;
	}
    }

#if GRID_DEBUG
    fprintf(stderr, "multiplot_set_grid: set %d x %d\n",  mp_rows, mp_cols);
#endif

    return err;
}

static int get_subplot_index (int i, int j)
{
    return (int) gretl_matrix_get(mp_layout, i, j) - 1;
}

/* write a layout matrix into a plot file, in a form that can
   be easily reconstituted via generate_matrix()
*/

static void output_layout_matrix (gretl_matrix *m, FILE *fp)
{
    double mij;
    int i, j;

    fputs("layout={", fp);
    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    mij = gretl_matrix_get(m, i, j);
	    fprintf(fp, "%d", (int) mij);
	    if (j < m->cols-1) {
		fputc(',', fp);
	    }
	}
	if (i < m->rows-1) {
	    fputc(';', fp);
	}
    }
    fputs("}\n", fp);
}

static void write_mp_spec_comment (int np, FILE *fp)
{
    fprintf(fp, "# grid_params: plots=%d, fontsize=%d, width=%d, height=%d, ",
	    np, mp_fontsize, mp_width, mp_height);
    if (mp_layout != NULL) {
	output_layout_matrix(mp_layout, fp);
    } else {
	fprintf(fp, "rows=%d, cols=%d\n", mp_rows, mp_cols);
    }
}

/* Write subplot buffer @buf to file, stripping out any
   "set term..." statements.
*/

static void filter_subplot (const char *buf, FILE *fp)
{
    size_t sz = 2048;
    char *line = calloc(sz, 1);

    bufgets_init(buf);
    while (safe_bufgets(&line, &sz, buf)) {
	if (strncmp(line, "set term", 8)) {
	    fputs(line, fp);
	}
    }
    bufgets_finalize(buf);
    free(line);
}

/* Write a multiplot specification to file, drawing on the
   strings array @a.

   @np indicates the number of included plots and @maxp the
   total number of subplots available (the length of the
   relevant array).
*/

static int output_multiplot_script (gretl_array *a,
				    int np, int maxp)
{
    const char *buf;
    int i, j, k, p;
    int err = 0;
    FILE *fp;

    /* insure against segfault */
    if (a == NULL) {
        fprintf(stderr, "output_multiplot_script: internal error!\n");
        return E_DATA;
    }

    fp = open_plot_input_file(PLOT_GRIDPLOT, 0, &err);
    if (err) {
        return err;
    }

    fputs("# literal lines = 1\n", fp);
    write_mp_spec_comment(np, fp);

    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n", mp_rows, mp_cols);
    gretl_push_c_numeric_locale();

    k = -1;
    p = 0;
    for (i=0; i<mp_rows; i++) {
	for (j=0; j<mp_cols; j++) {
	    if (mp_layout != NULL) {
		k = get_subplot_index(i, j);
	    } else {
		k++;
	    }
	    if (k < 0) {
		fputs("set multiplot next\n", fp);
	    } else {
		buf = NULL;
		if (k < maxp) {
		    buf = gretl_array_get_data(a, k);
		}
		if (buf != NULL) {
		    if (i + j > 0) {
			fputs("reset\n", fp);
		    }
		    fprintf(fp, "# subplot %d\n", ++p);
		    if (strstr(buf, "set term")) {
			filter_subplot(buf, fp);
		    } else {
			fputs(buf, fp);
		    }
		}
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);
    err = finalize_plot_input_file(fp);

    return err;
}

static int maybe_set_mp_layout (int *np)
{
    const char *s = get_optval_string(GRIDPLOT, OPT_L);
    const gretl_matrix *m = NULL;
    int err = 0;

    if (s != NULL) {
	m = get_matrix_by_name(s);
    }
    if (m != NULL) {
	int i, n = m->rows * m->cols;
        int nonzero = 0;
	double mi;

	for (i=0; i<n; i++) {
	    mi = m->val[i];
	    if (na(mi) || mi != floor(mi) ||
		mi < 0 || mi > *np) {
		gretl_errmsg_set(_("Invalid layout specification"));
		err = E_INVARG;
		break;
	    } else if (mi != 0) {
                nonzero++;
            }
	}
	if (!err) {
	    *np = nonzero;
	}
	if (!err) {
#if GRID_DEBUG
	    gretl_matrix_print(m, "m, in maybe_set_mp_layout");
#endif
	    set_mp_layout(m);
	}
    }

    return err;
}

/* respond to "end gpbuild" */

int gretl_multiplot_finalize (gretlopt opt)
{
    int err = 0;

    if (mp_array == NULL) {
	gretl_errmsg_set("end gpbuild: building not started");
	err = E_DATA;
    } else {
	gretl_multiplot_clear(0);
    }

    return err;
}

static int retrieve_plots_array (const char *argname,
				 gretl_array **pa,
				 int *pnp)
{
    gretl_array *a = NULL;
    int msg_set = 0;
    int err = 0;

    a = get_array_by_name(argname);
    if (a == NULL) {
	err = E_INVARG;
    } else {
	int i, n = gretl_array_get_length(a);
	const char *buf;

	/* check for embedded multiplots */
	for (i=0; i<n; i++) {
	    buf = gretl_array_get_data(a, i);
	    if (strstr(buf, "set multiplot")) {
		err = invalid_mp_error(GRIDPLOT);
		msg_set = 1;
		break;
	    }
	}
	if (!err) {
	    *pa = a;
	    *pnp = n;
	}
    }

    if (err && !msg_set) {
	gretl_errmsg_set("Didn't get a valid array of plot strings");
    }

    return err;
}

/* This supports production of a multiplot from the array of
   plot-specification strings identified by @param. Such an
   array may be produced by "gpbuild", or may be assembled
   manually.
*/

int gretl_multiplot_from_array (const char *param, gretlopt opt)
{
    gretl_array *a = NULL;
    int maxp, np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return E_DATA;
    }

    err = retrieve_plots_array(param, &a, &np);
    if (err) {
        return err;
    }

    maxp = np;

    if (!err) {
	/* pick up any options */
	gretlopt myopt = opt;

	set_multiplot_defaults();
	myopt &= ~OPT_S;
	if (myopt) {
	    err = set_multiplot_sizes(myopt);
	}
	if (myopt & OPT_L) {
	    err = maybe_set_mp_layout(&np);
	}
	if (!err && mp_layout == NULL) {
	    err = multiplot_set_grid(np);
	}
    }

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script(a, np, maxp);
    }

    return err;
}

int check_multiplot_options (int ci, gretlopt opt)
{
    int err = 0;

    if (ci == GPBUILD) {
	/* no options accepted */
	err = opt == OPT_NONE ? 0 : E_BADOPT;
    } else {
	/* gridplot: can't have both --output and --outbuf */
	err = incompatible_options(opt, OPT_U | OPT_b);
    }

    return err;
}

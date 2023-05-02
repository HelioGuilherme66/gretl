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

/* Support for the built-in "gridplot" command for producing
   multiple plots.
*/

#include "libgretl.h"
#include "uservar.h"
#include "gretl_multiplot.h"

typedef struct {
    int row;
    int col;
    gchar *buf;
} mplot_element;

typedef struct {
    gretlopt flag;
    int *target;
    int min;
    int max;
} mplot_option;

static GArray *multiplot;
static int mp_fontsize = 10;
static int mp_width = 800;
static int mp_height = 600;
static int mp_rows;
static int mp_cols;
static int mp_collecting;

static void set_multiplot_defaults (void)
{
    mp_fontsize = 10;
    mp_width = 800;
    mp_height = 600;
    mp_rows = 0;
    mp_cols = 0;
}

static void mplot_element_clear (mplot_element *element)
{
    g_free(element->buf);
}

int gretl_multiplot_active (void)
{
    return multiplot != NULL && mp_collecting;
}

void gretl_multiplot_destroy (void)
{
    if (multiplot != NULL) {
        g_array_unref(multiplot);
        multiplot = NULL;
    }
}

static const mplot_option mp_options[] = {
    { OPT_F, &mp_fontsize, 4, 24 },
    { OPT_W, &mp_width,    200, 2048 },
    { OPT_H, &mp_height,   200, 2048 },
    { OPT_R, &mp_rows,     1, 12 },
    { OPT_C, &mp_cols,     1, 12 }
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

int gretl_multiplot_start (gretlopt opt)
{
    int err = 0;

    if (multiplot == NULL) {
	if (opt == OPT_NONE) {
	    set_multiplot_defaults();
	} else {
	    err = set_multiplot_sizes(opt);
	}
	if (!err) {
	    multiplot = g_array_new(FALSE, FALSE, sizeof(mplot_element));
	    g_array_set_clear_func(multiplot, (GDestroyNotify) mplot_element_clear);
	    mp_collecting = 1;
	}
    } else {
        gretl_errmsg_set("multplot: cannot be nested");
        err = E_DATA;
    }

    return err;
}

int gretl_multiplot_add_plot (int row, int col, gchar *buf)
{
    if (multiplot != NULL && buf != NULL) {
        mplot_element element = {row, col, buf};

        g_array_append_val(multiplot, element);
        return 0;
    } else {
	gretl_errmsg_set("gretl_multiplot_add_plot: failed");
        return E_DATA;
    }
}

static int multiplot_set_grid (int n)
{
    int err = 0;

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

    return err;
}

static int output_multiplot_script (void)
{
    mplot_element *element;
    int i, err = 0;
    FILE *fp;

    fp = open_plot_input_file(PLOT_USER_MULTI, 0, &err);
    if (err) {
        return err;
    }

    fputs("# literal lines = 1\n", fp);
    fprintf(fp, "# grid_params: fontsize=%d, width=%d, height=%d, "
            "rows=%d, cols=%d, plots=%d\n", mp_fontsize, mp_width,
            mp_height, mp_rows, mp_cols, multiplot->len);
    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n", mp_rows, mp_cols);
    gretl_push_c_numeric_locale();
    for (i=0; i<multiplot->len; i++) {
        fprintf(fp, "# gridplot: subplot %d\n", i);
        if (i > 0) {
            fputs("reset\n", fp);
        }
        element = &g_array_index(multiplot, mplot_element, i);
        fputs(element->buf, fp);
    }
    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);
    err = finalize_plot_input_file(fp);

    return err;
}

int gretl_multiplot_finalize (gretlopt opt)
{
    int np, err = 0;

    if (multiplot == NULL) {
	gretl_errmsg_set("end multiplot: multiplot not started");
	return E_DATA;
    }

    mp_collecting = 0;
    np = multiplot->len;
    if (np > 0) {
	err = multiplot_set_grid(np);
    }

    if (np > 0 && !err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script();
    }

    gretl_multiplot_destroy();

    return err;
}

static void filter_multiplot_buffer (const char *buf,
				     int np, FILE *fp)
{
    char line[1024];

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	if (!strncmp(line, "# grid_params: ", 15)) {
	    /* substitute revised parameters in comment */
	    fprintf(fp, "# grid_params: fontsize=%d, width=%d, height=%d, "
		    "rows=%d, cols=%d, plots=%d\n", mp_fontsize, mp_width,
		    mp_height, mp_rows, mp_cols, np);
	} else if (!strncmp(line, "set multiplot", 13)) {
	    /* and revise multiplot layout spec */
	    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n",
		    mp_rows, mp_cols);
	} else {
	    /* simply transcribe */
	    fputs(line, fp);
	}
    }

    bufgets_finalize(buf);
}

static int revise_multiplot_script (const char *buf,
				    int filter,
				    int np)
{
    int err = 0;
    FILE *fp;

    fp = open_plot_input_file(PLOT_USER_MULTI, 0, &err);
    if (err) {
        return err;
    }

    if (filter) {
	filter_multiplot_buffer(buf, np, fp);
    } else {
	fputs(buf, fp);
    }

    err = finalize_plot_input_file(fp);

    return err;
}

static int retrieve_grid_params (const char *buf,
				 int *parms,
				 int *np)
{
    char line[256];
    int i = 0;

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && 1 < 10) {
        if (!strncmp(line, "# grid_params: ", 15)) {
            sscanf(line+15, "fontsize=%d, width=%d, height=%d, "
		   "rows=%d, cols=%d, plots=%d", &parms[0], &parms[1],
		   &parms[2], &parms[3], &parms[4], np);
            break;
        }
        i++;
    }

    bufgets_finalize(buf);

    return 0;
}

static void get_prior_plot_spec (gretlopt opt,
				 const char **pbuf,
				 gchar **pgbuf)
{
    const char *s;

    if (opt & OPT_i) {
	s = get_optval_string(GRIDPLOT, OPT_i);
	*pbuf = get_string_by_name(s);
    } else if (opt & OPT_I) {
	gboolean ok;

	s = get_optval_string(GRIDPLOT, OPT_I);
	ok = g_file_get_contents(s, pgbuf, NULL, NULL);
	if (ok) {
	    *pbuf = *pgbuf;
	}
    }
}

int gretl_multiplot_revise (gretlopt opt)
{
    const char *argname;
    const char *buf = NULL;
    gchar *gbuf = NULL;
    gretlopt myopt = opt;
    int params[5] = {0};
    int filter = 0;
    int i, np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return E_DATA;
    }

    err = incompatible_options(opt, OPT_i | OPT_I);
    if (!err) {
	get_prior_plot_spec(opt, &buf, &gbuf);
    }
    if (buf == NULL) {
        gretl_errmsg_set("Couldn't find an input specification");
        return E_DATA;
    }

    /* what sort of output? */
    argname = get_optval_string(GRIDPLOT, OPT_U);
    if (argname != NULL) {
        fprintf(stderr, "gretl_multiplot_revise: output '%s'\n",
                argname);
    }

    /* extract the dimensions recorded in @buf */
    retrieve_grid_params(buf, params, &np);

    /* let @opt override previous dimensions, if relevant */
    myopt &= ~OPT_i;
    myopt &= ~OPT_U;
    if (myopt != OPT_NONE) {
        err = set_multiplot_sizes(opt);
        if (!err) {
            err = multiplot_set_grid(np);
        }
    }

    if (!err) {
	/* Are any optional parameters changed? If, so, we'll
	   need to filter @buf, not just transcribe it.
	*/
	for (i=0; i<n_mp_options; i++) {
	    if (params[i] != *(mp_options[i].target)) {
		filter = 1;
		break;
	    }
	}
    }

#if 0
    fprintf(stderr, "*** here's the prior input (filter = %d) ***\n", filter);
    fputs(buf, stderr);
#endif

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
	revise_multiplot_script(buf, filter, np);
    }

    if (gbuf != NULL) {
	g_free(gbuf);
    }

    return err;
}

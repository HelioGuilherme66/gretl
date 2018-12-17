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

/* database.c for gretl */

#include "gretl.h"
#include "boxplots.h"
#include "database.h"
#include "datafiles.h"
#include "gretl_xml.h"
#include "gretl_www.h"
#include "gretl_untar.h"
#include "gretl_zip.h"
#include "gretl_string_table.h"
#include "csvdata.h"
#include "menustate.h"
#include "treeutils.h"
#include "textbuf.h"
#include "winstack.h"
#include "toolbar.h"
#include "dlgutils.h"
#include "fncall.h"
#include "dbread.h"
#include "fncall.h"
#include "varinfo.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>
#include <dirent.h>
#include <errno.h>

#if G_BYTE_ORDER == G_BIG_ENDIAN
# include <netinet/in.h>
#endif

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#define DB_SEARCH_DEBUG 0

/* private functions */
static GtkWidget *database_window (windata_t *vwin);
static int add_local_db_series_list (windata_t *vwin);
static int add_remote_db_series_list (windata_t *vwin, char *buf);
static int add_rats_db_series_list (windata_t *vwin);
static int add_pcgive_db_series_list (windata_t *vwin);
static dbwrapper *get_db_series_info (windata_t *vwin, int action);
static int *db_get_selection_list (windata_t *vwin);
static void gui_get_db_series (windata_t *vwin, int cmd);

enum db_data_actions {
    DB_DISPLAY,
    DB_GRAPH,
    DB_IMPORT
};

enum db_codebook_type {
    CB_NONE,
    CB_TEXT,
    CB_PDF
};

/* columns of individual database window */
enum {
    DBCOL_VARNAME,
    DBCOL_DESCRIP,
    DBCOL_OBSINFO,
    DBCOL_OFFSET
};

/* columns of list-of-databases window */
enum {
    COL_DBNAME,
    COL_DBINFO,
    COL_DBPATH
};

static int utf8_correct (char *orig)
{
    int err = 0;

    if (!g_utf8_validate(orig, -1, NULL)) {
	GError *gerr = NULL;
	gchar *conv;
	gsize wrote;

	conv = g_convert(orig, -1,
			 "UTF-8",
			 "ISO-8859-1",
			 NULL, &wrote, &gerr);

	if (gerr != NULL) {
	    errbox(gerr->message);
	    g_error_free(gerr);
	    strcpy(orig, "invalid string");
	    err = 1;
	} else {
	    strcpy(orig, conv);
	    g_free(conv);
	}
    }

    return err;
}

static void update_statusline (windata_t *vwin, const char *s)
{
    gchar *tmp = g_strdup_printf(_("Network status: %s"), s);

    gtk_label_set_text(GTK_LABEL(vwin->status), tmp);

    while (gtk_events_pending()) {
	gtk_main_iteration();
    }

    g_free(tmp);
}

static void set_time_series (DATASET *dset)
{
    if (dset->pd != 1 || strcmp(dset->stobs, "1")) {
	dset->structure = TIME_SERIES;
    }
}

void show_network_error (windata_t *vwin)
{
    const char *msg = gretl_errmsg_get();
    char *buf = NULL;

    if (*msg != '\0') {
	buf = gretl_strdup(msg);
    }

    if (buf != NULL) {
	size_t n = strlen(buf);

	if (buf[n-1] == '\n') {
	    buf[n-1] = '\0';
	}
	if (vwin != NULL) {
	    update_statusline(vwin, buf);
	} else {
	    errbox(buf);
	}
	free(buf);
    } else if (vwin != NULL) {
	update_statusline(vwin, _("Error retrieving data from server"));
    } else {
	errbox(_("Error retrieving data from server"));
    }
}

static int gui_get_remote_db_data (windata_t *vwin, SERIESINFO *sinfo,
				   double **Z)
{
    char *dbbase = vwin->fname;
    int err;

    update_statusline(vwin, _("Retrieving data..."));

    err = get_remote_db_data(dbbase, sinfo, Z);

    if (err) {
	show_network_error(vwin);
	return E_FOPEN;
    } else {
	update_statusline(vwin, "OK");
    }

    return err;
}

static void display_dbdata (DATASET *dbset)
{
    PRN *prn;
    int width = 36;

    if (bufopen(&prn)) {
	return;
    }

    if (dbset->v > 1) {
	width = 72;
    }

    printdata(NULL, NULL, dbset, OPT_O, prn);
    view_buffer(prn, width, 350, _("gretl: display database series"), PRINT,
		NULL);
}

static void graph_dbdata (DATASET *dbset)
{
    int *list;
    int err;

    list = gretl_consecutive_list_new(1, dbset->v - 1);
    if (list == NULL) {
	nomem();
	return;
    }

    if (dbset->structure == CROSS_SECTION) {
	err = boxplots(list, NULL, dbset, OPT_NONE);
    } else {
	err = gnuplot(list, NULL, dbset, OPT_G | OPT_O | OPT_T);
    }

    free(list);
    gui_graph_handler(err);
}

static int expand_data_dialog (int src_pd, int targ_pd, int *interpol,
			       GtkWidget *parent)
{
    int mult = targ_pd / src_pd;
    int resp;

    if ((targ_pd == 4 && src_pd == 1) ||
	(targ_pd == 12 && src_pd == 4)) {
	/* interpolation is an option */
	const char *opts[] = {
	    N_("Repeat the lower frequency values"),
	    N_("Interpolate higher frequency values")
	};

	resp = radio_dialog("gretl", _("Adding a lower frequency series to a\n"
				       "higher frequency dataset"),
			    opts, 2, 0, EXPAND, parent);
	if (resp == 1) {
	    *interpol = 1;
	    resp = GRETL_YES;
	} else if (resp == 0) {
	    resp = GRETL_YES;
	}
    } else {
	/* can only do expansion via replication */
	gchar *msg =
	    g_strdup_printf(_("Do you really want to add a lower frequency series\n"
			      "to a higher frequency dataset?\n\n"
			      "If you say 'yes' I will expand the source data by\n"
			      "repeating each value %d times.  In general, this is\n"
			      "not a valid thing to do."),
			    mult);
	resp = yes_no_dialog("gretl", msg, parent);
	g_free(msg);
    }

    return resp;
}

static int obs_overlap_check (int pd, const char *stobs,
			      const char *endobs,
			      const char *varname)
{
    int err;

    err = db_range_check(pd, stobs, endobs, varname, dataset);
    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int pd_conversion_check (int db_pd)
{
    int err;

    err = check_db_import_conversion(db_pd, dataset);

    if (err) {
	warnbox(_("Sorry, can't handle this frequency conversion"));
    }

    return err;
}

static const char *compact_method_string (CompactMethod method)
{
    if (method == COMPACT_SUM) {
	return "sum";
    } else if (method == COMPACT_SOP) {
	return "first";
    } else if (method == COMPACT_EOP) {
	return "last";
    } else if (method == COMPACT_SPREAD) {
	return "spread";
    } else {
	return NULL;
    }
}

static const char *trimmed_db_name (const char *fname)
{
    const char *s = fname;

    if (strstr(fname, gretl_binbase()) != NULL) {
	s = path_last_slash_const(fname);
	if (s != NULL) {
	    return s + 1;
	}
    }

    return NULL;
}

static void record_db_open_command (dbwrapper *dw)
{
    const char *dbname = (dw == NULL)? "dbnomics" : dw->fname;
    int dbtype = (dw == NULL)? GRETL_DBNOMICS : dw->dbtype;

    if (dbname != NULL) {
	const char *current_db = get_db_name();
	char *dbpath = NULL;

	if (dbtype == GRETL_PCGIVE_DB) {
	    dbpath = g_strdup_printf("%s.bn7", dbname);
	} else if (dbtype == GRETL_NATIVE_DB) {
	    dbpath = g_strdup_printf("%s.bin", dbname);
	} else {
	    dbpath = g_strdup(dbname);
	}

	/* record the "open" command if the database in
	   question is not already open */

	if (strcmp(current_db, dbpath)) {
	    int err = set_db_name(dbpath, dbtype, NULL);
	    int done = 0;

	    if (!err && dbtype == GRETL_NATIVE_DB) {
		const char *s = trimmed_db_name(dbpath);

		if (s != NULL) {
		    lib_command_sprintf("open %s", s);
		    done = 1;
		}
	    }

	    if (!err && !done) {
		if (dbtype == GRETL_NATIVE_DB_WWW) {
		    lib_command_sprintf("open %s --www", dbpath);
		} else if (strchr(dbpath, ' ') != NULL) {
		    lib_command_sprintf("open \"%s\"", dbpath);
		} else {
		    lib_command_sprintf("open %s", dbpath);
		}
	    }

	    if (!err) {
		record_command_verbatim();
	    }
	}

	g_free(dbpath);
    }
}

static char *dbnomics_id;

static void set_dbnomics_id (const char *s)
{
    dbnomics_id = g_strdup(s);
}

static void unset_dbnomics_id (void)
{
    free(dbnomics_id);
    dbnomics_id = NULL;
}

/* record successful importation in command log */

static void record_db_import (const char *vname,
			      int compact,
			      int interpol,
			      CompactMethod method)
{
    const char *cstr = NULL;

    if (dbnomics_id != NULL) {
	if (compact && (cstr = compact_method_string(method)) != NULL) {
	    lib_command_sprintf("data %s --name=%s --compact=%s",
				dbnomics_id, vname, cstr);
	} else if (interpol) {
	    lib_command_sprintf("data %s --name=%s --interpolate",
				dbnomics_id, vname);
	} else {
	    lib_command_sprintf("data %s --name=%s", dbnomics_id, vname);
	}
    } else {
	if (compact && (cstr = compact_method_string(method)) != NULL) {
	    lib_command_sprintf("data %s --compact=%s", vname,
				cstr);
	} else if (interpol) {
	    lib_command_sprintf("data %s --interpolate", vname);
	} else {
	    lib_command_sprintf("data %s", vname);
	}
    }

    record_command_verbatim();
}

static int handle_compact_spread (double **dbZ,
				  SERIESINFO *sinfo,
				  DATASET *dset,
				  DATASET *dbset)
{
    PRN *prn = NULL;
    int err = 0;

    if (bufopen(&prn)) {
	return E_ALLOC;
    }

    err = lib_spread_db_data(dbZ, sinfo, dset, dbset, prn);

    if (err) {
	char *buf = gretl_print_steal_buffer(prn);

	if (buf != NULL && *buf != '\0') {
	    errbox(buf);
	} else {
	    gui_errmsg(err);
	}
    }

    gretl_print_destroy(prn);

    return err;
}

static void maybe_retrieve_compact_method (int v, CompactMethod *pm)
{
    int m = series_get_compact_method(dataset, v);

    if (m != COMPACT_NONE) {
	*pm = m;
    }
}

static int
add_single_series_to_dataset (windata_t *vwin, DATASET *dbset)
{
    CompactMethod cmethod = COMPACT_AVG;
    int dbv, resp, overwrite = 0;
    int existing = 0;
    int compact = 0;
    int interpol = 0;
    int err = 0;

    err = pd_conversion_check(dbset->pd);
    if (!err) {
	err = obs_overlap_check(dbset->pd, dbset->stobs, dbset->endobs,
				dbset->varname[1]);
    }
    if (err) {
	return err;
    }

    record_db_open_command(NULL);

    /* is there a series of this name already in the dataset? */
    dbv = series_index(dataset, dbset->varname[1]);
    if (dbv < dataset->v) {
	existing = 1;
	maybe_retrieve_compact_method(dbv, &cmethod);
    }

    if (dbset->pd < dataset->pd) {
	/* the incoming series needs to be expanded */
	resp = expand_data_dialog(dbset->pd, dataset->pd, &interpol,
				  vwin->main);
	if (resp != GRETL_YES) {
	    return 0;
	}
    } else if (dbset->pd > dataset->pd) {
	/* the incoming series needs to be compacted */
	data_compact_dialog(dbset->pd, &dataset->pd, NULL,
			    &cmethod, NULL, vwin->main);
	if (cmethod == COMPACT_NONE) {
	    return 0; /* canceled */
	} else if (cmethod == COMPACT_SPREAD) {
	    return handle_compact_spread(NULL, NULL, dataset, dbset);
	}
	compact = 1;
    }

    if (existing) {
	resp = yes_no_dialog("gretl",
			     _("There is already a variable of this name\n"
			       "in the dataset.  OK to overwrite it?"),
			     vwin_toplevel(vwin));
	if (resp != GRETL_YES) {
	    return 0;
	}
	overwrite = 1;
    }

    if (!overwrite) {
	err = dataset_add_series(dataset, 1);
	if (err) {
	    nomem();
	    return err;
	}
    }

    err = transcribe_db_data(dataset, dbv, dbset->Z[1], dbset->pd,
			     dbset->n, dbset->stobs, cmethod,
			     interpol);

    if (!err) {
	const char *vlabel = series_get_label(dbset, 1);

	strcpy(dataset->varname[dbv], dbset->varname[1]);
	if (vlabel != NULL && *vlabel != '\0') {
	    series_set_label(dataset, dbv, vlabel);
	}
	record_db_import(dbset->varname[1], compact, interpol, cmethod);
    } else {
	if (!overwrite) {
	    dataset_drop_last_variables(dataset, 1);
	}
	gui_errmsg(err);
    }

    return err;
}

/* multiple series version of data adding function */

static int
add_db_series_to_dataset (windata_t *vwin, double **dbZ, dbwrapper *dw)
{
    SERIESINFO *sinfo;
    CompactMethod cmethod = COMPACT_AVG;
    int resp, warned = 0, chosen = 0;
    int i, err = 0;

    sinfo = &dw->sinfo[0];
    err = pd_conversion_check(sinfo->pd);
    if (err) {
	return err;
    }

    record_db_open_command(dw);

    for (i=0; i<dw->nv && !err; i++) {
	int v, dbv;
	int existing = 0;
	int overwrite = 0;
	int compact = 0;
	int interpol = 0;

	sinfo = &dw->sinfo[i];
	v = sinfo->v;

	if (obs_overlap_check(sinfo->pd, sinfo->stobs, sinfo->endobs,
			      sinfo->varname)) {
	    continue;
	}

	dbv = series_index(dataset, sinfo->varname);
	if (dbv < dataset->v) {
	    existing = 1;
	    maybe_retrieve_compact_method(dbv, &cmethod);
	}

	if (sinfo->pd < dataset->pd) {
	    /* the incoming series needs to be expanded */
	    if (!warned) {
		resp = expand_data_dialog(sinfo->pd, dataset->pd, &interpol,
					  vwin->main);
		if (resp != GRETL_YES) {
		    return 0;
		}
		warned = 1;
	    }
	} else if (sinfo->pd > dataset->pd) {
	    /* the incoming series needs to be compacted */
	    if (!chosen) {
		data_compact_dialog(sinfo->pd, &dataset->pd, NULL,
				    &cmethod, NULL, vwin->main);
		if (cmethod == COMPACT_NONE) {
		    /* canceled */
		    return 0;
		}
		chosen = 1;
	    }
	    compact = 1;
	    if (cmethod == COMPACT_SPREAD) {
		err = handle_compact_spread(dbZ, sinfo, dataset, NULL);
		if (err) {
		    break;
		} else {
		    record_db_import(sinfo->varname, compact, interpol, cmethod);
		    continue;
		}
	    }
	}

	if (existing) {
	    if (dw->nv == 1) {
		resp = yes_no_dialog("gretl",
				     _("There is already a variable of this name\n"
				       "in the dataset.  OK to overwrite it?"),
				     vwin_toplevel(vwin));
		if (resp != GRETL_YES) {
		    return 0;
		}
	    }
	    overwrite = 1;
	}

	if (!overwrite) {
	    err = dataset_add_series(dataset, 1);
	    if (err) {
		nomem();
		return err;
	    }
	}

	err = transcribe_db_data(dataset, dbv, dbZ[v], sinfo->pd,
				 sinfo->nobs, sinfo->stobs, cmethod,
				 interpol);
	if (err) {
	    gui_errmsg(err);
	    if (!overwrite) {
		dataset_drop_last_variables(dataset, 1);
	    }
	} else {
	    strcpy(dataset->varname[dbv], sinfo->varname);
	    series_set_label(dataset, dbv, sinfo->descrip);
	    record_db_import(sinfo->varname, compact, interpol, cmethod);
	}
    }

    return err;
}

static void add_dbdata (windata_t *vwin, DATASET *dbset,
			dbwrapper *dw, int *freeit)
{
    SERIESINFO *sinfo = NULL;
    int i, err = 0;

    if (data_status) {
	/* we already have data in gretl's workspace */
	if (dw != NULL) {
	    add_db_series_to_dataset(vwin, dbset->Z, dw);
	} else {
	    add_single_series_to_dataset(vwin, dbset);
	}
    } else {
	/* no data open: start new data set from db */
	destroy_dataset(dataset);
	dataset = dbset;
	*freeit = 0;

	if (dw != NULL) {
	    record_db_open_command(dw);
	    for (i=1; i<=dw->nv && !err; i++) {
		sinfo = &dw->sinfo[i-1];
		strcpy(dataset->varname[i], sinfo->varname);
		series_set_label(dataset, i, sinfo->descrip);
		lib_command_sprintf("data %s", sinfo->varname);
		record_command_verbatim();
	    }
	} else {
	    record_db_open_command(NULL);
	    record_db_import(dbset->varname[1], 0, 0, 0);
	}

	data_status |= (GUI_DATA | MODIFIED_DATA);
    }

    if (!err) {
	register_data(DATA_APPENDED);
    }
}

static void db_display_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_db_series(vwin, DB_DISPLAY);
}

static void db_graph_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_db_series(vwin, DB_GRAPH);
}

static void db_import_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_db_series(vwin, DB_IMPORT);
}

void sync_db_windows (void)
{
    const char *dname = get_db_name();

    if (*dname != '\0') {
	windata_t *vwin = get_browser_for_gretl_database(dname);

	if (vwin != NULL) {
	    add_local_db_series_list(vwin);
	}
    }
}

static void db_delete_callback (GtkWidget *w, windata_t *vwin)
{
    int *list = db_get_selection_list(vwin);
    gchar *query;
    int resp, err = 0;

    if (list == NULL) {
	return;
    }

    query = g_strdup_printf(_("Really delete the selected series\n"
			      "from the database '%s'?"),
			    gtk_window_get_title(GTK_WINDOW(vwin->main)));

    resp = yes_no_dialog("gretl", query, vwin_toplevel(vwin));

    g_free(query);

    if (resp == GRETL_YES) {
	err = db_delete_series_by_number(list, vwin->fname);
	if (err) {
	    gui_errmsg(err);
	} else {
	    /* revise window contents */
	    add_local_db_series_list(vwin);
	}
    }
}

/* this pointer to the source "vwin" for drag-and-drop is
   set by drag selection callback, do_db_drag()
*/
static windata_t *vwin_drag_src;

/* import series from a database window into main
   gretl workspace */

void drag_import_db_series (void)
{
    windata_t *vwin = vwin_drag_src;

    if (vwin != NULL) {
	gui_get_db_series(vwin, DB_IMPORT);
	vwin_drag_src = NULL;
    }
}

static int diffdate (double d1, double d0, int pd)
{
    double x;

    if (pd == 4 || pd == 12) {
	char s[16];
	int maj, min;
	int dmaj, dmin;

	gretl_push_c_numeric_locale();

	sprintf(s, "%g", d1);
	sscanf(s, "%d.%d", &dmaj, &dmin);

	sprintf(s, "%g", d0);
	sscanf(s, "%d.%d", &maj, &min);

	gretl_pop_c_numeric_locale();

	dmaj -= maj;
	dmin -= min;

	x = dmaj * pd + dmin;
    } else {
	x = d1 - d0;
    }

    return x;
}

static DATASET *dataset_from_dbwrapper (dbwrapper *dw)
{
    DATASET *dset = NULL;
    SERIESINFO *sinfo;
    char stobs[OBSLEN], endobs[OBSLEN];
    double xd, xdmax = 0, xdmin = DBL_MAX;
    int n0 = 0, nmax = 0;
    int i;

    /* Here we construct a dataset which can accommodate all
       the seleted series (in case there's more than one).
       Note that while multiple selection is enabled in a
       database window we have made it impossible to compose
       a selection that includes series of differing
       frequencies (i.e. they must all be quarterly, all
       annual, or whatever).
    */

    for (i=0; i<dw->nv; i++) {
	sinfo = &dw->sinfo[i];
	xd = get_date_x(sinfo->pd, sinfo->stobs);
	fprintf(stderr, "var %d: nobs=%d, pd=%d, stobs='%s', sd0=%g\n",
		i, sinfo->nobs, sinfo->pd, sinfo->stobs, xd);
	if (xd < xdmin) {
	    strcpy(stobs, sinfo->stobs);
	    xdmin = xd;
	    if (sinfo->nobs > n0) {
		n0 = sinfo->nobs;
	    }
	}
	if (xd > xdmax) {
	    strcpy(endobs, sinfo->endobs);
	    xdmax = xd;
	}
	if (sinfo->nobs > nmax) {
	    nmax = sinfo->nobs;
	}
	sinfo->v = i + 1;
	sinfo->t2 = sinfo->nobs - 1;
    }

    if (xdmax > xdmin) {
	int ni, dd;

	for (i=0; i<dw->nv; i++) {
	    sinfo = &dw->sinfo[i];
	    ni = sinfo->nobs;
	    xd = get_date_x(sinfo->pd, sinfo->stobs);
	    if (xd > xdmin) {
		dd = diffdate(xd, xdmin, sinfo->pd);
		ni += dd;
		sinfo->t1 = dd;
		sinfo->t2 += dd;
	    }
	    if (ni > nmax) {
		nmax = ni;
	    }
	}
    }

    fprintf(stderr, "min(sd0) = %g, stobs='%s', n = %d\n", xdmin,
	    stobs, nmax);

    dset = create_new_dataset(dw->nv + 1, nmax, 0);

    if (dset != NULL) {
	dset->pd = dw->sinfo[0].pd;
	strcpy(dset->stobs, stobs);
	strcpy(dset->endobs, endobs);
	colonize_obs(dset->stobs);
	colonize_obs(dset->endobs);
	dset->sd0 = xdmin;
	set_time_series(dset);
    }

    return dset;
}

static void gui_get_db_series (windata_t *vwin, int cmd)
{
    int dbcode = vwin->role;
    DATASET *dbset = NULL;
    dbwrapper *dw;
    int freeit = 1;
    int i, err = 0;

    /* build wrapper for selected row(s) of @vwin */
    dw = get_db_series_info(vwin, dbcode);
    if (dw == NULL) {
	return;
    }

    /* build "superset" dataset, allowing multiple selection */
    dbset = dataset_from_dbwrapper(dw);
    if (dbset == NULL) {
	dbwrapper_destroy(dw);
	nomem();
	return;
    }

    for (i=0; i<dw->nv; i++) {
	SERIESINFO *sinfo = &dw->sinfo[i];

	if (dbcode == NATIVE_SERIES) {
	    err = get_native_db_data(vwin->fname, sinfo, dbset->Z);
	} else if (dbcode == REMOTE_SERIES) {
	    err = gui_get_remote_db_data(vwin, sinfo, dbset->Z);
	} else if (dbcode == RATS_SERIES) {
	    err = get_rats_db_data(vwin->fname, sinfo, dbset->Z);
	} else if (dbcode == PCGIVE_SERIES) {
	    err = get_pcgive_db_data(vwin->fname, sinfo, dbset->Z);
	}

	if (cmd == DB_IMPORT && err == DB_MISSING_DATA) {
	    warnbox(_("Warning: series has missing observations"));
	}

	if (err && err != DB_MISSING_DATA && dbcode != REMOTE_SERIES) {
	    errbox(_("Couldn't access binary datafile"));
	    goto bailout;
	}

	strcpy(dbset->varname[i+1], sinfo->varname);
	series_set_label(dbset, i+1, sinfo->descrip);
    }

    if (cmd == DB_DISPLAY) {
	display_dbdata(dbset);
    } else if (cmd == DB_GRAPH) {
	graph_dbdata(dbset);
    } else if (cmd == DB_IMPORT) {
	add_dbdata(vwin, dbset, dw, &freeit);
    }

 bailout:

    if (freeit) {
	destroy_dataset(dbset);
    }

    dbwrapper_destroy(dw);
}

/* double-click callback */

void display_db_series (windata_t *vwin)
{
    gui_get_db_series(vwin, DB_DISPLAY);
}

static void db_view_codebook (GtkWidget *w, windata_t *vwin)
{
    gchar *cbname;

    if (vwin->flags & VWIN_CB_PDF) {
	cbname = g_strdup_printf("%s.pdf", vwin->fname);
	gretl_show_pdf(cbname, NULL);
	g_free(cbname);
    } else {
	cbname = g_strdup_printf("%s.cb", vwin->fname);
	view_file(cbname, 0, 0, 78, 350, VIEW_CODEBOOK);
    }

    g_free(cbname);
}

static void db_show_index (GtkWidget *w, windata_t *vwin)
{
    show_native_dbs();
}

static void build_db_content_popup (windata_t *vwin, int cb, int del)
{
    if (vwin->popup != NULL) {
	return;
    }

    vwin->popup = gtk_menu_new();

    add_popup_item(_("Display"), vwin->popup,
		   G_CALLBACK(db_display_series),
		   vwin);
    add_popup_item(_("Graph"), vwin->popup,
		   G_CALLBACK(db_graph_series),
		   vwin);
    add_popup_item(_("Import"), vwin->popup,
		   G_CALLBACK(db_import_series),
		   vwin);
    if (del) {
	add_popup_item(_("Delete"), vwin->popup,
		       G_CALLBACK(db_delete_callback),
		       vwin);
    }

    if (cb) {
	add_popup_item(_("Codebook"), vwin->popup,
		       G_CALLBACK(db_view_codebook),
		       vwin);
    }
}

enum {
    DEL_BTN = 1,
    CB_BTN,
    IDX_BTN
};

static GretlToolItem db_items[] = {
    { N_("Display values"), GTK_STOCK_MEDIA_PLAY, G_CALLBACK(db_display_series), 0 },
    { N_("Graph"),          GRETL_STOCK_TS,   G_CALLBACK(db_graph_series), 0 },
    { N_("Add to dataset"), GTK_STOCK_ADD,    G_CALLBACK(db_import_series), 0 },
    { N_("List databases"), GTK_STOCK_INDEX,  G_CALLBACK(db_show_index), IDX_BTN },
    { N_("Delete"),         GTK_STOCK_DELETE, G_CALLBACK(db_delete_callback), DEL_BTN },
    { N_("Codebook"),       GRETL_STOCK_BOOK, G_CALLBACK(db_view_codebook), CB_BTN }
};

static int n_db_items = G_N_ELEMENTS(db_items);

static void make_db_toolbar (windata_t *vwin, int cb, int del,
			     int index_button)
{
    GtkWidget *hbox;
    GretlToolItem *item;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new(NULL);

    for (i=0; i<n_db_items; i++) {
	item = &db_items[i];
	if (!del && item->flag == DEL_BTN) {
	    continue;
	} else if (!cb && item->flag == CB_BTN) {
	    continue;
	} else if (!index_button && item->flag == IDX_BTN) {
	    continue;
	}
	gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    vwin_add_winlist(vwin);
    vwin_add_finder(vwin);
    gtk_widget_show_all(hbox);
}

static int db_has_codebook (const char *fname)
{
    gchar *testname;
    int err, ret = CB_NONE;

    /* try first for *.cb (plain text) file */
    testname = g_strdup_printf("%s.cb", fname);
    err = gretl_test_fopen(testname, "rb");
    if (err == 0) {
	ret = CB_TEXT;
    } else {
	/* try for PDF documentation? */
	g_free(testname);
	testname = g_strdup_printf("%s.pdf", fname);
	err = gretl_test_fopen(testname, "rb");
	if (err == 0) {
	    ret = CB_PDF;
	}
    }

    g_free(testname);

    return ret;
}

static int db_is_writable (int action, const char *fname)
{
    int ret = 0;

    if (action == NATIVE_SERIES) {
	char testname[MAXLEN];
	int err;

	strcpy(testname, fname);
	strcat(testname, ".bin");
	err = gretl_write_access(testname);
	if (!err) {
	    ret = 1;
	}
    }

    return ret;
}

static gboolean
db_col_callback (GtkWidget *w, GdkEventMotion *event, gpointer p)
{
    GtkTreeViewColumn *col =
	gtk_tree_view_get_column(GTK_TREE_VIEW(w), 1);

    if (gtk_tree_view_column_get_max_width(col) > 0) {
	/* remove the width constraint */
	gtk_tree_view_column_set_max_width(col, -1);
    }

    return 0;
}

static void
maybe_adjust_descrip_column (windata_t *vwin)
{
    GtkTreeViewColumn *col;
    GdkWindow *window;
    gint w0, w1, lw, w1max;

    col = gtk_tree_view_get_column(GTK_TREE_VIEW(vwin->listbox), DBCOL_VARNAME);
    w0 = gtk_tree_view_column_get_width(col);

    col = gtk_tree_view_get_column(GTK_TREE_VIEW(vwin->listbox), DBCOL_DESCRIP);
    w1 = gtk_tree_view_column_get_width(col);

    window = gtk_widget_get_window(vwin->listbox);

#if GTK_MAJOR_VERSION >= 3
    lw = gdk_window_get_width(window);
#else
    gdk_drawable_get_size(window, &lw, NULL);
#endif

    w1max = lw - w0 - 140;

    if (w1 > w1max) {
	gtk_tree_view_column_set_max_width(col, w1max);
	g_signal_connect(vwin->listbox, "motion-notify-event",
			 G_CALLBACK(db_col_callback), NULL);
    }
}

static int
make_db_index_window (int action, char *fname, char *buf,
		      int index_button)
{
    GtkWidget *listbox;
    gchar *title;
    windata_t *vwin;
    int db_width = 700, db_height = 420;
    int cb = 0, del = 0;
    int err = 0;

    vwin = get_browser_for_database(fname);

    if (vwin != NULL) {
	gtk_window_present(GTK_WINDOW(vwin->main));
	return 0;
    }

    if (action == REMOTE_SERIES && buf == NULL) {
	return 1;
    }

    if (buf == NULL && path_last_slash(fname) != NULL) {
	title = path_last_slash(fname) + 1;
    } else {
	title = fname;
    }

    vwin = gretl_browser_new(action, title);
    if (vwin == NULL) {
	return 1;
    }

    db_width *= gui_scale;
    db_height *= gui_scale;
    gtk_window_set_default_size(GTK_WINDOW(vwin->main), db_width, db_height);

    if (action == NATIVE_SERIES || action == PCGIVE_SERIES) {
	strip_extension(fname);
    }

    strcpy(vwin->fname, fname);

    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    cb = db_has_codebook(fname);
    del = db_is_writable(action, fname);

    if (cb == CB_PDF) {
	vwin->flags |= VWIN_CB_PDF;
    }

    make_db_toolbar(vwin, cb, del, index_button);
    build_db_content_popup(vwin, cb, del);

    listbox = database_window(vwin);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), listbox, TRUE, TRUE, 0);

    if (action == REMOTE_SERIES) {
	GtkWidget *hbox;

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
	vwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    }

    if (action == NATIVE_SERIES) {
	err = add_local_db_series_list(vwin);
    } else if (action == REMOTE_SERIES) {
	err = add_remote_db_series_list(vwin, buf);
    } else if (action == RATS_SERIES) {
	err = add_rats_db_series_list(vwin);
    } else if (action == PCGIVE_SERIES) {
	err = add_pcgive_db_series_list(vwin);
    }

    if (err) {
	gtk_widget_destroy(vwin->main);
    } else {
	gtk_widget_show_all(vwin->main);
	maybe_adjust_descrip_column(vwin);
	listbox_select_first(vwin);
    }

    return err;
}

void open_rats_window (char *fname)
{
    make_db_index_window(RATS_SERIES, fname, NULL, 0);
}

void open_bn7_window (char *fname)
{
    make_db_index_window(PCGIVE_SERIES, fname, NULL, 0);
}

static int check_serinfo (char *str, char *sername, int *nobs)
{
    char stobs[OBSLEN], endobs[OBSLEN];
    char pdc = 0;
    int err = 0;

    *stobs = *endobs = '\0';
    *nobs = 0;

    if (!isalpha((unsigned char) *sername) ||
	sscanf(str, "%c %10s - %10s %*s = %d",
	       &pdc, stobs, endobs, nobs) != 4 ||
	!isdigit((unsigned char) *stobs) ||
	!isdigit((unsigned char) *endobs) ||
	(pdc != 'M' && pdc != 'A' && pdc != 'Q' && pdc != 'U' &&
	 pdc != 'D' && pdc != 'B' && pdc != 'S')) {
	errbox_printf(_("Database parse error at variable '%s'"), sername);
	fprintf(stderr, "%s: stobs='%s', endobs='%s', pdc='%c', nobs = %d\n",
		sername, stobs, endobs, pdc, *nobs);
	err = 1;
    }

    return err;
}

static char *start_trim (char *s)
{
    while (*s == ' ') s++;

    return s;
}

static void do_db_drag (GtkWidget *w, GdkDragContext *context,
			GtkSelectionData *sel, guint info, guint t,
			windata_t *vwin)
{
    vwin_drag_src = vwin;
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_INTEGER, 8,
			   (const guchar *) &vwin, sizeof vwin);
}

static void db_drag_connect (windata_t *vwin, int i)
{
    gtk_drag_source_set(vwin->listbox, GDK_BUTTON1_MASK,
			&gretl_drag_targets[i],
			1, GDK_ACTION_COPY);
    g_signal_connect(G_OBJECT(vwin->listbox), "drag-data-get",
		     G_CALLBACK(do_db_drag),
		     vwin);
}

#define DB_LINELEN 512

static int add_local_db_series_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;
    gchar *row[3];
    char sername[VNAMELEN];
    char line1[DB_LINELEN], line2[128];
    char dbidx[MAXLEN];
    FILE *fp;
    size_t n;
    int offset = 0;
    int err = 0;

    strcpy(dbidx, vwin->fname);
    strcat(dbidx, ".idx");
    fp = gretl_fopen(dbidx, "r");

    if (fp == NULL) {
	file_read_errbox(dbidx);
	return 1;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    while (fgets(line1, sizeof line1, fp) != NULL && !err) {
	int len, nobs = 0;

	if (*line1 == '#') {
	    continue;
	}

	len = strlen(line1);
	if (line1[len-1] != '\n') {
	    errbox_printf("Database index line is too long: max is %d characters",
			  DB_LINELEN - 1);
	    break;
	}

	err = utf8_correct(line1);

	tailstrip(line1);
	gretl_charsub(line1, '\t', ' ');

	if (gretl_scan_varname(line1, sername) != 1) {
	    break;
	}

	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	if (fgets(line2, sizeof line2, fp) == NULL) {
	    break;
	}

	tailstrip(line2);
	row[2] = line2;

	if (!err) {
	    err = check_serinfo(line2, sername, &nobs);
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter,
			   DBCOL_VARNAME, row[0],
			   DBCOL_DESCRIP, row[1],
			   DBCOL_OBSINFO, row[2],
			   DBCOL_OFFSET,  offset * sizeof(dbnumber),
			   -1);

	offset += nobs;
    }

    fclose(fp);
    db_drag_connect(vwin, GRETL_DBSERIES_PTR);

    return 0;
}

static int add_remote_db_series_list (windata_t *vwin, char *buf)
{
    GtkListStore *store;
    GtkTreeIter iter;
    gchar *row[3];
    char sername[VNAMELEN];
    char line1[256], line2[256];
    int offset = 0;
    int n, err = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    bufgets_init(buf);

    while (bufgets(line1, sizeof line1, buf) && !err) {
	int nobs = 0;

	if (*line1 == '#') {
	    continue;
	}
	tailstrip(line1);
	gretl_charsub(line1, '\t', ' ');
	err = utf8_correct(line1);
	if (gretl_scan_varname(line1, sername) != 1) {
	    break;
	}

	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	if (bufgets(line2, sizeof line2, buf) == NULL) {
	    break;
	}

	row[2] = tailstrip(line2);
	if (!err) {
	    err = check_serinfo(line2, sername, &nobs);
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set (store, &iter,
			    DBCOL_VARNAME, row[0],
			    DBCOL_DESCRIP, row[1],
			    DBCOL_OBSINFO, row[2],
			    DBCOL_OFFSET, nobs * sizeof(dbnumber),
			    -1);

	offset += nobs;
    }

    bufgets_finalize(buf);
    db_drag_connect(vwin, GRETL_DBSERIES_PTR);

    return 0;
}

static gchar *iso_comment_to_utf8 (const gchar *src, int *err)
{
    gchar *conv = NULL;

    if (!g_utf8_validate(src, -1, NULL)) {
	GError *gerr = NULL;
	gsize wrote;

	conv = g_convert(src, -1,
			 "UTF-8",
			 "ISO-8859-1",
			 NULL, &wrote, &gerr);

	if (gerr != NULL) {
	    if (err != NULL) {
		errbox(gerr->message);
		*err = 1;
	    }
	    g_error_free(gerr);
	}
    } else {
	conv = g_strdup(src);
    }

    return conv;
}

static gchar *format_obs_info (SERIESINFO *sinfo)
{
    int pdc = '1';

    if (sinfo->pd == 4) {
	pdc = 'Q';
    } else if (sinfo->pd == 12) {
	pdc = 'M';
    } else if (sinfo->pd == 5) {
	pdc = 'B';
    } else if (sinfo->pd == 6) {
	pdc = 'S';
    } else if (sinfo->pd == 7) {
	pdc = 'D';
    } else if (sinfo->pd == 52) {
	pdc = 'W';
    }

    return g_strdup_printf("%c  %s - %s  n = %d", pdc,
			   sinfo->stobs, sinfo->endobs,
			   sinfo->nobs);
}

static void insert_and_free_dbwrapper (dbwrapper *dw, GtkWidget *w)
{
    GtkTreeView *view = GTK_TREE_VIEW(w);
    GtkListStore *store;
    GtkTreeIter iter;
    int i, err = 0;
    int *perr = &err;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<dw->nv; i++) {
	gchar *obsinfo;

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, DBCOL_VARNAME,
			   dw->sinfo[i].varname, -1);

	if (*dw->sinfo[i].descrip != 0) {
	    gchar *comment;

	    comment = iso_comment_to_utf8(dw->sinfo[i].descrip, perr);
	    if (perr != NULL && *perr) {
		/* don't keep displaying error messages */
		perr = NULL;
	    }
	    if (comment != NULL) {
		gtk_list_store_set(store, &iter, DBCOL_DESCRIP,
				   comment, -1);
		g_free(comment);
	    }
	}

	obsinfo = format_obs_info(&dw->sinfo[i]);
	gtk_list_store_set(store, &iter, DBCOL_OBSINFO, obsinfo, -1);
	g_free(obsinfo);

	gtk_list_store_set(store, &iter, DBCOL_OFFSET, dw->sinfo[i].offset, -1);
    }

    dbwrapper_destroy(dw);
}

static int add_rats_db_series_list (windata_t *vwin)
{
    FILE *fp;
    dbwrapper *dw;

    fp = gretl_fopen(vwin->fname, "rb");
    if (fp == NULL) {
	file_read_errbox(vwin->fname);
	return 1;
    }

    /* extract catalog from RATS file */
    dw = read_rats_db(vwin->fname, fp);
    fclose(fp);

    if (dw == NULL) {
	gui_errmsg(1);
	return 1;
    }

    insert_and_free_dbwrapper(dw, vwin->listbox);
    vwin->active_var = 0;
    db_drag_connect(vwin, GRETL_DBSERIES_PTR);

    return 0;
}

static int add_pcgive_db_series_list (windata_t *vwin)
{
    char in7name[FILENAME_MAX];
    FILE *fp;
    dbwrapper *dw;

    *in7name = 0;
    strncat(in7name, vwin->fname, FILENAME_MAX - 5);
    strcat(in7name, ".in7");

    fp = gretl_fopen(in7name, "r");
    if (fp == NULL) {
	file_read_errbox(in7name);
	return 1;
    }

    /* extract catalog from PcGive file */
    dw = read_pcgive_db(vwin->fname, fp);
    fclose(fp);

    if (dw == NULL) {
	gui_errmsg(1);
	return 1;
    }

    insert_and_free_dbwrapper(dw, vwin->listbox);
    vwin->active_var = 0;
    db_drag_connect(vwin, GRETL_DBSERIES_PTR);

    return 0;
}

static GtkWidget *database_window (windata_t *vwin)
{
    const char *titles[] = {
	N_("Name"),
	N_("Description"),
	N_("Observations")
    };
    GType types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_INT
    };
    GtkWidget *box;

    /* FIXME column widths: we should ensure we show at
       least part of the observation-info column */

    box = gtk_vbox_new(FALSE, 0);
    vwin_add_list_box(vwin, GTK_BOX(box), 4, 1, types, titles, 0);
    g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
		     G_CALLBACK(popup_menu_handler),
		     vwin->popup);
    gtk_widget_show(box);

    return box;
}

static void
add_series_to_list (GtkTreeModel *model, GtkTreePath *path,
		    GtkTreeIter *iter, int *list)
{
    int row = tree_path_get_row_number(path);

    list[0] += 1;
    list[list[0]] = row;
}

static void db_fill_selection_list (windata_t *vwin, int *list)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    GtkTreeSelection *sel;

    sel = gtk_tree_view_get_selection(view);
    gtk_tree_selection_selected_foreach(sel,
					(GtkTreeSelectionForeachFunc)
					add_series_to_list, list);
}

static int *db_get_selection_list (windata_t *vwin)
{
    int n = vwin_selection_count(vwin, NULL);
    int *list = NULL;

    if (n == 0) {
	return NULL;
    }

    list = gretl_list_new(n);
    if (list == NULL) {
	nomem();
	return NULL;
    }

    if (n == 1) {
	list[1] = vwin->active_var;
    } else {
	list[0] = 0;
	db_fill_selection_list(vwin, list);
    }

    return list;
}

static int db_role_to_dbtype (int role)
{
    if (role == PCGIVE_SERIES) {
	return GRETL_PCGIVE_DB;
    } else if (role == REMOTE_SERIES) {
	return GRETL_NATIVE_DB_WWW;
    } else if (role == RATS_SERIES) {
	return GRETL_RATS_DB;
    } else {
	return GRETL_NATIVE_DB;
    }
}

static dbwrapper *get_db_series_info (windata_t *vwin, int action)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    int *rowlist = NULL;
    int i, nsel, row = 0;
    char stobs[OBSLEN], endobs[OBSLEN];
    char pdc;
    dbwrapper *dw;
    int err = 0;

    /* get list of selected rows in database window */
    rowlist = db_get_selection_list(vwin);
    if (rowlist == NULL) {
	return NULL;
    }

    /* count of selected rows */
    nsel = rowlist[0];

    /* construct wrapper with space for info on @nsel series */
    dw = dbwrapper_new(nsel, vwin->fname, db_role_to_dbtype(vwin->role));
    if (dw == NULL) {
	free(rowlist);
	return NULL;
    }

    dw->nv = nsel;

    /* loop across the selected rows, retrieve info on the
       associated db series and enter into the wrapper
    */

    for (i=0; i<rowlist[0]; i++) {
	SERIESINFO *sinfo = &dw->sinfo[i];
	gchar *tmp = NULL;

	row = rowlist[i+1];
	tree_view_get_int(view, row, DBCOL_OFFSET, &sinfo->offset);

	*sinfo->varname = '\0';
	tree_view_get_string(view, row, DBCOL_VARNAME, &tmp);
	strncat(sinfo->varname, tmp, VNAMELEN - 1);
	g_free(tmp);

	tmp = NULL;
	tree_view_get_string(view, row, DBCOL_DESCRIP, &tmp);
	if (tmp != NULL) {
	    series_info_set_description(sinfo, tmp);
	    g_free(tmp);
	}

	tmp = NULL;
	tree_view_get_string(view, row, DBCOL_OBSINFO, &tmp);
	if (sscanf(tmp, "%c %10s %*s %10s %*s %*s %d",
		   &pdc, stobs, endobs, &sinfo->nobs) != 4) {
	    errbox(_("Failed to parse series information"));
	    err = 1;
	    goto bailout;
	}
	g_free(tmp);

	sinfo->pd = 1;
	sinfo->undated = 0;

	if (pdc == 'M') {
	    sinfo->pd = 12;
	} else if (pdc == 'Q') {
	    sinfo->pd = 4;
	} else if (pdc == 'B') {
	    sinfo->pd = 5;
	} else if (pdc == 'S') {
	    sinfo->pd = 6;
	} else if (pdc == 'D') {
	    sinfo->pd = 7;
	} else if (pdc == 'W') {
	    sinfo->pd = 52;
	} else if (pdc == 'U') {
	    sinfo->undated = 1;
	}

	if (i > 0 && sinfo->pd != dw->sinfo[0].pd) {
	    /* this shouldn't happen */
	    errbox("Can't operate on series with different frequencies");
	    err = 1;
	    goto bailout;
	}

	if (strchr(stobs, '/')) {
	    /* daily data */
	    char *q = stobs;
	    char *p = strchr(stobs, '/');

	    if (p - q == 4) {
		strcpy(sinfo->stobs, q + 2);
	    }
	    q = endobs;
	    p = strchr(endobs, '/');
	    if (p && p - q == 4) {
		strcpy(sinfo->endobs, q + 2);
	    }
	} else {
	    sinfo->stobs[0] = 0;
	    sinfo->endobs[0] = 0;
	    strncat(sinfo->stobs, stobs, OBSLEN - 1);
	    strncat(sinfo->endobs, endobs, OBSLEN - 1);
	}
    }

 bailout:

    if (err) {
	dbwrapper_destroy(dw);
	dw = NULL;
    }

    free(rowlist);

    return dw;
}

/* the following two functions are used when a database has been
   specified on the gretl command line */

gboolean open_named_db_index (char *dbname)
{
    gboolean ret = FALSE;
    int action;
    FILE *fp;

    if (has_suffix(dbname, ".rat")) {
	action = RATS_SERIES;
    } else if (has_suffix(dbname, ".bn7")) {
	action = PCGIVE_SERIES;
    } else {
	action = NATIVE_SERIES;
    }

    if (action == NATIVE_SERIES && !strcmp(dbname, "dbnomics")) {
	warnbox("Sorry, this access to dbnomics not ready yet");
	return ret;
    }

    fp = gretl_fopen(dbname, "rb");

    if (fp == NULL && action == NATIVE_SERIES &&
	!has_suffix(dbname, ".bin")) {
	strcat(dbname, ".bin");
	fp = gretl_fopen(dbname, "rb");
    }

    if (fp == NULL) {
	file_read_errbox(dbname);
    } else {
	fclose(fp);
	make_db_index_window(action, dbname, NULL, 0);
	ret = TRUE;
    }

    return ret;
}

gboolean open_named_remote_db_index (char *dbname)
{
    gboolean ret = FALSE;
    char *getbuf = NULL;
    int err;

    err = retrieve_remote_db_index(dbname, &getbuf);

    if (err) {
	show_network_error(NULL);
    } else if (getbuf != NULL && !strncmp(getbuf, "Couldn't open", 13)) {
	errbox(getbuf);
    } else {
	err = make_db_index_window(REMOTE_SERIES, dbname, getbuf, 0);
	if (!err) {
	    ret = TRUE;
	}
    }

    free(getbuf);

    return ret;
}

void dbnomics_specific_series (GtkAction *action, gpointer data)
{
    char *datacode = NULL;
    int resp;

    if (data != NULL) {
	windata_t *vwin = (windata_t *) data;

	resp = dbnomics_dialog(&datacode, vwin->main);
    } else {
	resp = dbnomics_dialog(&datacode, NULL);
    }

    if (!canceled(resp)) {
	dbnomics_get_series_call(datacode);
    }

    free(datacode);
}

void open_dbnomics_series (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *scode = NULL, *path = NULL;
    gchar *datacode;

    path = g_object_get_data(G_OBJECT(vwin->listbox), "path");
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox),
			 vwin->active_var, COL_DBNAME, &scode);

    datacode = g_strdup_printf("%s/%s", path, scode);
    dbnomics_get_series_call(datacode);
    g_free(datacode);
}

static int dbn_general_search_results (const gchar *key,
				       gretl_array *a)
{
    int n = gretl_array_get_length(a);
    PRN *prn = NULL;
    int err, n_ok = 0;

    if (n == 0) {
	return 0;
    }

    err = bufopen(&prn);

    if (!err) {
	const char *pcode, *dcode, *name;
	gretl_bundle *b;
	int i;

	pprintf(prn, _("Results of DB.NOMICS search on '%s'\n"), key);
	pputs(prn, "Provider/Dataset : description\n\n");

	for (i=0; i<n; i++) {
	    b = gretl_array_get_bundle(a, i);
	    pcode = gretl_bundle_get_string(b, "provider_code", NULL);
	    dcode = gretl_bundle_get_string(b, "code", NULL);
	    name = gretl_bundle_get_string(b, "name", NULL);
	    if (pcode != NULL && dcode != NULL && name != NULL) {
		pprintf(prn, "%d <@dbn=\"%s/%s\"> : %s\n\n", i+1,
			pcode, dcode, name);
		n_ok++;
	    }
	}
	if (n_ok > 0) {
	    const char *title = "gretl: DB.NOMICS search";

	    view_buffer(prn, 78, 350, title, VIEW_DBSEARCH, NULL);
	} else {
	    gretl_print_destroy(prn);
	}
    }

    return n_ok;
}

static int dbn_dataset_search_results (const char *key,
				       const char *prov,
				       const char *dset,
				       int offset,
				       gretl_array *a,
				       windata_t *prev_vwin)
{
    int n = gretl_array_get_length(a);
    PRN *prn = NULL;
    int err, n_ok = 0;

    if (n == 0) {
	return 0;
    }

    err = bufopen(&prn);

    if (!err) {
	const char *scode, *name;
	gretl_bundle *b;
	int more = 0;
	int i, ntot;

	b = gretl_array_get_bundle(a, 0);
	ntot = gretl_bundle_get_int(b, "ntot", NULL);

	pprintf(prn, _("DB.NOMICS search on '%s' in dataset %s/%s\n"),
		key, prov, dset);
	pprintf(prn, "Matching series %d to %d of %d", offset + 1,
		offset + n, ntot);

	if (ntot > offset + n) {
	    pputs(prn, " [ <@dbn=\"_NEXT_\"> ]\n\n");
	    more = 1;
	} else {
	    pputs(prn, "\n\n");
	}

	for (i=0; i<n; i++) {
	    b = gretl_array_get_bundle(a, i);
	    scode = gretl_bundle_get_string(b, "code", NULL);
	    name = gretl_bundle_get_string(b, "name", NULL);
	    if (scode != NULL && name != NULL) {
		pprintf(prn, "%d <@dbn=\"%s;%s/%s/%s\"> : %s\n\n", i+1,
			scode, prov, dset, scode, name);
		n_ok++;
	    }
	}

	if (prev_vwin != NULL) {
	    /* we should replace the content of @prev_vwin */
	    char *buf = gretl_print_steal_buffer(prn);

	    /* update the offset record on prev_vwin->text */
	    widget_set_int(prev_vwin->text, "offset", offset);
	    /* clear and replace with new results */
	    textview_clear_text(prev_vwin->text);
	    textview_set_text_dbsearch(prev_vwin, buf);
	    gretl_print_destroy(prn);
	    return 1;
	} else if (n_ok > 0) {
	    /* we'll open a new window */
	    const char *title = "gretl: DB.NOMICS search";
	    windata_t *vwin;

	    vwin = view_buffer(prn, 78, 350, title, VIEW_DBSEARCH, NULL);
	    if (vwin != NULL && more) {
		widget_set_int(vwin->text, "offset", offset);
		g_object_set_data_full(G_OBJECT(vwin->text), "key",
				       g_strdup(key), g_free);
		g_object_set_data_full(G_OBJECT(vwin->text), "prov",
				       g_strdup(prov), g_free);
		g_object_set_data_full(G_OBJECT(vwin->text), "dset",
				       g_strdup(dset), g_free);
	    }
	} else {
	    gretl_print_destroy(prn);
	}
    }

    return n_ok;
}

#define SEARCH_CHUNK 80

/* The @key string is passed here when "all DB.NOMICS"
   is selected as the search space in the dbnomics
   providers window, or when "this database" is selected
   in a dbnomics dataset window.
*/

void dbnomics_search (gchar *key, windata_t *vwin)
{
    gretl_array *a = NULL;
    int n_found = 0;
    int err = 0;

    if (vwin->role == VIEW_DBSEARCH) {
	/* we're called in "next results" mode */
	const gchar *key, *prov, *dset;
	int offset;

	key = g_object_get_data(G_OBJECT(vwin->text), "key");
	prov = g_object_get_data(G_OBJECT(vwin->text), "prov");
	dset = g_object_get_data(G_OBJECT(vwin->text), "dset");
	offset = widget_get_int(vwin->text, "offset");

	offset += SEARCH_CHUNK;
	a = dbnomics_search_call(key, prov, dset, SEARCH_CHUNK, offset, &err);
	if (!err) {
	    n_found = dbn_dataset_search_results(key, prov, dset, offset, a, vwin);
	}
	key = NULL; /* don't free it! */
    } else if (vwin->role == DBNOMICS_DB) {
	/* searching "selected dataset" in datasets window */
	const gchar *prov = g_object_get_data(G_OBJECT(vwin->listbox),
					       "provider");
	gchar *dset = NULL;

	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox),
			     vwin->active_var, COL_DBNAME, &dset);
	a = dbnomics_search_call(key, prov, dset, SEARCH_CHUNK, 0, &err);
	if (!err) {
	    n_found = dbn_dataset_search_results(key, prov, dset, 0, a, NULL);
	}
	g_free(dset);
    } else if (vwin->role == DBNOMICS_SERIES) {
	/* searching from a particular dataset window */
	const gchar *path = g_object_get_data(G_OBJECT(vwin->listbox),
					      "path");
	const gchar *p = strchr(path, '/');
	const gchar *dset = p + 1;
	gchar *prov = g_strndup(path, p - path);

	a = dbnomics_search_call(key, prov, dset, SEARCH_CHUNK, 0, &err);
	if (!err) {
	    n_found = dbn_dataset_search_results(key, prov, dset, 0, a, NULL);
	}
	g_free(prov);
    } else {
	/* top-level search */
	a = dbnomics_search_call(key, NULL, NULL, SEARCH_CHUNK, 0, &err);
	if (!err) {
	    n_found = dbn_general_search_results(key, a);
	}
    }

    if (!err && n_found == 0) {
	warnbox(_("No matches were found"));
    }

    gretl_array_destroy(a);

    /* only if this arg is a GTK-allocated string */
    if (key != NULL) {
	g_free(key);
    }
}

void open_db_index (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *fname = NULL, *dbdir = NULL;
    char dbfile[MAXLEN];
    int action = NATIVE_SERIES;
    int idx = 0;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox),
			 vwin->active_var, COL_DBNAME, &fname);

    if (has_suffix(fname, ".rat")) {
	action = RATS_SERIES;
    }

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox),
			 vwin->active_var, (action == RATS_SERIES)? 1 : 2,
			 &dbdir);

    gretl_build_path(dbfile, dbdir, fname, NULL);
    g_free(fname);
    g_free(dbdir);

    if (action == NATIVE_SERIES) {
	GtkTreeModel *mod;

	mod = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
	idx = tree_model_count_rows(mod) > 1;
    }

    make_db_index_window(action, dbfile, NULL, idx);

    if (action == NATIVE_SERIES) {
	/* close the window from which this db was selected */
	gtk_widget_destroy(vwin->main);
    }
}

void open_remote_db_index (GtkWidget *w, gpointer data)
{
    GtkTreeIter iter;
    GtkTreeModel *model;
    GtkTreeSelection *sel;
    char *getbuf = NULL;
    gchar *fname = NULL;
    windata_t *vwin = (windata_t *) data;
    int err;

    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_selection_get_selected(sel, &model, &iter)) {
	return;
    }

    gtk_tree_model_get(model, &iter, 0, &fname, -1);
    if (fname == NULL || *fname == '\0') {
	g_free(fname);
	return;
    }

    update_statusline(vwin, _("Retrieving data..."));
    err = retrieve_remote_db_index(fname, &getbuf);

    if (err) {
	show_network_error(vwin);
    } else {
	update_statusline(vwin, "OK");
	make_db_index_window(REMOTE_SERIES, fname, getbuf, 0);
    }

    g_free(fname);
    free(getbuf);
}

void open_dbnomics_provider (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    GtkTreeIter iter;
    GtkTreeModel *model;
    GtkTreeSelection *sel;
    gchar *pname = NULL;

    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_selection_get_selected(sel, &model, &iter)) {
	return;
    }

    gtk_tree_model_get(model, &iter, 0, &pname, -1);
    if (pname != NULL && *pname != '\0') {
	display_files(DBNOMICS_DB, pname);
    }
    g_free(pname);
}

void open_dbnomics_dataset (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    GtkTreeIter iter;
    GtkTreeModel *model;
    GtkTreeSelection *sel;
    gchar *provider, *arg;
    gchar *dsname = NULL;

    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_selection_get_selected(sel, &model, &iter)) {
	return;
    }
    gtk_tree_model_get(model, &iter, 0, &dsname, -1);
    if (dsname == NULL || *dsname == '\0') {
	g_free(dsname);
	return;
    }
    provider = g_object_get_data(G_OBJECT(vwin->listbox), "provider");
    if (provider == NULL) {
	g_free(dsname);
	return;
    }
    arg = g_strdup_printf("%s/%s", provider, dsname);
    g_free(dsname);
    display_files(DBNOMICS_SERIES, arg);
    g_free(arg);
}

#define INFOLEN 100

static int parse_db_header (const char *buf, unsigned *idxlen,
			    unsigned *datalen, unsigned *cblen,
			    int *pdfdoc)
{
    char *p;
    int err = 0;

    *cblen = 0;

    /* length of index file (required) */
    if (sscanf(buf, "%u", idxlen) != 1) {
	err = 1;
    }

    /* length of data (required under "new" system) */
    if (!err) {
	p = strchr(buf, '\n');
	if (p == NULL) {
	    err = 1;
	} else if (sscanf(p + 1, "%u", datalen) != 1) {
	    err = 1;
	}
    }

    /* length of codebook (optional) */
    if (!err) {
	p = strchr(p + 1, '\n');
	if (p != NULL) {
	    int cbl;

	    if (sscanf(p + 1, "%u", &cbl) == 1) {
		*cblen = cbl;
		if (strstr(p, ".pdf")) {
		    *pdfdoc = 1;
		}
	    }
	}
    }

    return err;
}

static int ggz_extract (char *ggzname)
{
    gzFile fgz = NULL;
    FILE *fidx = NULL, *fbin = NULL, *fcbk = NULL;
    unsigned idxlen, datalen, cblen = 0;
    int bgot, bytesleft, pdfdoc = 0;
    char idxname[MAXLEN], binname[MAXLEN], cbname[MAXLEN];
    char gzbuf[GRETL_BUFSIZE];
#if G_BYTE_ORDER == G_BIG_ENDIAN
    netfloat nf;
    float val;
#endif
    int err = 0;

    switch_ext(idxname, ggzname, "idx");
    switch_ext(binname, ggzname, "bin");
    cbname[0] = '\0';

    fgz = gzopen(ggzname, "rb");
    if (fgz == NULL) {
	file_read_errbox(ggzname);
        return E_FOPEN;
    }

    fidx = gretl_fopen(idxname, "wb");
    if (fidx == NULL) {
	file_write_errbox(idxname);
	err = E_FOPEN;
	goto bailout;
    }

    fbin = gretl_fopen(binname, "wb");
    if (fbin == NULL) {
	file_write_errbox(binname);
	err = E_FOPEN;
	goto bailout;
    }

    memset(gzbuf, 0, GRETL_BUFSIZE);
    gzread(fgz, gzbuf, INFOLEN);

    if (parse_db_header(gzbuf, &idxlen, &datalen, &cblen, &pdfdoc)) {
	fputs("Error reading info buffer: failed to get byte counts\n",
	      stderr);
	fprintf(stderr, "bad infobuf:\n%s\n", gzbuf);
	err = 1;
	goto bailout;
    }

    if (cblen > 0) {
	switch_ext(cbname, ggzname, pdfdoc ? "pdf" : "cb");
	fcbk = gretl_fopen(cbname, "wb");
	if (fcbk == NULL) {
	    cblen = pdfdoc = 0;
	}
    }

    bytesleft = idxlen;
    while (bytesleft > 0) {
	memset(gzbuf, 0, GRETL_BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > GRETL_BUFSIZE)?
		      GRETL_BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fidx);
    }

    if (bytesleft > 0) {
	fputs("Error reading database info buffer\n", stderr);
	err = 1;
	goto bailout;
    }

    bytesleft = datalen;

    while (bytesleft > 0) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
        if ((bgot = gzread(fgz, gzbuf, sizeof(long) + sizeof(short))) > 0) {
	    /* read "netfloats" and write floats */
	    memcpy(&(nf.frac), gzbuf, sizeof(long));
	    memcpy(&(nf.exp), gzbuf + sizeof(long), sizeof(short));
	    val = retrieve_float(nf);
	    fwrite(&val, sizeof(float), 1, fbin);
	    bytesleft -= sizeof(dbnumber);
	} else break;
#else
	memset(gzbuf, 0, GRETL_BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > GRETL_BUFSIZE)?
		      GRETL_BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fbin);
#endif
    }

    if (bytesleft > 0) {
	fputs("Error reading database data\n", stderr);
	err = 1;
	goto bailout;
    }

    bytesleft = cblen;

    while (bytesleft > 0) {
	memset(gzbuf, 0, GRETL_BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > GRETL_BUFSIZE)?
		      GRETL_BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fcbk);
    }

    if (bytesleft > 0) {
	fputs("Error reading database codebook\n", stderr);
	err = 1;
    }

 bailout:

    if (fgz != NULL) gzclose(fgz);
    if (fidx != NULL) fclose(fidx);
    if (fbin != NULL) fclose(fbin);
    if (fcbk != NULL) fclose(fcbk);

    if (err) {
	/* clean up botched files */
	gretl_remove(idxname);
	gretl_remove(binname);
	if (cbname[0] != '\0') {
	    gretl_remove(cbname);
	}
    }

    gretl_remove(ggzname);

    return err;
}

static void offer_db_open (char *target, windata_t *vwin)
{
    int resp = yes_no_dialog("gretl",
			     _("Database installed.\n"
			       "Open it now?"),
			     vwin_toplevel(vwin));

    if (resp == GRETL_YES) {
	char dbpath[MAXLEN];

	strcpy(dbpath, target);
	strcpy(strrchr(dbpath, '.'), ".bin");
	open_named_db_index(dbpath);
    }
}

static int get_target_in_home (char *targ, int code,
			       const char *objname,
			       const char *ext)
{
#ifdef OS_OSX
    const char *savedir = gretl_app_support_dir();
#else
    const char *savedir = gretl_dotdir();
#endif
    int err = 0;

    if (savedir == NULL || *savedir == '\0') {
	err = E_FOPEN;
    } else {
	int subdir = 0;

	if (code == REMOTE_FUNC_FILES) {
	    sprintf(targ, "%sfunctions", savedir);
	    subdir = 1;
	} else if (code == REMOTE_DB) {
	    sprintf(targ, "%sdb", savedir);
	    subdir = 1;
	} else if (code == REMOTE_DATA_PKGS) {
	    sprintf(targ, "%sdata", savedir);
	    subdir = 1;
	} else {
	    sprintf(targ, "%s%s%s", savedir, objname, ext);
	}

	if (subdir) {
	    err = gretl_mkdir(targ);
	    if (!err) {
		strcat(targ, SLASHSTR);
		strcat(targ, objname);
		strcat(targ, ext);
	    }
	}
    }

    return err;
}

#ifndef OS_OSX

static void get_system_target (char *targ, int code,
			       const char *objname,
			       const char *ext)
{
    if (code == REMOTE_DB) {
	get_default_dir_for_action(targ, SAVE_REMOTE_DB);
    } else if (code == REMOTE_DATA_PKGS) {
	get_default_dir_for_action(targ, SAVE_DATA_PKG);
    } else if (code == REMOTE_FUNC_FILES) {
	get_default_dir_for_action(targ, SAVE_FUNCTIONS);
    }

    strcat(targ, objname);
    strcat(targ, ext);
}

#endif

enum {
    REAL_INSTALL,
    TMP_INSTALL
};

/* try to find a suitable path, for which the user has write
   permission, for installing a database, collection of
   data files, or function package
*/

static char *get_writable_target (int code, char *objname)
{
    const char *ext;
    char *targ;
    int done_home = 0;
    int err = 0;

#if 0
    fprintf(stderr, "get_writable_target, starting\n");
#endif

    targ = mymalloc(MAXLEN);
    if (targ == NULL) {
	return NULL;
    }

    *targ = '\0';

    if (code == REMOTE_DB) {
	ext = ".ggz";
    } else if (code == REMOTE_FUNC_FILES) {
	ext = ".gfn";
    } else {
	ext = ".tar.gz";
    }

#ifdef OS_OSX
    /* we prefer writing to ~/Library/Application Support
       rather than /Applications/Gretl.app
    */
    err = get_target_in_home(targ, code, objname, ext);
    done_home = 1;
#else
    get_system_target(targ, code, objname, ext);
#endif

    if (!err) {
	err = gretl_test_fopen(targ, "w");
	if (err == EACCES && !done_home) {
	    /* permissions problem: write to home dir instead */
	    err = get_target_in_home(targ, code, objname, ext);
	}
    }

    if (err) {
	file_write_errbox(targ);
	free(targ);
	targ = NULL;
    }

#if 0
    fprintf(stderr, "writable targ: '%s'\n", targ);
#endif

    return targ;
}

static int unpack_book_data (const char *fname)
{
    char *p, path[FILENAME_MAX];
    int err = 0;

    strcpy(path, fname);
    p = path_last_slash(path);
    if (p != NULL) {
	*p = '\0';
    }

    if (gretl_chdir(path) != 0) {
	err = E_FOPEN;
    }

    if (!err) {
	err = gretl_untar(fname);
    }

    return err;
}

static gchar *make_gfn_path (const char *pkgname)
{
    const char *fpath = gretl_function_package_path();

    return g_strdup_printf("%s%s%c%s.gfn", fpath,
			   pkgname, SLASH, pkgname);
}

#define STATUS_COLUMN  5
#define ZIPFILE_COLUMN 6

/* note : @vwin here is the source viewer window displaying the
   remote file (database, or datafiles package, or function package)
   that is to be installed onto the local machine.
*/

void install_file_from_server (GtkWidget *w, windata_t *vwin)
{
    gchar *objname = NULL;
    char *targ = NULL;
    gboolean zipfile = FALSE;
    int err = 0;

    /* note: addon files are handled separately, by the function
       install_addon_callback() in datafiles.c
    */

    if (vwin->role == REMOTE_DB) {
	/* database files */
	GtkTreeSelection *sel;
	GtkTreeModel *model;
	GtkTreeIter iter;

	sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
	if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
	    gtk_tree_model_get(model, &iter, 0, &objname, -1);
	}
    } else {
	/* datafiles or function package */
	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox),
			     vwin->active_var, 0, &objname);
	if (vwin->role == REMOTE_FUNC_FILES) {
	    tree_view_get_bool(GTK_TREE_VIEW(vwin->listbox),
			       vwin->active_var, ZIPFILE_COLUMN, &zipfile);
	}
    }

    if (objname == NULL || *objname == '\0') {
	g_free(objname);
	return;
    }

    if (!zipfile) {
	targ = get_writable_target(vwin->role, objname);
	if (targ == NULL) {
	    g_free(objname);
	    return;
	}
    }

    if (vwin->role == REMOTE_FUNC_FILES) {
	if (zipfile) {
	    const char *path = gretl_function_package_path();
	    gchar *basename = g_strdup_printf("%s.zip", objname);
	    gchar *fullname;

	    fullname = g_strdup_printf("%s%s", path, basename);
	    err = retrieve_remote_function_package(basename, fullname);
	    if (!err) {
		err = gretl_unzip_into(fullname, path);
		gretl_remove(fullname);
	    }
	    g_free(fullname);
	    g_free(basename);
	} else {
	    err = retrieve_remote_function_package(objname, targ);
	}
    } else if (vwin->role == REMOTE_DATA_PKGS) {
	gchar *tarname = g_strdup_printf("%s.tar.gz", objname);

	err = retrieve_remote_datafiles_package(tarname, targ);
	g_free(tarname);
    } else if (vwin->role == REMOTE_DB) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
	err = retrieve_remote_db(objname, targ, GRAB_NBO_DATA);
#else
	err = retrieve_remote_db(objname, targ, GRAB_DATA);
#endif
    }

    if (err) {
	show_network_error(NULL);
    } else {
	windata_t *local = get_local_viewer(vwin->role);

	if (vwin->role == REMOTE_FUNC_FILES) {
	    int notified = 0;

	    if (zipfile) {
		gchar *gfnpath = make_gfn_path(objname);

		notified = gui_function_pkg_query_register(gfnpath, vwin->main);
		g_free(gfnpath);
	    } else {
		notified = gui_function_pkg_query_register(targ, vwin->main);
	    }
	    if (!notified) {
		infobox(_("Installed"));
	    }
	    list_store_set_string(GTK_TREE_VIEW(vwin->listbox),
				  vwin->active_var, STATUS_COLUMN,
				  _("Up to date"));
	    if (local != NULL) {
		populate_filelist(local, NULL);
	    }
	} else if (vwin->role == REMOTE_DATA_PKGS) {
	    fprintf(stderr, "downloaded '%s'\n", targ);
	    err = unpack_book_data(targ);
	    remove(targ);
	    if (err) {
		errbox(_("Error unzipping compressed data"));
	    } else {
		infobox("Restart gretl to access this database");
	    }
	} else {
	    /* gretl-zipped database package */
	    fprintf(stderr, "downloaded '%s'\n", targ);
	    err = ggz_extract(targ);
	    if (err) {
		if (err != E_FOPEN) {
		    errbox(_("Error unzipping compressed data"));
		}
	    } else {
		tree_store_set_string(GTK_TREE_VIEW(vwin->listbox),
				      vwin->active_var, 2,
				      _("Up to date"));
		if (local != NULL) {
		    populate_filelist(local, NULL);
		} else {
		    offer_db_open(targ, vwin);
		}
	    }
	}
    }

    g_free(objname);
    free(targ);
}

/* called from within datafiles.c, when dragging a
   remote database or function package from its
   "on server" window to the associated local
   window */

void drag_file_from_server (guint info)
{
    windata_t *vwin = NULL;

    if (info == GRETL_REMOTE_DB_PTR ||
	info == GRETL_REMOTE_FNPKG_PTR) {
	vwin = vwin_drag_src;
	vwin_drag_src = NULL;
    }

    if (vwin != NULL) {
	install_file_from_server(NULL, vwin);
    }
}

/* Called when the "install" command is used to install a function
   package via console or script: try to sync the local and/or remote
   function-package windows, if they're open. Also present the
   package's menu-attachment option, if any.
*/

void maybe_update_pkgview (const char *filename,
			   const char *pkgname,
			   int zipfile,
			   GtkWidget *parent)
{
    windata_t *vwin;

    /* update local package browser? */
    vwin = get_browser_for_role(FUNC_FILES, NULL);
    if (vwin != NULL) {
	populate_filelist(vwin, NULL);
    }

    /* update remote package browser? */
    vwin = get_browser_for_role(REMOTE_FUNC_FILES, NULL);
    if (vwin != NULL && find_package_in_viewer(vwin, pkgname)) {
	list_store_set_string(GTK_TREE_VIEW(vwin->listbox),
			      vwin->active_var, STATUS_COLUMN,
			      _("Up to date"));
    }

    if (parent != NULL) {
	/* offer menu attachment if applicable */
	if (zipfile) {
	    gchar *gfnpath = make_gfn_path(pkgname);

	    gui_function_pkg_query_register(gfnpath, parent);
	    g_free(gfnpath);
	} else {
	    gui_function_pkg_query_register(filename, parent);
	}
    }
}

void pkg_info_from_server (GtkWidget *w, windata_t *vwin)
{
    static int idx;
    gchar *path, *objname = NULL;
    int zipfile = 0;
    int err = 0;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox),
			 vwin->active_var, 0, &objname);
    tree_view_get_bool(GTK_TREE_VIEW(vwin->listbox),
		       vwin->active_var, ZIPFILE_COLUMN, &zipfile);

    path = g_strdup_printf("%sdltmp.%d", gretl_dotdir(), idx++);

    if (zipfile) {
	gchar *zipname = g_strdup_printf("%s.zip", objname);

	err = retrieve_remote_gfn_content(zipname, path);
	g_free(zipname);
    } else {
	err = retrieve_remote_function_package(objname, path);
    }

    if (err) {
	show_network_error(NULL);
    } else {
	display_function_package_data(objname, path, VIEW_PKG_INFO);
    }

    g_free(objname);
    g_free(path);
}

static gchar *
real_get_db_description (const char *fullname, const char *binname,
			 const char *dbdir)
{
    FILE *fp;
    char idxname[FILENAME_MAX];
    char *p, *descrip = NULL;

    if (fullname == NULL) {
	gretl_build_path(idxname, dbdir, binname, NULL);
    } else {
	strcpy(idxname, fullname);
    }

    p = strrchr(idxname, '.');
    if (p != NULL) {
	strcpy(p, ".idx");
    }

    fp = gretl_fopen(idxname, "r");

    if (fp != NULL) {
	char tmp[DB_DESCRIP_LEN + 32] = {0};

	if (fgets(tmp, sizeof tmp, fp) != NULL) {
	    if (*tmp == '#' && strlen(tmp) > 2) {
		char *s = tmp + 2;

		tailstrip(s);
		utf8_correct(s);
		descrip = g_strdup(s);
	    }
	}
	fclose(fp);
    }

    return descrip;
}

gchar *get_db_description (const char *binname)
{
    return real_get_db_description(binname, NULL, NULL);
}

int write_db_description (const char *binname, const char *descrip)
{
    FILE *fnew, *fbak;
    char idxname[FILENAME_MAX];
    char idxtmp[FILENAME_MAX];
    char tmp[72];
    char *p;
    int err = 0;

    strcpy(idxname, binname);
    p = strrchr(idxname, '.');
    if (p != NULL) {
	strcpy(p, ".idx");
    }

    strcpy(idxtmp, idxname);
    p = strrchr(idxtmp, '.');
    if (p != NULL) {
	strcpy(p, ".idxtmp");
    }

    err = copyfile(idxname, idxtmp);

    if (!err) {
	fnew = gretl_fopen(idxname, "w");
	if (fnew == NULL) {
	    err = 1;
	} else {
	    fbak = gretl_fopen(idxtmp, "r");
	    if (fbak == NULL) {
		fclose(fnew);
		err = 1;
	    }
	}
    }

    if (!err) {
	const char *p = descrip;
	char line[256];

	while (isspace((unsigned char) *p)) {
	    p++;
	}
	sprintf(tmp, "# %.61s\n", p);
	tailstrip(tmp);
	fprintf(fnew, "%s\n", tmp);
	while (fgets(line, sizeof line, fbak)) {
	    if (*line != '#') {
		fputs(line, fnew);
	    }
	}
	fclose(fnew);
	fclose(fbak);
	gretl_remove(idxtmp);
    }

    return err;
}

static int
read_db_files_in_dir (GDir *dir, int dbtype, const char *path,
		      GtkListStore *store, GtkTreeIter *iter)
{
    const gchar *fname;
    gchar *name, *descrip;
    int len, ndb = 0;

    while ((fname = g_dir_read_name(dir)) != NULL) {
	len = strlen(fname);
	if (!g_ascii_strcasecmp(fname + len - 4, ".bin")) {
	    name = g_strndup(fname, len - 4);
	    descrip = real_get_db_description(NULL, fname, path);
	    if (name != NULL && descrip != NULL) {
#if DB_SEARCH_DEBUG
		fprintf(stderr, "  found '%s'\n", name);
#endif
		gtk_list_store_append(store, iter);
		gtk_list_store_set(store, iter,
				   COL_DBNAME, name,
				   COL_DBINFO, descrip,
				   COL_DBPATH, path, -1);
		ndb++;
	    }
	    g_free(name);
	    g_free(descrip);
	}
    }

    return ndb;
}

static void get_local_object_status (const char *fname,
				     int role,
				     const char **status,
				     time_t remtime)
{
    char fullname[MAXLEN];
    struct stat fbuf;
    char **dirs;
    SearchType stype;
    int found = 0;
    int i, n_dirs;
    int err = 0;

    if (role == REMOTE_DB) {
	stype = DB_SEARCH;
    } else if (role == REMOTE_DATA_PKGS) {
	stype = DATA_SEARCH;
    } else if (role == REMOTE_FUNC_FILES) {
	stype = FUNCS_SEARCH;
    } else {
	*status = N_("Unknown: access error");
	return;
    }

    dirs = get_plausible_search_dirs(stype, &n_dirs);

    for (i=0; i<n_dirs; i++) {
	gretl_build_path(fullname, dirs[i], fname, NULL);
	errno = 0;
	if (gretl_stat(fullname, &fbuf) == 0) {
	    found = 1;
	    break;
	} else if (errno != ENOENT) {
	    err = 1;
	    break;
	}
    }

    strings_array_free(dirs, n_dirs);

    if (found) {
	double dt;

	dt = difftime(remtime, fbuf.st_ctime);
	if (dt > 360) {
	    dt = difftime(remtime, fbuf.st_mtime);
	}
	if (dt > 360) {
	    *status = N_("Not up to date");
	} else {
	    *status = N_("Up to date");
	}
    } else if (err) {
	*status = N_("Unknown: access error");
    } else {
	*status = N_("Not installed");
    }
}

static int read_remote_filetime (char *line, char *fname,
				 time_t *date, char *tbuf)
{
    char month[4], hrs[9];
    int mday, yr, mon = 0;
    const char *months[] = {
	"Jan", "Feb", "Mar", "Apr",
	"May", "Jun", "Jul", "Aug",
	"Sep", "Oct", "Nov", "Dec"
    };
    int i;

    /* We're expecting a string of the form:

       "<bytes> <day> <mon> <mday> 00:00:00 <year> <filename>"

       e.g.

       "   2719 Foo Dec 17 00:00:00 2010     addlist.gfn"

       where <mon> is 3-letter month, <mday> is 2 digits,
       and <year> is 4-digit year; in this context we actually
       discard the first two fields, @bytes and @day.
    */

    if (sscanf(line, "%*s%*s%3s%2d%8s%4d%31s",
	       month, &mday, hrs, &yr, fname) != 5) {
	return 1;
    }

    for (i=0; i<12; i++) {
	if (!strcmp(month, months[i])) {
	    mon = i;
	}
    }

    if (date != NULL) {
	struct tm mytime;

	hrs[2] = 0;

	mytime.tm_sec = 0;
	mytime.tm_min = 0;
	mytime.tm_wday = 0;
	mytime.tm_yday = 0;
	mytime.tm_isdst = -1;
	mytime.tm_hour = atoi(hrs);
	mytime.tm_year = yr - 1900;
	mytime.tm_mday = mday;
	mytime.tm_mon = mon;

	*date = mktime(&mytime);
    }

    if (tbuf != NULL) {
	sprintf(tbuf, "%d-%02d-%02d", yr, mon + 1, mday);
    }

    return 0;
}

/* below: mechanism for tucking individual databases under a 'twisty',
   when there's a large number of databases from a single source.
*/

static char *get_source_string (char *src, const char *s)
{
    char *p;

    *src = 0;
    strncat(src, s, 95);

    p = strchr(src, '(');
    if (p == NULL) {
	p = strstr(src, "--");
    }
    if (p != NULL) {
	*p = '\0';
    }
    tailstrip(src);

    return src;
}

struct src_info {
    int start;
    int ndb;
};

static struct src_info *dbsrc;
static int n_src;

static int push_src_info (int start, int ndb)
{
    struct src_info *tmp;

    tmp = realloc(dbsrc, (n_src + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    dbsrc = tmp;
    dbsrc[n_src].start = start;
    dbsrc[n_src].ndb = ndb;
    n_src++;

    return 0;
}

static int get_ndbs (int lineno)
{
    int i;

    for (i=0; i<n_src; i++) {
	if (dbsrc[i].start == lineno) {
	    return dbsrc[i].ndb;
	}
    }

    return 1;
}

static void free_src_info (void)
{
    free(dbsrc);
    dbsrc = NULL;
    n_src = 0;
}

gint populate_remote_db_list (windata_t *vwin)
{
    GtkTreeStore *store;
    GtkTreeIter iter, child_iter;
    char *getbuf = NULL;
    char line[1024];
    char fname[32];
    char src[96], srcbak[96];
    const char *status;
    gchar *row[3];
    time_t remtime;
    int start, parent, kids;
    int i, ndb, err = 0;

    err = list_remote_dbs(&getbuf);
    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

    i = 0;
    ndb = start = 0;
    src[0] = srcbak[0] = '\0';

    /* first pass: figure "parentage" of databases */

    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	if (strstr(line, "idx")) {
	    continue;
	}
	if (read_remote_filetime(line, fname, &remtime, NULL)) {
	    continue;
	}
	if (bufgets(line, sizeof line, getbuf)) {
	    get_source_string(src, line + 2);
	    if (strcmp(src, srcbak)) {
		if (ndb > 3) {
		    push_src_info(start, ndb);
		}
		start = i;
		ndb = 1;
	    } else {
		ndb++;
	    }
	    strcpy(srcbak, src);
	}
	i++;
    }

    bufgets_finalize(getbuf);

    if (i == 0) {
	errbox(_("No database files found"));
	free_src_info();
	free(getbuf);
	return 1;
    }

    /* second pass: insert databases into tree view */

    store = GTK_TREE_STORE(gtk_tree_view_get_model
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_tree_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    i = 0;
    parent = kids = 0;

    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	if (strstr(line, "idx")) {
	    continue;
	}
	if (read_remote_filetime(line, fname, &remtime, NULL)) {
	    continue;
	}

	status = "";
	get_local_object_status(fname, vwin->role, &status, remtime);
	row[0] = strip_extension(fname);
	row[1] = NULL;
	row[2] = _(status);

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    utf8_correct(line);
	    row[1] = line + 2;
	    ndb = get_ndbs(i);
	    if (ndb > 1) {
		get_source_string(src, row[1]);
		parent = 1;
		kids = ndb;
	    }
	}
	if (parent) {
	    /* header for child databases */
	    gtk_tree_store_append(store, &iter, NULL);
	    gtk_tree_store_set(store, &iter, 0, "", 1, src, -1);
	    parent = 0;
	}
	if (kids > 0) {
	    /* insert child under heading */
	    gtk_tree_store_insert_before(store, &child_iter,
					 &iter, NULL);
	    gtk_tree_store_set(store, &child_iter, 0, row[0],
			       1, row[1], 2, row[2], -1);
	    kids--;
	} else {
	    /* insert at top level */
	    gtk_tree_store_append(store, &iter, NULL);
	    gtk_tree_store_set(store, &iter, 0, row[0], 1, row[1],
			       2, row[2], -1);
	}
	i++;
    }

    bufgets_finalize(getbuf);
    free(getbuf);
    free_src_info();

    if (!err) {
	db_drag_connect(vwin, GRETL_REMOTE_DB_PTR);
    }

    return err;
}

gint populate_dbnomics_provider_list (windata_t *vwin)
{
    gretl_array *A;
    GtkListStore *store;
    GtkTreeIter iter;
    int i, ndb = 0;
    int err = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* call dbnomics package to get listing */
    A = dbnomics_get_providers_call(&err);
    if (err) {
	return err;
    } else {
	gretl_bundle *b;
	char *code, *name;
	int n;

	n = gretl_array_get_length(A);
	for (i=0; i<n; i++) {
	    b = gretl_array_get_bundle(A, i);
	    code = (char *) gretl_bundle_get_string(b, "code", &err);
	    name = (char *) gretl_bundle_get_string(b, "name", &err);
	    if (!err) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter,
				   COL_DBNAME, code,
				   COL_DBINFO, name, -1);
		ndb++;
	    }
	}
	gretl_array_destroy(A);
    }

    if (ndb == 0) {
	errbox(_("No database files found"));
	err = 1;
    }

    return err;
}

struct dbn_pager {
    int offset;  /* offset to pass to dbnomics */
    int n;       /* number of datasets or series in window */
    int ntotal;  /* total, as reported by dbnomics */
    int chunk;   /* limit to number of items to grab */
};

#define DBN_DSETS_CHUNK 40
#define DBN_SERIES_CHUNK 100

static struct dbn_pager *dbn_pager_new (windata_t *vwin)
{
    struct dbn_pager *p = malloc(sizeof *p);

    if (p != NULL) {
	p->offset = 0;
	p->n = 0;
	p->ntotal = 0;
	if (vwin->role == DBNOMICS_DB) {
	    p->chunk = DBN_DSETS_CHUNK;
	} else {
	    p->chunk = DBN_SERIES_CHUNK;
	}
	vwin->data = p;
    }

    return p;
}

void dbnomics_pager_call (GtkWidget *button, windata_t *vwin)
{
    int action = widget_get_int(button, "action");
    struct dbn_pager *pgr = vwin->data;
    int oldoff = pgr->offset;
    int newoff;

    /* action: 1-based enumeration; see make_files_toolbar()
       int datafiles.c
    */

    if (action == 1) {
	/* first */
	pgr->offset = 0;
    } else if (action == 2) {
	/* previous */
	newoff = pgr->offset - pgr->chunk;
	pgr->offset = newoff < 0 ? 0 : newoff;
    } else if (action == 3) {
	/* next */
	int maxoff = pgr->ntotal - 1;

	newoff = pgr->offset + pgr->chunk;
	pgr->offset = newoff > maxoff ? maxoff : newoff;
    } else if (action == 4) {
	/* last */
	pgr->offset = pgr->ntotal - pgr->chunk + 1;
    }

    if (pgr->offset != oldoff) {
	if (vwin->role == DBNOMICS_DB) {
	    populate_dbnomics_dataset_list(vwin, NULL);
	} else {
	    populate_dbnomics_series_list(vwin, NULL);
	}
	listbox_select_first(vwin);
    }
}

static void set_dbn_pager_status (windata_t *vwin)
{
    GtkWidget *b1 = g_object_get_data(G_OBJECT(vwin->mbar), "first-button");
    GtkWidget *b2 = g_object_get_data(G_OBJECT(vwin->mbar), "prev-button");
    GtkWidget *b3 = g_object_get_data(G_OBJECT(vwin->mbar), "next-button");
    GtkWidget *b4 = g_object_get_data(G_OBJECT(vwin->mbar), "last-button");
    struct dbn_pager *pgr = vwin->data;
    int first = pgr->offset + 1;
    int last = pgr->offset + pgr->n;
    gchar *tmp;

    gtk_widget_set_sensitive(b1, pgr->offset > 0);
    gtk_widget_set_sensitive(b2, pgr->offset > 0);
    gtk_widget_set_sensitive(b3, last < pgr->ntotal);
    gtk_widget_set_sensitive(b4, last < pgr->ntotal);

    if (vwin->role == DBNOMICS_DB) {
	tmp = g_strdup_printf(_("showing datasets %d-%d of %d"),
			      first, last, pgr->ntotal);
    } else {
	tmp = g_strdup_printf(_("showing series %d-%d of %d"),
			      first, last, pgr->ntotal);
    }

    gtk_label_set_text(GTK_LABEL(vwin->status), tmp);
    while (gtk_events_pending()) {
	gtk_main_iteration();
    }
    g_free(tmp);
}

/* list the datasets available for a given provider */

gint populate_dbnomics_dataset_list (windata_t *vwin, gpointer p)
{
    gchar *provider = (gchar *) p;
    gretl_array *C, *N;
    gretl_bundle *b;
    char *code, *name;
    GtkListStore *store;
    GtkTreeIter iter;
    struct dbn_pager *pgr;
    int starting = 1;
    int i, imin, imax;
    int err = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (vwin->data == NULL) {
	/* starting: we don't have a pager yet */
	pgr = dbn_pager_new(vwin);
	b = dbnomics_dataset_list(provider, &err);
	if (err) {
	    return err;
	}
    } else {
	b = g_object_get_data(G_OBJECT(vwin->listbox), "dset-list");
	provider = g_object_get_data(G_OBJECT(vwin->listbox), "provider");
	pgr = vwin->data;
	starting = 0;
    }

    if (!err) {
	C = gretl_bundle_get_array(b, "codes", &err);
	N = gretl_bundle_get_array(b, "names", &err);
	if (!err && starting) {
	    pgr->ntotal = gretl_array_get_length(C);
	    if (pgr->ntotal == 0) {
		errbox(_("No datasets were found"));
		err = 1;
	    }
	}
    }

    if (err) {
	return err;
    }

    if (pgr->ntotal <= pgr->chunk) {
	/* just display everything */
	imin = 0;
	imax = pgr->ntotal;
    } else {
	/* display the first so many after offset */
	imin = pgr->offset;
	imax = imin + pgr->chunk;
	imax = imax > pgr->ntotal ? pgr->ntotal : imax;
    }

    pgr->n = imax - imin;

    for (i=imin; i<imax; i++) {
	code = gretl_array_get_element(C, i, NULL, &err);
	name = gretl_array_get_element(N, i, NULL, &err);
	if (!err) {
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter,
			       COL_DBNAME, code,
			       COL_DBINFO, name, -1);
	}
    }

    if (pgr->ntotal <= pgr->chunk) {
	/* no need to keep the dataset-list bundle */
	gretl_bundle_destroy(b);
    } else if (starting) {
	g_object_set_data_full(G_OBJECT(vwin->listbox), "dset-list", b,
			       (GDestroyNotify) gretl_bundle_destroy);
    }

    /* set and show status */
    set_dbn_pager_status(vwin);
    if (starting) {
	/* and make the provider name available downstream */
	g_object_set_data_full(G_OBJECT(vwin->listbox), "provider",
			       provider, g_free);
    }

    return err;
}

gint populate_dbnomics_series_list (windata_t *vwin, gpointer p)
{
    gchar *dsref = (gchar *) p;
    gretl_array *A = NULL;
    gretl_bundle *b;
    char *s, *prov, *dset;
    char *code, *name;
    GtkListStore *store;
    GtkTreeIter iter;
    struct dbn_pager *pgr;
    int starting = 1;
    int alen = 0;
    int i, err = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (vwin->data == NULL) {
	/* starting: we don't have a pager yet */
	pgr = dbn_pager_new(vwin);
    } else {
	dsref = g_object_get_data(G_OBJECT(vwin->listbox), "path");
	pgr = vwin->data;
	starting = 0;
    }

    s = strchr(dsref, '/');
    dset = g_strdup(s + 1);
    prov = g_strndup(dsref, s - dsref);

    /* Note: the length of the retrieved array, which we write into
       @alen, may be less (perhaps a lot less) than the "num_found"
       field that we access below, since the latter records the
       total number of series, regardless of the max number set on
       the dbnomics query.
    */

    A = dbnomics_probe_series(prov, dset, pgr->chunk, pgr->offset, &err);
    if (!err) {
	alen = gretl_array_get_length(A);
	if (alen == 0) {
	    errbox(_("No series were found"));
	    err = 1;
	}
    }

    if (err) {
	return err;
    }

    pgr->n = 0;
    for (i=0; i<alen; i++) {
	b = gretl_array_get_bundle(A, i);
	if (i == 0) {
	    pgr->ntotal = gretl_bundle_get_int(b, "num_found", NULL);
	}
	code = (char *) gretl_bundle_get_string(b, "code", &err);
	name = (char *) gretl_bundle_get_string(b, "name", &err);
	if (!err) {
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter,
			       COL_DBNAME, code,
			       COL_DBINFO, name, -1);
	    pgr->n += 1;
	}
    }

    gretl_array_destroy(A); /* we're done with this */

    if (pgr->n == 0) {
	errbox(_("No series were found"));
	err = 1;
    } else {
	/* set and show status */
	set_dbn_pager_status(vwin);
	if (starting) {
	    /* and make the dataset 'path' available downstream */
	    g_object_set_data_full(G_OBJECT(vwin->listbox), "path",
				   dsref, g_free);
	}
    }

    return err;
}

static int get_addon_info (xmlNodePtr node, xmlDocPtr doc, char **S)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    int err = 0;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "version")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &S[1]);
	} else if (!xmlStrcmp(cur->name, (XUC) "date")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &S[2]);
	} else if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &S[4]);
	}

	cur = cur->next;
    }

    if (S[1] == NULL || S[2] == NULL || S[4] == NULL) {
	err = E_DATA;
    } else {
	char *path = gretl_function_package_get_path(S[0], PKG_SUBDIR);

	if (path != NULL) {
	    S[3] = installed_addon_status_string(path, S[1]);
	    free(path);
	} else {
	    S[3] = gretl_strdup(_("Not installed"));
	}
    }

    return err;
}

gint populate_remote_addons_list (windata_t *vwin)
{
    xmlDocPtr doc;
    xmlNodePtr node = NULL;
    char *getbuf = NULL;
    int n = 0, err = 0;

    err = query_sourceforge("/addons-data/addons.xml", &getbuf);
    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

#if 0
    fprintf(stderr, "populate_remote_addons_list: got XML:\n");
    fprintf(stderr, "'%s'\n", getbuf);
#endif

    xmlKeepBlanksDefault(0);

    doc = xmlParseMemory(getbuf, strlen(getbuf));
    if (doc == NULL) {
	err = E_DATA;
    } else {
	node = xmlDocGetRootElement(doc);
	if (node == NULL || xmlStrcmp(node->name, (XUC) "gretl-addons")) {
	    err = E_DATA;
	}
    }

    if (!err) {
	GtkListStore *store;
	GtkTreeIter iter;
	char *S[5] = { NULL };
	int i;

	store = GTK_LIST_STORE(gtk_tree_view_get_model
			       (GTK_TREE_VIEW(vwin->listbox)));
	gtk_list_store_clear(store);
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	node = node->xmlChildrenNode;

	while (node != NULL && !err) {
	    if (!xmlStrcmp(node->name, (XUC) "gretl-addon")) {
		gretl_xml_get_prop_as_string(node, "name", &S[0]);
		if (S[0] == NULL) {
		    err = E_DATA;
		} else {
		    err = get_addon_info(node, doc, S);
		    if (!err) {
			gtk_list_store_append(store, &iter);
			gtk_list_store_set(store, &iter,
					   0, S[0], 1, S[1],
					   2, S[2], 3, S[3],
					   4, S[4], -1);
			for (i=0; i<5; i++) {
			    free(S[i]);
			    S[i] = NULL;
			}
			n++;
		    }
		}
	    }
	    node = node->next;
	}
    }

    if (err) {
	gui_errmsg(err);
    } else if (n == 0) {
	warnbox(_("No function packages found"));
	err = 1;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    free(getbuf);

    return err;
}

static void check_gfn_drag_connection (windata_t *vwin)
{
    int dc = widget_get_int(vwin->main, "drag-connected");

    if (!dc) {
	db_drag_connect(vwin, GRETL_REMOTE_FNPKG_PTR);
	widget_set_int(vwin->main, "drag-connected", 1);
    }
}

/* Fill a list box with name, version number, author,
   and short description of function packages, retrieved
   from server.
*/

gint populate_remote_func_list (windata_t *vwin, int filter)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char *getbuf = NULL;
    char line[1024];
    char fname[128];
    char basename[32];
    const char *status;
    time_t remtime;
    int n, err = 0;

    err = list_remote_function_packages(&getbuf, filter);
    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

#if 0
    fprintf(stderr, "getbuf: '%s'\n", getbuf);
#endif

    store = GTK_LIST_STORE(gtk_tree_view_get_model
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    n = 0;
    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	char *descrip = NULL;
	char *version = NULL;
	char *author = NULL;
	char datestr[12];
	gboolean zipfile;

	if (read_remote_filetime(line, fname, &remtime, datestr)) {
	    continue;
	}

	strcpy(basename, fname);
	strip_extension(basename);

	zipfile = has_suffix(fname, ".zip");
	if (zipfile) {
	    /* local status: look for PKG/PKG.gfn, not PKG.zip */
	    sprintf(fname, "%s%c%s.gfn", basename, SLASH, basename);
	}

	status = "";
	get_local_object_status(fname, vwin->role, &status, remtime);

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    utf8_correct(line);
	    descrip = gretl_strdup(line + 2);
	    maybe_ellipsize_string(descrip, 48);
	}

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    version = gretl_strdup(line + 2);
	}

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    author = gretl_strdup(line + 2);
	    maybe_ellipsize_string(author, 26);
	}

	if (descrip != NULL && version != NULL && author != NULL) {
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter,
			       0, basename,
			       1, version,
			       2, datestr,
			       3, author,
			       4, descrip,
			       5, _(status),
			       6, zipfile,
			       -1);
	    n++;
	}

	free(descrip);
	free(version);
	free(author);
    }

    bufgets_finalize(getbuf);
    free(getbuf);

    if (n == 0) {
	warnbox(_("No function packages found"));
	err = 1;
    }

    if (!err) {
	check_gfn_drag_connection(vwin);
    }

    return err;
}

/* Fill a list box with names and short descriptions
   of data file packages, retrieved from server.
*/

gint populate_remote_data_pkg_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char *getbuf = NULL;
    char line[256];
    char fname[32];
    char tstr[16];
    char *descrip;
    int n, err = 0;

    err = list_remote_data_packages(&getbuf);

    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    n = 0;
    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	descrip = NULL;
	*tstr = '\0';

	if (read_remote_filetime(line, fname, NULL, tstr)) {
	    continue;
	}

	strip_extension(fname);

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    utf8_correct(line);
	    descrip = gretl_strdup(line + 2);
	}

	if (descrip == NULL) {
	    descrip = gretl_strdup("");
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter,
			   0, fname,
			   1, descrip,
			   2, tstr,
			   -1);

	free(descrip);
	n++;
    }

    bufgets_finalize(getbuf);
    free(getbuf);

    if (n == 0) {
	warnbox(_("No data packages found"));
	err = 1;
    }

    return err;
}

static int maybe_replace_db_path (const char *name,
				  GtkTreeIter *iter,
				  GtkTreeIter *iprev,
				  GtkTreeModel *mod)
{
    gchar *db, *dbprev;
    gchar *path, *pathprev;
    struct stat buf, bufprev;
    int err = 0, ret = 0;

    fprintf(stderr, "databases: found a duplicate of %s\n", name);

    gtk_tree_model_get(mod, iter, 2, &path, -1);
    gtk_tree_model_get(mod, iprev, 2, &pathprev, -1);
    /* formulate full names of the two .bin files */
    db = g_strdup_printf("%s%c%s.bin", path, SLASH, name);
    dbprev = g_strdup_printf("%s%c%s.bin", pathprev, SLASH, name);

    err = gretl_stat(db, &buf);
    if (!err) {
	err = gretl_stat(dbprev, &bufprev);
    }
    if (!err && buf.st_mtime > bufprev.st_mtime) {
	/* @db is newer than @dbprev, so replace path */
	fprintf(stderr, " using newer version in %s\n", path);
	fprintf(stderr, " masking version in %s\n", pathprev);
	gtk_list_store_set(GTK_LIST_STORE(mod), iprev, 2, path, -1);
	ret = 1;
    } else {
	fprintf(stderr, " keeping version in %s\n", pathprev);
	fprintf(stderr, " ignoring version in %s\n", path);
    }

    g_free(path);
    g_free(pathprev);
    g_free(db);
    g_free(dbprev);

    return ret;
}

/* Purge any duplicates from the list of database files to
   display -- in case of duplicates we keep the newer file
   as assessed by stat's st_mtime.
*/

static void maybe_prune_db_list (GtkTreeView *tview,
				 int *pndb)
{
    char **S;
    GtkTreeModel *mod;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkTreeIter *icpy;
    int ndb, i = 0;

    mod = gtk_tree_view_get_model(tview);
    store = GTK_LIST_STORE(mod);
    if (!gtk_tree_model_get_iter_first(mod, &iter)) {
	return;
    }

    ndb = *pndb;
    S = strings_array_new(ndb);
    icpy = malloc(ndb * sizeof *icpy);
    if (S == NULL || icpy == NULL) {
	return;
    }

    while (1) {
	int j, drop = 0;

	gtk_tree_model_get(mod, &iter, 0, &S[i], -1);
	icpy[i] = iter;
	for (j=0; j<i; j++) {
	    if (S[j] != NULL && !strcmp(S[i], S[j])) {
		/* found a duplicate */
		drop = 1;
		*pndb -= 1;
		maybe_replace_db_path(S[j], &iter, &icpy[j], mod);
		gtk_list_store_remove(store, &iter);
		iter = icpy[i-1]; /* back up one row */
		g_free(S[i]);
		S[i] = NULL;
		break;
	    }
	}
	if (!gtk_tree_model_iter_next(mod, &iter)) {
	    break;
	}
	if (!drop) {
	    i++;
	}
    }

    for (i=0; i<ndb; i++) {
	g_free(S[i]);
    }
    free(S);
    free(icpy);
}

static int prep_dbnomics_series (gretl_bundle *b,
				 DATASET *dbset)
{
    gretl_array *A;
    gretl_matrix *v;
    const char *id;
    int T, err = 0;

    T = gretl_bundle_get_int(b, "T", &err);
    A = gretl_bundle_get_array(b, "period", &err);
    v = gretl_bundle_get_matrix(b, "value", &err);
    id = gretl_bundle_get_string(b, "series_code", &err);

    if (!err && (T <= 0 || A == NULL || v == NULL)) {
	err = E_DATA;
    }

    if (!err) {
	char **S = gretl_array_get_strings(A, &T);
	gchar *fname;
	FILE *fp;
	int t;

	fname = g_strdup_printf("%sdnomics_tmp.txt", gretl_dotdir());
	fp = gretl_fopen(fname, "w");
	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    gretl_push_c_numeric_locale();
	    fputs("obs dbnomics_data\n", fp);
	    for (t=0; t<T; t++) {
		if (na(v->val[t])) {
		    fprintf(fp, "%s NA\n", S[t]);
		} else {
		    fprintf(fp, "%s %.12g\n", S[t], v->val[t]);
		}
	    }
	    gretl_pop_c_numeric_locale();
	    fclose(fp);
	    err = import_csv(fname, dbset, OPT_NONE, NULL);
	    if (!err && id != NULL) {
		series_set_display_name(dbset, 1, id);
	    }
	    gretl_remove(fname);
	}
	g_free(fname);
    }

    return err;
}

int add_dbnomics_data (windata_t *vwin)
{
    gretl_bundle *b = vwin->data;
    DATASET *dbset = NULL;
    int err;

    dbset = datainfo_new();
    if (dbset == NULL) {
	nomem();
	return E_ALLOC;
    } else {
	int freeit = 1;

	err = prep_dbnomics_series(b, dbset);
	if (!err) {
	    char vname[VNAMELEN];
	    const char *S[4];
	    gchar *full_id = NULL;
	    char *descrip = NULL;
	    int cancel = 0;

	    /* construct a default name for the series */
	    *vname = '\0';
	    S[0] = gretl_bundle_get_string(b, "series_code", &err);
	    if (!err) {
		normalize_join_colname(vname, S[0], 0, 0);
	    }
	    /* construct its description */
	    S[1] = gretl_bundle_get_string(b, "provider_code", &err);
	    S[2] = gretl_bundle_get_string(b, "dataset_code", &err);
	    S[3] = gretl_bundle_get_string(b, "series_name", &err);
	    if (!err) {
		full_id = g_strdup_printf("%s/%s/%s", S[1], S[2], S[0]);
		descrip = g_strdup_printf("%s: %s", full_id, S[3]);
	    }
	    name_new_series_dialog(vname, descrip, vwin, &cancel);
	    if (!cancel) {
		set_dbnomics_id(full_id);
		strcpy(dbset->varname[1], vname);
		series_set_label(dbset, 1, descrip);
		series_set_display_name(dbset, 1, "");
		add_dbdata(vwin, dbset, NULL, &freeit);
		unset_dbnomics_id();
	    }
	    g_free(descrip);
	    g_free(full_id);
	}
	if (freeit) {
	    destroy_dataset(dbset);
	}
    }

    return err;
}

int show_dbnomics_data (windata_t *vwin, int plot)
{
    gretl_bundle *b = vwin->data;
    DATASET dbset = {0};
    int err;

    err = prep_dbnomics_series(b, &dbset);
    if (!err) {
	if (plot) {
	    graph_dbdata(&dbset);
	} else {
	    display_dbdata(&dbset);
	}
    }

    free_Z(&dbset);
    clear_datainfo(&dbset, CLEAR_FULL);

    return err;
}

gint populate_dbfilelist (windata_t *vwin, int *pndb)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char **dirnames;
    GDir *dir;
    int i, n_dirs;
    int nf, ndb = 0;
    int err = 0;

#if DB_SEARCH_DEBUG
    fprintf(stderr, "populate_dbfilelist...\n");
#endif

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    dirnames = get_plausible_search_dirs(DB_SEARCH, &n_dirs);

    for (i=0; i<n_dirs; i++) {
	dir = gretl_opendir(dirnames[i]);
	if (dir != NULL) {
	    nf = read_db_files_in_dir(dir, vwin->role, dirnames[i], store, &iter);
#if DB_SEARCH_DEBUG
	    fprintf(stderr, " found %d db files in '%s'\n", nf, dirnames[i]);
#endif
	    ndb += nf;
	    g_dir_close(dir);
	}
    }

    strings_array_free(dirnames, n_dirs);

    if (ndb == 0) {
	errbox(_("No database files found"));
	err = 1;
    } else {
	maybe_prune_db_list(GTK_TREE_VIEW(vwin->listbox), &ndb);
	presort_treelist(vwin);
    }

    if (pndb != NULL) {
	*pndb = ndb;
    }

    return err;
}

void set_db_dir_callback (windata_t *vwin, char *path)
{
    GDir *dir = gretl_opendir(path);
    int ndb = 0;

    if (dir != NULL) {
	GtkListStore *store;
	GtkTreeIter iter;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	ndb = read_db_files_in_dir(dir, vwin->role, path, store, &iter);
	g_dir_close(dir);
    }

    if (ndb == 0) {
	warnbox(_("No database files found"));
    }
}

static void set_compact_info_from_default (int method)
{
    int i;

    for (i=1; i<dataset->v; i++) {
	if (series_get_compact_method(dataset, i) == COMPACT_NONE) {
	    series_set_compact_method(dataset, i, method);
	}
    }
}

void do_compact_data_set (void)
{
    CompactMethod method = COMPACT_AVG;
    int err, newpd = 0, monstart = 1;
    int repday = 0;
    int *pmonstart = NULL;

    if (maybe_restore_full_data(COMPACT)) {
	return;
    }

    if (dated_seven_day_data(dataset)) {
	pmonstart = &monstart;
    }

    data_compact_dialog(dataset->pd, &newpd, pmonstart,
			&method, &repday, mdata->main);

    if (method == COMPACT_NONE) {
	/* the user cancelled */
	fprintf(stderr, "canceled!\n");
	return;
    }

    err = compact_data_set(dataset, newpd, method, monstart, repday);

    if (err) {
	gui_errmsg(err);
    } else {
	const char *mstr = compact_method_string(method);

	if (mstr != NULL) {
	    lib_command_sprintf("dataset compact %d %s", newpd, mstr);
	} else {
	    lib_command_sprintf("dataset compact %d", newpd);
	}
	record_command_verbatim();

	mark_dataset_as_modified();
	if (method == COMPACT_SPREAD) {
	    populate_varlist();
	}
	set_compact_info_from_default(method);
    }
}

void do_expand_data_set (void)
{
    int err, newpd, interpol = 1;

    if (maybe_restore_full_data(EXPAND)) {
	return;
    }

    /* supported: annual to quarterly or quarterly to monthly */
    newpd = (dataset->pd == 1)? 4 : 12;

    data_expand_dialog(dataset->pd, &interpol, mdata->main);

    if (interpol < 0) {
	/* canceled */
	return;
    }

    err = expand_data_set(dataset, newpd, interpol);

    if (err) {
	gui_errmsg(err);
    } else {
	mark_dataset_as_modified();
    }
}

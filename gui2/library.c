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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* library.c for gretl -- main interface to libgretl functions */

#include "gretl.h"
#include "var.h"
#include "textbuf.h"
#include "gpt_control.h"
#include "graph_page.h"
#include "console.h"
#include "system.h"
#include "gretl_restrict.h"
#include "gretl_func.h"
#include "modelspec.h"
#include "forecast.h"
#include "dbwrite.h"
#include "menustate.h"
#include "dlgutils.h"
#include "ssheet.h"
#include "lib_private.h"
#include "cmd_private.h"
#include "libset.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_panel.h"

#ifdef G_OS_WIN32 
# include <io.h>
# include "gretlwin32.h"
#else
# include <unistd.h>
# include <sys/stat.h>
#endif

#include "session.h"
#include "selector.h"
#include "boxplots.h"
#include "series_view.h"
#include "objectsave.h"
#include "datafiles.h"
#include "model_table.h"
#include "cmdstack.h"
#include "filelists.h"

#define CMD_DEBUG 0

#ifdef HAVE_TRAMO
extern char tramo[];
extern char tramodir[];
#endif

/* private functions */
static void update_model_tests (windata_t *vwin);
static int finish_genr (MODEL *pmod, dialog_t *dlg);
#ifndef G_OS_WIN32
static int get_terminal (char *s);
#endif

const char *CANTDO = N_("Can't do this: no model has been estimated yet\n");

/* file scope state variables */
static CMD cmd;
static char cmdline[MAXLINE];
static int replay;
static MODELSPEC *modelspec;
static gretl_equation_system *sys;
static gretl_restriction_set *rset;
static int original_n;

char *get_lib_cmdline (void)
{
    return cmdline;
}

CMD *get_lib_cmd (void)
{
    return &cmd;
}

void lib_cmd_destroy_context (void)
{
    gretl_cmd_destroy_context(&cmd);
}

void lib_modelspec_free (void)
{
    if (modelspec != NULL) {
	free_modelspec(modelspec);
	modelspec = NULL;
    }
}

void set_original_n (int n)
{
    original_n = n;
}

int get_original_n (void)
{
    return original_n;
} 

void library_command_init (void)
{
    gretl_cmd_init(&cmd);
}

void library_command_free (void)
{
    gretl_cmd_free(&cmd);
}

int replaying (void)
{
    return replay;
}

void set_replay_on (void)
{
    replay = 1;
}

void set_replay_off (void)
{
    replay = 0;
}

void register_graph (void)
{
    const char *msg;

    gnuplot_show_png(gretl_plotfile(), NULL, 0);

    msg = get_gretl_errmsg();
    if (msg != NULL && *msg != '\0') {
	errbox(msg);
    }
}

static void gui_graph_handler (int err)
{
    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	const char *msg = get_gretl_errmsg();

	if (*msg) {
	    errbox(msg);
	} else {
	    errbox(_("gnuplot command failed"));
	}
    } else {
	register_graph();
    }
}

static void launch_gnuplot_interactive (void)
{
# ifdef G_OS_WIN32
    gchar *gpline;

    gpline = g_strdup_printf("\"%s\" \"%s\" -", paths.gnuplot,
			     gretl_plotfile());
    create_child_process(gpline, NULL);
    g_free(gpline);
# else
    char term[8];
    char plotfile[MAXLEN];

    strcpy(plotfile, gretl_plotfile());

    if (get_terminal(term)) {
	return;
    } else {
	GError *error = NULL;
	gchar *argv[] = { 
	    term, "+sb", "+ls", 
	    "-geometry", "40x4",
	    "-title", "gnuplot: type q to quit",
	    "-e", paths.gnuplot, plotfile, "-",
	    NULL 
	};
	int ok;

	ok = g_spawn_async(NULL, /* working dir */
			   argv,
			   NULL, /* env */
			   G_SPAWN_SEARCH_PATH,
			   NULL, /* child_setup */
			   NULL, /* user_data */
			   NULL, /* child_pid ptr */
			   &error);
	if (!ok) {
	    errbox(error->message);
	    g_error_free(error);
	}
    }
# endif
}

void set_sample_label_special (void)
{
    char labeltxt[80];

    sprintf(labeltxt, _("Undated: Full range n = %d; current sample"
	    " n = %d"), get_full_length_n(), datainfo->n);
    gtk_label_set_text(GTK_LABEL(mdata->status), labeltxt);

    time_series_menu_state(FALSE);
    panel_menu_state(FALSE);
    ts_or_panel_menu_state(FALSE);
}

int gretl_command_sprintf (const char *template, ...)
{
    va_list args;
    int len;

    memset(cmdline, 0, MAXLINE);

    va_start(args, template);
    len = vsprintf(cmdline, template, args);
    va_end(args);

#if 0
    fprintf(stderr, "gretl_command_sprintf: cmdline = '%s'\n", cmdline);
#endif

    return len;
}

int gretl_command_strcpy (const char *s)
{
    memset(cmdline, 0, MAXLINE);
    strcpy(cmdline, s);
    
    return 0;
}

static int gretl_command_strcat (const char *s)
{
    strcat(cmdline, s);

    return 0;
}

int user_fopen (const char *fname, char *fullname, PRN **pprn)
{
    int err = 0;

    strcpy(fullname, paths.userdir);
    strcat(fullname, fname);

    *pprn = gretl_print_new_with_filename(fullname);

    if (*pprn == NULL) {
	errbox(_("Couldn't open file for writing"));
	err = 1;
    }

    return err;
}

gint bufopen (PRN **pprn)
{
    *pprn = gretl_print_new(GRETL_PRINT_BUFFER);

    if (*pprn == NULL) {
	nomem();
	return 1;
    }

    return 0;
}

static int freq_error (FreqDist *freq, PRN *prn)
{
    int err = 0;

    if (freq == NULL) {
	if (prn == NULL) {
	    nomem();
	} else {
	    pputs(prn, _("Out of memory!"));
	    pputc(prn, '\n');
	}
	err = 1;
    } else if (get_gretl_errno()) {
	if (prn == NULL) {
	    gui_errmsg(get_gretl_errno());
	} else {
	    errmsg(get_gretl_errno(), prn);
	}
	free_freq(freq);
	err = 1;
    }

    return err;
}

static void maybe_quote_filename (char *s, char *cmd)
{
    size_t len = strlen(cmd);

    if (strlen(s) > len + 1) {
	char *p = s + len + 1;

	if (*p == '"' || *p == '\'') {
	    return;
	}
	
	if (strchr(p, ' ')) {
	    char tmp[MAXLEN];

	    *tmp = 0;
	    strcpy(tmp, p);
	    sprintf(s, "%s \"%s\"", cmd, tmp);
	}
    }
}

static int cmd_init (char *s)
{
    PRN *echo;
    int err = 0;

    /* note "cmd.*" elements are filled out already, if
       check_specific_command() has been called on the 
       command string */

#if CMD_DEBUG
    fprintf(stderr, "cmd_init: got cmdstr: '%s'\n", s);
    fprintf(stderr, "cmd.word: '%s'\n", cmd.word);
    fprintf(stderr, "cmd.param: '%s'\n", cmd.param);
    fprintf(stderr, "cmd.opt: %d\n", (int) cmd.opt);
#endif

    if (cmd.ci == OPEN || cmd.ci == RUN) {
	maybe_quote_filename(s, cmd.word);
    }

    /* arrange to have the command recorded on a stack */
    if (bufopen(&echo)) {
	err = 1;
    } else {
	const char *buf;

	echo_cmd(&cmd, datainfo, s, 0, echo);
	buf = gretl_print_get_buffer(echo);
#if CMD_DEBUG
	fprintf(stderr, "from echo_cmd: buf='%s'\n", buf);
#endif
	err = add_command_to_stack(buf);
	gretl_print_destroy(echo);
    }

    return err;
}

static int lib_cmd_init (void)
{
    return cmd_init(cmdline);
}

/* checks command line s for validity, but does not
   record the command */

int check_specific_command (char *s)
{
    int err;

#if CMD_DEBUG
    fprintf(stderr, "check_specific_command: s = '%s'\n", s);
#endif

    /* "cmd" is global */
    err = parse_command_line(s, &cmd, &Z, datainfo); 
    if (err) {
	gui_errmsg(err);
    } else {
	/* At this point we're not just replaying 
	   saved session commands. */
	set_replay_off();
    }

    return err;
}

static int check_lib_command (void)
{
    return check_specific_command(cmdline);
}

int check_and_record_command (void)
{
    return (check_specific_command(cmdline) || cmd_init(cmdline));
}

/* checks command for errors, and if OK returns an allocated
   copy of the command list */

int *command_list_from_string (char *s)
{
    CMD mycmd;
    int *list = NULL;
    int err;

    err = gretl_cmd_init(&mycmd);

    if (!err) {
	err = parse_command_line(s, &mycmd, &Z, datainfo);
	if (!err) {
	    list = gretl_list_copy(mycmd.list);
	}
    }

    if (err) {
	gui_errmsg(err);
    } 

    gretl_cmd_free(&mycmd);

    return list;    
}

/* Model estimated via console or script: unlike a gui model, which is
   kept in memory so long as its window is open, these models are
   immediately discarded.  So if we want to be able to refer back to
   them later we need to record their specification */

static gint stack_script_modelspec (MODEL *pmod)
{
    int err;

    attach_subsample_to_model(models[0], datainfo);
    err = modelspec_save(models[0], &modelspec);

    return err;
}

void add_mahalanobis_data (windata_t *vwin)
{
    MahalDist *md = (MahalDist *) vwin->data;
    const double *dx;
    const int *mlist;
    char *liststr;
    char vname[VNAMELEN];
    int v, t;

    if (md == NULL) {
	errbox(_("Error adding variables"));
	return;
    }

    dx = mahal_dist_get_distances(md);
    mlist = mahal_dist_get_varlist(md);
    if (dx == NULL || mlist == NULL) {
	errbox(_("Error adding variables"));
	return;
    }

    if (dataset_add_series(1, &Z, datainfo)) {
	nomem();
	return;
    }

    v = datainfo->v - 1;

    strcpy(vname, "mdist");
    make_varname_unique(vname, 0, datainfo);
    strcpy(datainfo->varname[v], vname);
    sprintf(VARLABEL(datainfo, v), _("Mahalanobis distances"));

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return;
    }

    for (t=0; t<datainfo->n; t++) {
	Z[v][t] = dx[t];
    }

    liststr = gretl_list_to_string(mlist);
    gretl_command_sprintf("mahal %s --save", liststr);
    free(liststr);
    check_and_record_command();
}

void add_pca_data (windata_t *vwin)
{
    int err, oldv = datainfo->v;
    gretlopt oflag = OPT_D;
    VMatrix *corrmat = (VMatrix *) vwin->data;

    err = call_pca_plugin(corrmat, &Z, datainfo, &oflag, NULL);

    if (err) {
	gui_errmsg(err);
	return;
    }

    if (datainfo->v > oldv) {
	/* if data were added, register the command */
	if (oflag == OPT_O || oflag == OPT_A) {
	    char *liststr = gretl_list_to_string(corrmat->list);
	    
	    if (liststr != NULL) {
		const char *flagstr = print_flags(oflag, PCA); 

		gretl_command_sprintf("pca %s%s", liststr, flagstr);
		check_and_record_command();
		free(liststr);
	    }
	}
    }
}

static void make_fcast_save_name (char *vname, const char *s)
{
    strcpy(vname, s); 
    gretl_trunc(vname, 5);
    if (strlen(vname) < 5) {
	strcat(vname, "_hat");
    } else {
	strcat(vname, "hat");
    }
    make_varname_unique(vname, 0, datainfo);
}

void add_fcast_data (windata_t *vwin)
{
    char stobs[OBSLEN], endobs[OBSLEN];
    FITRESID *fr = (FITRESID *) vwin->data;
    char vname[VNAMELEN];
    int v, t;

    if (dataset_add_series(1, &Z, datainfo)) {
	nomem();
	return;
    }

    v = datainfo->v - 1;

    make_fcast_save_name(vname, fr->depvar);
    strcpy(datainfo->varname[v], vname);
    sprintf(VARLABEL(datainfo, v), _("forecast of %s"), fr->depvar);

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return;
    }

    for (t=0; t<datainfo->n; t++) {
	Z[v][t] = fr->fitted[t];
    }

    ntodate(stobs, fr->t1, datainfo);
    ntodate(endobs, fr->t2, datainfo);

    gretl_command_sprintf("fcast %s %s %s", stobs, endobs, datainfo->varname[v]);
    model_command_init(fr->model_ID);

    /* nothing else need be done, since we're called by
       add_data_callback() */
}

static const char *selected_varname (void)
{
    return datainfo->varname[mdata_active_var()];
}

void do_menu_op (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    char title[48];
    char *liststr = NULL;
    int err = 0;
    gpointer obj = NULL;
    gretlopt opt = OPT_NONE;
    gint hsize = 78, vsize = 380;

    strcpy(title, "gretl: ");

    if (action == CORR || action == SUMMARY ||
	action == PCA || action == MAHAL) {
	liststr = main_window_selection_as_string();
	if (liststr == NULL) return;
    }

    switch (action) {
    case CORR:
	gretl_command_sprintf("corr%s", liststr);
	strcat(title, _("correlation matrix"));
	action = CORR;
	break;
    case PCA:
	gretl_command_sprintf("pca%s", liststr);
	strcat(title, _("principal components"));
	break;
    case MAHAL:
	gretl_command_sprintf("mahal%s", liststr);
	hsize = 60;
	strcat(title, _("Mahalanobis distances"));
	break;
    case FREQ:
	gretl_command_sprintf("freq %s", selected_varname());
	strcat(title, _("frequency distribution"));
	vsize = 340;
	break;
    case RUNS:
	gretl_command_sprintf("runs %s", selected_varname());
	strcat(title, _("runs test"));
	vsize = 200;
	break;
    case SUMMARY:
	gretl_command_sprintf("summary%s", liststr);
	strcat(title, _("summary statistics"));
	action = SUMMARY;
	break;
    case VAR_SUMMARY:
	gretl_command_sprintf("summary %s", selected_varname());
	strcat(title, _("summary stats: "));
	strcat(title, selected_varname());
	vsize = 300;
	break;
    default:
	break;
    }

    if (liststr != NULL) {
	free(liststr);
    }

    /* check the command and initialize output buffer */
    if (check_and_record_command() || bufopen(&prn)) {
	return;
    }

    /* execute the command */
    switch (action) {

    case CORR:
	obj = corrlist(cmd.list, (const double **) Z, datainfo);
	if (obj == NULL) {
	    errbox(_("Failed to generate correlation matrix"));
	    gretl_print_destroy(prn);
	    return;
	} 
	matrix_print_corr(obj, datainfo, prn);
	break;

    case FREQ:
	err = freqdist(cmd.list[1], (const double **) Z, datainfo,
		       0, OPT_NONE, prn);
	break;

    case RUNS:
	err = runs_test(cmd.list[1], (const double **) Z, datainfo, prn);
	break;

    case PCA:
	obj = corrlist(cmd.list, (const double **) Z, datainfo);
	if (obj == NULL) {
	    errbox(_("Failed to generate correlation matrix"));
	    gretl_print_destroy(prn);
	    return;
	} else {
	    err = call_pca_plugin((VMatrix *) obj, &Z, datainfo, 
				  NULL, prn);
	}
	break;

    case MAHAL:
	if (cmd.list[0] <= 4) {
	    opt = OPT_V;
	}
	obj = get_mahal_distances(cmd.list, &Z, datainfo, opt, prn);
	if (obj == NULL) {
	    errbox(_("Command failed"));
	    gretl_print_destroy(prn);
	    return;
	}
	break;

    case SUMMARY:
    case VAR_SUMMARY:	
	obj = summary(cmd.list, (const double **) Z, datainfo, prn);
	if (obj == NULL) {
	    errbox(_("Failed to generate summary statistics"));
	    gretl_print_destroy(prn);
	    return;
	}	    
	print_summary(obj, datainfo, prn);
	break;
    }

    if (err) {
	gui_errmsg(err);
    } 


    view_buffer(prn, hsize, vsize, title, action, obj);
}

int do_coint (selector *sr)
{
    const char *buf = selector_list(sr);
    int action = selector_code(sr);
    GRETL_VAR *jvar = NULL;
    PRN *prn;
    int err = 0, order = 0;

    if (buf == NULL) return 1;

    cmd.opt = selector_get_opts(sr);

    if (action == COINT) {
	gretl_command_sprintf("coint %s%s", buf, print_flags(cmd.opt, action));
    } else {
	gretl_command_sprintf("coint2 %s%s", buf, print_flags(cmd.opt, action));
    }	

    if (check_and_record_command() || bufopen(&prn)) {
	return 1;
    }

    order = atoi(cmd.param);
    if (!order) {
	errbox(_("Couldn't read cointegration order"));
	gretl_print_destroy(prn);
	return 1;
    }

    if (action == COINT) {
	err = coint(order, cmd.list, &Z, datainfo, cmd.opt, prn);
    } else {
	jvar = johansen_test(order, cmd.list, &Z, datainfo, cmd.opt, prn);
	if ((err = jvar->err)) {
	    gretl_VAR_free(jvar);
	}
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return err;
    } 

    view_buffer(prn, 78, 400, _("gretl: cointegration test"), 
		action, (action == COINT2)? jvar : NULL);

    return 0;
}

static int ok_obs_in_series (int varno)
{
    int t, t1, t2;

    for (t=datainfo->t1; t<datainfo->t2; t++) {
	if (!na(Z[varno][t])) break;
    }

    t1 = t;

    for (t=datainfo->t2; t>=datainfo->t1; t--) {
	if (!na(Z[varno][t])) break;
    }

    t2 = t;

    return t2 - t1 + 1;
}

void unit_root_test (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    const char *adf_opts[] = {
	N_("test without constant"),
	N_("with constant"),
	N_("with constant and trend"),
	N_("with constant, trend and trend squared"),
	N_("include seasonal dummies"),
	N_("show regression results"),
	N_("test down from maximum lag order"),
	N_("use level of variable"),
	N_("use first difference of variable")
    };
    const char *kpss_opts[] = {
	N_("include a trend"),
	N_("show regression results"),
	N_("use level of variable"),
	N_("use first difference of variable")
    };

    const char *adf_title = N_("gretl: ADF test");
    const char *kpss_title = N_("gretl: KPSS test");
    const char *adf_spintext = N_("Lag order for ADF test:");
    const char *kpss_spintext = N_("Lag order for KPSS test:");
    const char *title, *spintext, **opts;

    /* save the user's settings, per session */
    static int adf_active[] = { 0, 1, 1, 1, 0, 0, 0 };
    static int kpss_active[] = { 1, 0 };
    static int order = 1;

    int difference = 0;
    int v = mdata_active_var();

    int okT, omax, err;

    if (order < 0) {
	order = -order;
    }

    if (action == ADF) {
	title = adf_title;
	spintext = adf_spintext;
	opts = adf_opts;

	okT = ok_obs_in_series(v);
	omax = okT / 2;
    } else {
	title = kpss_title;
	spintext = kpss_spintext;
	opts = kpss_opts;

	okT = ok_obs_in_series(v);	
	omax = okT / 2;
	order = 4.0 * pow(okT / 100.0, 0.25);
    }

    if (order > omax) {
	order = omax;
    }  

    if (action == ADF && datainfo->pd == 1) {
	adf_active[4] = -1;
    }

    err = checks_dialog(_(title), opts,
			(action == ADF)? 7 : 2,
			(action == ADF)? adf_active : kpss_active,
			2, &difference,
			&order, _(spintext),
			0, omax, action);
    if (err < 0) {
	return;
    }

    if (action == ADF) {
	if (adf_active[0] == 0 &&
	    adf_active[1] == 0 &&
	    adf_active[2] == 0 &&
	    adf_active[3] == 0) {
	    return;
	}
    }

    gretl_command_sprintf("%s %d %s", (action == ADF)? "adf" : "kpss", order, 
			  selected_varname());

    if (action == ADF) {
	if (adf_active[0]) gretl_command_strcat(" --nc");
	if (adf_active[1]) gretl_command_strcat(" --c");
	if (adf_active[2]) gretl_command_strcat(" --ct");
	if (adf_active[3]) gretl_command_strcat(" --ctt");
	if (adf_active[4] > 0) gretl_command_strcat(" --seasonals");
	if (adf_active[5]) gretl_command_strcat(" --verbose");
	if (adf_active[6]) order = -order; /* auto-trim the lag order */
    } else {
	if (kpss_active[0]) gretl_command_strcat(" --trend");
	if (kpss_active[1]) gretl_command_strcat(" --verbose");
    } 

    if (difference) {
	gretl_command_strcat(" --difference");
    }

    if (check_and_record_command() || bufopen(&prn)) {
	return;
    }

    if (action == ADF) {
	err = adf_test(order, cmd.list[1], &Z, datainfo, cmd.opt, prn);
    } else {
	err = kpss_test(order, cmd.list[1], &Z, datainfo, cmd.opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 350, title, action, NULL);
    }    
}

int do_spearman (selector *sr)
{
    const char *buf = selector_list(sr);
    PRN *prn;
    char title[64];
    gint err;

    if (buf == NULL) return 1;
    
    gretl_command_sprintf("spearman%s --verbose", buf);

    if (check_and_record_command() || bufopen(&prn)) {
	return 1;
    }

    err = spearman(cmd.list, (const double **) Z, datainfo, OPT_V, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        return 1;
    }

    strcpy(title, "gretl: ");
    strcat(title, _("rank correlation"));

    view_buffer(prn, 78, 400, title, SPEARMAN, NULL); 

    return 0;
}

int do_two_var_test (selector *sr)
{
    int action = selector_code(sr);
    const char *buf = selector_list(sr);
    PRN *prn;
    char title[64];
    int err = 0;

    if (buf == NULL) return 1;

    strcpy(title, "gretl: ");

    if (action == MEANTEST) {
	gretl_command_sprintf("meantest %s", buf);
	strcat(title, _("means test"));
    } else if (action == MEANTEST2) {
	gretl_command_sprintf("meantest %s --unequal-vars", buf);
	strcat(title, _("means test"));
    } else if (action == VARTEST) {
	gretl_command_sprintf("vartest %s", buf);
	strcat(title, _("variances test"));
    } else {
	dummy_call();
	return 1;
    }

    if (check_and_record_command() || bufopen(&prn)) {
	return 1;
    }

    if (action == MEANTEST) {
	err = means_test(cmd.list, (const double **) Z, datainfo, OPT_NONE, prn);
    } else if (action == MEANTEST2) {
	err = means_test(cmd.list, (const double **) Z, datainfo, OPT_O, prn);
    } else if (action == VARTEST) {
	err = vars_test(cmd.list, (const double **) Z, datainfo, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 300, title, action, NULL); 
    }

    return err;
}

void open_info (gpointer data, guint edit, GtkWidget *widget)
{
    if (datainfo->descrip == NULL) {
	if (yes_no_dialog(_("gretl: add info"), 
			  _("The data file contains no informative comments.\n"
			    "Would you like to add some now?"), 
			  0) == GRETL_YES) {
	    edit_header(NULL, 0, NULL);
	}
    } else {
	char *buf = g_strdup(datainfo->descrip);
	PRN *prn;
	
	if (buf != NULL) {
	    prn = gretl_print_new_with_buffer(buf);
	    view_buffer(prn, 80, 400, _("gretl: data info"), INFO, NULL);
	}
    }
}

void gui_errmsg (const int errcode)
{
    const char *msg = get_gretl_errmsg();
    char errtext[MAXLEN];

    *errtext = 0;

    if (*msg != '\0') {
	errbox(msg);
    } else {
	msg = get_errmsg(errcode, errtext, NULL);
	if (msg != NULL) {
	    errbox(msg);
	} else {
	    errbox(_("Unspecified error"));
	}
    }
}

/* OPT_M  drop all obs with missing data values 
   OPT_O  sample using dummy variable
   OPT_R  sample using boolean expression
   OPT_N  random sub-sample
   OPT_C  replace current restriction
*/

int bool_subsample (gretlopt opt)
{
    PRN *prn;
    const char *smplmsg;
    int err;

    if (bufopen(&prn)) {
	return 1;
    }

    if (opt & OPT_M) {
	err = restrict_sample(NULL, NULL, &Z, &datainfo, opt, prn);
    } else {
	err = restrict_sample(cmdline, NULL, &Z, &datainfo, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	goto alldone;
    } 

    smplmsg = gretl_print_get_buffer(prn);
    if (smplmsg != NULL && *smplmsg != 0) {
	infobox(smplmsg);
	goto alldone;
    }

    if (dataset_is_panel(datainfo) || dataset_is_time_series(datainfo)) {
	set_sample_label(datainfo);
    } else {
	/* special for undated data */
	set_sample_label_special();
    }

    restore_sample_state(TRUE);

    if (opt & OPT_M) {
	infobox(_("Sample now includes only complete observations"));
    } 

 alldone:

    gretl_print_destroy(prn);

    return err;
}

void drop_all_missing (gpointer data, guint opt, GtkWidget *w)
{
    int err = bool_subsample(OPT_M);

    if (!err) {
	gretl_command_strcpy("smpl --no-missing");
	check_and_record_command();
    }
}

void do_samplebool (GtkWidget *widget, dialog_t *dlg)
{
    const gchar *buf = edit_dialog_get_text(dlg);
    gretlopt opt;
    int err;

    if (buf == NULL) return;

    opt = edit_dialog_get_opt(dlg);

    if (opt & OPT_C) { 
	gretl_command_sprintf("smpl %s --restrict --replace", buf); 
    } else {
	gretl_command_sprintf("smpl %s --restrict", buf);
    }

    if (check_and_record_command()) {
	return;
    }

    err = bool_subsample(opt | OPT_R);

    if (!err) {
	close_dialog(dlg);
    }
}

int do_set_sample (void)
{
    return set_sample(cmdline, (const double **) Z, datainfo);
}

void count_missing (void)
{
    PRN *prn;

    if (bufopen(&prn)) return;
    if (count_missing_values(&Z, datainfo, prn)) {
	view_buffer(prn, 78, 300, _("gretl: missing values info"), 
		    SMPL, NULL);
    } else {
	infobox(_("No missing data values"));
	gretl_print_destroy(prn);
    }
}

void do_add_markers (const char *fname) 
{
    if (add_obs_markers_from_file(datainfo, fname)) { 
	errbox(_("Failed to add case markers"));
    } else {
	mark_dataset_as_modified();
	add_remove_markers_state(TRUE);
    }
}

void do_remove_markers (gpointer data, guint u, GtkWidget *w) 
{
    dataset_destroy_obs_markers(datainfo);
    mark_dataset_as_modified();
    add_remove_markers_state(FALSE);
}

int out_of_sample_info (int add_ok, int *t2)
{
    const char *can_add = 
	N_("There are no observations available for forecasting\n"
	   "out of sample.  You can add some observations now\n"
	   "if you wish.");
    int err = 0;

    if (add_ok) {
	int n = add_obs_dialog(can_add, 0);

	if (n < 0) {
	    err = 1;
	} else if (n > 0) {
	    set_original_n(datainfo->n);
	    err = dataset_add_observations(n, &Z, datainfo, OPT_A);
	    if (err) {
		gui_errmsg(err);
	    } else {
		mark_dataset_as_modified();
		drop_obs_state(TRUE);
		*t2 += n;
	    }
	} 
    } else if (!expert) {
	infobox(_("There are no observations available for forecasting\n"
		  "out of sample.  If you wish, you can add observations\n"
		  "(Data menu, Edit data), or you can shorten the sample\n"
		  "range over which the model is estimated (Sample menu)."));
    }

    return err;
}

void do_forecast (gpointer data, guint u, GtkWidget *w) 
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    char startobs[OBSLEN], endobs[OBSLEN];
    int t2, t1 = 0;
    int dyn_ok, add_obs_ok;
    int premax, pre_n = 0;
    int t1min = 0;
    int dt2 = datainfo->n - 1;
    int st2 = datainfo->n - 1;
    gretlopt opt = OPT_NONE;
    FITRESID *fr;
    PRN *prn;
    int resp, err;

    /* try to figure which options might be applicable */
    forecast_options_for_model(pmod, (const double **) Z, datainfo, &dyn_ok, 
			       &add_obs_ok, &dt2, &st2);

    if (dyn_ok) {
	t2 = dt2;
    } else {
	t2 = st2;
    }

    /* if no out-of-sample obs are available in case of time-
       series data, alert the user */
    if (t2 <= pmod->t2 && dataset_is_time_series(datainfo)) {
	err = out_of_sample_info(add_obs_ok, &t2);
	if (err) {
	    return;
	}
    }

    /* max number of pre-forecast obs in "best case" */
    premax = datainfo->n - 1;

    /* if there are spare obs available, default to an
       out-of-sample forecast */
    if (t2 > pmod->t2) {
	t1 = pmod->t2 + 1;
	pre_n = pmod->t2 / 2;
	if (pre_n > 100) {
	    pre_n = 100;
	}
	if (pmod->ci == GARCH) {
	    /* force out-of-sample fcast */
	    t1min = t1;
	}
    } else {
	pre_n = 0;
    }
    
    resp = forecast_dialog(t1min, t2, &t1,
			   0, t2, &t2,
			   0, premax, &pre_n,
			   dyn_ok);
    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt = OPT_D;
    } else if (resp == 2) {
	opt = OPT_S;
    }

    ntodate(startobs, t1, datainfo);
    ntodate(endobs, t2, datainfo);

    gretl_command_sprintf("fcasterr %s %s%s", startobs, endobs,
			  print_flags(opt, FCASTERR));

    if (check_and_record_command() || bufopen(&prn)) {
	return;
    }

    fr = get_forecast(pmod, t1 - pre_n, t1, t2, &Z, datainfo, opt);

    if (fr == NULL) {
	errbox(_("Failed to generate fitted values"));
	gretl_print_destroy(prn);
    } else if (fr->err) {
	gui_errmsg(fr->err);
	free_fit_resid(fr);
    } else {
	gretlopt popt = (LIMDEP(pmod->ci))? OPT_NONE : OPT_P;
	int width = 78;
	
	err = text_print_forecast(fr, &Z, datainfo, popt, prn);
	if (!err && popt == OPT_P) {
	    register_graph();
	}
	if (fr->sderr == NULL) {
	    width = 50;
	}
	view_buffer(prn, width, 400, _("gretl: forecasts"), FCASTERR, fr);
    }
}

int do_coeff_sum (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    PRN *prn;
    char title[48];
    MODEL *pmod;
    gint err;

    if (buf == NULL) {
	return 0;
    }

    gretl_command_sprintf("coeffsum %s", buf);

    if (check_lib_command() || bufopen(&prn)) {
	return 1;
    }

    pmod = vwin->data;
    err = sum_test(cmd.list, pmod, &Z, datainfo, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        return 1;
    }

    model_command_init(pmod->ID);

    strcpy(title, "gretl: ");
    strcat(title, _("Sum of coefficients"));

    view_buffer(prn, 78, 200, title, COEFFSUM, NULL); 

    return 0;
}

int do_add_omit (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    PRN *prn;
    char title[48];
    MODEL *orig, *pmod;
    gint err;

    if (buf == NULL) {
	return 1;
    }

    orig = vwin->data;

    if (selector_code(sr) == ADD) {
        gretl_command_sprintf("addto %d %s", orig->ID, buf);
    } else {
        gretl_command_sprintf("omitfrom %d %s", orig->ID, buf);
    }

    if (check_lib_command() || bufopen(&prn)) {
	return 1;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	gretl_print_destroy(prn);
	return 0;
    }

    if (selector_code(sr) == ADD) { 
        err = add_test(cmd.list, orig, pmod, &Z, datainfo, OPT_S, prn);
    } else {
        err = omit_test(cmd.list, orig, pmod, &Z, datainfo, OPT_S, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
	gretl_model_free(pmod);
        return err;
    }

    update_model_tests(vwin);

    if (lib_cmd_init()) {
	errbox(_("Error saving model information"));
	return 0;
    }

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);

    gretl_object_ref(pmod, GRETL_OBJ_EQN);

    sprintf(title, _("gretl: model %d"), pmod->ID);
    view_model(prn, pmod, 78, 420, title);

    return 0;
}

int do_VAR_omit (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    int *omitlist;
    GRETL_VAR *orig;
    GRETL_VAR *var = NULL;
    PRN *prn;
    gint err;

    if (buf == NULL) {
	return 1;
    }

    orig = vwin->data;

    if (bufopen(&prn)) {
	return 1;
    }

    omitlist = gretl_list_from_string(buf);
    if (omitlist == NULL) {
	err = E_ALLOC;
    } else {
	var = gretl_VAR_omit_test(omitlist, orig, &Z, datainfo, prn, &err);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
		    VAR, var);
    }

    return err;
}

int do_confidence_region (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    MODEL *pmod;

    char *mask = NULL;
    char iname[VNAMELEN];
    char jname[VNAMELEN];
    gretl_matrix *V = NULL;
    int v[2];
    double b[2];
    double t, kF;
    int err;

    if (buf == NULL || sscanf(buf, "%d %d", &v[0], &v[1]) != 2) {
	return 0;
    }

    pmod = (MODEL *) vwin->data;
    if (pmod == NULL) {
	return 0;
    }

    mask = calloc(pmod->ncoeff, 1);
    if (mask == NULL) {
	return 0;
    }

    mask[v[0]] = mask[v[1]] = 1;

    V = gretl_vcv_matrix_from_model(pmod, mask);
    if (V == NULL) {
	free(mask);
	return 0;
    }

    b[0] = pmod->coeff[v[0]];
    b[1] = pmod->coeff[v[1]];

    t = tcrit95(pmod->dfd);
    kF = 2.0 * f_critval(.05, 2, pmod->dfd);

    gretl_model_get_param_name(pmod, datainfo, v[0], iname);
    gretl_model_get_param_name(pmod, datainfo, v[1], jname);

    err = confidence_ellipse_plot(V, b, t, kF, iname, jname);
    gui_graph_handler(err);

    gretl_matrix_free(V);
    free(mask);

    return 0;
}

static void print_test_to_window (const MODEL *pmod, GtkWidget *w)
{
    if (w != NULL) {
	GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
	GtkTextIter iter;
	const char *txt;
	PRN *prn;

	if (bufopen(&prn)) return;

	gretl_model_print_last_test(pmod, prn);
	txt = gretl_print_get_buffer(prn);

	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_insert(buf, &iter, txt, -1);
	gretl_print_destroy(prn);
    }
}

static void update_model_tests (windata_t *vwin)
{
    MODEL *pmod = (MODEL *) vwin->data;

    if (pmod->ntests > vwin->n_model_tests) {
	print_test_to_window(pmod, vwin->w);
	vwin->n_model_tests += 1;
    }
}

void do_lmtest (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    PRN *prn;
    char title[64];
    int err = 0;

    if (bufopen(&prn)) return;

    strcpy(title, _("gretl: LM test "));

    if (action == LMTEST_WHITE) {
	gretl_command_strcpy("lmtest --white");
	err = whites_test(pmod, &Z, datainfo, OPT_S, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcat(title, _("(heteroskedasticity)"));
	}
    } else if (action == LMTEST_GROUPWISE) {
	gretl_command_strcpy("lmtest --panel");
	err = groupwise_hetero_test(pmod, &Z, datainfo, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcpy(title, _("gretl: groupwise heteroskedasticity"));
	}
    } else {
	int aux = (action == LMTEST_SQUARES)? AUX_SQ : AUX_LOG;

	if (action == LMTEST_SQUARES) { 
	    gretl_command_strcpy("lmtest --squares");
	} else {
	    gretl_command_strcpy("lmtest --logs");
	}
	clear_model(models[0]);
	err = nonlinearity_test(pmod, &Z, datainfo, aux, OPT_S, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcat(title, _("(non-linearity)"));
	} 
    }

    if (!err) {
	update_model_tests(vwin);
	model_command_init(pmod->ID);
	view_buffer(prn, 78, 400, title, LMTEST, NULL); 
    }
}

void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    PRN *prn;
    gretlopt opt = OPT_NONE;
    int err;

    if (bufopen(&prn)) {
	return;
    }	
	
    err = panel_diagnostics(pmod, &Z, datainfo, opt, prn);
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 400, _("gretl: panel model diagnostics"), 
		    PANEL, NULL);
    }
}

static void set_model_id_on_window (GtkWidget *w, int ID)
{
    g_object_set_data(G_OBJECT(w), "model_ID", 
		      GINT_TO_POINTER(ID));
}

static int get_model_id_from_window (GtkWidget *w)
{
    return GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "model_ID"));
}

static int make_and_display_graph (void)
{
    if (gnuplot_make_graph()) {
	errbox(_("gnuplot command failed"));
	return 1;
    } 

    register_graph();

    return 0;
}

void add_leverage_data (windata_t *vwin)
{
    void *handle;
    unsigned char (*leverage_data_dialog) (void);
    gretl_matrix *m = (gretl_matrix *) vwin->data;
    unsigned char opt;
    int err;

    if (m == NULL) return;

    leverage_data_dialog = gui_get_plugin_function("leverage_data_dialog",
						   &handle);
    if (leverage_data_dialog == NULL) return;

    opt = leverage_data_dialog();
    close_plugin(handle);

    if (opt == 0) return;

    err = add_leverage_values_to_dataset(&Z, datainfo, m, opt);

    if (err) {
	gui_errmsg(err);
    } else {
	int ID = get_model_id_from_window(vwin->dialog);

	gretl_command_strcpy("leverage --save");
	model_command_init(ID);
    }
}

void do_leverage (gpointer data, guint u, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, double ***, 
				     DATAINFO *, PRN *, int);
    PRN *prn;
    gretl_matrix *m;

    model_leverage = gui_get_plugin_function("model_leverage", 
					     &handle);
    if (model_leverage == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    m = (*model_leverage)(pmod, &Z, datainfo, prn, 1);
    close_plugin(handle);

    if (m != NULL) {
	windata_t *levwin;

	levwin = view_buffer(prn, 78, 400, _("gretl: leverage and influence"), 
			     LEVERAGE, m); 
	set_model_id_on_window(levwin->dialog, pmod->ID);

	make_and_display_graph();

	gretl_command_strcpy("leverage");
	model_command_init(pmod->ID);
    } else {
	errbox(_("Command failed"));
    }
}

void do_vif (gpointer data, guint u, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    int (*print_vifs) (MODEL *, double ***, DATAINFO *, PRN *);
    void *handle;
    int err;
    PRN *prn;

    print_vifs = gui_get_plugin_function("print_vifs", &handle);
    if (print_vifs == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    err = (*print_vifs)(pmod, &Z, datainfo, prn);
    close_plugin(handle);

    if (!err) {
	windata_t *vifwin;

	vifwin = view_buffer(prn, 78, 400, _("gretl: collinearity"), 
			     PRINT, NULL); 

	gretl_command_strcpy("vif");
	model_command_init(pmod->ID);
    } else {
	errbox(_("Command failed"));
    }
}

static int reject_scalar (int vnum)
{
    if (var_is_scalar(datainfo, vnum)) {
	errbox(_("variable %s is a scalar"), datainfo->varname[vnum]);
	return 1;
    }

    return 0;
}

void do_gini (gpointer data, guint u, GtkWidget *w)
{
    gretlopt opt = OPT_NONE;
    PRN *prn;
    int v = mdata_active_var();
    int err;

    if (reject_scalar(v)) {
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    err = gini(v, &Z, datainfo, opt, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("Gini coefficient"));

	view_buffer(prn, 78, 200, title, PRINT, NULL);
	g_free(title);
	register_graph();
    } 
}

void do_kernel (gpointer data, guint u, GtkWidget *w)
{
    void *handle;
    int (*kernel_density) (int, const double **, const DATAINFO *,
			   double, gretlopt);
    gretlopt opt = OPT_NONE;
    double bw = 1.0;
    int v = mdata_active_var();
    int err;

    if (reject_scalar(v)) {
	return;
    }

    err = density_dialog(v, &bw);
    if (err < 0) {
	return;
    }

    if (err > 0) {
	opt |= OPT_O;
    }

    kernel_density = gui_get_plugin_function("kernel_density", 
					     &handle);
    if (kernel_density == NULL) {
	return;
    }

    err = (*kernel_density)(v, (const double **) Z, 
			    datainfo, bw, opt);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } else {
	make_and_display_graph();
    } 
}

void do_chow_cusum (gpointer data, guint action, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    PRN *prn;
    int err = 0;

    if (pmod->ci != OLS) {
	errbox(_("This test only implemented for OLS models"));
	return;
    }

    if (action == CHOW) {
	char brkstr[OBSLEN];
	int resp, brk = (pmod->t2 - pmod->t1) / 2;

	set_window_busy(vwin);
	resp = get_obs_dialog(_("gretl: Chow test"), 
			      _("Observation at which to split the sample:"),
			      NULL, NULL, 
			      pmod->t1 + 1, pmod->t2 - 1, &brk,
			      0, 0, NULL);
	unset_window_busy(vwin);

	if (resp < 0) {
	    return;
	}

	ntodate(brkstr, brk, datainfo);
	gretl_command_sprintf("chow %s", brkstr);
    } else if (action == QLRTEST) {
	gretl_command_strcpy("qlrtest");
    } else {
	gretl_command_strcpy("cusum");
    }

    if (bufopen(&prn)) {
	return;
    }

    if (action == CHOW || action == QLRTEST) {
	err = chow_test(cmdline, pmod, &Z, datainfo, OPT_S, prn);
    } else {
	err = cusum_test(pmod, &Z, datainfo, OPT_S, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	if (action == CUSUM) {
	    register_graph();
	}

	update_model_tests(vwin);
	model_command_init(pmod->ID);

	view_buffer(prn, 78, 400, (action == CHOW)?
		    _("gretl: Chow test output") : 
		    (action == QLRTEST)?
		    _("gretl: QLR test output") : 
		    _("gretl: CUSUM test output"),
		    action, NULL);
    }
}

void do_reset (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    PRN *prn;
    char title[40];
    int err = 0;

    if (bufopen(&prn)) return;

    strcpy(title, _("gretl: RESET test"));

    err = reset_test(pmod, &Z, datainfo, OPT_S, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	update_model_tests(vwin);
	gretl_command_strcpy("reset");
	model_command_init(pmod->ID);
	view_buffer(prn, 78, 400, title, RESET, NULL); 
    }
}

void do_autocorr (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    PRN *prn;
    char title[64];
    int order, err = 0;

    order = default_lag_order(datainfo);

    set_window_busy(vwin);
    err = spin_dialog(_("gretl: autocorrelation"), 
		      &order, _("Lag order for test:"),
		      1, datainfo->n / 2, LMTEST);
    unset_window_busy(vwin);

    if (err < 0) {
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    strcpy(title, _("gretl: LM test (autocorrelation)"));

    if (dataset_is_panel(datainfo)) {
	err = panel_autocorr_test(pmod, order, Z, datainfo,
				  OPT_S, prn);
    } else {
	err = autocorr_test(pmod, order, &Z, datainfo, OPT_S, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	update_model_tests(vwin);
	gretl_command_sprintf("lmtest --autocorr %d", order);
	model_command_init(pmod->ID);
	view_buffer(prn, 78, 400, title, LMTEST, NULL); 
    }
}

void do_arch (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = vwin->data;
    PRN *prn;
    char tmpstr[26];
    int smpl_t1 = datainfo->t1;
    int smpl_t2 = datainfo->t2;
    int i, order;
    int err = 0;

    order = default_lag_order(datainfo);

    set_window_busy(vwin);
    err = spin_dialog(_("gretl: ARCH test"),
		      &order, _("Lag order for ARCH test:"),
		      1, datainfo->n / 2, ARCH);
    unset_window_busy(vwin);
    if (err < 0) {
	return;
    }

    gretl_command_sprintf("arch %d ", order);

    for (i=1; i<=pmod->list[0]; i++) {
	sprintf(tmpstr, "%d ", pmod->list[i]);
	gretl_command_strcat(tmpstr);
    }

    if (check_and_record_command()) {
	return;
    }

    order = atoi(cmd.param);
    if (!order) {
	errbox(_("Couldn't read ARCH order"));
	return;
    }

    if (bufopen(&prn)) return;

    clear_model(models[1]);

    /* temporarily reimpose the sample range in effect
       when pmod was estimated */
    impose_model_smpl(pmod, datainfo);

    /* FIXME opt? */
    *models[1] = arch_test(pmod, order, &Z, datainfo, OPT_S, prn);
    if ((err = (models[1])->errcode)) { 
	gui_errmsg(err);
    } else {
	update_model_tests(vwin);
    }

    clear_model(models[1]);

    /* restore original sample */
    datainfo->t1 = smpl_t1;
    datainfo->t2 = smpl_t2;

    if (err) {
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 400, _("gretl: ARCH test"), ARCH, NULL);
    }
}

static int model_error (MODEL *pmod)
{
    int err = 0;

    if (pmod->errcode) {
	if (pmod->errcode != E_CANCEL) {
	    gui_errmsg(pmod->errcode);
	}
	gretl_model_free(pmod);
	err = 1;
    }

    return err;
}

static int model_output (MODEL *pmod, PRN *prn)
{
    int err = 0;

    if (model_error(pmod)) {
	err = 1;
    } else {
	printmodel(pmod, datainfo, OPT_NONE, prn);
    }

    return err;
}

static gint check_model_cmd (void)
{
    int err = parse_command_line(cmdline, &cmd, &Z, datainfo); 

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int 
record_model_commands_from_buf (const gchar *buf, const MODEL *pmod,
				int got_start, int got_end)
{
    bufgets_init(buf);

    if (!got_start) {
	gretl_command_strcpy("restrict");
	model_command_init(pmod->ID);
    }

    gretl_cmd_set_context(&cmd, RESTRICT);

    while (bufgets(cmdline, MAXLINE, buf)) {
	if (string_is_blank(cmdline)) {
	    continue;
	}
	top_n_tail(cmdline);
	model_command_init(pmod->ID);
    }

    gretl_cmd_destroy_context(&cmd);

    if (!got_end) {
	gretl_command_strcpy("end restrict");
	model_command_init(pmod->ID);
    }

    return 0;
}

void do_restrict (GtkWidget *widget, dialog_t *dlg)
{
    MODEL *pmod = NULL;
    gretl_equation_system *sys = NULL;
    GRETL_VAR *vecm = NULL;

    gchar *buf;
    PRN *prn;
    char title[64], bufline[MAXLINE];
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    gretl_restriction_set *my_rset = NULL;
    int got_start_line = 0, got_end_line = 0;
    int height = 300;
    int err = 0;

    if (vwin->role == VIEW_MODEL) {
	pmod = (MODEL *) vwin->data;
    } else if (vwin->role == VECM) {
	vecm = (GRETL_VAR *) vwin->data;
    } else if (vwin->role == SYSTEM) {
	sys = (gretl_equation_system *) vwin->data;
    }

    if (pmod == NULL && vecm == NULL && sys == NULL) {
	close_dialog(dlg);
	return;
    }

    buf = edit_dialog_special_get_text(dlg);
    if (buf == NULL) return;

    bufgets_init(buf);

    while (bufgets(bufline, MAXLINE, buf) && !err) {
	if (string_is_blank(bufline)) {
	    continue;
	}

	top_n_tail(bufline);

	if (!strcmp(bufline, "end restrict")) {
	    got_end_line = 1;
	    break;
	} else if (!strncmp(bufline, "restrict", 8)) {
	    got_start_line = 1;
	}

	if (my_rset == NULL) {
	    if (pmod != NULL) {
		my_rset = restriction_set_start(bufline, pmod, datainfo,
						OPT_NONE);
	    } else if (sys != NULL) {
		my_rset = cross_restriction_set_start(bufline, sys);
	    } else {
		my_rset = var_restriction_set_start(bufline, vecm);
	    }
	    if (my_rset == NULL) {
 		err = 1;
		gui_errmsg(err);
	    }
	} else {
	    err = restriction_set_parse_line(my_rset, bufline);
	    if (err) {
		gui_errmsg(err);
	    }
	}
    }

    if (err) {
	g_free(buf);
	return;
    }

    close_dialog(dlg);

    if (bufopen(&prn)) return; 

    err = gretl_restriction_set_finalize(my_rset, (const double **) Z, 
					 datainfo, prn);

    if (err) {
	errmsg(err, prn);
    } else {
	if (pmod != NULL) {
	    record_model_commands_from_buf(buf, pmod, got_start_line,
					   got_end_line);
	} else if (sys != NULL) {
	    gretl_equation_system_estimate(sys, &Z, datainfo, OPT_NONE, prn);
	    height = 450;
	} else if (vecm != NULL) {
	    height = 450;
	}
    }

    g_free(buf);

    strcpy(title, "gretl: ");
    strcat(title, _("linear restrictions"));

    view_buffer(prn, 78, height, title, PRINT, NULL);
}

static int 
record_sys_commands_from_buf (const gchar *buf, const char *startline, 
			      int got_end_line)
{
    char bufline[MAXLINE];    

    bufgets_init(buf);

    while (bufgets(bufline, MAXLINE, buf)) {
	if (string_is_blank(bufline)) {
	    continue;
	}
	if (!strncmp(bufline, "system", 6)) {
	    add_command_to_stack(startline);
	} else {
	    top_n_tail(bufline);
	    add_command_to_stack(bufline);
	}
    }

    if (!got_end_line) {
	add_command_to_stack("end system");
    }

    return 0;
}

static void maybe_grab_system_name (const char *s, char *name)
{
    s = strstr(s, "name=");
    if (s != NULL) {
	s += 5;
	if (*s == '"') {
	    sscanf(s + 1, "%31[^\"]", name);
	} else {
	    sscanf(s, "%31s", name);
	}
    }
}

void do_eqn_system (GtkWidget *widget, dialog_t *dlg)
{
    gretl_equation_system *my_sys = NULL;
    gchar *buf;
    PRN *prn;
    char sysname[32];
    char bufline[MAXLINE];
    int *slist = NULL;
    char *startline = NULL;
    int got_end_line = 0;
    int method, err = 0;

    buf = edit_dialog_special_get_text(dlg);
    if (buf == NULL) {
	return;
    }

    method = edit_dialog_get_opt(dlg);

    bufgets_init(buf);
    *sysname = 0;

    while (bufgets(bufline, MAXLINE, buf) && !err) {
	if (string_is_blank(bufline) || *bufline == '#') {
	    continue;
	}

	top_n_tail(bufline);

	if (!strcmp(bufline, "end system")) {
	    got_end_line = 1;
	    break;
	}	    

	if (!strncmp(bufline, "system", 6)) {
	    maybe_grab_system_name(bufline, sysname);
	    continue;
	} 

	if (my_sys == NULL) {
	    startline = g_strdup_printf("system method=%s", 
					system_method_short_string(method));
	    my_sys = system_start(startline, OPT_NONE); /* FIXME opt? */
	    if (my_sys == NULL) {
		fprintf(stderr, "do_eqn_system: sys is NULL\n");
 		err = 1;
	    }
	}

	if (err) {
	    gui_errmsg(err);
	    break;
	}

	if (!strncmp(bufline, "equation", 8)) {
	    slist = command_list_from_string(bufline);
	    if (slist == NULL) {
		/* error message should be emitted by selector in this
		   case? */
		err = 1;
	    } else {
		err = gretl_equation_system_append(my_sys, slist);
		free(slist);
		if (err) {
		    /* note: sys is destroyed on error */
		    gui_errmsg(err);
		} 
	    }
	} else {
	    err = system_parse_line(my_sys, bufline, datainfo);
	    if (err) {
		/* sys is destroyed on error */
		gui_errmsg(err);
	    }
	} 
    }

    if (err) {
	g_free(buf);
	return;
    }

    close_dialog(dlg);

    if (bufopen(&prn)) {
	g_free(buf);
	return; 
    }

    err = gretl_equation_system_finalize(my_sys, &Z, datainfo, prn);
    if (err) {
	errmsg(err, prn);
    } else {
	record_sys_commands_from_buf(buf, startline, got_end_line);
	if (*sysname != 0) {
	    my_sys->name = g_strdup(sysname);
	}
    }

    g_free(buf);
    g_free(startline);

    view_buffer(prn, 78, 450, 
		(my_sys->name != NULL)? my_sys->name: 
		_("gretl: simultaneous equations system"), 
		SYSTEM, my_sys);
}

void do_saved_eqn_system (GtkWidget *widget, dialog_t *dlg)
{
    gretl_equation_system *my_sys;
    PRN *prn;
    int err = 0;

    my_sys = (gretl_equation_system *) edit_dialog_get_data(dlg);
    if (my_sys == NULL) {
	return;
    }

    my_sys->method = edit_dialog_get_opt(dlg);

    close_dialog(dlg);

    if (bufopen(&prn)) {
	return; 
    }

    err = gretl_equation_system_estimate(my_sys, &Z, datainfo,
					 OPT_NONE, prn);
    if (err) {
	errmsg(err, prn);
    } 

    /* ref count? */

    view_buffer(prn, 78, 450, my_sys->name, SYSTEM, my_sys);
}

static int do_nls_genr (void)
{
    int err;

    if (check_and_record_command()) {
	err = 1;
    } else {
	err = finish_genr(NULL, NULL);
    }

    return err;
}

static int is_genr_line (char *s)
{
    if (!strncmp(s, "genr", 4) ||
	!strncmp(s, "series", 6) ||
	!strncmp(s, "scalar", 6)) {
	return 1;
    } else if (!strncmp(s, "param ", 6) && strchr(s, '=')) {
	gchar *tmp = g_strdup_printf("genr %s", s + 6);
	
	strcpy(s, tmp);
	g_free(tmp);
	return 1;
    } else {
	return 0;
    }
}

static void real_do_nonlinear_model (dialog_t *dlg, int ci)
{
    gchar *buf = edit_dialog_special_get_text(dlg);
    gretlopt opt = edit_dialog_get_opt(dlg);
    char realline[MAXLINE];
    char bufline[MAXLINE];
    char title[26];
    int err = 0, started = 0;
    MODEL *pmod = NULL;
    const char *cstr;
    const char *endstr;
    PRN *prn;

    if (buf == NULL) {
	return;
    }
    
    if (ci == NLS) {
	cstr = "nls";
	endstr = "end nls";
    } else {
	cstr = "mle";
	endstr = "end mle";
    }

    bufgets_init(buf);
    *realline = 0;

    while (bufgets(bufline, MAXLINE-1, buf) && !err) {
	int len, cont = 0;

	if (string_is_blank(bufline)) {
	    *realline = 0;
	    continue;
	}

	cont = top_n_tail(bufline);

	len = strlen(bufline) + strlen(realline);
	if (len > MAXLINE - 1) {
	    errbox("command line is too long (maximum is %d characters)", 
		   MAXLINE - 1);
	    err = 1;
	    break;
	}

	strcat(realline, bufline);

	if (cont) {
	    continue;
	}

	if (started && !strncmp(realline, endstr, 7)) {
	    break;
	}

	if (!started && is_genr_line(realline)) {
	    gretl_command_strcpy(realline);
	    err = do_nls_genr();
	    *realline = 0;
	    continue;
	}

	if (!started && strncmp(realline, cstr, 3)) {
	    char tmp[MAXLINE];
	    
	    strcpy(tmp, realline);
	    strcpy(realline, cstr);
	    strcat(realline, " ");
	    strcat(realline, tmp);
	} 

	err = nls_parse_line(ci, realline, (const double **) Z, datainfo, NULL);
	started = 1;

	if (err) {
	    gui_errmsg(err);
	} else {
	    gretl_command_strcpy(realline);
	    err = lib_cmd_init();
	}

	*realline = 0;
    }

    g_free(buf);

    if (err) {
	return;
    }

    /* if the user didn't give "end XXX", supply it */
    if (strncmp(bufline, endstr, 7)) {
	gretl_command_strcpy(endstr);
	lib_cmd_init();
    }    

    if (bufopen(&prn)) return;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	return;
    }

    *pmod = nls(&Z, datainfo, opt, prn);
    err = model_output(pmod, prn);

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    close_dialog(dlg);

    sprintf(title, _("gretl: model %d"), pmod->ID);

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);

    gretl_object_ref(pmod, GRETL_OBJ_EQN);
    
    view_model(prn, pmod, 78, 420, title); 
}

void do_nls_model (GtkWidget *widget, dialog_t *dlg)
{
    real_do_nonlinear_model(dlg, NLS);
}

void do_mle_model (GtkWidget *widget, dialog_t *dlg)
{
    real_do_nonlinear_model(dlg, MLE);
}

static int logistic_model_get_lmax (CMD *cmd)
{
    double ymax, lmax;
    int err;

    err = logistic_ymax_lmax(Z[cmd->list[1]], datainfo, &ymax, &lmax);

    if (!err) {
	lmax_dialog(&lmax, ymax);
	if (na(lmax)) {
	    err = 1;
	} else if (lmax == 0.0) {
	    /* canceled */
	    err = -1;
	} else {
	    free(cmd->param);
	    cmd->param = g_strdup_printf("ymax=%g", lmax);
	    gretl_command_strcat(" ");
	    gretl_command_strcat(cmd->param);
	}
    }

    return err;
}

int do_model (selector *sr) 
{
    const char *buf;
    PRN *prn;
    MODEL *pmod;
    char title[26], estimator[9];
    int action;
    double rho;
    int err = 0;

    if (selector_error(sr)) {
	return 1;
    }

    buf = selector_list(sr);
    if (buf == NULL) {
	return 1;
    }

    action = selector_code(sr);
    strcpy(estimator, gretl_command_word(action));

    cmd.opt = selector_get_opts(sr);

    gretl_command_sprintf("%s %s%s", estimator, buf, 
			  print_flags(cmd.opt, action));

#if 0
    fprintf(stderr, "do_model: cmdline = '%s'\n", cmdline);
#endif

    if (check_model_cmd() || bufopen(&prn)) {
	return 1;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	return 1;
    }

    switch (action) {

    case CORC:
    case HILU:
    case PWE: 
	rho = estimate_rho(cmd.list, &Z, datainfo, action, 
			   &err, (cmd.opt | OPT_P), prn);
	if (err) {
	    gui_errmsg(err);
	    break;
	}
	*pmod = ar1_lsq(cmd.list, &Z, datainfo, action, OPT_NONE, rho);
	err = model_output(pmod, prn);
	if (action == HILU) {
	    register_graph();
	}
	break;

    case OLS:
    case WLS:
    case HCCM:
	*pmod = lsq(cmd.list, &Z, datainfo, action, cmd.opt);
	err = model_output(pmod, prn);
	break;

    case PANEL:
	*pmod = panel_model(cmd.list, &Z, datainfo, cmd.opt, prn);
	err = model_output(pmod, prn);
	break;

    case HSK:
	*pmod = hsk_func(cmd.list, &Z, datainfo);
	err = model_output(pmod, prn);
	break;

    case TSLS:
	*pmod = tsls_func(cmd.list, TSLS, &Z, datainfo, cmd.opt);
	err = model_output(pmod, prn);
	break;

    case AR:
	*pmod = ar_func(cmd.list, &Z, datainfo, OPT_NONE, prn);
	err = model_error(pmod);
	break;

    case LOGIT:
    case PROBIT:
	*pmod = logit_probit(cmd.list, &Z, datainfo, action, cmd.opt);
	err = model_output(pmod, prn);
	break;

    case TOBIT:
	*pmod = tobit_model(cmd.list, &Z, datainfo, 
			    (cmd.opt & OPT_V)? prn : NULL); 
	err = model_output(pmod, prn);
	break;

    case POISSON:
	*pmod = poisson_model(cmd.list, &Z, datainfo,
			      (cmd.opt & OPT_V)? prn : NULL);
	err = model_output(pmod, prn);
	break;

    case ARMA:
	*pmod = arma(cmd.list, (const double **) Z, datainfo,
		     cmd.opt, prn);
	err = model_output(pmod, prn);
	break;

    case GARCH:
	*pmod = garch(cmd.list, &Z, datainfo, cmd.opt, prn); 
	err = model_output(pmod, prn);
	break;

    case LOGISTIC:
	err = logistic_model_get_lmax(&cmd);
	if (err < 0) {
	    return 1;
	} else if (err) {
	    gui_errmsg(err);
	    break;
	} else {
	    *pmod = logistic_model(cmd.list, &Z, datainfo, cmd.param);
	    err = model_output(pmod, prn);
	}
	break;	

    case LAD:
	*pmod = lad(cmd.list, &Z, datainfo);
	err = model_output(pmod, prn);
	break;	

    case MPOLS:
	*pmod = mp_ols(cmd.list, (const double **) Z, datainfo);
	err = model_output(pmod, prn);
	break;	

    default:
	errbox(_("Sorry, not implemented yet!"));
	break;
    }

    if (err) {
	if (action == ARMA && (cmd.opt & OPT_V) && !(cmd.opt & OPT_X)) {
	    /* non-convergence info? */
	    view_buffer(prn, 78, 400, _("gretl: ARMA"), PRINT, NULL);
	} else if (action == GARCH && (cmd.opt & OPT_V)) {
	    /* ditto */
	    view_buffer(prn, 78, 400, _("gretl: GARCH"), PRINT, NULL);
	} else {
	    gretl_print_destroy(prn);
	}
	return err;
    }

    if (lib_cmd_init()) {
	errbox(_("Error saving model information"));
	return 0;
    }

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);

    gretl_object_ref(pmod, GRETL_OBJ_EQN);
    
    sprintf(title, _("gretl: model %d"), pmod->ID);
    view_model(prn, pmod, 78, 420, title); 

    return 0;
}

int do_vector_model (selector *sr) 
{
    GRETL_VAR *var;
    char estimator[9];
    const char *buf;
    PRN *prn;
    int order, action;
    int err = 0;

    if (selector_error(sr)) {
	return 1;
    }

    buf = selector_list(sr);
    if (buf == NULL) {
	return 1;
    }

    cmd.opt = selector_get_opts(sr);
    action = selector_code(sr);

    if (action == VLAGSEL) {
	cmd.opt |= OPT_L;
	action = VAR;
    }

    strcpy(estimator, gretl_command_word(action));

    gretl_command_sprintf("%s %s%s", estimator, buf, 
			  print_flags(cmd.opt, action));

#if 0
    fprintf(stderr, "do_vector_model: cmdline = '%s'\n", cmdline);
#endif

    if (check_model_cmd() || bufopen(&prn)) {
	return 1;
    }

    sscanf(buf, "%d", &order);
    if (order > var_max_order(cmd.list, datainfo)) {
	errbox(_("Insufficient degrees of freedom for regression"));
	gretl_print_destroy(prn);
	return 1;
    }    

    if (action == VAR && !(cmd.opt & OPT_L)) {
	/* regular VAR, not VAR lag selection */
	var = gretl_VAR(order, cmd.list, &Z, datainfo, cmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
			VAR, var);
	}
    } else if (action == VAR) {
	/* VAR lag selection */
	gretl_VAR(order, cmd.list, &Z, datainfo, cmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 72, 350, _("gretl: VAR lag selection"), 
			PRINT, NULL);
	}	
    } else if (action == VECM) {
	/* Vector Error Correction Model */
	var = vecm(order, atoi(cmd.extra), cmd.list, &Z, datainfo, cmd.opt, 
		   prn, &err);
	if (!err) {
	    view_buffer(prn, 78, 450, _("gretl: VECM"), VECM, var);
	}
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    }

    return err;
}

void do_graph_model (GPT_SPEC *spec)
{
    char *buf;
    PRN *prn;
    MODEL *pmod = NULL;
    char title[26];
    int err = 0;

    if (spec == NULL || spec->reglist == NULL) {
	return;
    }

    buf = gretl_list_to_string(spec->reglist);
    if (buf == NULL) {
	return;
    }

    gretl_command_sprintf("ols%s", buf);
    free(buf);

    if (check_model_cmd() || bufopen(&prn)) {
	return;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	return;
    }

    *pmod = lsq(cmd.list, &Z, datainfo, OLS, cmd.opt);
    err = model_output(pmod, prn);

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (lib_cmd_init()) {
	errbox(_("Error saving model information"));
	return;
    }

    attach_subsample_to_model(pmod, datainfo);

    gretl_object_ref(pmod, GRETL_OBJ_EQN);
    
    sprintf(title, _("gretl: model %d"), pmod->ID);
    view_model(prn, pmod, 78, 420, title);     
}

void do_minibuf (GtkWidget *widget, dialog_t *dlg) 
{
    const gchar *buf = edit_dialog_get_text(dlg);
    int oldv = datainfo->v;
    int err;

    if (buf == NULL) return;

    gretl_command_sprintf("%s", buf);

    if (dlg != NULL) {
	close_dialog(dlg);
    }

    console_record_sample(datainfo);

    err = gui_exec_line(cmdline, NULL, CONSOLE_EXEC, NULL);
    if (err) {
	gui_errmsg(err);
    }

    /* update variable listing in main window if needed */
    if (datainfo->v != oldv || !strncmp(cmdline, "rename", 6)) {
	populate_varlist();
    }

    /* update sample info and options if needed */
    if (console_sample_changed(datainfo)) {
	set_sample_label(datainfo);
    }
}

void do_genr (GtkWidget *widget, dialog_t *dlg) 
{
    const gchar *buf = edit_dialog_get_text(dlg);
    char test[8];

    if (buf == NULL) return;

    while (isspace((unsigned char) *buf)) buf++;
    sscanf(buf, "%7s", test);
    if (!strcmp(test, "series")) {
	gretl_command_sprintf("%s", buf);
    } else {
	gretl_command_sprintf("genr %s", buf);
    }

    if (check_and_record_command()) {
	return;
    }

    finish_genr(NULL, dlg);
}

void do_model_genr (GtkWidget *widget, dialog_t *dlg) 
{
    const gchar *buf;
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    MODEL *pmod = vwin->data;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    gretl_command_sprintf("genr %s", buf);

    if (model_command_init(pmod->ID)) {
	return;
    }

    finish_genr(pmod, dlg);
}

static void do_seed (guint32 newseed)
{
    guint32 oldseed = get_gretl_random_seed();

    if (newseed == oldseed) {
	return;
    }
	
    gretl_command_sprintf("set seed %u", newseed); 
    if (check_and_record_command()) {
	return;
    }

    gretl_rand_set_seed(newseed);
}

static void real_do_random (dialog_t *dlg, int action) 
{
    char vname[VNAMELEN];
    const gchar *buf;
    double f1, f2;
    double *s;
    int v;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    if (action == RANDOM_CHISQ || action == RANDOM_ST) {
	if (sscanf(buf, "%15s %d", vname, &v) != 2) {
	    errbox(_("Specification is malformed\n"
		     "Should be like \"foo 5\""));
	    return;
	}
    } else if (sscanf(buf, "%15s %lf %lf", vname, &f1, &f2) != 3) {
	if (action == RANDOM_NORMAL) {
	    errbox(_("Specification is malformed\n"
		     "Should be like \"foo 1 2.5\""));
	} else {
	    errbox(_("Specification is malformed\n"
		     "Should be like \"foo 0 10\""));
	}
	return;
    }

    if (action == RANDOM_NORMAL && f2 <= 0.0) {
	errbox(_("Can't have a negative standard deviation!"));
	return;
    } else if (action == RANDOM_UNIFORM && f1 >= f2) {
	errbox(_("Range is non-positive!"));
	return;
    } else if ((action == RANDOM_CHISQ || action == RANDOM_ST) 
	       && v < 1) {
	errbox(_("The degrees of freedom must be positive"));
	return;
    }	

    if (validate_varname(vname)) {
	return;
    }

    s = (double *) edit_dialog_get_data(dlg);
    if (s != NULL) {
	do_seed((guint32) *s);
    }

    if (action == RANDOM_NORMAL) {
	if (f1 != 0. || f2 != 1.) {
	    gretl_command_sprintf("genr %s = %g * normal() + %g", 
				  vname, f2, f1);
	} else {
	    gretl_command_sprintf("genr %s = normal()", vname); 
	}
    } else if (action == RANDOM_UNIFORM) {
	if (f1 != 0. || f2 != 1.) {
	    gretl_command_sprintf("genr %s = %g + (uniform() * %g)", 
				  vname, f1, (f2 - f1));
	} else {
	    gretl_command_sprintf("genr %s = uniform()", vname); 
	}
    } else if (action == RANDOM_CHISQ) {
	gretl_command_sprintf("genr %s = chisq(%d)", vname, v); 
    } else if (action == RANDOM_ST) {
	gretl_command_sprintf("genr %s = student(%d)", vname, v); 
    }

    if (check_and_record_command()) {
	return;
    }

    finish_genr(NULL, dlg);
}

void do_random_uniform (GtkWidget *widget, dialog_t *dlg) 
{
    real_do_random(dlg, RANDOM_UNIFORM);
}

void do_random_normal (GtkWidget *widget, dialog_t *dlg) 
{
    real_do_random(dlg, RANDOM_NORMAL);
}

void do_random_chisq (GtkWidget *widget, dialog_t *dlg) 
{
    real_do_random(dlg, RANDOM_CHISQ);
}

void do_random_st (GtkWidget *widget, dialog_t *dlg) 
{
    real_do_random(dlg, RANDOM_ST);
}

static int finish_genr (MODEL *pmod, dialog_t *dlg)
{
    int err = 0;

    err = generate(cmdline, &Z, datainfo, OPT_NONE, NULL); 

    if (err) {
	gui_errmsg(err);
	delete_last_command();
    } else {
	if (dlg != NULL) {
	    close_dialog(dlg);
	}
	populate_varlist();
	mark_dataset_as_modified();
    }

    return err;
}

static int real_do_setmiss (double missval, int varno) 
{
    int i, t, count = 0;
    int start = 1, end = datainfo->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	if (var_is_scalar(datainfo, i)) {
	    continue;
	}
	for (t=0; t<datainfo->n; t++) {
	    if (Z[i][t] == missval) {
		Z[i][t] = NADBL;
		count++;
	    }
	}	
    }

    return count;
}

void do_global_setmiss (GtkWidget *widget, dialog_t *dlg)
{
    const gchar *buf;
    double missval;
    int count, err;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    if ((err = check_atof(buf))) {
	gui_errmsg(err);
	return;
    }

    missval = atof(buf);
    count = real_do_setmiss(missval, 0);

    close_dialog(dlg);

    if (count) {
	infobox(_("Set %d values to \"missing\""), count);
	mark_dataset_as_modified();
    } else {
	errbox(_("Didn't find any matching observations"));
    }	
}

void do_variable_setmiss (GtkWidget *widget, dialog_t *dlg)
{
    const gchar *buf;
    double missval;
    int v = mdata_active_var();
    int count, err;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    if (var_is_scalar(datainfo, v)) {
	close_dialog(dlg);
	errbox(_("This variable is a scalar"));
	return;
    }

    if ((err = check_atof(buf))) {
	gui_errmsg(err);
	return;
    }    

    missval = atof(buf);
    count = real_do_setmiss(missval, v);

    close_dialog(dlg);

    if (count) {
	infobox(_("Set %d observations to \"missing\""), count);
	mark_dataset_as_modified();
    } else {
	errbox(_("Didn't find any matching observations"));
    }
}

int do_rename_variable (int v, const char *newname, int full)
{
    int err = 0;

    gretl_command_sprintf("rename %d %s", v, newname);

    if (full) {
	err = check_and_record_command();
    } else {
	err = check_lib_command();
    }

    if (!err) {
	strcpy(datainfo->varname[v], newname);
    }

    return err;
}

int record_varlabel_change (int v)
{
    gretl_command_sprintf("label %s -d \"%s\" -n \"%s\"", 
			  datainfo->varname[v],
			  VARLABEL(datainfo, v), 
			  DISPLAYNAME(datainfo, v));

    return check_and_record_command();
}

static void normal_test (MODEL *pmod, FreqDist *freq)
{
    ModelTest *test = model_test_new(GRETL_TEST_NORMAL);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_NORMAL_CHISQ);
	model_test_set_dfn(test, 2);
	model_test_set_value(test, freq->test);
	model_test_set_pvalue(test, chisq_cdf_comp(freq->test, 2));
	maybe_add_test_to_model(pmod, test);
    }
}

void do_resid_freq (gpointer data, guint action, GtkWidget *widget)
{
    FreqDist *freq;
    PRN *prn;
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***rZ;
    DATAINFO *rinfo;

    if (bufopen(&prn)) return;
    
    if (pmod->dataset != NULL) {
	rZ = &pmod->dataset->Z;
	rinfo = pmod->dataset->dinfo;
    } else {
	rZ = &Z;
	rinfo = datainfo;
    }

    if (genr_fit_resid(pmod, rZ, rinfo, GENR_RESID, 1)) {
	nomem();
	return;
    }

    freq = get_freq(rinfo->v - 1, (const double **) *rZ, rinfo, 
		    pmod->ncoeff, OPT_NONE);

    dataset_drop_last_variables(1, rZ, rinfo);

    if (freq_error(freq, NULL)) {
	gretl_print_destroy(prn);
	return;
    }
    
    normal_test(pmod, freq);
    update_model_tests(vwin);

    gretl_command_strcpy("testuhat");

    if (model_command_init(pmod->ID)) {
	return;
    }
 
    print_freq(freq, prn);

    view_buffer(prn, 78, 300, _("gretl: residual dist."), TESTUHAT,
		NULL);

    /* show the graph too */
    if (plot_freq(freq, DIST_NORMAL) == 0) {
	register_graph();
    }

    free_freq(freq);
}

static int 
series_has_negative_vals (const double *x)
{
    int t;

    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (x[t] < 0.0) {
	    return 1;
	}
    }

    return 0;
}

void do_freqplot (gpointer data, guint dist, GtkWidget *widget)
{
    FreqDist *freq;
    gretlopt opt = (dist == DIST_GAMMA)? OPT_O : OPT_NONE;
    int v = mdata_active_var();

    gretl_command_sprintf("freq %s%s", datainfo->varname[v],
			  (dist == DIST_GAMMA)? " --gamma" : "");

    if (check_and_record_command()) {
	return;
    }

    freq = get_freq(v, (const double **) Z, datainfo, 1, opt);

    if (!freq_error(freq, NULL)) { 
	if (dist == DIST_GAMMA && series_has_negative_vals(Z[v])) {
	    errbox(_("Data contain negative values: gamma distribution not "
		   "appropriate"));
	} else {
	    if (plot_freq(freq, dist)) {
		errbox(_("gnuplot command failed"));
	    } else {
		register_graph();
	    }
	}
	free_freq(freq);
    }
}

#if defined(HAVE_TRAMO) || defined (HAVE_X12A)

void do_tramo_x12a (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    int v = mdata_active_var();
    int graph = 0, oldv = datainfo->v;
    gchar *databuf;
    void *handle;
    int (*write_tx_data) (char *, int, 
			  double ***, DATAINFO *, int *,
			  const char *, const char *, char *);
    PRN *prn;
    char fname[MAXLEN] = {0};
    char errtext[MAXLEN];
    char *prog = NULL, *workdir = NULL;

    if (opt == TRAMO) {
#ifdef HAVE_TRAMO
	prog = tramo;
	workdir = tramodir;
#else
	return;
#endif
    } else {
#ifdef HAVE_X12A
	prog = paths.x12a;
	workdir = paths.x12adir;
#else
	return;
#endif
    }

    if (var_is_scalar(datainfo, v)) {
	errbox(_("Can't do this analysis on a scalar"));
	return;
    }

    if (opt != TRAMO) {
	/* we'll let tramo handle annual data */
	if (datainfo->pd == 1 || !dataset_is_time_series(datainfo)) {
	    errbox(_("This analysis is applicable only to seasonal time series"));
	    return;
	}
    }

    write_tx_data = gui_get_plugin_function("write_tx_data", 
					    &handle);
    if (write_tx_data == NULL) {
	return;
    }

    *errtext = 0;

    err = write_tx_data (fname, v, &Z, datainfo, 
			 &graph, prog, workdir, errtext);
    
    close_plugin(handle);

    if (err) {
	if (*errtext != 0) {
	    errbox(errtext);
	} else {
	    errbox((opt == TRAMO)? _("TRAMO command failed") : 
		   _("X-12-ARIMA command failed"));
	}
	return;
    } else if (*fname == '\0') {
	return;
    }


    g_file_get_contents(fname, &databuf, NULL, NULL);
    if (databuf == NULL) {
	errbox((opt == TRAMO)? _("TRAMO command failed") : 
	       _("X-12-ARIMA command failed"));
	return;
    }

    prn = gretl_print_new_with_buffer(databuf);

    view_buffer(prn, (opt == TRAMO)? 106 : 84, 500, 
		(opt == TRAMO)? _("gretl: TRAMO analysis") :
		_("gretl: X-12-ARIMA analysis"),
		opt, NULL);

    if (graph) {
	make_and_display_graph();
    }

    if (datainfo->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }

}

#endif /* HAVE_TRAMO || HAVE_X12A */

void do_range_mean (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    int v = mdata_active_var();
    void *handle;
    int (*range_mean_graph) (int, const double **, 
			     const DATAINFO *, PRN *);
    PRN *prn;

    if (reject_scalar(v)) {
	return;
    }

    range_mean_graph = gui_get_plugin_function("range_mean_graph", 
					       &handle);
    if (range_mean_graph == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = range_mean_graph(v, (const double **) Z, 
			   datainfo, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: range-mean statistics"), RMPLOT, 
		NULL);
}

void do_hurst (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    int v = mdata_active_var();
    void *handle;
    int (*hurst_exponent) (int, const double **, 
			   const DATAINFO *, PRN *);
    PRN *prn;

    if (reject_scalar(v)) {
	return;
    }

    hurst_exponent = gui_get_plugin_function("hurst_exponent", 
					     &handle);
    if (hurst_exponent == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = hurst_exponent(v, (const double **) Z,
			 datainfo, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: Hurst exponent"), HURST, 
		NULL);
}

enum {
    SELECTED_VAR,
    MODEL_VAR
};

static void real_do_corrgm (double ***pZ, DATAINFO *pdinfo, int code)
{
    char title[64];
    int order, err = 0;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    PRN *prn;

    strcpy(title, "gretl: ");
    strcat(title, _("correlogram"));

    order = auto_acf_order(pdinfo->pd, T);

    err = spin_dialog(title, &order, _("Maximum lag:"),
		      1, T - 1, CORRGM);
    if (err < 0) {
	return;
    }    

    if (bufopen(&prn)) return;

    if (code == SELECTED_VAR) {
	gretl_command_sprintf("corrgm %s %d", selected_varname(), order);
	if (check_and_record_command()) {
	    gretl_print_destroy(prn);
	    return;
	}
	err = corrgram(cmd.list[1], order, 0, pZ, pdinfo, prn, OPT_NONE);
    } else {
	err = corrgram(pdinfo->v - 1, order, 0, pZ, pdinfo, prn, OPT_R);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    }

    register_graph();

    view_buffer(prn, 78, 360, title, CORRGM, NULL);
}

void do_corrgm (gpointer data, guint u, GtkWidget *widget)
{
    real_do_corrgm(&Z, datainfo, SELECTED_VAR);
}

void residual_correlogram (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv;
    double ***gZ;
    DATAINFO *ginfo;

    origv = (pmod->dataset != NULL)? 
	pmod->dataset->dinfo->v : datainfo->v;

    /* add residuals to data set temporarily */
    if (add_fit_resid(pmod, GENR_RESID, 1)) return;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }    

    real_do_corrgm(gZ, ginfo, MODEL_VAR);

    dataset_drop_last_variables(ginfo->v - origv, gZ, ginfo);    
}

static void 
real_do_pergm (guint bartlett, double ***pZ, DATAINFO *pdinfo, int code)
{
    gint err;
    PRN *prn;

    if (bufopen(&prn)) return;

    if (code == SELECTED_VAR) {
	if (bartlett) {
	    gretl_command_sprintf("pergm %s --bartlett", selected_varname());
	} else {
	    gretl_command_sprintf("pergm %s", selected_varname());
	}
	if (check_and_record_command()) {
	    gretl_print_destroy(prn);
	    return;
	}
	err = periodogram(cmd.list[1], pZ, pdinfo, cmd.opt, prn);
    } else {
	gretlopt opt = OPT_R;
	if (bartlett) {
	    opt |= OPT_O;
	}
	err = periodogram(pdinfo->v - 1, pZ, pdinfo, opt, prn);
    }

    if (err) {
	gretl_errmsg_set(_("Periodogram command failed"));
	gui_errmsg(1);
	gretl_print_destroy(prn);
	return;
    }

    register_graph();

    view_buffer(prn, 60, 400, _("gretl: periodogram"), PERGM, 
		NULL);
}

void do_pergm (gpointer data, guint opt, GtkWidget *widget)
{
    real_do_pergm(opt, &Z, datainfo, SELECTED_VAR);
}

void residual_periodogram (gpointer data, guint opt, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv;
    double ***gZ;
    DATAINFO *ginfo;

    origv = (pmod->dataset != NULL)? 
	pmod->dataset->dinfo->v : datainfo->v;

    /* add residuals to data set temporarily */
    if (add_fit_resid(pmod, GENR_RESID, 1)) return;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }    

    real_do_pergm(1, gZ, ginfo, MODEL_VAR);

    dataset_drop_last_variables(ginfo->v - origv, gZ, ginfo); 
}

void do_coeff_intervals (gpointer data, guint u, GtkWidget *w)
{
    PRN *prn;
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    CoeffIntervals *cf;

    if (bufopen(&prn)) return;

    cf = gretl_model_get_coeff_intervals(pmod, datainfo);

    if (cf != NULL) {
	text_print_model_confints(cf, prn);
	view_buffer(prn, 78, 300, 
		    _("gretl: coefficient confidence intervals"), 
		    COEFFINT, cf);
    }
}

void do_outcovmx (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    VMatrix *vcv = NULL;

    if (Z == NULL || datainfo == NULL) {
	errbox(_("Data set is gone"));
	return;
    }

    if (bufopen(&prn)) return;

    vcv = gretl_model_get_vcv(pmod, datainfo);

    if (vcv == NULL) {
	errbox(_("Error generating covariance matrix"));
    } else {
	text_print_vmatrix(vcv, prn);
	view_buffer(prn, 80, 300, _("gretl: coefficient covariances"), 
		    COVAR, vcv);
    }
}

void add_dummies (gpointer data, guint u, GtkWidget *widget)
{
    gretlopt opt = OPT_NONE;
    gint err;

    if (u == TS_DUMMIES) {
	gretl_command_strcpy("genr dummy");
    } else if (dataset_is_panel(datainfo)) {
	if (u == PANEL_UNIT_DUMMIES) {
	    gretl_command_strcpy("genr unitdum");
	} else {
	    gretl_command_strcpy("genr timedum");
	    opt = OPT_T;
	}
    } else {
	errbox(_("Data set is not recognized as a panel.\n"
		 "Please use \"Sample/Set frequency, startobs\"."));
	return;
    }

    if (check_and_record_command()) {
	return;
    }

    if (u == TS_DUMMIES) {
	err = dummy(&Z, datainfo, 0) == 0;
    } else {
	err = panel_dummies(&Z, datainfo, opt);
    } 

    if (err) {
	gui_errmsg(err);
    } else {
	populate_varlist();
    }
}

void add_index (gpointer data, guint tm, GtkWidget *widget)
{
    gretl_command_strcpy((tm)? "genr time" : "genr index");

    if (check_and_record_command()) {
	return;
    }

    if (genrtime(&Z, datainfo, tm)) {
	errbox((tm)? _("Error generating time trend") :
	       _("Error generating index variable"));
    } else {
	populate_varlist();
    }
}

void do_add_obs (gpointer data, guint u, GtkWidget *widget)
{
    int n = add_obs_dialog(NULL, 1);
    int err = 0;

    if (n > 0) {
	err = dataset_add_observations(n, &Z, datainfo, OPT_A);
	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_dataset_as_modified();
	}
    }
}

void do_remove_obs (gpointer data, guint u, GtkWidget *widget)
{
    int drop = 0;

    if (complex_subsampled()) {
	errbox(_("The data set is currently sub-sampled.\n"));
	drop_obs_state(FALSE);
    } else {
	drop = datainfo->n - get_original_n();
    }

    if (drop > 0) {
	gchar *msg;
	int resp;

	msg = g_strdup_printf(_("Really delete the last %d observations?"),
			      drop);
	resp = yes_no_dialog(_("gretl: drop observations"), msg, 0);
	g_free(msg);

	if (resp == GRETL_YES) {
	    int err = dataset_drop_observations(drop, &Z, datainfo);

	    if (err) {
		gui_errmsg(err);
	    } else {
		mark_dataset_as_modified();
	    }
	    drop_obs_state(FALSE);
	}
    } else {
	errbox(_("There are no extra observations to drop"));
	drop_obs_state(FALSE);
    }
}

void add_logs_etc (gpointer data, guint action, GtkWidget *widget)
{
    char *liststr;
    int order = 0;
    int err = 0;

    liststr = main_window_selection_as_string();
    if (liststr == NULL) {
	return;
    }

    if (action == LAGS) {
	int resp;

	order = default_lag_order(datainfo);

	resp = spin_dialog(_("gretl: generate lags"), 
			   &order, _("Number of lags to create:"), 
			   1, datainfo->n - 1, 0);
	if (resp < 0) {
	    free(liststr);
	    return;
	}

	if (order > 0) {
	    gretl_command_sprintf("lags %d ;%s", order, liststr);
	} else {
	    gretl_command_sprintf("lags%s", liststr);
	}
    } else {
	gretl_command_sprintf("%s%s", gretl_command_word(action), liststr);
    }

    free(liststr);

    if (check_and_record_command()) {
	return;
    }

    if (action == LAGS) {
	err = list_laggenr(&cmd.list, order, &Z, datainfo);
    } else if (action == LOGS) {
	err = list_loggenr(cmd.list, &Z, datainfo);
    } else if (action == SQUARE) {
	err = list_xpxgenr(&cmd.list, &Z, datainfo, OPT_NONE);
    } else if (action == DIFF || action == LDIFF || action == SDIFF) {
	err = list_diffgenr(cmd.list, action, &Z, datainfo);
    } else if (action == DUMMIFY) {
	int i;

	for (i=1; i<=cmd.list[0]; i++) {
	    if (!var_is_discrete(datainfo, cmd.list[i])) {
		err++; 
	    }
	}
	if (err < cmd.list[0]) {
	    err = list_dumgenr(&cmd.list, &Z, datainfo);
	} else {
	    errbox(_("No discrete variables were selected"));
	    return;
	}
    }

    if (err) {
	errbox(_("Error adding variables"));
    } else {
	populate_varlist();
    }
}

/* 
   add_fit_resid: if undo = 1, don't bother with the label, don't
   update the var display in the main window, and don't add to command
   log.
*/

int add_fit_resid (MODEL *pmod, int code, int undo)
{
    int err;

    if (pmod->dataset != NULL) {
	if (!undo) {
	    return 1;
	} 
	err = genr_fit_resid(pmod, 
			     &pmod->dataset->Z, 
			     pmod->dataset->dinfo, 
			     code, undo);
    } else {
	err = genr_fit_resid(pmod, &Z, datainfo, code, undo);
    }

    if (err) {
	nomem();
	return 1;
    }

    if (!undo) {
	int v = datainfo->v - 1;

	/* give the user a chance to choose a different name */
	varinfo_dialog(v, 0);

	if (*datainfo->varname[v] == '\0') {
	    /* the user canceled */
	    dataset_drop_last_variables(1, &Z, datainfo);
	    return 0;
	}	

	populate_varlist();

	if (code == GENR_RESID) {
	    gretl_command_sprintf("genr %s = $uhat", datainfo->varname[v]);
	} else if (code == GENR_FITTED) {
	    gretl_command_sprintf("genr %s = $yhat", datainfo->varname[v]);
	} else if (code == GENR_RESID2) {
	    gretl_command_sprintf("genr %s = $uhat*$uhat", datainfo->varname[v]);
	} else if (code == GENR_H) {
	    gretl_command_sprintf("genr %s = $h", datainfo->varname[v]);
	} else if (code == GENR_AHAT) {
	    gretl_command_sprintf("genr %s = $ahat", datainfo->varname[v]);
	}

	model_command_init(pmod->ID);
	mark_dataset_as_modified();
    }

    return 0;
}

int add_system_resid (gpointer data, int eqnum, int ci)
{
    windata_t *vwin = (windata_t *) data;
    int err, v;

    if (ci == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;

	err = gretl_VAR_add_resids_to_dataset(var, eqnum,
					      &Z, datainfo);
    } else {
	const char *sysname = (const char *) vwin->data;

	err = gretl_system_add_resids_to_dataset(sysname, eqnum,
						 &Z, datainfo);
    }	

    if (err) {
	nomem();
	return 1;
    }

    v = datainfo->v - 1;

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return 0;
    }    

    populate_varlist();
    mark_dataset_as_modified();

    return 0;
}

void add_model_stat (MODEL *pmod, int which)
{
    char vname[VNAMELEN], vlabel[MAXLABEL];
    char statname[8];
    int i, n;

    if (dataset_add_scalar(&Z, datainfo)) {
	nomem();
	return;
    }

    i = datainfo->v - 1;
    n = datainfo->n;

    switch (which) {
    case ESS:
	sprintf(vname, "ess_%d", pmod->ID);
	sprintf(vlabel, _("error sum of squares from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->ess;
	strcpy(statname, "$ess");
	break;
    case R2:
	sprintf(vname, "r2_%d", pmod->ID);
	sprintf(vlabel, _("R-squared from model %d"), pmod->ID);
	Z[i][0] = pmod->rsq;
	strcpy(statname, "$rsq");
	break;
    case TR2:
	sprintf(vname, "trsq%d", pmod->ID);
	sprintf(vlabel, _("T*R-squared from model %d"), pmod->ID);
	Z[i][0] = pmod->nobs * pmod->rsq;
	strcpy(statname, "$trsq");
	break;
    case DF:
	sprintf(vname, "df_%d", pmod->ID);
	sprintf(vlabel, _("degrees of freedom from model %d"), 
		pmod->ID);
	Z[i][0] = (double) pmod->dfd;
	strcpy(statname, "$df");
	break;
    case SIGMA:
	sprintf(vname, "sgma_%d", pmod->ID);
	sprintf(vlabel, _("std err of residuals from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->sigma;
	strcpy(statname, "$sigma");
	break;
    case LNL:
	sprintf(vname, "lnl_%d", pmod->ID);
	sprintf(vlabel, _("log likelihood from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->lnL;
	strcpy(statname, "$lnl");
	break;	
    case AIC:
	sprintf(vname, "aic_%d", pmod->ID);
	sprintf(vlabel, _("Akaike Information Criterion from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->criterion[C_AIC];
	strcpy(statname, "$aic");
	break;
    case BIC:
	sprintf(vname, "bic_%d", pmod->ID);
	sprintf(vlabel, _("Bayesian Information Criterion from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->criterion[C_BIC];
	strcpy(statname, "$bic");
	break;
    case HQC:
	sprintf(vname, "hqc_%d", pmod->ID);
	sprintf(vlabel, _("Hannan-Quinn Information Criterion from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->criterion[C_HQC];
	strcpy(statname, "$hqc");
	break;
    }

    strcpy(datainfo->varname[i], make_varname_unique(vname, i, datainfo));
    strcpy(VARLABEL(datainfo, i), vlabel);

    /* give the user a chance to choose a different name */
    varinfo_dialog(i, 0);

    if (*datainfo->varname[i] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return;
    }

    gretl_command_sprintf("genr %s = %s", datainfo->varname[i], statname);

    populate_varlist();
    model_command_init(pmod->ID);

    /* note: since this is a scalar, which will not be saved by
       default on File/Save data, we will not mark the data set
       as "modified" here. (FIXME saving scalars?) */
}

void resid_plot (gpointer data, guint xvar, GtkWidget *widget)
{
    GnuplotFlags flags = 0;
    int err, origv, ts, plotlist[4], lines[1];
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    int pdum = vwin->active_var; 
    int yno, uhatno;
    double ***gZ;
    DATAINFO *ginfo;

    /* special case: GARCH model (show fitted variance) */
    if (pmod->ci == GARCH && xvar == 0) {
	err = garch_resid_plot(pmod, datainfo);
	if (err) {
	    errbox(_("gnuplot command failed"));
	} else {
	    register_graph();
	}
	return;
    }

    origv = (pmod->dataset != NULL)? 
	pmod->dataset->dinfo->v : datainfo->v;

    /* add residuals to data set temporarily */
    if (add_fit_resid(pmod, GENR_RESID, 1)) {
	return;
    }

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }    

    flags = GP_GUI | GP_RESIDS;
    if (pdum) {
	flags |= GP_DUMMY;
    }

    ts = dataset_is_time_series(ginfo);
    uhatno = ginfo->v - 1; /* residual: last var added */

    plotlist[0] = 1;
    plotlist[1] = uhatno; 

    strcpy(ginfo->varname[uhatno], _("residual"));
    yno = gretl_model_get_depvar(pmod);
    sprintf(VARLABEL(ginfo, uhatno), "residual for %s", 
	    ginfo->varname[yno]);

    if (xvar) { 
	/* plot against specified xvar */
	plotlist[0] = 2;
	plotlist[2] = xvar;
	lines[0] = 0;
    } else {    
	/* plot against obs index or time */
	flags |= GP_IDX;
	lines[0] = (ts)? 1 : 0;
    } 

    /* plot separated by dummy variable? */
    if (pdum) {
	plotlist[0] += 1;
	plotlist[plotlist[0]] = pdum;
    }

    /* generate graph */
    err = gnuplot(plotlist, lines, NULL, gZ, ginfo, &plot_count, flags);

    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
    
    dataset_drop_last_variables(ginfo->v - origv, gZ, ginfo);
}

void fit_actual_plot (gpointer data, guint xvar, GtkWidget *widget)
{
    GnuplotFlags flags = GP_GUI | GP_FA;
    int err, origv, plotlist[4], lines[2] = {0};
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***gZ;
    DATAINFO *ginfo;
    char *formula;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }

    formula = gretl_model_get_fitted_formula(pmod, xvar, (const double **) *gZ,
					     ginfo);

    if (formula != NULL) {
	/* fitted value can be represented as a formula: if feasible,
	   produces a better-looking graph */
	plotlist[0] = 3;
	plotlist[1] = 0; /* placeholder entry */
	plotlist[2] = gretl_model_get_depvar(pmod);
	plotlist[3] = xvar;
	err = gnuplot(plotlist, lines, formula, gZ, ginfo,
		      &plot_count, flags);
	if (err) {
	    errbox(_("gnuplot command failed"));
	} else {
	    register_graph();
	}
	free(formula);
	return;
    }

    origv = ginfo->v;

    /* add fitted values to data set temporarily */
    if (add_fit_resid(pmod, GENR_FITTED, 1)) {
	return;
    }

    plotlist[0] = 3;
    plotlist[1] = ginfo->v - 1; /* last var added (fitted vals) */

    /* depvar from regression */
    plotlist[2] = gretl_model_get_depvar(pmod);

    if (xvar) { 
	/* plot against specified xvar */
	plotlist[3] = xvar;
	/* is it a simple regression? */
	if ((pmod->ifc && pmod->list[0] == 3) || pmod->list[0] == 2) {
	    lines[0] = 1;
	} else {
	    lines[0] = 0;
	}
	lines[1] = 0;
    } else { 
	/* plot against obs */
	int ts = dataset_is_time_series(ginfo);

	plotlist[0] -= 1;
	flags |= GP_IDX;
	lines[0] = (ts)? 1 : 0; 
	lines[1] = (ts)? 1 : 0;
    } 

    err = gnuplot(plotlist, lines, NULL, gZ, ginfo,
		  &plot_count, flags);

    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }

    dataset_drop_last_variables(ginfo->v - origv, gZ, ginfo);
}

void fit_actual_splot (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***gZ;
    DATAINFO *ginfo;
    int list[4];
    int err;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }    

    /* Y, X, Z */

    list[0] = 3;
    list[1] = pmod->list[4];
    list[2] = pmod->list[3];
    list[3] = pmod->list[1];

    err = gnuplot_3d(list, NULL, gZ, ginfo,
		     &plot_count, GP_GUI | GP_FA);

    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	launch_gnuplot_interactive();
    }
}

#define MAXDISPLAY 4096
/* max number of observations for which we expect to be able to 
   use the buffer approach for displaying data, as opposed to
   disk file */

void display_selected (gpointer data, guint action, GtkWidget *widget)
{
    char *liststr; 
    PRN *prn = NULL;
    int *list = NULL;
    int n = datainfo->t2 - datainfo->t1 + 1;

    liststr = main_window_selection_as_string();
    if (liststr == NULL || *liststr == '\0') {
	return;
    }

    list = gretl_list_from_string(liststr);
    free(liststr);

    if (list == NULL) {
	return;
    }

    /* special case: showing only one series */
    if (list[0] == 1) {
	display_var();
	free(list);
	return;
    }

    if (list[0] * n > MAXDISPLAY) { 
	/* use disk file */
	char fname[MAXLEN];

	if (user_fopen("data_display_tmp", fname, &prn)) {
	    return;
	}

	printdata(list, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	free(list);
	view_file(fname, 0, 1, 78, 350, VIEW_DATA);
    } else { 
	/* use buffer */
	multi_series_view *mview = NULL;
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	err = printdata(list, (const double **) Z, datainfo, OPT_O, prn);
	if (err) {
	    nomem();
	    gretl_print_destroy(prn);
	    free(list);
	    return;
	}
	if (get_printdata_blocks() == 1) {
	    mview = multi_series_view_new(list);
	}
	view_buffer(prn, 78, 350, _("gretl: display data"), PRINT, mview);
    }
}

void display_fit_resid (gpointer data, guint code, GtkWidget *widget)
{
    PRN *prn;
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    FITRESID *fr;

    if (bufopen(&prn)) return;

    fr = get_fit_resid(pmod, (const double **) Z, datainfo);

    if (fr == NULL) {
	errbox(_("Failed to generate fitted values"));
	gretl_print_destroy(prn);
    } else {
	text_print_fit_resid(fr, datainfo, prn);
	view_buffer(prn, 78, 350, _("gretl: display data"), FCAST, fr);  
    }  
}

/* Before deleting specified variables, check that they are not
   required by any saved models; also, don't delete variables 
   whose deletion would result in the renumbering of variables
   used in saved models.
*/

static int maybe_prune_delete_list (int *list)
{
    int vsave = 0, pruned = 0;
    int i, vmax;

    /* check open model windows */
    vmax = highest_numbered_variable_in_winstack();
    if (vmax > vsave) {
	vsave = vmax;
    }

    /* check models saved as icons */
    vmax = highest_numbered_variable_in_session();
    if (vmax > vsave) {
	vsave = vmax;
    }
    
    for (i=1; i<=list[0]; i++) {
	if (list[i] <= vsave) {
	    gretl_list_delete_at_pos(list, i);
	    i--;
	    pruned = 1;
	}
    }

    return pruned;
}

void delete_selected_vars (int id)
{
    int err, renumber, pruned = 0;
    char *liststr = NULL;
    char *msg;

    if (dataset_locked()) {
	return;
    }

    if (complex_subsampled()) {
	errbox(_("Can't delete a variable when in sub-sample"
		 " mode\n"));
	return;
    }

    if (id > 0) {
	/* delete single specified var */
	int testlist[2];

	testlist[0] = 1;
	testlist[1] = id;

	if (maybe_prune_delete_list(testlist)) {
	    errbox(_("Cannot delete %s; variable is in use"), 
		   datainfo->varname[id]);
	    return;
	} else {
	    msg = g_strdup_printf(_("Really delete %s?"), datainfo->varname[id]);
	}
    } else {
	/* delete list of vars */
	liststr = main_window_selection_as_string();
	if (liststr == NULL) {
	    return;
	}
	msg = g_strdup_printf(_("Really delete %s?"), liststr);
    }

    if (yes_no_dialog(_("gretl: delete"), msg, 0) != GRETL_YES) {
	g_free(msg);
	if (liststr != NULL) {
	    free(liststr);
	}
	return;
    }

    g_free(msg);

    if (id > 0) {
	gretl_command_sprintf("delete %d", id);
    } else {
	gretl_command_sprintf("delete%s", liststr);
	free(liststr);  
    } 

    if (check_and_record_command()) {
	return;
    }

    if (id == 0) {
	pruned = maybe_prune_delete_list(cmd.list);
    }

    if (cmd.list[0] == 0) {
	errbox(_("Cannot delete the specified variables"));
	return;
    } else if (pruned) {
	errbox(_("Cannot delete all of the specified variables"));
    }

    err = dataset_drop_listed_variables(cmd.list, &Z, datainfo, &renumber);

    if (err) {
	nomem();
    } else {
	refresh_data();
	if (renumber) {
	    infobox(_("Take note: variables have been renumbered"));
	}
	maybe_clear_selector(cmd.list);
	mark_dataset_as_modified();
    }
}

static void do_stacked_ts_plot (int varnum)
{
    int list[2];
    int err;
    
    list[0] = 1;
    list[1] = varnum;

    err = gretl_panel_ts_plot(list, (const double **) Z, datainfo);

    gui_graph_handler(err);
}

void do_graph_var (int varnum)
{
    int err, lines[1];

    if (varnum <= 0) return;

    if (datainfo->structure == STACKED_TIME_SERIES &&
	datainfo->n / datainfo->pd < 10 &&
	balanced_panel(datainfo)) {
	do_stacked_ts_plot(varnum);
	return;
    }

    if (!dataset_is_time_series(datainfo) &&
	datainfo->structure != STACKED_TIME_SERIES) {
	do_freqplot(NULL, 0, NULL);
	return;
    }

    gretl_command_sprintf("gnuplot %s --time-series", 
			  datainfo->varname[varnum]);

    if (check_and_record_command()) {
	return;
    }

    lines[0] = 1;
    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, GP_GUI | GP_IDX);

    gui_graph_handler(err);
}

void ts_plot_var (gpointer data, guint opt, GtkWidget *widget)
{
    do_graph_var(mdata_active_var());
}

void do_boxplot_var (int varnum)
{
    if (varnum < 0) return;

    gretl_command_sprintf("boxplot %s", datainfo->varname[varnum]);

    if (check_and_record_command()) {
	return;
    }

    if (boxplots(cmd.list, NULL, &Z, datainfo, 0)) {
	errbox (_("boxplot command failed"));
    }
}

int do_scatters (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    GnuplotFlags flags = 0;
    int err; 

    if (buf == NULL) return 1;

    if (opt & OPT_L) {
	gretl_command_sprintf("scatters %s --with-lines", buf);
	flags |= GP_LINES;
    } else {
	gretl_command_sprintf("scatters %s", buf);
    }

    if (check_and_record_command()) {
	return 1;
    }

    err = multi_scatters(cmd.list, (const double **) Z, datainfo, NULL, flags);

    if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }

    return 0;
}

void do_box_graph (GtkWidget *widget, dialog_t *dlg)
{
    int action = edit_dialog_get_action(dlg);
    const char *buf = edit_dialog_get_text(dlg);
    int err;

    if (buf == NULL) return;

    if (strchr(buf, '(')) {
	err = boolean_boxplots(buf, &Z, datainfo, (action == GR_NBOX));
    } else {
	gretl_command_sprintf("boxplot %s%s", 
			      (action == GR_NBOX)? "--notches " : "", buf);

	if (check_and_record_command()) {
	    return;
	}
	err = boxplots(cmd.list, NULL, &Z, datainfo, (action == GR_NBOX));
    }

    if (err) {
	errbox(_("boxplot command failed"));
    } else {
	close_dialog(dlg);
    }
}

/* X, Y scatter with separation by dummy (factor) */

int do_dummy_graph (selector *sr)
{
    const char *buf = selector_list(sr);
    gint err, lines[1] = {0}; 

    if (buf == NULL) return 1;

    gretl_command_sprintf("gnuplot %s --dummy", buf);

    if (check_and_record_command()) {
	return 1;
    }

    if (cmd.list[0] != 3 || 
	!gretl_isdummy(datainfo->t1, datainfo->t2, Z[cmd.list[3]])) {
	errbox(_("You must supply three variables, the last\nof which "
	       "is a dummy variable (values 1 or 0)"));
	return 1;
    }

    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, GP_GUI | GP_DUMMY);

    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }

    return 0;
}

int do_graph_from_selector (selector *sr)
{
    GnuplotFlags flags = GP_GUI;
    const char *buf = selector_list(sr);
    gint i, err, *lines = NULL;
    gint imp = (selector_code(sr) == GR_IMP);

    if (buf == NULL) return 1;

    gretl_command_sprintf("gnuplot %s%s", buf, 
			 (imp)? " --with-impulses" : "");

    if (selector_code(sr) == GR_PLOT) { 
        gretl_command_strcat(" --time-series");
	flags |= GP_IDX;
    }

    if (check_and_record_command()) {
	return 1;
    }

    if (imp) {
	flags |= GP_IMPULSES;
    } else {
	lines = mymalloc((cmd.list[0] - 1) * sizeof *lines);
	if (lines == NULL) {
	    return 0;
	}
	for (i=0; i<cmd.list[0]-1 ; i++) {
	    if (selector_code(sr) == GR_PLOT) {
		lines[i] = 1;
	    } else {
		lines[i] = 0;
	    }
	}
    }

    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, flags);

    gui_graph_handler(err);

    if (lines != NULL) {
	free(lines);
    }

    return 0;
}

#ifndef G_OS_WIN32

#include <signal.h>
#include <errno.h>

static int executable_exists (const char *fname)
{
    struct stat sbuf;
    int ok = 0;

    if (stat(fname, &sbuf) == 0) {
	if (sbuf.st_mode & S_IXOTH) {
	    ok = 1;
	} else if (getuid() == sbuf.st_uid &&
		   (sbuf.st_mode & S_IXUSR)) {
	    ok = 1;
	} else if (getgid() == sbuf.st_gid &&
		   (sbuf.st_mode & S_IXGRP)) {
	    ok = 1;
	}
    }

    return ok;
}

static int mywhich (const char *prog)
{
    char test[FILENAME_MAX];
    char *path, *p;
    gchar *mypath;
    int i, ret = 0;

    path = getenv("PATH");
    if (path == NULL) return 0;

    mypath = g_strdup(path);
    if (mypath == NULL) return 0;

    for (i=0; ; i++) {
	p = strtok((i)? NULL : mypath, ":");
	if (p == NULL) break;
	sprintf(test, "%s/%s", p, prog);
	if (executable_exists(test)) {
	    ret = 1;
	    break;
	}
    }

    g_free(mypath);

    return ret;
}

static int get_terminal (char *s)
{
    if (mywhich("xterm")) {
	strcpy(s, "xterm");
	return 0;
    }

    if (mywhich("rxvt")) {
	strcpy(s, "rxvt");
	return 0;
    }

    errbox(_("Couldn't find xterm or rxvt"));
    return 1;
}

#endif /* not G_OS_WIN32 */

int do_splot_from_selector (selector *sr)
{
    char line[MAXLINE];
    const char *buf = selector_list(sr);
    int err;

    if (buf == NULL) return 1;

    sprintf(line, "gnuplot %s", buf);
    if (check_specific_command(line) || cmd.list[0] != 3) {
	return 1;
    }

    err = gnuplot_3d(cmd.list, NULL, &Z, datainfo,
		     &plot_count, GP_GUI);

    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	launch_gnuplot_interactive();
    }

    return 0;
}

static int list_position (int v, const int *list)
{
    int i;

    for (i=list[0]; i>=1; i--) {
	if (v == list[i]) {
	    return i;
	}
    }

    return 0;
}

static int maybe_reorder_list (char *liststr)
{
    const char *query = _("X-axis variable");
    int *list = gretl_list_from_string(liststr);

    if (list == NULL) {
	return 1;
    } else {
	int xvar = select_var_from_list(list, query);

	if (xvar < 0) {
	    /* the user cancelled */
	    return 1;
	}

	if (xvar != list[list[0]]) {
	    int tmp = list[list[0]];
	    int pos = list_position(xvar, list);
	    int i;

	    list[list[0]] = xvar;
	    list[pos] = tmp;
	    *liststr = '\0';
	    for (i=1; i<=list[0]; i++) {
		char numstr[8];

		sprintf(numstr, " %d", list[i]);
		strcat(liststr, numstr);
	    }
	}
	free(list);
    }

    return 0;
}

void plot_from_selection (gpointer data, guint action, GtkWidget *widget)
{
    GnuplotFlags flags = GP_GUI;
    char *liststr;
    int *lines = NULL;
    int llen = 1;
    gint i, err;

    liststr = main_window_selection_as_string();
    if (liststr == NULL || *liststr == 0) {
	return;
    }

    if (action == GR_XY) {
	err = maybe_reorder_list(liststr);
	if (err) return;
    } else if (action == GR_PLOT) {
	flags |= GP_IDX;
    }

    gretl_command_sprintf("gnuplot%s%s", liststr, 
			  (action == GR_PLOT)? " --time-series" : "");
    free(liststr);

    if (check_and_record_command()) {
	return;
    }

    if (cmd.list[0] - 1 > llen) {
	llen = cmd.list[0] - 1;
    }

    lines = mymalloc(llen * sizeof *lines);
    if (lines == NULL) {
	return;
    }

    for (i=0; i<llen; i++) {
	lines[i] = (action == GR_PLOT);
    }

    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, flags);

    gui_graph_handler(err);

    free(lines);
}

void display_var (void)
{
    int list[2];
    PRN *prn;
    windata_t *vwin;
    int height = 400;
    int vec = 1;
    int n = datainfo->t2 - datainfo->t1 + 1;

    list[0] = 1;
    list[1] = mdata_active_var();

    if (var_is_scalar(datainfo, list[1])) {
	vec = 0;
	height = 140;
    }

    if (n > MAXDISPLAY) { 
	/* use disk file */
	char fname[MAXLEN];

	if (user_fopen("data_display_tmp", fname, &prn)) {
	    return;
	}

	printdata(list, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 28, height, VIEW_DATA);
    } else { 
	/* use buffer */
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	err = printdata(list, (const double **) Z, datainfo, OPT_O, prn);
	if (err) {
	    nomem();
	    gretl_print_destroy(prn);
	    return;
	}
	vwin = view_buffer(prn, 36, height, 
			   datainfo->varname[list[1]], 
			   (vec)? VIEW_SERIES : VIEW_SCALAR, 
			   NULL);
	series_view_connect(vwin, list[1]);
    }
}

#define PGRAB
#undef SCRIPT_TO_FILE

void do_run_script (gpointer data, guint code, GtkWidget *w)
{
    PRN *prn;
    char *runfile = NULL;
#ifdef SCRIPT_TO_FILE
    char fname[MAXLEN];
#endif
    int err;

#ifdef SCRIPT_TO_FILE
    if (user_fopen("gretl_output_tmp", fname, &prn)) return;
#else
    if (bufopen(&prn)) {
	return ;
    }
#endif

    if (code == SCRIPT_EXEC) {
	runfile = scriptfile;
	gretl_command_sprintf("run %s", scriptfile);
	check_and_record_command();
    } else if (code == SESSION_EXEC) {
	runfile = cmdfile;
    }

    if (data != NULL) { 
	/* get commands from file view buffer */
#ifdef PGRAB
	GdkCursor *plswait;
#endif
	windata_t *vwin = (windata_t *) data;
	gchar *buf;

	buf = textview_get_text(vwin->w);

	if (buf == NULL || *buf == '\0') {
	    errbox("No commands to execute");
	    gretl_print_destroy(prn);
	    if (buf != NULL) {
		g_free(buf);
	    }
	    return;
	}

#ifdef PGRAB
	plswait = gdk_cursor_new(GDK_WATCH);
	gdk_pointer_grab(vwin->dialog->window, TRUE,
			 GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK |
			 GDK_BUTTON_RELEASE_MASK,
			 NULL, plswait,
			 GDK_CURRENT_TIME); 
	gdk_cursor_destroy(plswait);
#endif

	err = execute_script(NULL, buf, prn, code);
	g_free(buf);

#ifdef PGRAB
	gdk_pointer_ungrab(GDK_CURRENT_TIME);
#endif
    } else {
	/* get commands from file */
	err = execute_script(runfile, NULL, prn, code);
    }

#ifdef SCRIPT_TO_FILE
    gretl_print_destroy(prn);
#endif

    if (err == -1) return;

    refresh_data();

#ifdef SCRIPT_TO_FILE
    view_file(fname, 1, 1, 78, 450, SCRIPT_OUT);
#else
    view_buffer(prn, 78, 450, NULL, SCRIPT_OUT, NULL);
#endif

    /* re-establish command echo */
    set_gretl_echo(1);
}

void do_open_script (void)
{
    int n = strlen(paths.scriptdir);
    FILE *fp;

    fp = fopen(tryfile, "r");
    if (fp == NULL) {
	errbox(_("Couldn't open %s"), tryfile);
	delete_from_filelist(FILE_LIST_SESSION, tryfile);
	delete_from_filelist(FILE_LIST_SCRIPT, tryfile);
	return;
    } else {
	fclose(fp);
    }
	
    strcpy(scriptfile, tryfile);

    mkfilelist(FILE_LIST_SCRIPT, scriptfile);

    if (strncmp(scriptfile, paths.scriptdir, n)) { 
	view_file(scriptfile, 1, 0, 78, 370, EDIT_SCRIPT);
    } else {
	view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT);
    }
}

void do_new_script (gpointer data, guint u, GtkWidget *widget) 
{
    char temp[MAXLEN];
    FILE *fp;

    sprintf(temp, "%sscript_tmp", paths.userdir);
    fp = gretl_tempfile_open(temp);
    if (fp == NULL) {
	return;
    }

    fclose(fp);
    strcpy(scriptfile, temp); /* ?? */
    view_file(scriptfile, 1, 1, 78, 370, EDIT_SCRIPT);
}

void maybe_display_string_table (void)
{
    static int s_table_waiting;

    if (gretl_string_table_written() || s_table_waiting) {
	char stname[MAXLEN];

	if (mdata == NULL) {
	    s_table_waiting = 1;
	    return;
	} 

	s_table_waiting = 0;
	build_path(stname, paths.userdir, "string_table.txt", NULL);
	view_file(stname, 0, 0, 78, 350, VIEW_FILE);
    }
}

void do_open_csv_box (char *fname, int code, int append)
{
    int err;
    PRN *prn;
    char buf[32];
    char dtype[8];

    if (bufopen(&prn)) return;

    if (code == OPEN_BOX) {
	err = import_box(&Z, &datainfo, fname, prn);
	strcpy(dtype, "BOX");
    } else if (code == OPEN_OCTAVE) {
	err = import_octave(&Z, &datainfo, fname, prn);
	strcpy(dtype, "Octave");
    } else {
	err = import_csv(&Z, &datainfo, fname, prn); 
	strcpy(dtype, "CSV");
    }

    sprintf(buf, _("gretl: import %s data"), dtype);

    view_buffer(prn, 78, 350, buf, IMPORT, NULL); 

    if (err) return;

    maybe_display_string_table();
    data_status |= IMPORT_DATA;

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	strcpy(paths.datfile, fname);
	register_data(fname, NULL, 1);
    }
}

int dataset_is_subsampled (void)
{
    int ret = 0;

    if (mdata->ifac != NULL) {
	GtkWidget *w = gtk_item_factory_get_item(mdata->ifac, 
						 "/Sample/Restore full range");

	if (w != NULL && GTK_IS_WIDGET(w) && GTK_WIDGET_IS_SENSITIVE(w)) {
	    ret = 1;
	}
    }

    return ret;
}

int dataset_is_restricted (void)
{
    /* Should we indicate "restricted" if t1 and t2 are reset, or only
       if a sub-sampling mask is in place?  For now we'll go with the
       broader option.
    */

    return dataset_is_subsampled();
}

int maybe_restore_full_data (int action)
{
    if (dataset_is_subsampled()) {
	int resp = GRETL_CANCEL;

	if (action == SAVE_DATA) {
	    resp = yes_no_dialog(_("gretl: save data"), 
				 _("The data set is currently sub-sampled.\n"
				   "Would you like to restore the full range?"), 1);
	} else if (action == COMPACT) {
	    resp = yes_no_dialog(_("gretl: Compact data"), 
				 _("The data set is currently sub-sampled.\n"
				   "You must restore the full range before compacting.\n"
				   "Restore the full range now?"), 1);
	} else if (action == EXPAND) {
	    resp = yes_no_dialog(_("gretl: Expand data"), 
				 _("The data set is currently sub-sampled.\n"
				   "You must restore the full range before expanding.\n"
				   "Restore the full range now?"), 1);
	}

	if (resp == GRETL_YES) {
	    gui_restore_sample();
	} else if (resp == GRETL_CANCEL || resp < 0 || 
		   action == COMPACT || action == EXPAND) {
	    return 1;
	}
    } 

    return 0;
}

void gui_transpose_data (gpointer p, guint u, GtkWidget *w)
{
    int i, resp;

    for (i=1; i<datainfo->v; i++) {
	if (var_is_scalar(datainfo, i)) {
	    errbox(_("Dataset contains scalars, can't transpose"));
	    return;
	}
    }

    resp = yes_no_dialog(_("gretl: transpose data"), 
			 _("Transposing means that each variable becomes interpreted\n"
			   "as an observation, and each observation as a variable.\n"
			   "Do you want to proceed?"), 0);

    if (resp == GRETL_YES) {
	int err = transpose_data(&Z, datainfo);
    
	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_dataset_as_modified();
	    populate_varlist();
	    infobox(_("Data transposed"));
	}
    }
}

static int db_write_response (const char *savename, const int *list)
{
    gchar *msg;
    int resp, ret = 0;

    msg = g_strdup_printf("%s\n%s", get_gretl_errmsg(),
			  _("OK to overwrite?"));

    resp = yes_no_dialog("gretl", msg, 0);
    if (resp == GRETL_NO) {
	ret = 1;
    } else {
	ret = write_db_data(savename, list, OPT_F,
			    (const double **) Z, datainfo);
    }

    g_free(msg);  

    return ret;
}

#define DATA_EXPORT(o) (o & (OPT_M | OPT_R | OPT_G | OPT_A | OPT_C | OPT_D))

#define WRITING_DB(o) (o & OPT_D)

int do_store (char *savename, gretlopt opt, int overwrite)
{
    gchar *tmp = NULL;
    FILE *fp;
    int showlist = 1;
    int err = 0;

    /* if the data set is sub-sampled, give a chance to rebuild
       the full data range before saving */
    if (maybe_restore_full_data(SAVE_DATA)) {
	goto store_get_out;
    }

    /* "storelist" is a global */
    if (storelist == NULL) {
	showlist = 0;
    }

    if (opt != OPT_NONE) { 
	/* not a bog-standard native save */
	const char *flagstr = print_flags(opt, STORE);

	tmp = g_strdup_printf("store '%s' %s%s", savename, 
			      (showlist)? storelist : "", flagstr);
    } else if (has_suffix(savename, ".dat")) { 
	/* saving in "traditional" mode as ".dat" */
	tmp = g_strdup_printf("store '%s' %s -t", savename, 
			      (showlist)? storelist : "");
	opt = OPT_T;
    } else {
	/* standard data save */
	tmp = g_strdup_printf("store '%s' %s", savename, 
			      (showlist)? storelist : ""); 
    }

    if (!overwrite) {
	fp = gretl_fopen(savename, "rb");
	if (fp != NULL) {
	    fclose(fp);
	    if (yes_no_dialog(_("gretl: save data"), 
			      _("There is already a data file of this name.\n"
				"OK to overwrite it?"), 
			      0) == GRETL_NO) {
		goto store_get_out;
	    }
	}
    }

    err = check_specific_command(tmp);
    if (err) goto store_get_out;

    if (!WRITING_DB(opt)) {
	err = cmd_init(tmp);
	if (err) goto store_get_out;
    }

    /* back up existing datafile if need be (not for databases) */
    if (!WRITING_DB(opt)) {
	if ((fp = gretl_fopen(savename, "rb")) && fgetc(fp) != EOF &&
	    fclose(fp) == 0) {
	    tmp = g_strdup_printf("%s~", savename);
	    if (copyfile(savename, tmp)) {
		err = 1;
		goto store_get_out;
	    }
	}
    } 

    /* actually write the data to file */
    err = write_data(savename, cmd.list, (const double **) Z, datainfo, 
		     opt, &paths);

    if (err) {
	if (WRITING_DB(opt) && err == E_DB_DUP) {
	    err = db_write_response(savename, cmd.list);
	    if (err) {
		goto store_get_out;
	    }
	} else {
	    errbox(_("Write of data file failed\n%s"), get_gretl_errmsg());
	    err = 1;
	    goto store_get_out;
	} 
    }   

    /* record that data have been saved, etc. */
    if (!DATA_EXPORT(opt)) {
	mkfilelist(FILE_LIST_DATA, savename);
	if (paths.datfile != savename) {
	    strcpy(paths.datfile, savename);
	}
	data_status = (HAVE_DATA | USER_DATA);
	if (is_gzipped(paths.datfile)) {
	    data_status |= GZIPPED_DATA;
	} 
	edit_info_state(TRUE);
	set_sample_label(datainfo);	
    }

    /* tell the user */
    if (WRITING_DB(opt)) {
	database_description_dialog(savename);
    } 

 store_get_out:

    if (storelist != NULL) {
	free(storelist);
	storelist = NULL;
    }

    g_free(tmp);

    return err;
}

#ifdef G_OS_WIN32

static int get_latex_path (char *latex_path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, latex, NULL, MAXLEN, latex_path, &p);

    return (ret == 0);
}

#else

static int spawn_latex (char *texsrc)
{
    GError *error = NULL;
    gchar *errout = NULL, *sout = NULL;
    gchar *argv[] = {
	latex,
	"\\batchmode",
	"\\input",
	texsrc,
	NULL
    };
    int ok, status;
    int ret = LATEX_OK;

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (paths.userdir, /* working dir */
		       argv,
		       NULL,    /* envp */
		       G_SPAWN_SEARCH_PATH,
		       NULL,    /* child_setup */
		       NULL,    /* user_data */
		       &sout,   /* standard output */
		       &errout, /* standard error */
		       &status, /* exit status */
		       &error);

    if (!ok) {
	errbox(error->message);
	g_error_free(error);
	ret = LATEX_EXEC_FAILED;
    } else if (errout && *errout) {
	errbox(errout);
	ret = LATEX_ERROR;
    } else if (status != 0) {
	gchar *errmsg;

	errmsg = g_strdup_printf("%s\n%s", 
				 _("Failed to process TeX file"),
				 sout);
	errbox(errmsg);
	g_free(errmsg);
	ret = LATEX_ERROR;
    }

    if (errout != NULL) g_free(errout);
    if (sout != NULL) g_free(sout);

    return ret;
}

#endif /* !G_OS_WIN32 */

int latex_compile (char *texshort)
{
#ifdef G_OS_WIN32
    static char latex_path[MAXLEN];
    char tmp[MAXLEN];
#endif
    int err = LATEX_OK;

#ifdef G_OS_WIN32
    if (*latex_path == 0 && get_latex_path(latex_path)) {
	DWORD dw = GetLastError();
	win_show_error(dw);
	return LATEX_EXEC_FAILED;
    }

    sprintf(tmp, "\"%s\" %s", latex_path, texshort);
    if (winfork(tmp, paths.userdir, SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE)) {
	return LATEX_EXEC_FAILED;
    }
#else
    err = spawn_latex(texshort);
#endif /* G_OS_WIN32 */

    return err;
}

#ifdef OSX_BUILD

#include <Carbon/Carbon.h>

int osx_open_file (const char *path)
{
    FSRef r;
    int err;
    
    err = FSPathMakeRef(path, &r, NULL);
    if (!err) {
	err = LSOpenFSRef(&r, NULL);
    }

    return err;
}

int osx_open_url (const char *url)
{
    CFStringRef s;
    CFURLRef u;
    int err;
    
    s = CFStringCreateWithBytes(NULL, url, strlen(url), 
                                kCFStringEncodingASCII, 
				0);
    if (s == NULL) {
        err = 1;
    } else {
        u = CFURLCreateWithString(NULL, s, NULL);
        if (u == NULL) {
	    err = 1;
        } else {
	    err = LSOpenCFURLRef(u, NULL);
	    CFRelease(u);
	}
	CFRelease(s);
    }

    return err;
}

#endif /* OSX_BUILD */

static void view_or_save_latex (PRN *bprn, const char *fname, int saveit)
{
    char texfile[MAXLEN], texbase[MAXLEN], tmp[MAXLEN];
    int dot, err = LATEX_OK;
    char *texshort = NULL;
    const char *buf;
    PRN *fprn;

    *texfile = 0;

    if (fname != NULL) {
	strcpy(texfile, fname);
    } else {
	sprintf(texfile, "%swindow.tex", paths.userdir);
    } 

    fprn = gretl_print_new_with_filename(texfile);
    if (fprn == NULL) {
	errbox(_("Couldn't write to %s"), texfile);
	return;
    }

    gretl_tex_preamble(fprn, tex_eqn_format(bprn));
    buf = gretl_print_get_buffer(bprn);
    pputs(fprn, buf);
    pputs(fprn, "\n\\end{document}\n");

    gretl_print_destroy(fprn);
	
    if (saveit) {
	return;
    }

    dot = dotpos(texfile);
    *texbase = 0;
    strncat(texbase, texfile, dot);

    texshort = strrchr(texbase, SLASH) + 1;
    if (texshort == NULL) {
	errbox(_("Failed to process TeX file"));
	return;
    } 

    err = latex_compile(texshort);

    if (err == LATEX_OK) {
#if defined(G_OS_WIN32)
	if (!strncmp(latex, "pdf", 3)) {
	    sprintf(tmp, "%s.pdf", texbase);
	    if ((int) ShellExecute(NULL, "open", tmp, NULL, NULL, SW_SHOW) <= 32) {
		DWORD dw = GetLastError();
		win_show_error(dw);
	    }
	} else {
	    sprintf(tmp, "\"%s\" \"%s.dvi\"", viewdvi, texbase);
	    if (WinExec(tmp, SW_SHOWNORMAL) < 32) {
		DWORD dw = GetLastError();
		win_show_error(dw);
	    }
	}
#elif defined(OSX_BUILD)
	if (!strncmp(latex, "pdf", 3)) {
	    sprintf(tmp, "%s.pdf", texbase);
	} else {
	    sprintf(tmp, "%s.dvi", texbase);
	}
	if (osx_open_file(tmp)) {
	    errbox(_("Couldn't open %s"), tmp);
	}
#else
	if (!strncmp(latex, "pdf", 3)) {
	    sprintf(tmp, "%s.pdf", texbase);
	    gretl_fork(viewpdf, tmp);
	} else {
	    gretl_fork(viewdvi, texbase);
	}
#endif
    }

#ifdef KILL_DVI_FILE
    sleep(2); /* let forked xdvi get the DVI file */
    sprintf(tmp, "%s.dvi", texbase);
    remove(tmp);
#endif

    sprintf(tmp, "%s.log", texbase);
    if (err == LATEX_ERROR) {
	view_file(tmp, 0, 1, 78, 350, VIEW_FILE);
    } else {
	remove(texfile);
	remove(tmp);
    }

    sprintf(tmp, "%s.aux", texbase);
    remove(tmp);
}

void view_latex (PRN *prn)
{
    view_or_save_latex(prn, NULL, 0);
}

void save_latex (PRN *prn, const char *fname)
{
    if (prn != NULL) {
	view_or_save_latex(prn, fname, 1);
    } else {
	save_graph_page(fname);
    }
}

#if 0
static const char *exec_string (int i)
{
    switch (i) {
    case CONSOLE_EXEC: return "CONSOLE_EXEC";
    case SCRIPT_EXEC: return "SCRIPT_EXEC";
    case SESSION_EXEC: return "SESSION_EXEC";
    case SAVE_SESSION_EXEC: return "SAVE_SESSION_EXEC";
    default: return "Unknown";
    }
}
#endif

static int ok_script_file (const char *runfile)
{
    FILE *fp;
    char myline[32];
    int content = 0;

    fp = gretl_fopen(runfile, "r");
    if (fp == NULL) {
	errbox(_("Couldn't open script"));
	return 0;
    }

    /* check that the file has something in it */
    while (fgets(myline, sizeof myline, fp)) {
	const char *p = myline;

	while (*p) {
	    if (!isspace(*p)) {
		content = 1;
		break;
	    }
	    p++;
	}
	if (content) break;
    }

    fclose(fp);

    if (!content) {
	errbox(_("No commands to execute"));
	return 0;
    }

    return 1;
}

static void output_line (const char *line, PRN *prn) 
{
    int n = strlen(line);

    if ((line[0] == '(' && line[1] == '*') ||
	(line[n-1] == ')' && line[n-2] == '*')) {
	pprintf(prn, "\n%s\n", line);
    } else if (line[0] == '#') {
	if (gretl_compiling_loop()) {
	    pprintf(prn, "> %s\n", line);;
	} else {
	    pprintf(prn, "%s\n", line);
	}
    } else if (!string_is_blank(line)) {
	safe_print_line(line, prn);
    }
}

/* run commands from runfile or buf, output to prn */

int execute_script (const char *runfile, const char *buf,
		    PRN *prn, int exec_code)
{
    FILE *fb = NULL;
    char line[MAXLINE] = {0};
    char tmp[MAXLINE] = {0};
    int including = (exec_code & INCLUDE_EXEC);
    int exec_err = 0;

#if 0
    debug_print_model_info(models[0], "Start of execute_script, models[0]");
#endif

    if (runfile != NULL) { 
	/* we'll get commands from file */
	if (!ok_script_file(runfile)) {
	    return -1;
	}
	fb = gretl_fopen(runfile, "r");
    } else { 
	/* no runfile, commands from buffer */
	if (buf == NULL || *buf == '\0') {
	    errbox(_("No commands to execute"));
	    return -1;	
	}
	bufgets_init(buf);
    }

    /* reset model count to 0 if starting/saving session (?) */
    if (exec_code == SESSION_EXEC) {
	reset_model_count();
    }

    if (!including) {
	gui_script_logo(prn);
    }

    *cmd.word = '\0';

    while (strcmp(cmd.word, "quit")) {
	if (gretl_execute_loop()) { 
	    exec_err = gretl_loop_exec(line, &Z, &datainfo, models, prn);
	    if (exec_err) {
		goto endwhile;
	    }
	} else { 
	    char *gotline = NULL;

	    *line = '\0';

	    if (gretl_executing_function()) {
		gotline = gretl_function_get_line(line, MAXLINE,
						  &Z, &datainfo, &exec_err);
	    } else if (fb != NULL) {
		gotline = fgets(line, MAXLINE, fb);
	    } else {
		gotline = bufgets(line, MAXLINE, buf);
	    }

	    if (gotline == NULL) {
		/* done reading */
		goto endwhile;
	    }

	    while (top_n_tail(line) && !exec_err) {
		/* handle backslash-continued lines */
		*tmp = '\0';

		if (gretl_executing_function()) {
		    gretl_function_get_line(tmp, MAXLINE, &Z, &datainfo, &exec_err);
		} else if (fb != NULL) {
		    fgets(tmp, MAXLINE, fb);
		} else {
		    bufgets(tmp, MAXLINE, buf); 
		}

		if (!exec_err && *tmp != '\0') {
		    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
			pprintf(prn, _("Maximum length of command line "
				       "(%d bytes) exceeded\n"), MAXLINE);
			exec_err = 1;
		    } else {
			strcat(line, tmp);
			compress_spaces(line);
		    }
		}		
	    }

	    if (!exec_err) {
		if (!strncmp(line, "noecho", 6)) {
		    set_gretl_echo(0);
		} else if (!strncmp(line, "(* saved objects:", 17)) { 
		    strcpy(line, "quit"); 
		} else if (gretl_echo_on() && !including) {
		    output_line(line, prn);
		}
		strcpy(tmp, line);
		exec_err = gui_exec_line(line, prn, exec_code, runfile);
	    }

	    if (exec_err) {
		pprintf(prn, _("\nError executing script: halting\n"));
		pprintf(prn, "> %s\n", tmp);
		goto endwhile;
	    }
	} /* end non-loop command processor */
    } /* end while command != quit */

 endwhile:

    if (fb != NULL) {
	fclose(fb);
    }

    refresh_data();

    return exec_err;
}

static int script_model_test (int test_ci, int model_id, PRN *prn)
{
    int m, mc;
    const char *no_gui_test = 
	N_("Sorry, can't do this.\nTo operate on a model estimated "
	   "via the graphical interface, please use the\nmenu items in "
	   "the model window.\n");

    if (model_id != 0) {
	m = modelspec_index_from_model_id(modelspec, model_id);
    } else {
	m = modelspec_last_index(modelspec);
    }  

#ifdef MSPEC_DEBUG
    fprintf(stderr, "model_test_start: test_ci=%d, model_id=%d, m=%d\n",
	    test_ci, model_id, m);
#endif

    mc = get_model_count();

    if (m < 0) { 
	/* reference model not found */
	if (mc == 0) {
	    pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	} else if (model_id == 0) {
	    /* requested "the last model" */
	    pputs(prn, _(no_gui_test));
	} else if (model_id > mc) {
	    /* requested specific, out-of-range model */
	    pprintf(prn, _("Can't do this: there is no model %d\n"), model_id);
	} else {
	    /* requested specific, in-range model, but it's a gui model */
	    pputs(prn, _(no_gui_test));
	}
	return 1;
    }
     
    if (!command_ok_for_model(test_ci, 
			      model_ci_from_modelspec(modelspec, m))) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
	return 1;
    }			      

    if (model_sample_issue(NULL, modelspec, m, datainfo)) {
	pputs(prn, _("Can't do: the current data set is different from "
		     "the one on which\nthe reference model was estimated\n"));
	return 1;
    }

    return 0;
}

static void do_autofit_plot (PRN *prn)
{
    int lines[1];
    int plotlist[3];
    int err;

    plotlist[0] = 2;
    plotlist[1] = gretl_model_get_depvar(models[0]);
    plotlist[2] = varindex(datainfo, "autofit");

    lines[0] = 1; 
    err = gnuplot(plotlist, lines, NULL, &Z, datainfo,
		  &plot_count, OPT_T); 

    if (err) {
	pprintf(prn, _("gnuplot command failed\n"));
    } else {
	register_graph();
    }
}

int gui_exec_line (char *line, PRN *prn, int exec_code, const char *myname) 
{
    int lines[1];
    int dbdata = 0, alt_model = 0;
    int script_code = exec_code;
    double rho;
    char runfile[MAXLEN], datfile[MAXLEN];
    char linecopy[1024];
    char texfile[MAXLEN];
    GnuplotFlags plotflags = 0;
    gretlopt testopt = OPT_NONE;
    MODEL tmpmod;
    GRETL_VAR *var = NULL;
    PRN *outprn = NULL;
    int k, err = 0;

    static int in_comment;

#if CMD_DEBUG
    fprintf(stderr, "gui_exec_line: exec_code = %d\n",
	    exec_code);
#endif

    if (string_is_blank(line)) {
	return 0;
    }

    if (gretl_compiling_function()) {
#if CMD_DEBUG
	fprintf(stderr, "gui_exec_line: compiling function\n");
#endif
	err = gretl_function_append_line(line);
	if (err) {
	    errmsg(err, prn);
	}
	return err;
    }     

    if (!in_comment) {
	/* catch requests relating to saved objects, which are not
	   really "commands" as such */
	k = saved_object_action(line, prn);
	if (k == 1) return 0;   /* action was OK */
	if (k == -1) return 1;  /* action was faulty */

	/* are we ready for this? */
	if (!data_status && !cmd.ignore && !ready_for_command(line)) {
	    pprintf(prn, _("You must open a data file first\n"));
	    return 1;
	}
    }

    *linecopy = 0;
    strncat(linecopy, line, sizeof linecopy - 1);

    /* if we're stacking commands for a loop, parse "lightly" */
    if (gretl_compiling_loop()) { 
	err = get_command_index(line, &cmd, datainfo);
    } else {
	err = parse_command_line(line, &cmd, &Z, datainfo);
    }

#if CMD_DEBUG
    fprintf(stderr, "gui_exec_line: '%s'\n cmd.ci = %d\n", line, cmd.ci);
#endif

    if (err) {
        errmsg(err, prn);
        return 1;
    }

    /* are we in a multi-line comment block? */
    in_comment = cmd.ignore;

    if (cmd.ci < 0) {
	return 0; /* nothing there, or a comment */
    }

    if (sys != NULL && cmd.ci != END && cmd.ci != EQUATION &&
	cmd.ci != SYSTEM) {
	pprintf(prn, _("Command '%s' ignored; not valid within "
		       "equation system\n"), line);
	gretl_equation_system_destroy(sys);
	sys = NULL;
	return 1;
    }

    if (cmd.ci == LOOP && exec_code == CONSOLE_EXEC) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }

    if (cmd.ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd.ci)) {
            pprintf(prn, _("Sorry, this command is not available in loop mode\n"));
            return 1;
        }
	err = gretl_loop_append_line(line, cmd.ci, cmd.opt, &Z, datainfo);
	if (err) {
	    print_gretl_errmsg(prn);
	    return 1;
	} 
	return 0;
    } 

    /* Attach outprn to a specific buffer, if wanted */
    if (*cmd.savename != '\0' && TEXTSAVE_OK(cmd.ci)) {
	if (bufopen(&outprn)) return 1;
    } else {
	outprn = prn;
    }

    check_for_loop_only_options(cmd.ci, cmd.opt, prn);

    switch (cmd.ci) {

    case ADDOBS:
    case ADF: 
    case COINT2:
    case COINT: 
    case CORR: 
    case CRITERIA: 
    case CRITICAL: 
    case DATA:
    case DIFF: 
    case DISCRETE:
    case DUMMIFY:
    case ESTIMATE:
    case FNCALL:
    case FUNC:
    case FUNCERR:
    case GRAPH: 
    case PLOT: 
    case HURST: 
    case INFO: 
    case KPSS:
    case LABELS: 
    case LAGS: 
    case LDIFF: 
    case LOGS:
    case MAHAL:
    case MATRIX:
    case MEANTEST: 
    case MULTIPLY: 
    case OUTFILE: 
    case PCA:
    case PRINT: 
    case REMEMBER:
    case RENAME:
    case RHODIFF:
    case RMPLOT: 
    case RUNS: 
    case SDIFF:
    case SETINFO:
    case SHELL:
    case SPEARMAN: 
    case SQUARE: 
    case STORE:
    case SUMMARY:
    case TRANSPOSE:
    case VARLIST:
    case VARTEST: 
    case XTAB:
	err = simple_commands(&cmd, line, &Z, datainfo, outprn);
	if (err) {
	    errmsg(err, prn);
	} else if (cmd.ci == DATA) {
	    register_data(NULL, NULL, 0);
	}
	break;

    case ADD:
    case OMIT:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
    plain_add_omit:
	clear_model(models[1]);
	if (cmd.ci == ADD || cmd.ci == ADDTO) {
	    err = add_test(cmd.list, models[0], models[1], 
			   &Z, datainfo, cmd.opt, outprn);
	} else {
	    err = omit_test(cmd.list, models[0], models[1],
			    &Z, datainfo, cmd.opt, outprn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	} else {
	    /* for command-line use, we keep a stack of 
	       two models, and recycle the places */
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(models[0], models[1]);
	    }
	    clear_model(models[1]);
	}
	break;	

    case ADDTO:
    case OMITFROM:
	k = atoi(cmd.param);
	if ((err = script_model_test(cmd.ci, k, prn))) {
	    break;
	}
	if (k == (models[0])->ID) {
	    goto plain_add_omit;
	}
	err = re_estimate(modelspec_get_command_by_id(modelspec, k), 
			  &tmpmod, &Z, datainfo);
	if (err) {
	    pprintf(prn, _("Failed to reconstruct model %d\n"), k);
	    break;
	} 
	clear_model(models[1]);
	tmpmod.ID = k;
	if (cmd.ci == ADDTO) {
	    err = add_test(cmd.list, &tmpmod, models[1], 
			   &Z, datainfo, cmd.opt, outprn);
	} else {
	    err = omit_test(cmd.list, &tmpmod, models[1],
			    &Z, datainfo, cmd.opt, outprn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	    break;
	} else {
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(models[0], models[1]);
	    }
	    clear_model(models[1]);
	}
	clear_model(&tmpmod);
	break;

    case AR:
	clear_model(models[0]);
	*models[0] = ar_func(cmd.list, &Z, datainfo, cmd.opt, outprn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	}
	break;

    case ARCH:
	clear_model(models[1]);
	*models[1] = arch_model(cmd.list, cmd.order, &Z, datainfo, 
				cmd.opt, outprn);
	if ((err = (models[1])->errcode)) {
	    errmsg(err, prn);
	}
	if ((models[1])->ci == ARCH) {
	    alt_model = 1;
	    swap_models(models[0], models[1]);
	}
	clear_model(models[1]);
	break;

    case ARMA:
	clear_model(models[0]);
	*models[0] = arma(cmd.list, (const double **) Z, datainfo,
			  cmd.opt, outprn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	} else {	
	    printmodel(models[0], datainfo, cmd.opt, outprn);
	}
	break;

    case BXPLOT:
	if (cmd.nolist) { 
	    err = boolean_boxplots(line, &Z, datainfo, (cmd.opt != 0));
	} else {
	    err = boxplots(cmd.list, NULL, &Z, datainfo, (cmd.opt != 0));
	}
	break;

    case CHOW:
    case CUSUM:
    case QLRTEST:
    case RESET:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	if (cmd.ci == CHOW || cmd.ci == QLRTEST) {
	    err = chow_test(line, models[0], &Z, datainfo, testopt, outprn);
	} else if (cmd.ci == CUSUM) {
	    err = cusum_test(models[0], &Z, datainfo, testopt, outprn);
	} else {
	    err = reset_test(models[0], &Z, datainfo, testopt, outprn);
	}
	if (err) {
	    errmsg(err, prn);
	} 
	break;

    case COEFFSUM:
    case VIF:
        if ((err = script_model_test(cmd.ci, 0, prn))) break;
	if (cmd.ci == COEFFSUM) {
	    err = sum_test(cmd.list, models[0], &Z, datainfo, outprn);
	} else {
	    err = vif_test(models[0], &Z, datainfo, outprn);
	}
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, &Z, datainfo, cmd.ci,
			   &err, cmd.opt, outprn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	clear_model(models[0]);
	*models[0] = ar1_lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, rho);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	} else {
	    printmodel(models[0], datainfo, cmd.opt, outprn);
	}
	break;

    case CORRGM:
	k = atoi(cmd.param);
	err = corrgram(cmd.list[1], k, 0, &Z, datainfo, outprn, OPT_A);
	if (err) {
	    pprintf(prn, _("Failed to generate correlogram\n"));
	}
	break;

    case DELEET:
	if (dataset_locked()) {
	    break;
	}
	if (complex_subsampled()) {
	    pputs(prn, _("Can't delete a variable when in sub-sample"
		    " mode\n"));
	    break;
	}
	maybe_prune_delete_list(cmd.list);
	if (cmd.list[0] == 0) {
	    err = 1;
	} else {
	    err = dataset_drop_listed_variables(cmd.list, &Z, datainfo, &k);
	}
	if (err) {
	    pputs(prn, _("Failed to shrink the data set"));
	} else {
	    if (k) {
		pputs(prn, _("Take note: variables have been renumbered"));
		pputc(prn, '\n');
	    }
	    maybe_clear_selector(cmd.list);
	    maybe_list_vars(datainfo, prn);
	}
	break;

    case END:
	if (!strcmp(cmd.param, "system")) {
	    err = gretl_equation_system_finalize(sys, &Z, datainfo, outprn);
	    if (err) {
		errmsg(err, prn);
	    }
	    sys = NULL;
	} else if (!strcmp(cmd.param, "mle") || !strcmp(cmd.param, "nls")) {
	    clear_model(models[0]);
	    *models[0] = nls(&Z, datainfo, cmd.opt, outprn);
	    if ((err = (models[0])->errcode)) {
		errmsg(err, prn);
	    } else {
		alt_model = 1;
		printmodel(models[0], datainfo, cmd.opt, outprn);
	    }
	} else if (!strcmp(cmd.param, "restrict")) {
	    err = gretl_restriction_set_finalize(rset, (const double **) Z, 
						 datainfo, prn);
	    if (err) {
		errmsg(err, prn);
	    }
	    rset = NULL;
	} else {
	    err = 1;
	}
	break;

    case BREAK:
    case ENDLOOP:
	pprintf(prn, _("You can't end a loop here, "
		       "you haven't started one\n"));
	err = 1;
	break;

    case EQUATION:
	/* one equation within a system */
	err = gretl_equation_system_append(sys, cmd.list);
	if (err) {
	    sys = NULL;
	    errmsg(err, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	if ((models[0])->errcode == E_NAN) {
	    pprintf(prn, _("Couldn't format model\n"));
	    break;
	}
	if ((err = script_model_test(cmd.ci, 0, prn))) {
	    break;
	}
	strcpy(texfile, cmd.param);
	err = texprint(models[0], datainfo, texfile, 
		       (cmd.ci == EQNPRINT)? (cmd.opt | OPT_E) :
		       cmd.opt);
	if (err) {
	    pprintf(prn, _("Couldn't open tex file for writing\n"));
	} else {
	    pprintf(prn, _("Model printed to %s\n"), texfile);
	}
	break;

    case FCAST:
    case FIT:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	if (cmd.ci == FIT) {
	    err = add_forecast("fcast autofit", models[0], &Z, datainfo, cmd.opt);
	} else {
	    err = add_forecast(line, models[0], &Z, datainfo, cmd.opt);
	}
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (cmd.ci == FIT) {
		pprintf(prn, _("Retrieved fitted values as \"autofit\"\n"));
	    }
	    maybe_list_vars(datainfo, prn);
	    if (cmd.ci == FIT && exec_code == CONSOLE_EXEC && 
		dataset_is_time_series(datainfo)) {
		do_autofit_plot(prn);
	    }
	}
	break;

    case FCASTERR:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = display_forecast(line, models[0], &Z, datainfo, 
			       cmd.opt, outprn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case FREQ:
	err = freqdist(cmd.list[1], (const double **) Z, 
		       datainfo, (exec_code == CONSOLE_EXEC),
		       cmd.opt, prn);
	if (!err && exec_code == CONSOLE_EXEC) {
	    register_graph();
	}
	break;

    case GENR:
	err = generate(line, &Z, datainfo, cmd.opt, prn);
	if (err) {
	    errmsg(err, prn);
	} 
	break;

    case GNUPLOT:
    case SCATTERS:
	plotflags = gp_flags((exec_code == SCRIPT_EXEC), cmd.opt);
	if (cmd.ci == GNUPLOT) {
	    if ((cmd.opt & OPT_M) || (cmd.opt & OPT_Z) || (cmd.opt & OPT_S)) { 
		err = gnuplot(cmd.list, NULL, cmd.param, &Z, datainfo,
			      &plot_count, plotflags); 
	    } else {
		lines[0] = (cmd.opt != 0);
		err = gnuplot(cmd.list, lines, cmd.param, 
			      &Z, datainfo, &plot_count, plotflags);
	    }
	} else {
	    err = multi_scatters(cmd.list, (const double **) Z, datainfo, 
				 &plot_count, plotflags);
	}

	if (err) {
	    pputs(prn, (cmd.ci == GNUPLOT)? 
		  _("gnuplot command failed\n") :
		  _("scatters command failed\n"));
	} else {
	    if (exec_code == CONSOLE_EXEC && *cmd.savename == '\0') {
		register_graph();
	    } else if (exec_code == SCRIPT_EXEC) {
		pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	    }
	    err = maybe_save_graph(&cmd, gretl_plotfile(),
				   GRETL_OBJ_GRAPH, prn);
	}
	break;

    case HAUSMAN:
	err = script_model_test(cmd.ci, 0, prn);
	if (!err) {
	    err = panel_hausman_test(models[0], &Z, datainfo, cmd.opt, outprn);
	}
	break;

    case HELP:
	help(cmd.param, paths.cli_helpfile, prn);
	break;

    case IMPORT:
	if (dataset_locked()) {
	    break;
	}
        err = getopenfile(line, datfile, &paths, 0, 0);
        if (err) {
            pprintf(prn, _("import command is malformed\n"));
            break;
        }
	if (data_status & HAVE_DATA) {
	    close_session();
	}
        if (cmd.opt & OPT_B) {
            err = import_box(&Z, &datainfo, datfile, prn);
	} else if (cmd.opt & OPT_O) {
	    err = import_octave(&Z, &datainfo, datfile, prn);
        } else {
            err = import_csv(&Z, &datainfo, datfile, prn);
	}
        if (!err) { 
	    maybe_display_string_table();
	    data_status |= IMPORT_DATA;
	    register_data(datfile, NULL, 1);
            print_smpl(datainfo, 0, prn);
            varlist(datainfo, prn);
            pprintf(prn, _("You should now use the \"print\" command "
			   "to verify the data\n"));
            pprintf(prn, _("If they are OK, use the  \"store\" command "
			   "to save them in gretl format\n"));
        }
        break;

    case OPEN:
    case APPEND:
	if (dataset_locked()) {
	    break;
	}
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    errbox(_("'open' command is malformed"));
	    break;
	}
#if CMD_DEBUG
	fprintf(stderr, "OPEN in gui_exec_line, datfile='%s'\n", datfile);
#endif
	k = detect_filetype(datfile, &paths, prn);
	dbdata = (k == GRETL_NATIVE_DB || k == GRETL_RATS_DB);

	if (cmd.ci == OPEN && (data_status & HAVE_DATA) && !dbdata) {
	    close_session();
	}

	if (k == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, datfile, prn);
	} else if (k == GRETL_OCTAVE) {
	    err = import_octave(&Z, &datainfo, datfile, prn);
	} else if (k == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else if (k == GRETL_XML_DATA) {
	    err = gretl_read_gdt(&Z, &datainfo, datfile, &paths, data_status, prn, 1);
	} else if (WORKSHEET_IMPORT(k)) {
	    err = import_other(&Z, &datainfo, k, datfile, prn);
	} else if (dbdata) {
	    err = set_db_name(datfile, k, &paths, prn);
	} else {
	    err = gretl_get_data(&Z, &datainfo, datfile, &paths, data_status, prn);
	}
	if (err) {
	    gui_errmsg(err);
	    break;
	}
	if (!dbdata && cmd.ci != APPEND) {
	    strncpy(paths.datfile, datfile, MAXLEN-1);
	}
	if (k == GRETL_CSV_DATA || k == GRETL_BOX_DATA || dbdata) {
	    data_status |= IMPORT_DATA;
	    maybe_display_string_table();
	}
	if (datainfo->v > 0 && !dbdata) {
	    if (cmd.ci == APPEND) {
		register_data(NULL, NULL, 0);
	    } else {
		register_data(paths.datfile, NULL, 0);
	    }
	    varlist(datainfo, prn);
	}
	*paths.currdir = '\0'; 
	break;

    case LEVERAGE:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = leverage_test(models[0], &Z, datainfo, cmd.opt, outprn);
	if (err > 1) {
	    errmsg(err, prn);
	} else if (cmd.opt & OPT_S) {
	    maybe_list_vars(datainfo, prn);
	}
	break;

    case LMTEST:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = lmtest_driver(cmd.param, models[0], &Z, datainfo, 
			    cmd.opt, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case GARCH:
    case HSK:
    case LAD:
    case LOGISTIC:
    case LOGIT:
    case POISSON:
    case PROBIT:
    case TOBIT:
    case TSLS:
	clear_model(models[0]);
	if (cmd.ci == LOGIT || cmd.ci == PROBIT) {
	    *models[0] = logit_probit(cmd.list, &Z, datainfo, cmd.ci, cmd.opt);
	} else if (cmd.ci == LOGISTIC) {
	    *models[0] = logistic_model(cmd.list, &Z, datainfo, cmd.param);
	} else if (cmd.ci == TOBIT) {
	    *models[0] = tobit_model(cmd.list, &Z, datainfo,
				     (cmd.opt & OPT_V)? outprn : NULL);
	} else if (cmd.ci == POISSON) {
	    *models[0] = poisson_model(cmd.list, &Z, datainfo,
				       (cmd.opt & OPT_V)? outprn : NULL);
	} else if (cmd.ci == TSLS) {
	    *models[0] = tsls_func(cmd.list, TSLS, &Z, datainfo, cmd.opt);
	} else if (cmd.ci == HSK) {
	    *models[0] = hsk_func(cmd.list, &Z, datainfo);
	} else if (cmd.ci == LAD) {
	    *models[0] = lad(cmd.list, &Z, datainfo);
	} else if (cmd.ci == GARCH) {
	    *models[0] = garch(cmd.list, &Z, datainfo, cmd.opt, outprn);
	} else {
	    /* can't happen */
	    err = 1;
	    break;
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	} else {
	    printmodel(models[0], datainfo, cmd.opt, outprn);
	}
	break;

    case MODELTAB:
	err = modeltab_parse_line(line, models[0], prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case MLE:
    case NLS:
	err = nls_parse_line(cmd.ci, line, (const double **) Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	} else {
	    gretl_cmd_set_context(&cmd, cmd.ci);
	}
	break;

    case NULLDATA:
	if (dataset_locked()) {
	    break;
	}
	k = gretl_int_from_string(cmd.param, (const double **) Z, 
				  datainfo, &err);
	if (!err && k < 2) {
	    err = 1;
	}
	if (err) {
	    pputs(prn, _("Data series length count missing or invalid\n"));
	    break;
	}
	if (data_status & HAVE_DATA) {
	    close_session();
	}
	err = open_nulldata(&Z, datainfo, data_status, k, prn);
	if (err) { 
	    pprintf(prn, _("Failed to create empty data set\n"));
	} else {
	    register_data(NULL, NULL, 0);
	}
	break;

    case OLS:
    case WLS:
    case HCCM:
    case PANEL:
	clear_model(models[0]);
	if (cmd.ci == PANEL) {
	    *models[0] = panel_model(cmd.list, &Z, datainfo, cmd.opt, prn);
	} else {
	    *models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt);
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn); 
	} else {
	    printmodel(models[0], datainfo, cmd.opt, outprn);
	}
	break;

#ifdef ENABLE_GMP
    case MPOLS:
	clear_model(models[0]);
	*models[0] = mp_ols(cmd.list, (const double **) Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn); 
	} else {
	    printmodel(models[0], datainfo, cmd.opt, outprn);
	}
	break;
#endif

    case PERGM:
	err = periodogram(cmd.list[1], &Z, datainfo, cmd.opt | OPT_N, outprn);
	if (err) pprintf(prn, _("Failed to generate periodogram\n"));
	break;

    case PRINTF:
	err = do_printf(line, &Z, datainfo, prn);
	break;

    case PVALUE:
	if (strcmp(line, "pvalue") == 0) {
	    help("pvalue", paths.cmd_helpfile, prn);	    
	} else {
	    err = (batch_pvalue(line, (const double **) Z, datainfo, 
				outprn, OPT_NONE) == NADBL);
	}
	break;

    case QUIT:
	if (exec_code == CONSOLE_EXEC) {
	   pprintf(prn, _("Please use the Close button to exit\n")); 
	} else {
	    pprintf(prn, _("Script done\n"));
	} 
	break;

    case RUN:
    case INCLUDE:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pprintf(prn, _("Run command failed\n"));
	    break;
	}
	if (cmd.ci == INCLUDE && gretl_is_xml_file(runfile)) {
	    err = load_user_function_file(runfile);
	    if (err) {
		pputs(prn, _("Error reading function definitions\n"));
	    }
	    break;
	}
	if (myname != NULL && strcmp(runfile, myname) == 0) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    return 1;
	}
	if (exec_code == CONSOLE_EXEC) {
	    script_code = SCRIPT_EXEC;
	}
	if (cmd.ci == INCLUDE) {
	    pprintf(prn, _("%s opened OK\n"), runfile);
	    script_code |= INCLUDE_EXEC;
	}
	err = execute_script(runfile, NULL, prn, script_code);
	break;

    case SET:
	err = execute_set_line(line, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case SETOBS:
	err = set_obs(line, &Z, datainfo, cmd.opt);
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (datainfo->n > 0) {
		set_sample_label(datainfo);
		print_smpl(datainfo, 0, prn);
	    } else {
		pprintf(prn, _("setting data frequency = %d\n"), datainfo->pd);
	    }
	}
	break;	

    case SETMISS:
        set_miss(cmd.list, cmd.param, Z, datainfo, prn);
        break;

    case SMPL:
	k = 0;
	if (cmd.opt == OPT_F) {
	    gui_restore_sample();
	    k = 1;
	} else if (cmd.opt) {
	    err = restrict_sample(line, cmd.list, &Z, &datainfo, 
				  cmd.opt, prn);
	} else { 
	    err = set_sample(line, (const double **) Z, datainfo);
	}

	if (err) {
	    errmsg(err, prn);
	} else {
	    print_smpl(datainfo, get_full_length_n(), prn);
	    if (cmd.opt && cmd.opt != OPT_F) { /* FIXME? */
		set_sample_label_special();
	    } else {
		set_sample_label(datainfo);
	    }
	    if (!k) {
		restore_sample_state(TRUE);
	    }
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model or system */
	if (rset == NULL) {
	    if (*cmd.param == '\0') {
		/* if param is non-blank, we're restricting a named system */
		if ((err = script_model_test(cmd.ci, 0, prn))) break;
	    }	
	    rset = restriction_set_start(line, models[0], datainfo, cmd.opt);
	    if (rset == NULL) {
		err = 1;
		errmsg(err, prn);
	    } else {
		gretl_cmd_set_context(&cmd, RESTRICT);
	    }
	} else {
	    err = restriction_set_parse_line(rset, line);
	    if (err) {
		errmsg(err, prn);
		rset = NULL;
	    }	
	}
	break;

    case SYSTEM:
	/* system of equations */
	if (sys == NULL) {
	    sys = system_start(line, cmd.opt);
	    if (sys == NULL) {
		err = 1;
		errmsg(err, prn);
	    } else {
		gretl_cmd_set_context(&cmd, SYSTEM);
		maybe_save_system(&cmd, sys, prn);
	    }
	} else {
	    err = system_parse_line(sys, line, datainfo);
	    if (err) {
		errmsg(err, prn);
		sys = NULL;
	    }
	}
	break;

    case TESTUHAT:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	if (genr_fit_resid(models[0], &Z, datainfo, GENR_RESID, 1)) {
	    pputs(prn, _("Out of memory!"));
	    pputc(prn, '\n');
	    err = 1;
	} else {
	    FreqDist *freq; 
	 
	    freq = get_freq(datainfo->v - 1, (const double **) Z, datainfo, 
			    (models[0])->ncoeff, OPT_NONE);
	    dataset_drop_last_variables(1, &Z, datainfo);
	    if (!(err = freq_error(freq, prn))) {
		print_freq(freq, outprn); 
		free_freq(freq);
	    }
	}
	break;

    case VAR:
	var = gretl_VAR(cmd.order, cmd.list, &Z, datainfo, cmd.opt, 
			outprn, &err);
	if (var != NULL) {
	    err = maybe_save_var(&cmd, &var, prn);
	}
	break;

    case VECM:
	var = vecm(cmd.order, atoi(cmd.extra), cmd.list, &Z, datainfo, 
		   cmd.opt, outprn, &err);
	if (var != NULL) {
	    err = maybe_save_var(&cmd, &var, prn);
	}
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		       "in libgretl\n"), cmd.word);
	break;
    } /* end of command switch */

    /* clean up in case a user function bombed */
    if (err && gretl_executing_function()) {
	gretl_function_stop_on_error(&Z, &datainfo, prn);
    }    

    /* log the specific command? */
    if (exec_code == CONSOLE_EXEC && !err) {
	cmd_init(line);
    }

    /* save specific output (text) buffer? */
    if (outprn != NULL && outprn != prn) {
	if (!err) {
	    err = save_text_buffer(outprn, cmd.savename, prn);
	} else {
	    gretl_print_destroy(outprn);
	}
	outprn = NULL;
    }

    if (!err && (is_model_cmd(cmd.word) || alt_model)
	&& !is_quiet_model_test(cmd.ci, cmd.opt)) {
	if (is_model_cmd(cmd.word)) {
	    stack_script_modelspec(models[0]);
	}
	maybe_save_model(&cmd, models[0], prn);
    }

    return (err != 0);
}


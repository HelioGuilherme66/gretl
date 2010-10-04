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

/* TRAMO/SEATS, X-12-ARIMA plugin for gretl */

#include "libgretl.h"
#include "estim_private.h"

#include <gtk/gtk.h>
#include "tramo_x12a.h"

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#ifdef WIN32
# include <windows.h>
#else
# include <signal.h>
#endif

enum prog_codes {
    TRAMO_SEATS,
    TRAMO_ONLY,
    X12A
};

const char *x12a_save_strings[] = {
    "d11", /* seasonally adjusted */
    "d12", /* trend/cycle */
    "d13", /* irregular */
    NULL
};

static int tramo_got_irfin;

const char *tramo_save_strings[] = {
    "safin.t", /* final seasonally adjusted series */
    "trfin.t", /* final trend */
    "irfin.t", /* final irregular factor (component) */
    "irreg.t", /* irregular component (logs) */
    NULL
};

const char *tx_descrip_formats[] = {
    N_("seasonally adjusted %s"),
    N_("trend/cycle for %s"),
    N_("irregular component of %s")
};

const char *default_mdl = {
    "# ARIMA specifications that will be tried\n"
    "(0 1 1)(0 1 1) X\n"
    "(0 1 2)(0 1 1) X\n"
    "(2 1 0)(0 1 1) X\n"
    "(0 2 2)(0 1 1) X\n"
    "(2 1 2)(0 1 1)\n"
};  

#ifndef WIN32

#define SP_DEBUG 0

static int glib_spawn (const char *workdir, const char *fmt, ...)
{
    GError *gerr = NULL;
    gchar *sout = NULL;
    gchar *serr = NULL;
    gchar *argv[8];
    char *s;
    va_list ap;
    int i, ok, nargs;
    int status = 0;
    int err = 0;

    argv[0] = g_strdup(fmt);
    argv[1] = NULL;
    i = nargs = 1;

    va_start(ap, fmt);

    while ((s = va_arg(ap, char *))) {
	argv[i] = g_strdup(s);
	argv[++i] = NULL;
    }

    va_end(ap);

    nargs = i;

#if SP_DEBUG
    fputs("spawning the following:\n", stderr);
    for (i=0; i<nargs; i++) {
	fprintf(stderr, " argv[%d] = '%s'\n", i, argv[i]);
    }
#endif

    gretl_error_clear();

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (workdir,
		       argv,
		       NULL,
		       G_SPAWN_SEARCH_PATH,
		       NULL,
		       NULL,
		       &sout,
		       &serr,
		       &status,
		       &gerr);

    if (!ok) {
	gretl_errmsg_set(gerr->message);
	fprintf(stderr, "spawn failed: '%s'\n", gerr->message);
	g_error_free(gerr);
	err = E_EXTERNAL;
    } else if (status != 0) {
	if (sout && *sout) {
	    gretl_errmsg_set(sout);
	    fprintf(stderr, "spawn: status = %d: '%s'\n", status, sout);
	} else {
	    gretl_errmsg_set(_("Command failed"));
	    fprintf(stderr, "spawn: status = %d\n", status);
	}
	err = E_DATA;
    } else if (serr && *serr) {
	fprintf(stderr, "stderr: '%s'\n", serr);
    }

    if (serr != NULL) g_free(serr);
    if (sout != NULL) g_free(sout);

    for (i=0; i<nargs; i++) {
	if (err) {
	    if (i == 0) {
		fputc(' ', stderr);
	    }
	    fprintf(stderr, "%s ", argv[i]);
	    if (i == nargs - 1) {
		fputc('\n', stderr);
	    }
	}
	free(argv[i]);
    }

    return err;
}

#endif /* !WIN32 */

static void toggle_outliers (GtkToggleButton *b, tx_request *request)
{
    request->outliers = gtk_toggle_button_get_active(b);
}

static void toggle_trading_days (GtkToggleButton *b, tx_request *request)
{
    request->trdays = gtk_toggle_button_get_active(b);
}

static void set_logtrans (GtkButton *b, tx_request *request)
{
    gpointer p = g_object_get_data(G_OBJECT(b), "transval");

    request->logtrans = GPOINTER_TO_INT(p);
}

static void show_x12a_options (tx_request *request, GtkBox *vbox)
{
    GtkWidget *tmp, *b[3];
    GSList *group;

    tmp = gtk_check_button_new_with_label(_("Detect and correct for outliers"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), request->outliers);
    g_signal_connect(GTK_TOGGLE_BUTTON(tmp), "toggled",
		     G_CALLBACK(toggle_outliers), request);

    tmp = gtk_check_button_new_with_label(_("Correct for trading days effect"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), request->trdays);
    g_signal_connect(GTK_TOGGLE_BUTTON(tmp), "toggled",
		     G_CALLBACK(toggle_trading_days), request);

    tmp = gtk_hseparator_new();
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);

    b[0] = gtk_radio_button_new_with_label(NULL, _("Log transformation"));
    gtk_widget_show(b[0]);
    gtk_box_pack_start(vbox, b[0], FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[0]), "toggled",
		     G_CALLBACK(set_logtrans), request);
    g_object_set_data(G_OBJECT(b[0]), "transval", GINT_TO_POINTER(1));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b[0]));
    b[1] = gtk_radio_button_new_with_label(group, _("No log transformation"));
    gtk_widget_show(b[1]);
    gtk_box_pack_start(vbox, b[1], FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[1]), "toggled",
		     G_CALLBACK(set_logtrans), request);
    g_object_set_data(G_OBJECT(b[1]), "transval", GINT_TO_POINTER(2));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b[1]));
    b[2] = gtk_radio_button_new_with_label(group, _("Automatic"));
    gtk_widget_show(b[2]);
    gtk_box_pack_start(vbox, b[2], FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[2]), "toggled",
		     G_CALLBACK(set_logtrans), request);
    g_object_set_data(G_OBJECT(b[2]), "transval", GINT_TO_POINTER(3));

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b[request->logtrans - 1]), 
				 TRUE); 

    tmp = gtk_hseparator_new();
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);

    tmp = gtk_label_new(_("Save data"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);

    tmp = gtk_check_button_new_with_label(_("Seasonally adjusted series"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    request->opts[TX_SA].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), 
				 (*request->popt & OPT_A)? TRUE : FALSE);

    tmp = gtk_check_button_new_with_label(_("Trend/cycle"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    request->opts[TX_TR].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), 
				 (*request->popt & OPT_B)? TRUE : FALSE);

    tmp = gtk_check_button_new_with_label(_("Irregular"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    request->opts[TX_IR].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), 
				 (*request->popt & OPT_C)? TRUE : FALSE);

    tmp = gtk_hseparator_new();
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);
    
    tmp = gtk_check_button_new_with_label(_("Generate graph"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    request->opts[TRIGRAPH].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), 
				 (*request->popt & OPT_G)? TRUE : FALSE);

    tmp = gtk_check_button_new_with_label(_("Show full output"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    request->opts[TEXTOUT].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp),
				 (*request->popt & OPT_Q)? FALSE : TRUE);
}

static int tx_dialog (tx_request *request)
{
    GtkWidget *hbox, *vbox;
    gint i, ret = 0;

    for (i=0; i<TX_MAXOPT; i++) {
	request->opts[i].check = NULL;
    }

    request->dialog = 
	gtk_dialog_new_with_buttons ((request->prog == TRAMO_SEATS)?
				     "TRAMO/SEATS" : "X-12-ARIMA",
				     NULL,
				     GTK_DIALOG_MODAL | 
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     GTK_STOCK_CANCEL,
				     GTK_RESPONSE_REJECT,
				     GTK_STOCK_OK,
				     GTK_RESPONSE_ACCEPT,
				     NULL);

    vbox = gtk_vbox_new(FALSE, 0);    

    if (request->prog == TRAMO_SEATS) {
	gtk_dialog_set_has_separator(GTK_DIALOG(request->dialog), FALSE);
	show_tramo_options(request, vbox);
    } else {
	show_x12a_options(request, GTK_BOX(vbox));
    }

    gtk_widget_show(vbox);
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(request->dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    ret = gtk_dialog_run(GTK_DIALOG(request->dialog));

    return (ret == GTK_RESPONSE_ACCEPT)? 1 : 0;
}

static void get_seats_command (char *seats, const char *tramo)
{
    char *p;

    strcpy(seats, tramo);
    p = strrchr(seats, SLASH);
    if (p != NULL) {
	strcpy(p + 1, "seats");
    } else {
	strcpy(seats, "seats");
    }
}

/* try to avoid collision of date and graph key, which
   defaults to right top
*/

static void set_keypos (const double *x, int t1, int t2,
			FILE *fp)
{
    int T = t2 - t1 + 1;

    if (T <= 12) {
	if (x[t2] > x[t1]) {
	    fputs("set key left top\n", fp);
	}
    } else {
	double m1, m2;
	int r = T / 6;

	m1 = gretl_mean(t1, t1 + r, x);
	m2 = gretl_mean(t2 - r, t2, x);

	if (m2 > m1) {
	    fputs("set key left top\n", fp);
	}
    }
}

static int graph_series (const double **Z, const DATAINFO *pdinfo, 
			 tx_request *req)
{
    FILE *fp = NULL;
    const double *obs;
    int v_sa = TX_SA + 1;
    int v_tr = TX_TR + 1;
    int v_ir = TX_IR + 1;
    double irbar, irmax;
    int sub1 = 0;
    double f1;
    char title[32];
    int t, err = 0;

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    fp = get_plot_input_stream(PLOT_TRI_GRAPH, &err);
    if (err) {
	return err;
    }

    gretl_push_c_numeric_locale();

    if (pdinfo->pd == 4) {
	if ((pdinfo->t2 - pdinfo->t1) / 4 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 4\n", fp);
	}
    } else if (pdinfo->pd == 12) {
	if ((pdinfo->t2 - pdinfo->t1) / 12 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 12\n", fp);
	}
    }

    if (req->seasonal_ok) {
	f1 = 0.33;
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.32\n", fp);
    } else {
	f1 = 0.5;
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
	tramo_got_irfin = 0; /* I _think_ this may be right */
    }

    if (req->prog == TRAMO_SEATS && tramo_got_irfin) {
	/* need to divide by 100? */
	irmax = 10.0;
    } else {
	irmax = 0.5;
    }

    irbar = gretl_mean(pdinfo->t1, pdinfo->t2, Z[v_ir]);
    if (irbar > irmax) {
	sub1 = 1;
    }

    /* irregular component */
    if (sub1) {
	sprintf(title, "%s - 1", _("irregular"));
    } else {
	sprintf(title, "%s", _("irregular"));
    }

    fprintf(fp, "set bars 0\n"
	    "set origin 0.0,0.0\n"
	    "set xzeroaxis\n"
	    "plot '-' using 1:%s title '%s' w impulses\n",
	    (sub1)? "($2-1.0)" : "2", title);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	double yt = Z[v_ir][t];

	if (req->prog == TRAMO_SEATS && tramo_got_irfin) {
	    yt /= 100.0;
	}

	fprintf(fp, "%.10g %.10g\n", obs[t], yt);
    }
    fputs("e\n", fp);

    set_keypos(Z[0], pdinfo->t1, pdinfo->t2, fp);

    /* actual (in var 0) vs trend/cycle */

    fprintf(fp, "set origin 0.0,%.2f\n"
	    "plot '-' using 1:2 title '%s' w l, \\\n"
	    " '-' using 1:2 title '%s' w l\n",
	    f1, pdinfo->varname[0], _("trend/cycle"));

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	fprintf(fp, "%.10g %.10g\n", obs[t], Z[0][t]);
    }
    fputs("e , \\\n", fp);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	fprintf(fp, "%.10g %.10g\n", obs[t], Z[v_tr][t]);
    }
    fputs("e\n", fp);

    if (req->seasonal_ok) {
	/* actual vs seasonally adjusted */
	fprintf(fp, "set origin 0.0,0.66\n"
		"plot '-' using 1:2 title '%s' w l, \\\n"
		" '-' using 1:2 title '%s' w l\n",
		pdinfo->varname[0], _("adjusted"));

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    fprintf(fp, "%.10g %.10g\n", obs[t], Z[0][t]);
	}
	fputs("e\n", fp);

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    fprintf(fp, "%.10g %.10g\n", obs[t], Z[v_sa][t]);
	}
	fputs("e\n", fp);
    }

    fputs("set nomultiplot\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static void copy_variable (double **targZ, DATAINFO *targinfo, int targv,
			   double **srcZ, DATAINFO *srcinfo, int srcv)
{
    int t;

    for (t=0; t<targinfo->n; t++) {
	targZ[targv][t] = srcZ[srcv][t];
    }

    strcpy(targinfo->varname[targv], srcinfo->varname[srcv]);
    strcpy(VARLABEL(targinfo, targv), VARLABEL(srcinfo, srcv));
}

static void clear_tramo_files (const char *path, const char *vname)
{
    char fname[MAXLEN];
    int i;

    for (i=0; tramo_save_strings[i] != NULL; i++) {
	sprintf(fname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
		tramo_save_strings[i]);
	gretl_remove(fname);
    }

    sprintf(fname, "%s%coutput%c%s.out", path, SLASH, SLASH, vname);
    gretl_remove(fname);
}

static void clear_x12a_files (const char *path, const char *vname)
{
    char fname[MAXLEN];
    int i;

    for (i=0; x12a_save_strings[i] != NULL; i++) {
	sprintf(fname, "%s%c%s.%s", path, SLASH, vname,
		x12a_save_strings[i]);
	gretl_remove(fname);
    }    

    sprintf(fname, "%s%c%s.out", path, SLASH, vname);
    gretl_remove(fname);

    sprintf(fname, "%s%c%s.err", path, SLASH, vname);
    gretl_remove(fname);
}

static int add_series_from_file (const char *path, int src,
				 double **Z, DATAINFO *pdinfo,
				 int targv, tx_request *request,
				 char *errmsg)
{
    FILE *fp;
    char line[128], sfname[MAXLEN];
    char varname[VNAMELEN], date[8];
    double x;
    int d, yr, per, err = 0;
    int t;

    if (request->prog == TRAMO_SEATS) {
	tramo_got_irfin = 1;
	sprintf(sfname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
		tramo_save_strings[src]);
    } else {
	char *p;

	strcpy(sfname, path);
	p = strrchr(sfname, '.');
	if (p != NULL) {
	    strcpy(p + 1, x12a_save_strings[src]);
	}
    }

    fp = gretl_fopen(sfname, "r");

    if (fp == NULL) {
	/* couldn't open the file we wanted */
	int gotit = 0;

	/* This is a bit of a pest: under some configurations, tramo/seats
	   outputs a series "irfin"; sometimes that is not created, but
	   we do get an "irreg".  So if we can't find the one, try looking
	   for the other.  Also, the seasonally adjusted series "safin"
	   is not always available.
	*/
	if (request->prog == TRAMO_SEATS) {
	    if (src == TX_IR) { 
		/* try "irreg" */
		sprintf(sfname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
			tramo_save_strings[src + 1]);
		fp = gretl_fopen(sfname, "r");
		if (fp != NULL) {
		    gotit = 1;
		}
		tramo_got_irfin = 0;
	    } else if (src == TX_SA) {
		/* scrub all use of seasonal series */
		request->seasonal_ok = 0;
		if (request->opts[src].save) {
		    request->opts[src].save = 0;
		    request->savevars -= 1;
		}
		return 0;
	    }
	}

	if (!gotit) {
	    fprintf(stderr, "Couldn't open %s\n", sfname);
	    sprintf(errmsg, _("Couldn't open %s"), sfname);
	    return 1;
	}
    }

    /* formulate name of new variable to add */
    if (request->prog == TRAMO_SEATS) {
	sprintf(varname, "%.5s_%.2s", pdinfo->varname[0], 
		tramo_save_strings[src]);
    } else {
	sprintf(varname, "%.4s_%s", pdinfo->varname[0], 
		x12a_save_strings[src]);
    }

    /* copy varname and label into place */
    strcpy(pdinfo->varname[targv], varname);
    sprintf(VARLABEL(pdinfo, targv), _(tx_descrip_formats[src]), pdinfo->varname[0]);

    if (request->prog == TRAMO_SEATS) {
	strcat(VARLABEL(pdinfo, targv), " (TRAMO/SEATS)");
    } else {
	strcat(VARLABEL(pdinfo, targv), " (X-12-ARIMA)");
    }	

    for (t=0; t<pdinfo->n; t++) {
	Z[targv][t] = NADBL;
    }

    gretl_push_c_numeric_locale();

    if (request->prog == TRAMO_SEATS) {
	int i = 0;

	t = pdinfo->t1;
	while (fgets(line, 127, fp)) {
	    i++;
	    if (i >= 7 && sscanf(line, " %lf", &x) == 1) {
		if (t >= pdinfo->n) {
		    fprintf(stderr, "t = %d >= pdinfo->n = %d\n", t, pdinfo->n);
		    err = 1;
		    break;
		}		
		Z[targv][t++] = x;
	    }
	}
    } else {
	/* grab the data from the x12arima file */
	while (fgets(line, 127, fp)) {
	    if (*line == 'd' || *line == '-') {
		continue;
	    }
	    if (sscanf(line, "%d %lf", &d, &x) != 2) {
		err = 1; 
		break;
	    }
	    yr = d / 100;
	    per = d % 100;
	    sprintf(date, "%d.%d", yr, per);
	    t = dateton(date, pdinfo);
	    if (t < 0 || t >= pdinfo->n) {
		err = 1;
		break;
	    }
	    Z[targv][t] = x;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

static int grab_deseasonal_series (double *y, const DATAINFO *pdinfo,
				   int prog, const char *path)
{
    FILE *fp;
    char line[128], sfname[MAXLEN], date[8];
    double yt;
    int d, yr, per, err = 0;
    int t;

    if (prog == TRAMO_SEATS) {
	sprintf(sfname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
		tramo_save_strings[TX_SA]);
    } else {
	char *p;

	strcpy(sfname, path);
	p = strrchr(sfname, '.');
	if (p != NULL) {
	    strcpy(p + 1, x12a_save_strings[TX_SA]);
	}
    }

    fp = gretl_fopen(sfname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    if (prog == TRAMO_SEATS) {
	int i = 0;

	t = pdinfo->t1;
	while (fgets(line, 127, fp)) {
	    i++;
	    if (i >= 7 && sscanf(line, " %lf", &yt) == 1) {
		if (t >= pdinfo->n) {
		    fprintf(stderr, "t = %d >= pdinfo->n = %d\n", t, pdinfo->n);
		    err = E_DATA;
		    break;
		}		
		y[t++] = yt;
	    }
	}
    } else {
	/* grab the data from the x12arima file */
	while (fgets(line, 127, fp)) {
	    if (*line == 'd' || *line == '-') {
		continue;
	    }
	    if (sscanf(line, "%d %lf", &d, &yt) != 2) {
		err = 1; 
		break;
	    }
	    yr = d / 100;
	    per = d % 100;
	    sprintf(date, "%d.%d", yr, per);
	    t = dateton(date, pdinfo);
	    if (t < 0 || t >= pdinfo->n) {
		err = E_DATA;
		break;
	    }
	    y[t] = yt;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

static void request_opts_init (tx_request *request, const DATAINFO *pdinfo)
{
    int i;

    request->savevars = 0;
    request->logtrans = 3; /* x12a: automatic logs or not */
    request->outliers = 1; /* x12a: detect outliers */
    request->trdays = (pdinfo->pd == 12); /* x12a: trading days */

    for (i=0; i<TX_MAXOPT; i++) {
	request->opts[i].save = 0;
    }

    request->seasonal_ok = 1;
}

static void set_opts (tx_request *request)
{
    int i;

    request->savevars = 0;

    *request->popt &= ~(OPT_A | OPT_B | OPT_C);

    for (i=0; i<TX_MAXOPT; i++) {
	if (request->opts[i].check != NULL && 
	    gtk_toggle_button_get_active
	    (GTK_TOGGLE_BUTTON(request->opts[i].check))) {
	    request->opts[i].save = 1;
	    if (i < TRIGRAPH) {
		request->savevars++;
		if (i == 0) {
		    *request->popt |= OPT_A;
		} else if (i == 1) {
		    *request->popt |= OPT_B;
		} else if (i == 2) {
		    *request->popt |= OPT_C;
		}
	    }
	} else {
	    request->opts[i].save = 0;
	} 
    }
}

static void cancel_savevars (tx_request *request)
{
    int i;

    request->savevars = 0;

    for (i=0; i<TX_MAXOPT; i++) {
	request->opts[i].save = 0;
    } 
}

static int write_tramo_file (const char *fname, 
			     const double *y, 
			     const char *vname, 
			     const DATAINFO *pdinfo,
			     tx_request *request) 
{
    int startyr, startper;
    int T = pdinfo->t2 - pdinfo->t1 + 1; 
    char *p, tmp[8];
    double x;
    FILE *fp;
    int t;

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    gretl_push_c_numeric_locale();

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    p = strchr(tmp, '.');
    if (p != NULL) {
	startper = atoi(p + 1);
    } else {
	startper = 1;
    }

    fprintf(fp, "%s\n", vname);
    fprintf(fp, "%d %d %d %d\n", T, startyr, startper, pdinfo->pd);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (t && t % pdinfo->pd == 0) fputc('\n', fp);
	if (na(y[t])) {
	    fputs("-99999 ", fp);
	} else {
	    fprintf(fp, "%g ", y[t]);
	}
    }
    fputc('\n', fp);

    if (request == NULL) {
	fputs("$INPUT rsa=3,out=2,$END\n", fp);
    } else if (print_tramo_options(request, fp) == 0) {
	/* not running SEATS */
	request->prog = TRAMO_ONLY; 
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static int write_spc_file (const char *fname, const double *y,
			   const char *vname,
			   const DATAINFO *pdinfo, 
			   const int *savelist,
			   int logtrans,
			   int outliers,
			   int trdays)
{
    int startyr, startper;
    char *p, tmp[8];
    double x;
    FILE *fp;
    int i, t;

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    gretl_push_c_numeric_locale();

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    p = strchr(tmp, '.');
    if (p != NULL) {
	startper = atoi(p + 1);
    } else {
	startper = 1;
    }

    fprintf(fp, "series{\n period=%d\n title=\"%s\"\n", pdinfo->pd, vname);
    fprintf(fp, " start=%d.%d\n", startyr, startper);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(y[t])) {
	    fputs(" missingcode=-99999\n", fp);
	    break;
	}
    }

    fputs(" data=(\n", fp);

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(y[t])) {
	    fputs("-99999 ", fp);
	} else {
	    fprintf(fp, "%g ", y[t]);
	}
	if ((i + 1) % 7 == 0) {
	    fputc('\n', fp);
	}
	i++;
    }
    fputs(" )\n}\n", fp);

    if (logtrans == 1) {
	fputs("transform{function=log}\n", fp);
    } else if (logtrans == 2) {
	fputs("transform{function=none}\n", fp);
    } else {
	fputs("transform{function=auto}\n", fp);
    }

    if (trdays) {
	fputs("regression{variables = td}\n", fp);
    }

    if (outliers) {
	fputs("outlier{}\n", fp);
    }

    fputs("automdl{}\n", fp); 

    fputs("x11{", fp);

    if (savelist[0] > 0) {
	if (savelist[0] == 1) {
	    fprintf(fp, " save=%s ", x12a_save_strings[savelist[1]]); 
	} else {
	    fputs(" save=( ", fp);
	    for (i=1; i<=savelist[0]; i++) {
		fprintf(fp, "%s ", x12a_save_strings[savelist[i]]);
	    }
	    fputs(") ", fp);
	}
    }

    fputs("}\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static void form_savelist (int *list, tx_request *request)
{
    int i, j = 1;

    list[0] = 0;

    for (i=0; i<TRIGRAPH; i++) {
	if (request->opts[TRIGRAPH].save || request->opts[i].save) {
	    list[0] += 1;
	    list[j++] = i;
	}
    }
}

static void copy_basic_data_info (DATAINFO *targ, DATAINFO *src)
{
    targ->sd0 = src->sd0;
    strcpy(targ->stobs, src->stobs); 
    targ->t1 = src->t1;
    targ->t2 = src->t2;
    targ->pd = src->pd;
    targ->structure = src->structure;
}

static int save_vars_to_dataset (double ***pZ, DATAINFO *pdinfo,
				 double **tmpZ, DATAINFO *tmpinfo,
				 int *varlist, tx_request *request,
				 char *errmsg)
{
    int i, v, j, addvars = 0;

    /* how many vars are wanted, and how many are new? */
    for (i=1; i<=varlist[0]; i++) {
	if (request->opts[varlist[i]].save && 
	    series_index(pdinfo, tmpinfo->varname[i]) == pdinfo->v) {
	    addvars++;
	}
    }

    if (addvars > 0 && dataset_add_series(addvars, pZ, pdinfo)) {
	return E_ALLOC;
    }

    j = pdinfo->v - addvars;

    for (i=1; i<=varlist[0]; i++) {
	if (request->opts[varlist[i]].save) {
	    v = series_index(pdinfo, tmpinfo->varname[i]);
	    if (v < pdinfo->v) {
		copy_variable(*pZ, pdinfo, v, tmpZ, tmpinfo, i);
	    } else {
		copy_variable(*pZ, pdinfo, j++, tmpZ, tmpinfo, i);
	    }
	}
    }

    return 0;
}

#ifdef WIN32

static int helper_spawn (const char *path, const char *vname,
			 const char *workdir, int prog)
{
    char *cmd = NULL;
    int err = 0;

    if (prog == TRAMO_ONLY) {
	cmd = g_strdup_printf("\"%s\" -i %s -k serie", path, vname);
    } else if (prog == TRAMO_SEATS) {
	cmd = g_strdup_printf("\"%s\" -OF %s", path, vname);
    } else if (prog == X12A) {
	cmd = g_strdup_printf("\"%s\" %s -r -p -q", path, vname);
    } else {
	return E_EXTERNAL;
    }

    if (cmd == NULL) {
	err = E_ALLOC;
    } else {
	err = win_run_sync(cmd, workdir);
	g_free(cmd);
    }

    return err;
}

#else

static int helper_spawn (const char *path, const char *vname,
			 const char *workdir, int prog)
{
    int err;

    if (prog == TRAMO_ONLY) {
	err = glib_spawn(workdir, path, "-i", vname, "-k", "serie", NULL);
    } else if (prog == TRAMO_SEATS) {
	err = glib_spawn(workdir, path, "-OF", vname, NULL);
    } else if (prog == X12A) {
	err = glib_spawn(workdir, path, vname, "-r", "-p", "-q", NULL);
    } else {
	err = E_EXTERNAL;
    }

    return err;
}

#endif

/* make a default x12a.mdl file if it doesn't already exist */

static int check_x12a_model_file (const char *workdir, char *fname)
{
    FILE *fp;
    int err = 0;

    sprintf(fname, "%s%cx12a.mdl", workdir, SLASH);
    fp = gretl_fopen(fname, "r");

    if (fp != NULL) {
	fclose(fp); /* assume we're OK */
    } else {
	fp = gretl_fopen(fname, "w");
	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    fprintf(fp, "%s", default_mdl);
	    fclose(fp);
	}
    } 

    return err;
}

static int check_sample_bound (int prog, const DATAINFO *pdinfo,
			       char *errmsg)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;

    if (prog == TRAMO_SEATS && T > 600) {
	strcpy(errmsg, _("TRAMO can't handle more than 600 observations.\n"
			 "Please select a smaller sample."));
	return E_EXTERNAL;
    } else if (prog == X12A) {
	int pdmax = get_x12a_maxpd();

	if (T > 50 * pdmax) {
	    sprintf(errmsg, _("X-12-ARIMA can't handle more than %d observations.\n"
			      "Please select a smaller sample."), 50 * pdmax);
	    return E_EXTERNAL;
	}
    }

    return 0;
}

int write_tx_data (char *fname, int varnum, 
		   double ***pZ, DATAINFO *pdinfo, gretlopt *opt, 
		   int tramo, int *graph_ok, char *errmsg)
{
    const char *exepath;
    const char *workdir;
    char vname[VNAMELEN];
    int savelist[4];
    tx_request request;
    double **tmpZ;
    DATAINFO *tmpinfo;
    int i, doit;
    int err = 0;

    *errmsg = '\0';

    if (tramo) {
	request.prog = TRAMO_SEATS;
	exepath = gretl_tramo();
	workdir = gretl_tramo_dir();
    } else {
	request.prog = X12A;
	exepath = gretl_x12_arima();
	workdir = gretl_x12_arima_dir();
    }	
	
    request_opts_init(&request, pdinfo);

    err = check_sample_bound(request.prog, pdinfo, errmsg);
    if (err) {
	return err;
    }

    request.pd = pdinfo->pd;
    request.popt = opt;

    /* show dialog and get option settings */
    doit = tx_dialog(&request); 
    if (!doit) {
	gtk_widget_destroy(request.dialog);
	return 0;
    }
    set_opts(&request);
    gtk_widget_destroy(request.dialog);

#if 0
    if (request.prog == TRAMO_SEATS) {
	print_tramo_options(&request, stderr);
	return 1;
    }
#endif

    /* create little temporary dataset */
    tmpinfo = create_auxiliary_dataset(&tmpZ, 4, pdinfo->n);
    if (tmpinfo == NULL) {
	return E_ALLOC;
    }

    copy_basic_data_info(tmpinfo, pdinfo);

    if (request.prog == X12A) { 
	err = check_x12a_model_file(workdir, fname);
	if (err) {
	    goto bailout;
	}
    } 

    strcpy(vname, pdinfo->varname[varnum]);
    form_savelist(savelist, &request);

    if (request.prog == X12A) { 
	/* write out the .spc file for x12a */
	sprintf(fname, "%s%c%s.spc", workdir, SLASH, vname);
	write_spc_file(fname, (*pZ)[varnum], vname, pdinfo, savelist,
		       request.logtrans, request.outliers, request.trdays);
    } else { 
	/* TRAMO, possibly plus SEATS */
	lower(vname);
	gretl_trunc(vname, 8);
	sprintf(fname, "%s%c%s", workdir, SLASH, vname);
	/* next line: this also sets request->prog = TRAMO_ONLY if
	   SEATS is not to be run */
	write_tramo_file(fname, (*pZ)[varnum], vname, pdinfo, &request);
	if (request.prog == TRAMO_ONLY) {
	    cancel_savevars(&request); /* FIXME later */
	    savelist[0] = 0;
	}
    }

    /* now run the program(s): we try to ensure that any
       old output files get deleted first 
    */

    if (request.prog == X12A) {
	clear_x12a_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, X12A);
    } else { 
	char seats[MAXLEN];

	clear_tramo_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, TRAMO_ONLY);

	if (!err && request.prog == TRAMO_SEATS) {
	    get_seats_command(seats, exepath);
	    err = helper_spawn(seats, vname, workdir, TRAMO_SEATS);
	}
    }

    if (err == E_EXTERNAL) {
	/* fatal: couldn't run program */
	*fname = '\0';
    } else if (err) {
	if (request.prog == X12A) {
	    sprintf(fname, "%s%c%s.err", workdir, SLASH, vname);
	} else {
	    sprintf(fname, "%s%coutput%c%s.out", workdir, SLASH, SLASH, vname);
	}
    } else {
	if (request.prog == X12A) {
	    sprintf(fname, "%s%c%s.out", workdir, SLASH, vname); 
	} else {
	    sprintf(fname, "%s%coutput%c%s.out", workdir, SLASH, SLASH, vname);
	    if (request.prog == TRAMO_ONLY) {
		/* no graph offered */
		request.opts[TRIGRAPH].save = 0;
		*graph_ok = 0;
	    }
	} 

	/* save vars locally if needed; graph if wanted */
	if (savelist[0] > 0) {
	    const char *path = (request.prog == X12A)? fname : workdir;

	    copy_variable(tmpZ, tmpinfo, 0, *pZ, pdinfo, varnum);

	    for (i=1; i<=savelist[0]; i++) {
		err = add_series_from_file(path, savelist[i], tmpZ, tmpinfo,
					   i, &request, errmsg);
		if (err) {
		    fprintf(stderr, "i = %d: add_series_from_file() failed\n", i);
		    if (request.prog == X12A) {
			/* switch to X12A error file */
			sprintf(fname, "%s%c%s.err", workdir, SLASH, vname);
		    }
		    break;
		} 
	    }

	    if (!err) {
		if (request.opts[TRIGRAPH].save) {
		    err = graph_series((const double **) tmpZ, tmpinfo, &request);
		    if (err) {
			fprintf(stderr, "graph_series() failed\n");
		    } else {
			*opt |= OPT_G;
		    }
		} else {
		    *opt &= ~OPT_G;
		}
	    }

	    if (request.prog == X12A) {
		if (request.opts[TEXTOUT].save) {
		    *opt &= ~OPT_Q;
		} else {
		    *opt |= OPT_Q;
		}
	    }
	}

	/* now save the local vars to main dataset, if wanted */
	if (!err && request.savevars > 0) {
	    err = save_vars_to_dataset(pZ, pdinfo, tmpZ, tmpinfo, savelist, 
				       &request, errmsg);
	}
    }

 bailout:

    destroy_dataset(tmpZ, tmpinfo);

    return err;
}

int adjust_series (const double *x, double *y, const DATAINFO *pdinfo, 
		   int tramo)
{
    int prog = (tramo)? TRAMO_SEATS : X12A;
    int savelist[2] = {1, TX_SA};
    const char *vname = "x";
    const char *exepath;
    const char *workdir;
    char fname[MAXLEN];
    int err = 0;

    if (prog == X12A) { 
	exepath = gretl_x12_arima();
	workdir = gretl_x12_arima_dir();
    } else { 
	exepath = gretl_tramo();
	workdir = gretl_tramo_dir();
    }    

    if (prog == X12A) { 
	err = check_x12a_model_file(workdir, fname);
	if (err) {
	    return err;
	}
    } 

    if (prog == X12A) { 
	sprintf(fname, "%s%c%s.spc", workdir, SLASH, vname);
	write_spc_file(fname, x, vname, pdinfo, savelist, 2, 0,
		       pdinfo->pd == 12);
    } else { 
	sprintf(fname, "%s%c%s", workdir, SLASH, vname);
	write_tramo_file(fname, x, vname, pdinfo, NULL); 
    }

    if (prog == X12A) {
	clear_x12a_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, X12A);
    } else { 
	char seats[MAXLEN];

	clear_tramo_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, TRAMO_ONLY);
	if (!err) {
	    get_seats_command(seats, exepath);
	    err = helper_spawn(seats, vname, workdir, TRAMO_SEATS);
	}
    }

    if (!err) {
	const char *path = (prog == X12A)? fname : workdir;

	err = grab_deseasonal_series(y, pdinfo, prog, path);
    }

    return err;
}

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

/* dialogs.c for gretl */

#include "gretl.h"

extern const char *version_string;

#include "selector.h"

extern GtkWidget *active_edit_id;
extern GtkWidget *active_edit_name;
extern GtkWidget *active_edit_text;

extern void show_spreadsheet (DATAINFO *pdinfo);

GtkWidget *open_dialog;
int session_saved;

/* ........................................................... */

void random_dialog (gpointer data, guint code, GtkWidget *widget) 
{
    if (code == GENR_UNIFORM) {
	edit_dialog (_("gretl: uniform variable"), 
		     _("Enter name for variable, and\n"
		     "minimum and maximum values:"), 
		     "unif 0 100",  
		     _("Apply"), do_random, NULL, 
		     _("  Cancel  "), NULL, NULL, GENR_UNIFORM, GENR);
    } else if (code == GENR_NORMAL) {
	edit_dialog (_("gretl: normal variable"), 
		     _("Enter name, mean and standard deviation:"), 
		     "norm 0 1", 
		     _("Apply"), do_random, NULL, 
		     _("  Cancel  "), NULL, NULL, GENR_NORMAL, GENR);
    }
}

/* ........................................................... */

static void fix_obsstr (char *str)
{
    char pt = get_local_decpoint();
    char *p;

    p = strchr(str, ':');
    if (p != NULL) {
	*p = pt;
    }

    if (pt == ',' && (p = strchr(str, '.'))) {
	*p = pt;
    }

    if (pt == '.' && (p = strchr(str, ','))) {
	*p = pt;
    }
}

/* ........................................................... */

static void prep_spreadsheet (GtkWidget *widget, dialog_t *data)
{
    gchar *edttext;
    char dataspec[32];
    char *test, stobs[9], endobs[9], firstvar[9];
    double sd0, ed0;

    edttext = gtk_entry_get_text (GTK_ENTRY (data->edit));
    strncpy(dataspec, edttext, 31);
    if (dataspec[0] == '\0') return;

    /* check validity of dataspec */
    if (sscanf(dataspec, "%8s %8s %8s", stobs, endobs, firstvar) != 3) {
	errbox(_("Insufficient dataset information supplied"));
	return;
    }

    /* daily data: special */
    if (datainfo->pd == 5 || datainfo->pd == 7) {
	int err = 0;
	sd0 = (double) get_epoch_day(stobs); 
	ed0 = (double) get_epoch_day(endobs);

	if (sd0 < 0) {
	    err = 1;
	    sprintf(errtext, _("Invalid starting observation '%s'"), stobs);
	}
	if (!err && ed0 < 0) {
	    err = 1;
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	}
	if (err) {
	    errbox(errtext);
	    return;
	}
    } else { /* not daily data */
	fix_obsstr(stobs);
	fix_obsstr(endobs);

	sd0 = strtod(stobs, &test);
	if (strcmp(stobs, test) == 0 || test[0] != '\0' || sd0 < 0) {
	    sprintf(errtext, _("Invalid starting observation '%s'"), stobs);
	    errbox(errtext);
	    return;
	}
	ed0 = strtod(endobs, &test);
	if (strcmp(endobs, test) == 0 || test[0] != '\0' || ed0 < 0) {
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	    errbox(errtext);
	    return;
	}
    }

    if (sd0 > ed0) {
	sprintf(errtext, _("Empty data range '%s - %s'"), stobs, endobs);
	errbox(errtext);
	return;
    }

    colonize_obs(stobs);
    colonize_obs(endobs);

    if (datainfo->pd == 999) { /* panel */
	char unit[8], period[8];

	/* try to infer structure from ending obs */
	if (sscanf(endobs, "%[^:]:%s", unit, period) == 2) { 
	    datainfo->pd = atoi(period);
	    fprintf(stderr, _("Setting data frequency = %d\n"), datainfo->pd);
	} else {
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	    errbox(errtext);
	    return;	    
	}
    }    

    if (datainfo->pd == 1) {
	size_t i, n;
	
	n = strlen(stobs);
	for (i=0; i<n; i++) {
	    if (!isdigit((unsigned char) stobs[i])) {
		sprintf(errtext, _("Invalid starting observation '%s'\n"
			"for data frequency 1"), stobs);
		errbox(errtext);
		return;
	    }
	}
	n = strlen(endobs);
	for (i=0; i<n; i++) {
	    if (!isdigit((unsigned char) endobs[i])) {
		sprintf(errtext, _("Invalid ending observation '%s'\n"
			"for data frequency 1"), endobs);
		errbox(errtext);
		return;
	    }
	}	
    } 
    else if (datainfo->pd != 5 && datainfo->pd != 7) { 
	char year[8], subper[8];

	if (sscanf(stobs, "%[^:]:%s", year, subper) != 2 ||
	    strlen(year) > 4 || atoi(subper) > datainfo->pd ||
	    (datainfo->pd < 10 && strlen(subper) != 1) ||
	    (datainfo->pd >= 10 && strlen(subper) != 2)) {
	    sprintf(errtext, _("Invalid starting observation '%s'\n"
		    "for data frequency %d"), stobs, datainfo->pd);
	    errbox(errtext);
	    return;
	}
	if (sscanf(endobs, "%[^:]:%s", year, subper) != 2 ||
	    strlen(year) > 4 || atoi(subper) > datainfo->pd ||
	    (datainfo->pd < 10 && strlen(subper) != 1) ||
	    (datainfo->pd >= 10 && strlen(subper) != 2)) {
	    sprintf(errtext, _("Invalid ending observation '%s'\n"
		    "for data frequency %d"), endobs, datainfo->pd);
	    errbox(errtext);
	    return;
	}	    
    }

    gtk_widget_destroy(data->dialog); 

    strcpy(datainfo->stobs, stobs);
    strcpy(datainfo->endobs, endobs);
    datainfo->sd0 = sd0;
    datainfo->n = -1;
    datainfo->n = dateton(datainfo->endobs, datainfo) + 1; 

    if (datainfo->n <= 0) {
	errbox("Got zero-length data series");
	return;
    }

    datainfo->v = 2;
    start_new_Z(&Z, datainfo, 0);
    datainfo->markers = 0;

    strcpy(datainfo->varname[1], firstvar);

    show_spreadsheet(datainfo);
}

/* ........................................................... */

void newdata_dialog (gpointer data, guint pd_code, GtkWidget *widget) 
{
    windata_t *wdata = NULL;
    gchar *obsstr = NULL;

    if (pd_code == 0) {
	datainfo->time_series = 0;
	datainfo->pd = 1;
    } else {
	datainfo->time_series = TIME_SERIES;
	datainfo->pd = pd_code;
    }

    switch (pd_code) {
    case 0:
	datainfo->pd = 1;
	obsstr = g_strdup_printf("1 50 %s", _("newvar"));
	break;
    case 1:
	obsstr = g_strdup_printf("1950 2001 %s", _("newvar"));
	break;
    case 4:
	obsstr = g_strdup_printf("1950:1 2001:4 %s", _("newvar"));
	break;
    case 5:
	obsstr = g_strdup_printf("99/01/18 01/03/31 %s", _("newvar"));
	break;
    case 7:
	obsstr = g_strdup_printf("99/01/18 01/03/31 %s", _("newvar"));
	break;
    case 12:
	obsstr = g_strdup_printf("1950:01 2001:12 %s", _("newvar"));
	break;
    case 24:
	obsstr = g_strdup_printf("0:01 0:24 %s", _("newvar"));
	break;
    case 52:
	obsstr = g_strdup_printf("1950:01 2001:52 %s", _("newvar"));
	break;
    }
    edit_dialog (_("gretl: create data set"), 
		 _("Enter start and end obs for new data set\n"
		 "and name of first var to add:"), 
		 obsstr, 
		 _("Apply"), prep_spreadsheet, wdata, 
		 _(" Cancel "), NULL, NULL, 0, 0);
    g_free(obsstr);
}

/* ........................................................... */

void start_panel_dialog (gpointer data, guint u, GtkWidget *widget) 
{
    windata_t *wdata = NULL;

    datainfo->pd = 999;

    edit_dialog (_("gretl: create panel data set"), 
		 _("Enter starting and ending observations and\n"
		 "the name of the first variable to add.\n"
		 "The example below is suitable for 20 units\n"
		 "observed over 10 periods"), 
		 "1.01 10.20 newvar", 
		 _("Apply"), prep_spreadsheet, wdata, 
		 _(" Cancel "), NULL, NULL, 0, 0);
}

/* ........................................................... */

void addvars_dialog (gpointer data, guint add_code, GtkWidget *widget)
{
    simple_selection (_("gretl: data transformations"), 
		      _("Apply"), add_logs_etc, add_code, NULL);    
}

/* ........................................................... */

void destroy_dialog_data (GtkWidget *w, gpointer data) 
{
    GList *list;
    dialog_t *ddata = (dialog_t *) data;

    gtk_main_quit();
    list = ddata->all_buttons;
    while (list != NULL) {
	if (list->data != ddata && list->data) g_free (list->data);
	list = list->next;
    }
    g_list_free (ddata->all_buttons);
    g_free (ddata);
    open_dialog = NULL;
    if (active_edit_id) active_edit_id = NULL;
    if (active_edit_name) active_edit_name = NULL;
    if (active_edit_text) active_edit_text = NULL;
}

/* ........................................................... */

static void dialog_table_setup (dialog_t *dlg, int hsize)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_widget_set_usize(sw, hsize, 200);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg->dialog)->vbox), 
		       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
				    GTK_POLICY_AUTOMATIC,
				    GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(sw), dlg->edit); 
    gtk_widget_show(dlg->edit);
    gtk_widget_show(sw);
}

/* ........................................................... */

static GtkWidget *text_edit_new (int *hsize)
{
    GtkWidget *tbuf;

    tbuf = gtk_text_new(NULL, NULL);

    gtk_text_set_editable(GTK_TEXT(tbuf), TRUE);
    gtk_text_set_word_wrap(GTK_TEXT(tbuf), FALSE);
    *hsize *= gdk_char_width(fixed_font, 'W');
    *hsize += 48;

    return tbuf;
}

/* ........................................................... */

void edit_dialog (char *diagtxt, char *infotxt, char *deftext, 
		  char *oktxt, void (*okfunc)(), void *okptr,
		  char *canceltxt, void (*cancelfunc)(), 
		  void *cancelptr, guint cmdcode, guint varclick)
{
    dialog_t *d, *cancel_d;
    GtkWidget *tempwid;

    if (open_dialog != NULL) {
	gdk_window_raise(open_dialog->window);
	return;
    }

    if ((d = mymalloc(sizeof *d)) == NULL)
 	return;
    d->data = okptr;

    if ((cancel_d = mymalloc(sizeof *cancel_d)) == NULL)
	return;
    cancel_d->data = cancelptr;

    cancel_d->all_buttons = d->all_buttons = NULL;
    d->dialog = cancel_d->dialog = gtk_dialog_new();
    open_dialog = d->dialog;
    d->code = cmdcode;

    gtk_window_set_title (GTK_WINDOW (d->dialog), diagtxt);
    gtk_window_set_policy (GTK_WINDOW (d->dialog), FALSE, FALSE, FALSE);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (d->dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (d->dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (d->dialog), "destroy", 
			GTK_SIGNAL_FUNC (destroy_dialog_data), 
			cancel_d);

    if (cmdcode == NLS) {
	int hsize = 64;
	gchar *lbl;

	lbl = g_strdup_printf("%s\n%s", infotxt,
			      _("(Please refer to Help for guidance)"));
	tempwid = gtk_label_new (lbl);
	gtk_label_set_justify(GTK_LABEL(tempwid), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    tempwid, TRUE, TRUE, 10);
	gtk_widget_show (tempwid);
	g_free(lbl);

	d->edit = cancel_d->edit = text_edit_new (&hsize);
	dialog_table_setup(d, hsize);	
    } else {
	tempwid = gtk_label_new (infotxt);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show (tempwid);

	d->edit = cancel_d->edit = gtk_entry_new ();
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    d->edit, TRUE, TRUE, FALSE);

	/* make the Enter key do the business */
	if (okfunc) 
	    gtk_signal_connect (GTK_OBJECT (d->edit), "activate", 
				GTK_SIGNAL_FUNC (okfunc), (gpointer) d);

	gtk_entry_set_visibility (GTK_ENTRY (d->edit), TRUE);
	if (deftext) {
	    gtk_entry_set_text (GTK_ENTRY (d->edit), deftext);
	    gtk_entry_select_region (GTK_ENTRY (d->edit), 0, strlen (deftext));
	}
	gtk_widget_show (d->edit);
    }

    if (varclick == VARCLICK_INSERT_ID)
	active_edit_id = d->edit; 
    else if (varclick == VARCLICK_INSERT_NAME)
	active_edit_name = d->edit;
    else if (varclick == VARCLICK_INSERT_TEXT)
	active_edit_text = d->edit;

    gtk_widget_grab_focus (d->edit);

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (oktxt);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    if (okfunc) {
	gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			    GTK_SIGNAL_FUNC (okfunc), (gpointer) d);
    }

    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create the "Cancel" button */
    tempwid = gtk_button_new_with_label (canceltxt);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    if (cancelfunc) 
	gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			    GTK_SIGNAL_FUNC(cancelfunc), (gpointer) cancel_d);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
    gtk_widget_show (tempwid);

    /* Create a "Help" button if wanted */
    if (cmdcode && cmdcode != PRINT) {
	tempwid = gtk_button_new_with_label (_("Help"));
	GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			    GTK_SIGNAL_FUNC (context_help), 
			    GINT_TO_POINTER (cmdcode));
	gtk_widget_show (tempwid);
    }

    d->all_buttons = g_list_append (d->all_buttons, d);
    d->all_buttons = g_list_append (d->all_buttons, cancel_d);
    cancel_d->all_buttons = d->all_buttons;

    gtk_widget_show (d->dialog); 
    gtk_main();
} 

#ifdef USE_GNOME

void about_dialog (gpointer data) 
{
    GtkWidget* dlg;
    char const *authors[] = {
	"Allin Cottrell",
	NULL
    };
    const gchar *blurb = N_("An econometrics program for the gnome desktop "
			    "issued under the GNU General Public License.  "
			    "http://gretl.sourceforge.net/");
    gchar *comment = NULL;

#ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	comment = g_strconcat(_(blurb), " ", _("translator_credits"),
			      NULL);
    }
#endif 

    dlg = gnome_about_new("gretl", version_string,
			  "(c) 2000-2003 Allin Cottrell", 
			  authors, (comment != NULL)? comment : _(blurb),
			  gnome_pixmap_file("gretl-logo.xpm") 
			  );

    if (comment != NULL) g_free(comment);

    gnome_dialog_set_parent(GNOME_DIALOG(dlg), GTK_WINDOW(mdata->w));
    gtk_widget_show(dlg);
}

#else /* plain GTK version of About dialog follows */

static int open_xpm (char *filename, GtkWidget *parent, GdkPixmap **pixmap, 
		     GdkBitmap **mask) 
{
    char exfile[MAXLEN];
    GtkStyle *style;

    if (*filename == '\0') return 1;
    strcpy(exfile, paths.gretldir);
    if (exfile[strlen(exfile) - 2] != SLASH)
	strcat(exfile, SLASHSTR);
    strcat(exfile, filename);

    style = gtk_widget_get_style (parent);
    *pixmap = gdk_pixmap_create_from_xpm (parent->window, 
					  mask, 
					  &style->bg[GTK_STATE_NORMAL], 
					  exfile);
    if (*pixmap == NULL) return 0;
    else return 1;
}

void about_dialog (gpointer data) 
{
    GtkWidget *tempwid, *notebook, *box, *label, *view, *vscroll;
    GdkPixmap *logo_pixmap;
    GdkBitmap *logo_mask;
    char *tempstr, *no_gpl, buf[MAXSTR];
    const gchar *tr_credit = "";
    GtkWidget *dialog;
    FILE *fd;

    no_gpl = 
	g_strdup_printf (_("Cannot find the license agreement file COPYING. "
			 "Please make sure it's in %s"), 
			 paths.gretldir);
    dialog = gtk_dialog_new ();
    gtk_window_set_title (GTK_WINDOW (dialog), _("About gretl"));
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
    gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
    gtk_signal_connect_object (GTK_OBJECT (dialog), 
			       "delete_event", GTK_SIGNAL_FUNC 
			       (gtk_widget_destroy), GTK_OBJECT (dialog));
    gtk_widget_realize (dialog);
      
    notebook = gtk_notebook_new ();
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			notebook, TRUE, TRUE, 0);
    gtk_widget_show (notebook);
   
    box = gtk_vbox_new (TRUE, 5);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);
   
    if (open_xpm ("gretl-logo.xpm", mdata->w, &logo_pixmap, &logo_mask)) {
	tempwid = gtk_pixmap_new (logo_pixmap, logo_mask);
	gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
	gtk_widget_show (tempwid);
    }

#ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	tr_credit = _("translator_credits");
    }
#endif  

    tempstr = g_strdup_printf ("gretl, version %s\n"
			       "Copyright (C) 2000-2001 Allin Cottrell "
			       "<cottrell@wfu.edu>\nHomepage: "
			       "http://gretl.sourceforge.net/\n%s",
			       version_string, tr_credit);
    tempwid = gtk_label_new (tempstr);
    g_free (tempstr);
    gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
    gtk_widget_show (tempwid);
   
    label = gtk_label_new (_("About"));
    gtk_widget_show (label);
   
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    box = gtk_vbox_new (FALSE, 5);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    tempwid = gtk_table_new (1, 2, FALSE);
    gtk_box_pack_start (GTK_BOX (box), tempwid, TRUE, TRUE, 0);
    gtk_widget_show (tempwid);

    view = gtk_text_new (NULL, NULL);
    gtk_text_set_editable (GTK_TEXT (view), FALSE);
    gtk_text_set_word_wrap (GTK_TEXT (view), TRUE);
    gtk_table_attach (GTK_TABLE (tempwid), view, 0, 1, 0, 1,
		      GTK_FILL | GTK_EXPAND, GTK_FILL | 
		      GTK_EXPAND | GTK_SHRINK, 0, 0);
    gtk_widget_show (view);

    vscroll = gtk_vscrollbar_new (GTK_TEXT (view)->vadj);
    gtk_table_attach (GTK_TABLE (tempwid), vscroll, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0);
    gtk_widget_show (vscroll);

    label = gtk_label_new (_("License Agreement"));
    gtk_widget_show (label);
   
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    tempwid = gtk_button_new_with_label (_("  Close  "));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			tempwid, FALSE, FALSE, 0);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tempstr = g_strdup_printf("%s/COPYING", paths.gretldir);
    if ((fd = fopen (tempstr, "r")) == NULL) {
	gtk_text_insert (GTK_TEXT (view), NULL, NULL, NULL, 
			 no_gpl, strlen (no_gpl));
	gtk_widget_show (dialog);
	g_free (tempstr);
	return;
    }
    g_free(tempstr);
   
    memset (buf, 0, sizeof (buf));
    while (fread (buf, 1, sizeof (buf) - 1, fd)) {
	gtk_text_insert (GTK_TEXT (view), 
			 fixed_font, NULL, NULL, buf, strlen (buf));
	memset (buf, 0, sizeof (buf));
    }
    fclose (fd);
    gtk_widget_show(dialog);
    g_free(no_gpl);
}         
#endif /* not GNOME */

/* ........................................................... */

void menu_exit_check (GtkWidget *w, gpointer data)
{
    int ret = exit_check(w, NULL, data);

    if (ret == FALSE) gtk_main_quit();
}

/* ........................................................... */

int work_done (void)
     /* See whether user has done any work, to determine whether or
	not to offer the option of saving commands/output.  Merely
	running a script, or opening a data file, or a few other
	trivial actions, do not count as "work done". */
{
    FILE *fp;
    char line[MAXLEN];
    int work = 0;
    
    fp = fopen(cmdfile, "r");
    if (fp == NULL) return -1;
    while (fgets(line, MAXLEN-1, fp)) {
	if (strlen(line) > 2 && 
	    strncmp(line, "run ", 4) &&
	    strncmp(line, "open", 4) &&
	    strncmp(line, "help", 4) &&
	    strncmp(line, "impo", 4) &&
	    strncmp(line, "info", 4) &&
	    strncmp(line, "labe", 4) &&
	    strncmp(line, "list", 4) &&
	    strncmp(line, "quit", 4)) {
	    work = 1;
	    break;
	}
    }
    fclose(fp);
    return work;
}

/* ......................................................... */

static void save_data_callback (void)
{
    file_save(NULL, SAVE_DATA, NULL);
    if (data_status & MODIFIED_DATA)
	data_status ^= MODIFIED_DATA;
    /* FIXME: need to do more here */
}

#ifdef USE_GNOME
/* ......................................................... */

int yes_no_dialog (char *title, char *message, int cancel)
{
    GtkWidget *dialog, *label;
    int button;

    if (cancel)
	dialog = gnome_dialog_new (title,
				   GNOME_STOCK_BUTTON_YES,
				   GNOME_STOCK_BUTTON_NO,
				   GNOME_STOCK_BUTTON_CANCEL,
				   NULL);
    else
	dialog = gnome_dialog_new (title,
				   GNOME_STOCK_BUTTON_YES,
				   GNOME_STOCK_BUTTON_NO,
				   NULL);

    gnome_dialog_set_parent (GNOME_DIALOG (dialog), 
			     GTK_WINDOW(mdata->w));

    label = gtk_label_new (message);
    gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (GNOME_DIALOG (dialog)->vbox), label, 
			TRUE, TRUE, 0);

    button = gnome_dialog_run_and_close (GNOME_DIALOG (dialog));

    if (button == 0) return GRETL_YES;
    if (button == 1) return GRETL_NO;
    if (button == 2) return GRETL_CANCEL;

    return GRETL_CANCEL;
}

#else /* not USE_GNOME */

struct yes_no_data {
    GtkWidget *dialog;
    int *ret;
    int button;
};

static void yes_no_callback (GtkWidget *w, gpointer data)
{
    struct yes_no_data *mydata = data;

    *(mydata->ret) = mydata->button;
    gtk_main_quit();
    gtk_widget_destroy(mydata->dialog);
}

/* ......................................................... */

gint yes_no_dialog (char *title, char *msg, int cancel)
{
   GtkWidget *tempwid, *dialog;
   int ret;
   struct yes_no_data yesdata, nodata, canceldata;

   dialog = gtk_dialog_new();

   yesdata.dialog = nodata.dialog = canceldata.dialog 
       = dialog;
   yesdata.ret = nodata.ret = canceldata.ret = &ret; 
   yesdata.button = GRETL_YES;
   nodata.button = GRETL_NO;
   canceldata.button = GRETL_CANCEL;
   
   gtk_grab_add (dialog);
   gtk_window_set_title (GTK_WINDOW (dialog), title);
   gtk_window_set_policy (GTK_WINDOW (dialog), FALSE, FALSE, FALSE);
   gtk_container_border_width 
       (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), 10);
   gtk_container_border_width 
       (GTK_CONTAINER (GTK_DIALOG (dialog)->action_area), 5);
   gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
   gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->action_area), 15);
   gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
   gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

   tempwid = gtk_label_new (msg);
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), tempwid, 
		       TRUE, TRUE, FALSE);
   gtk_widget_show(tempwid);

   /* "Yes" button */
   tempwid = gtk_button_new_with_label (_("Yes"));
   GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		       tempwid, TRUE, TRUE, TRUE);  
   gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC (yes_no_callback), &yesdata);
   gtk_widget_grab_default (tempwid);
   gtk_widget_show (tempwid);

   /* "No" button */
   tempwid = gtk_button_new_with_label (_("No"));
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		       tempwid, TRUE, TRUE, TRUE); 
   gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC (yes_no_callback), &nodata);
   gtk_widget_show (tempwid);

   /* Cancel button -- if wanted */
   if (cancel) {
       tempwid = gtk_button_new_with_label (_("Cancel"));
       gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			   tempwid, TRUE, TRUE, TRUE); 
       gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			   GTK_SIGNAL_FUNC (yes_no_callback), &canceldata);
       gtk_widget_show (tempwid);
   }

   gtk_widget_show (dialog);
   gtk_main();
   return ret;
}

#endif /* plain GTK */

/* ........................................................... */

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data) 
{
    char fname[MAXLEN];
    int button;
    extern int replay; /* lib.c */

    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");
    dump_cmd_stack(fname);

    /* FIXME: should make both save_session_callback() and
       save_data_callback() blocking functions */

    if (!expert && !replay && 
	(session_changed(0) || (work_done() && !session_saved))) {
	button = yes_no_dialog ("gretl", 		      
				_("Do you want to save the commands and\n"
				"output from this gretl session?"), 1);
	if (button == GRETL_YES) {
	    save_session_callback(NULL, 1, NULL);
	    return TRUE; /* bodge */
	}
	/* button -1 = wm close */
	else if (button == GRETL_CANCEL || button == -1) return TRUE;
	/* else button = GRETL_NO: so fall through */
    }

    if (data_status & MODIFIED_DATA) {
	button = yes_no_dialog ("gretl", 
				_("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (button == GRETL_YES) {
	    save_data_callback();
	    return TRUE; 
	}
	else if (button == GRETL_CANCEL || button == -1) return TRUE;
    }    

    write_rc();
    return FALSE;
}

typedef struct {
    GtkWidget *space_button;
    GtkWidget *point_button;
    gint delim;
    gint decpoint;
} csv_stuff;

#ifdef ENABLE_NLS
static void set_dec (GtkWidget *w, gpointer p)
{
    gint i;
    csv_stuff *csv = (csv_stuff *) p;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
	csv->decpoint = i;
	if (csv->decpoint == ',' && csv->delim == ',') {
	    csv->delim = ' ';
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (csv->space_button), 
					  TRUE);
	}
    }
}
#endif

static void set_delim (GtkWidget *w, gpointer p)
{
    gint i;
    csv_stuff *csv = (csv_stuff *) p;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
	csv->delim = i;
	if (csv->point_button != NULL && 
	    csv->delim == ',' && csv->decpoint == ',') {
	    csv->decpoint = '.';
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (csv->point_button), 
					  TRUE);
	}
    }
}

static void really_set_csv_stuff (GtkWidget *w, gpointer p)
{
    csv_stuff *stuff = (csv_stuff *) p;

    datainfo->delim = stuff->delim;
    datainfo->decpoint = stuff->decpoint;
}

static void destroy_delim_dialog (GtkWidget *w, gint *p)
{
    free(p);
    gtk_main_quit();
}

void delimiter_dialog (void)
{
    GtkWidget *dialog, *tempwid, *button;
    GSList *group;
    csv_stuff *csvptr = NULL;

    csvptr = mymalloc(sizeof *csvptr);
    if (csvptr == NULL) return;
    csvptr->delim = datainfo->delim;
    csvptr->decpoint = '.';
    csvptr->point_button = NULL;

    dialog = gtk_dialog_new();

    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: data delimiter"));
    gtk_window_set_policy (GTK_WINDOW (dialog), FALSE, FALSE, FALSE);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (dialog), "destroy", 
			destroy_delim_dialog, csvptr);

    tempwid = gtk_label_new (_("separator for data columns:"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_widget_show(tempwid);

    /* comma separator */
    button = gtk_radio_button_new_with_label (NULL, _("comma (,)"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (csvptr->delim == ',')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(','));
    gtk_widget_show (button);

    /* space separator */
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("space"));
    csvptr->space_button = button;
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (csvptr->delim == ' ')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(' '));    
    gtk_widget_show (button);

    /* tab separator */
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("tab"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (csvptr->delim == '\t')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER('\t'));    
    gtk_widget_show (button);

#ifdef ENABLE_NLS
    if (',' == get_local_decpoint()) {
	GSList *decgroup;

	tempwid = gtk_hseparator_new();
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show(tempwid);
	
	tempwid = gtk_label_new (_("decimal point character:"));
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show(tempwid);
 
	/* period decpoint */
	button = gtk_radio_button_new_with_label (NULL, _("period (.)"));
	csvptr->point_button = button;
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    button, TRUE, TRUE, FALSE);
	if (csvptr->decpoint == '.')
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
	gtk_signal_connect(GTK_OBJECT(button), "clicked",
			   GTK_SIGNAL_FUNC(set_dec), csvptr);
	gtk_object_set_data(GTK_OBJECT(button), "action", 
			    GINT_TO_POINTER('.'));
	gtk_widget_show (button);

	/* comma decpoint */
	decgroup = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
	button = gtk_radio_button_new_with_label(decgroup, _("comma (,)"));
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    button, TRUE, TRUE, FALSE);
	if (csvptr->decpoint == ',')
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
	gtk_signal_connect(GTK_OBJECT(button), "clicked",
			   GTK_SIGNAL_FUNC(set_dec), csvptr);
	gtk_object_set_data(GTK_OBJECT(button), "action", 
			    GINT_TO_POINTER(','));    
	gtk_widget_show (button);
    }
#endif

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (_("OK"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(really_set_csv_stuff), csvptr);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    gtk_widget_show (dialog);

    gtk_main();
}

struct varinfo_settings {
    GtkWidget *dlg;
    GtkWidget *name_entry;
    GtkWidget *label_entry;
    GtkWidget *display_name_entry;
    int varnum;
};

static void really_set_variable_info (GtkWidget *w, 
				      struct varinfo_settings *vset)
{
    const char *edttext;
    int v = vset->varnum;
    int changed = 0;
    int gui_changed = 0;

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->name_entry));
    if (validate_varname(edttext)) {
	return;
    } else {
	if (strcmp(datainfo->varname[v], edttext)) {
	    strcpy(datainfo->varname[v], edttext);
	    gui_changed = 1;
	}
    }

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->label_entry));
    if (strcmp(VARLABEL(datainfo, v), edttext)) {
	*VARLABEL(datainfo, v) = 0;
	strncat(VARLABEL(datainfo, v), edttext, MAXLABEL - 1);
	changed = 1;
	gui_changed = 1;
    }

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->display_name_entry));
    if (strcmp(DISPLAYNAME(datainfo, v), edttext)) {
	*DISPLAYNAME(datainfo, v) = 0;
	strncat(DISPLAYNAME(datainfo, v), edttext, MAXDISP - 1);
	changed = 1;
    }

    if (changed) {
	sprintf(line, "label %s -d \"%s\" -n \"%s\"", datainfo->varname[v],
		VARLABEL(datainfo, v), DISPLAYNAME(datainfo, v));
	verify_and_record_command(line);
    }

    if (gui_changed)
	populate_main_varlist();

    gtk_widget_destroy(vset->dlg);
}

static void free_vsettings (GtkWidget *w, 
			    struct varinfo_settings *vset)
{
    free(vset);
}

void varinfo_dialog (int varnum)
{
    GtkWidget *tempwid, *hbox;
    struct varinfo_settings *vset;

    vset = mymalloc(sizeof *vset);
    if (vset == NULL) return;

    vset->varnum = varnum;
    vset->dlg = gtk_dialog_new();

    gtk_signal_connect (GTK_OBJECT(vset->dlg), "destroy", 
			GTK_SIGNAL_FUNC(free_vsettings), vset);

    gtk_window_set_title(GTK_WINDOW(vset->dlg), _("gretl: variable attributes"));
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (vset->dlg)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (vset->dlg)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (vset->dlg)->vbox), 5);
    gtk_window_set_position (GTK_WINDOW (vset->dlg), GTK_WIN_POS_MOUSE);

    /* read/set name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("name of variable:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

    vset->name_entry = gtk_entry_new_with_max_length(8);
    gtk_entry_set_text(GTK_ENTRY(vset->name_entry), 
		       datainfo->varname[varnum]);
    gtk_box_pack_start(GTK_BOX(hbox), 
		       vset->name_entry, FALSE, FALSE, 0);
    gtk_widget_show(vset->name_entry); 

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox); 

    /* read/set descriptive string */
    hbox = gtk_hbox_new(FALSE, 0);
    tempwid = gtk_label_new (_("description:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    hbox = gtk_hbox_new(FALSE, 0);
    vset->label_entry = gtk_entry_new_with_max_length(MAXLABEL-1);
    gtk_entry_set_text(GTK_ENTRY(vset->label_entry), 
		       VARLABEL(datainfo, varnum));
    gtk_box_pack_start(GTK_BOX(hbox), vset->label_entry, TRUE, TRUE, 0);
    gtk_widget_show(vset->label_entry);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox); 

    /* read/set display name */
    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("display name (shown in graphs):"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

    vset->display_name_entry = gtk_entry_new_with_max_length(MAXDISP-1);
    gtk_entry_set_text(GTK_ENTRY(vset->display_name_entry), 
		       DISPLAYNAME(datainfo, varnum));
    gtk_box_pack_start(GTK_BOX(hbox), 
		       vset->display_name_entry, FALSE, FALSE, 0);
    gtk_widget_show(vset->display_name_entry); 

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox); 

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (_("OK"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (vset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(really_set_variable_info), vset);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(delete_widget), 
			vset->dlg);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* And a Cancel button */
    tempwid = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(vset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(delete_widget), vset->dlg);
    gtk_widget_show (tempwid);

    gtk_widget_show (vset->dlg);
}

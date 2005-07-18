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

/* gpt_dialog.c for gretl -- GUI gnuplot controller dialog */

#include "gretl.h"
#include "gpt_control.h"
#include "session.h"
#include "dlgutils.h"

#include "../pixmaps/mouse.xpm"

struct gpt_titles_t {
    char *description; /* How the field will show up in the options dialog */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
};

struct gpt_range_t {
    gint ID;
    GtkWidget *isauto;
    GtkWidget *min;
    GtkWidget *max;
};

GtkWidget *linetitle[MAX_PLOT_LINES];
GtkWidget *stylecombo[MAX_PLOT_LINES];
GtkWidget *yaxiscombo[MAX_PLOT_LINES];
GtkWidget *linescale[MAX_PLOT_LINES];

GtkWidget *labeltext[MAX_PLOT_LABELS];
GtkWidget *labeljust[MAX_PLOT_LABELS];
GtkWidget *labelpos[MAX_PLOT_LABELS];

static GtkWidget *gpt_control;
static GtkWidget *keycombo;
static GtkWidget *termcombo;
static GtkWidget *fitline_check;
static GtkWidget *border_check;
static GtkWidget *y2_check;
static GtkWidget *ttfcombo;
static GtkWidget *ttfspin;

GtkWidget *filesavebutton;

#define MAX_AXES 3

struct gpt_range_t axis_range[MAX_AXES];

#define NTITLES 4
#define PLOT_LABEL_POS_LEN 32

struct gpt_titles_t gpt_titles[] = {
    { N_("Title of plot"),  0, NULL },
    { N_("Title for axis"), 1, NULL },
    { N_("Title for axis"), 2, NULL },
    { N_("Title for axis"), 3, NULL },
};

#define frequency_plot(s)       (s->code == PLOT_FREQ_SIMPLE || \
			         s->code == PLOT_FREQ_NORMAL || \
			         s->code == PLOT_FREQ_GAMMA || \
                                 s->code == PLOT_KERNEL)

#define no_edit_lines(s) (s->code == PLOT_FORECAST || \
                          s->code == PLOT_FREQ_SIMPLE || \
                          s->code == PLOT_FREQ_NORMAL || \
                          s->code == PLOT_FREQ_GAMMA || \
                          s->code == PLOT_GARCH || \
                          s->code == PLOT_KERNEL || \
                          s->code == PLOT_CUSUM)

static const char *get_font_filename (const char *showname);


static void close_plot_controller (GtkWidget *widget, gpointer data) 
{
    GPT_SPEC *spec = (GPT_SPEC *) data;
    png_plot *plot = (png_plot *) spec->ptr;

    gpt_control = NULL;

    if (plot != NULL) { /* PNG plot window open */
	plot_remove_controller(plot);
    } else {
	free_plotspec(spec); 
    }
} 

static void flip_manual_range (GtkWidget *widget, gpointer data)
{
    gint axis = GPOINTER_TO_INT(data);

    if (GTK_TOGGLE_BUTTON (axis_range[axis].isauto)->active) {
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].min), FALSE);
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].max), FALSE);
    } else {
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].min), TRUE);
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].max), TRUE);
    }
}

/* Take text from a gtkentry and write to gnuplot spec string.
   Under gtk2, the entries will be in utf-8, and have to be converted
   to the locale for use with gnuplot.
*/

static void entry_to_gp_string (GtkWidget *w, char *targ, size_t n)
{
    const gchar *wstr;
    
    *targ = '\0';
    g_return_if_fail(GTK_IS_ENTRY(w));
    wstr = gtk_entry_get_text(GTK_ENTRY(w));

    if (wstr != NULL && *wstr != '\0') {
#ifndef OLD_GTK
	gchar *trstr = force_locale_from_utf8(wstr);

	if (trstr != NULL) {
	    strncat(targ, trstr, n-1);
	    g_free(trstr);
	}
#else
	strncat(targ, wstr, n-1);
#endif /* OLD_GTK */
    }
}

/* take a double (which might be NA) and format it for
   a gtkentry widget */

static void double_to_gp_entry (double x, GtkWidget *w)
{
    if (w != NULL && GTK_IS_ENTRY(w)) {
	gchar *numstr;

	if (na(x)) {
	    numstr = g_strdup("*");
	} else {
	    numstr = g_strdup_printf("%g", x);
	}
	gtk_entry_set_text(GTK_ENTRY(w), numstr);
	g_free(numstr);
    }
}

/* read a double from a gtkentry, with error checking */

static double entry_to_gp_double (GtkWidget *w)
{
    double ret = NADBL;

    if (w != NULL && GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));

	if (s != NULL && *s != '\0') {
#ifdef ENABLE_NLS
	    gchar *tmp = g_strdup(s);

	    charsub(tmp, ',', '.');
	    setlocale(LC_NUMERIC, "C");
#else
	    const gchar *tmp = s;
#endif
	    if (check_atof(tmp)) {
		errbox(get_gretl_errmsg());
	    } else {
		ret = atof(tmp);
	    }
#ifdef ENABLE_NLS
	    setlocale(LC_NUMERIC, "");
	    g_free(tmp);
#endif
	}
    }

    return ret;
}

/* Take text from a gnuplot spec string and put it into a gtkentry.
   Under gtk2, we have to ensure that the text is put into utf-8.
*/

static void gp_string_to_entry (GtkWidget *w, const char *str)
{
#ifndef OLD_GTK
# ifdef ENABLE_NLS
    int l2 = use_latin_2();
# else
    int l2 = 0;
# endif /* ENABLE_NLS */
    gchar *trstr;

    if (l2) {
	char lstr[MAXTITLE];
	
	sprint_html_to_l2(lstr, str);
	trstr = my_locale_to_utf8(lstr);
    } else {
	trstr = my_locale_to_utf8(str);
    }

    if (trstr != NULL) {
	gtk_entry_set_text(GTK_ENTRY(w), trstr);
	g_free(trstr);
    }    
#else 
    gtk_entry_set_text(GTK_ENTRY(w), str);
#endif /* OLD_GTK */
}

static int
get_label_pos_from_entry (GtkWidget *w, double *pos)
{
    int err = 0;

    if (GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
	int chk;
    
	chk = sscanf(s, "%lf %lf", &pos[0], &pos[1]);
	if (chk != 2) {
	    errbox(_("Invalid label position, must be X Y"));
	    gtk_editable_select_region(GTK_EDITABLE(w), 0, strlen(s));
	    pos[0] = pos[1] = NADBL;
	    err = 1;
	} 
    } else {
	err = 1;
    }

    return err;
}

static int validate_range (double *r)
{
    int err = 0;

    if (na(r[0]) || na(r[1])) {
	r[0] = r[1] = NADBL;
	err = 1;
    } else if (r[1] <= r[0]) {
	r[0] = r[1] = NADBL;
	err = 1;
    }

    return err;
}

static void apply_gpt_changes (GtkWidget *widget, GPT_SPEC *spec) 
{
    const gchar *yaxis;
    int supress_y2 = 0;
    int i, k, save = 0;
    int err = 0;

    /* entry_to_gp_string translates from utf-8 to the locale, if
       using NLS */

    if (widget == filesavebutton) {
	entry_to_gp_string(GTK_COMBO(termcombo)->entry, spec->termtype, 
			   sizeof spec->termtype);
	if (strcmp(spec->termtype, "screen")) save = 1;
    }
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].widget != NULL) {
	    entry_to_gp_string(gpt_titles[i].widget, spec->titles[i], 
			       sizeof spec->titles[0]);
	}
    }

    entry_to_gp_string(GTK_COMBO(keycombo)->entry, spec->keyspec, 
		       sizeof spec->keyspec);

    spec->flags &= ~GPTSPEC_Y2AXIS;

    if (y2_check != NULL) {
	if (GTK_TOGGLE_BUTTON(y2_check)->active) {
	    supress_y2 = 1;
	} 
    } 

    if (!no_edit_lines(spec)) {    
	for (i=0; i<spec->n_lines; i++) {
	    spec->lines[i].yaxis = 1;
	    if (supress_y2) {
		continue;
	    }
	    yaxis = 
		gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(yaxiscombo[i])->entry));
	    if (yaxis != NULL && *yaxis && !strcmp(yaxis, "right"))
		spec->lines[i].yaxis = 2;	
	    if (spec->lines[i].yaxis == 2) {
		spec->flags |= GPTSPEC_Y2AXIS;
	    }
	}
    }

    if (spec->code == PLOT_REGULAR) {
	k = (spec->flags & GPTSPEC_Y2AXIS)? 3 : 2;
	for (i=0; i<k; i++) {
	    if (axis_range[i].isauto != NULL) {
		if (GTK_TOGGLE_BUTTON(axis_range[i].isauto)->active) {
		    spec->range[i][0] = NADBL;
		    spec->range[i][1] = NADBL;
		} else {
		    spec->range[i][0] = entry_to_gp_double(axis_range[i].min);
		    spec->range[i][1] = entry_to_gp_double(axis_range[i].max);
		    err = validate_range(spec->range[i]);
		}
	    }
	}
    }

    if (!err && !no_edit_lines(spec)) {   
	for (i=0; i<spec->n_lines; i++) {
	    entry_to_gp_string(GTK_COMBO(stylecombo[i])->entry, 
			       spec->lines[i].style, 
			       sizeof spec->lines[0].style);
	    entry_to_gp_string(linetitle[i], 
			       spec->lines[i].title, 
			       sizeof spec->lines[0].title);
	    entry_to_gp_string(linescale[i], 
			       spec->lines[i].scale, 
			       sizeof spec->lines[0].scale);
	}
    }

    for (i=0; i<MAX_PLOT_LABELS && !err; i++) {
#ifdef OLD_GTK
	GtkWidget *active_item;
	int opt;
#endif

	entry_to_gp_string(labeltext[i], spec->labels[i].text, 
			   sizeof spec->labels[0].text);
	if (string_is_blank(spec->labels[i].text)) {
	    continue;
	}
	err = get_label_pos_from_entry(labelpos[i], spec->labels[i].pos);
	if (err) {
	    break;
	}
#ifdef OLD_GTK
	active_item = GTK_OPTION_MENU(labeljust[i])->menu_item;
	opt = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(active_item), 
						  "option"));
	spec->labels[i].just = opt;
#else
	spec->labels[i].just = 
	    gtk_option_menu_get_history(GTK_OPTION_MENU(labeljust[i]));
#endif
    } 

    if (!err && border_check != NULL) {
	if (GTK_TOGGLE_BUTTON(border_check)->active) {
	    spec->flags &= ~GPTSPEC_BORDER_HIDDEN;
	} else {
	    spec->flags |= GPTSPEC_BORDER_HIDDEN;
	}
    } 

    if (!err && fitline_check != NULL) {
	if (GTK_TOGGLE_BUTTON(fitline_check)->active) {
	    spec->flags |= GPTSPEC_OLS_HIDDEN;
	} else {
	    spec->flags &= ~GPTSPEC_OLS_HIDDEN;
	}
    }

    if (!err && ttfcombo != NULL && ttfspin != NULL) {
	const gchar *tmp = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry));
#ifdef OLD_GTK
	int ptsize = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ttfspin));
#else
	int ptsize = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(ttfspin));
#endif

	if (tmp != NULL && *tmp != '\0') {
	    const char *fname = get_font_filename(tmp);
	    char pngfont[32];

	    if (fname != NULL && ptsize > 5 && ptsize < 25) {
		sprintf(pngfont, "%s %d", fname, ptsize);
	    } else {
		*pngfont = '\0';
	    }
	    set_gretl_png_font(pngfont, &paths);
	}
    }

    if (!err) {
	if (save) { /* do something other than a screen graph? */
	    file_selector(_("Save gnuplot graph"), SAVE_GNUPLOT, spec);
	} else { 
	    png_plot *plot = (png_plot *) spec->ptr;

	    set_plot_has_y2_axis(plot, spec->flags & GPTSPEC_Y2AXIS);
	    redisplay_edited_png(plot);
	}
	session_changed(1);
    }
}

static void set_keyspec_sensitivity (GPT_SPEC *spec)
{
    if (!no_edit_lines(spec)) {
	int i; 
	const char *p;

	for (i=0; i<spec->n_lines; i++) {
	    p = gtk_entry_get_text(GTK_ENTRY(linetitle[i]));
	    if (p != NULL && *p != 0) {
		gtk_widget_set_sensitive(keycombo, TRUE);
		return;
	    }
	}
    }
    gtk_widget_set_sensitive(keycombo, FALSE);
}

#define TAB_MAIN_COLS 3

struct font_info {
    const char *fname;
    const char *showname;
};

static struct font_info ttf_fonts[] = {
    { "arial", "Arial", },
    { "georgia", "Georgia", },
#ifndef G_OS_WIN32
    { "luxirr", "Luxi Serif" },
    { "luxisr", "Luxi Sans" },
    { "Vera", "Vera" },
#endif
    { "tahoma", "Tahoma" },
    { "trebuc", "Trebuchet" },
    { "verdana", "Verdana" }
};

static const char *get_font_filename (const char *showname)
{
    int i, nfonts = sizeof ttf_fonts / sizeof ttf_fonts[0];

    for (i=0; i<nfonts; i++) {
	if (!strcmp(ttf_fonts[i].showname, showname)) {
	    return ttf_fonts[i].fname;
	}
    }
    return NULL;
}

#ifndef G_OS_WIN32

static int font_is_ok (const char *fname)
{
    char cmd[64];
    int err;

    sprintf(cmd, "set term png font %s 10", fname);
    err = gnuplot_test_command(cmd);

    return err == 0;
}

#endif

static struct font_info *get_gnuplot_ttf_list (int *nf)
{
    static struct font_info *retlist = NULL;
    static int goodfonts = -1;

    if (goodfonts >= 0) {
	*nf = goodfonts;
    } else {
	int i, j, nfonts = sizeof ttf_fonts / sizeof ttf_fonts[0];

	retlist = malloc(nfonts * sizeof *retlist);
	if (retlist == NULL) return NULL;

	j = 0;
	for (i=0; i<nfonts; i++) {
#ifdef G_OS_WIN32
	    retlist[j++] = ttf_fonts[i];
#else
	    if (font_is_ok(ttf_fonts[i].fname)) {
		retlist[j++] = ttf_fonts[i];
	    }
#endif
	}
	goodfonts = j;
	*nf = goodfonts;
    }

    return retlist;
}

static int font_match (const char *ttfname, const char *pngfont)
{
    return !strncmp(ttfname, pngfont, strlen(ttfname));
}

static int get_point_size (const char *font)
{
    int pts;

    if (sscanf(font, "%*s %d\n", &pts) == 1) {
	return pts;
    } else {
	return 10;
    }
}

static void strip_lr (gchar *txt)
{
    gchar test[16];
    gchar *p;

    sprintf(test, "(%s)", I_("left"));
    p = strstr(txt, test);
    if (p != NULL) {
	*p = '\0';
    } else {
	sprintf(test, "(%s)", I_("right"));
	p = strstr(txt, test);
	if (p != NULL) {
	   *p = '\0';
	}
    } 
}

static void toggle_axis_selection (GtkWidget *w, GPT_SPEC *spec)
{
    int no_y2 = GTK_TOGGLE_BUTTON(w)->active;
    int i;

    for (i=0; i<spec->n_lines; i++) {
	if (i >= MAX_PLOT_LINES) {
	    break;
	}
	gtk_widget_set_sensitive(yaxiscombo[i], !no_y2);
    }
}

static void gpt_tab_main (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl;
    int i, tbl_len;
    GList *keypos_list = NULL;

    gchar *keypos[] = {
	"left top",
	"right top",
	"left bottom",
	"right bottom",
	"outside",
	"none",
	NULL
    };

    for (i=0; keypos[i] != NULL; i++) {
	keypos_list = g_list_append(keypos_list, keypos[i]);
    }
   
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);
    
    label = gtk_label_new(_("Main"));
    gtk_widget_show(label);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, TAB_MAIN_COLS, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].tab == 0) {
	    GtkWidget *entry;

	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, TAB_MAIN_COLS);
	    label = gtk_label_new(_(gpt_titles[i].description));
	    gtk_table_attach_defaults(GTK_TABLE (tbl), 
				      label, 0, 1, tbl_len-1, tbl_len);
	    gtk_widget_show(label);

	    entry = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      entry, 1, TAB_MAIN_COLS, 
				      tbl_len-1, tbl_len);
				      
            if (spec->titles[i] != NULL && *spec->titles[i] != '\0') {
		gp_string_to_entry(entry, spec->titles[i]);
            }		

	    g_signal_connect(G_OBJECT(entry), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     spec);

	    gtk_widget_show(entry);
	    gpt_titles[i].widget = entry;
	}
    }

    /* specify position of plot key or legend */
    tbl_len++;
    label = gtk_label_new(_("key position"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(label);

    keycombo = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      keycombo, 1, TAB_MAIN_COLS, tbl_len-1, tbl_len);
    gtk_combo_set_popdown_strings(GTK_COMBO(keycombo), keypos_list); 
    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(keycombo)->entry), spec->keyspec);
    gtk_widget_show (keycombo);	

    /* give option of removing top & right border */
    if (!(spec->flags & GPTSPEC_Y2AXIS)) { 
	y2_check = NULL;
	tbl_len++;
	border_check = gtk_check_button_new_with_label(_("Show full border"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  border_check, 0, TAB_MAIN_COLS, 
				  tbl_len-1, tbl_len);
	if (!(spec->flags & GPTSPEC_BORDER_HIDDEN)) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(border_check),
					 TRUE);
	}	
	gtk_widget_show(border_check);
    } else {
	border_check = NULL;
	tbl_len++;
	y2_check = gtk_check_button_new_with_label(_("Use only one y axis"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  y2_check, 0, TAB_MAIN_COLS, 
				  tbl_len-1, tbl_len);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(y2_check),
				     FALSE);
	g_signal_connect(G_OBJECT(y2_check), "clicked", 
			 G_CALLBACK(toggle_axis_selection), spec);
	gtk_widget_show(y2_check);
    }

    /* give option of removing an auto-fitted line */
    if (spec->flags & GPTSPEC_AUTO_OLS) { 
	tbl_len++;
	fitline_check = gtk_check_button_new_with_label(_("Hide fitted line"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  fitline_check, 0, TAB_MAIN_COLS, 
				  tbl_len-1, tbl_len);
	if (spec->flags & GPTSPEC_OLS_HIDDEN) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(fitline_check),
					 TRUE);
	}	
	gtk_widget_show(fitline_check);
    } else {
	fitline_check = NULL;
    }

    /* set TT font (if gnuplot uses libgd and freetype) */
    if (gnuplot_has_ttf()) {
	GtkWidget *ebox, *hsep;
#ifdef OLD_GTK
	GtkObject *adj;
#endif
	GList *fontnames = NULL;
	struct font_info *ttflist;
	const char *default_font = NULL;
	int nfonts;

	ttflist = get_gnuplot_ttf_list(&nfonts);

	for (i=0; i<nfonts; i++) {
	    fontnames = g_list_append(fontnames, (gpointer) ttflist[i].showname);
	    if (font_match(ttflist[i].fname, gretl_png_font())) {
		default_font = ttflist[i].showname;
	    }
	}

	fontnames = g_list_append(fontnames, _("None"));
	if (default_font == NULL) {
	    default_font = _("None");
	}

	/* first a separator */
	tbl_len++;
	hsep = gtk_hseparator_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
				  tbl_len-1, tbl_len);  
	gtk_widget_show(hsep);

	tbl_len++;
	ebox = gtk_event_box_new();
	label = gtk_label_new(_("TrueType font"));
	gtk_container_add(GTK_CONTAINER(ebox), label);
	gtk_table_attach_defaults(GTK_TABLE (tbl), ebox, 0, 1, 
				  tbl_len-1, tbl_len);
	gtk_widget_show(label);
	gtk_widget_show(ebox);

	ttfcombo = gtk_combo_new();
	gtk_entry_set_max_length(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry), 15);
	gtk_table_attach_defaults(GTK_TABLE(tbl), ttfcombo, 1, 2, 
				  tbl_len-1, tbl_len);
	gtk_combo_set_popdown_strings(GTK_COMBO(ttfcombo), fontnames); 
	gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry), default_font);
#ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(GTK_COMBO(ttfcombo)->entry), "activate", 
			   GTK_SIGNAL_FUNC(apply_gpt_changes), 
			   spec);
#else
	gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry), 15);
	g_signal_connect(G_OBJECT(GTK_COMBO(ttfcombo)->entry), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
#endif
	gtk_widget_show (ttfcombo);

#ifdef OLD_GTK
	adj = gtk_adjustment_new(10, 6, 24, 1, 1, 1);
	ttfspin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
#else
	ttfspin = gtk_spin_button_new_with_range(6, 24, 1);
#endif
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(ttfspin), 
				  get_point_size(gretl_png_font()));
	gtk_table_attach_defaults(GTK_TABLE(tbl), ttfspin, 2, 3, 
				  tbl_len-1, tbl_len);
	gtk_widget_show(ttfspin);
    } else {
	ttfcombo = NULL;
	ttfspin = NULL;
    }

    if (gnuplot_has_specified_colors()) { 
	GtkWidget *hsep;
	int colmax;

	tbl_len++;
	hsep = gtk_hseparator_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
				  tbl_len-1, tbl_len);  
	gtk_widget_show(hsep);

	if (frequency_plot(spec)) {
	    colmax = 1;
	} else {
	    colmax = COLOR_MAX;
	}

	for (i=0; i<colmax; i++) {
	    GtkWidget *button, *hbox;
	    char labstr[16];

	    if (frequency_plot(spec)) {
		i = COLOR_MAX;
	    }

	    tbl_len++;
	    hbox = gtk_hbox_new(FALSE, 2);

	    if (i == COLOR_MAX) {
		sprintf(labstr, _("Fill color"));
	    } else {
		sprintf(labstr, _("Color %d"), i + 1);
	    }

	    label = gtk_label_new(labstr);
	    gtk_container_add(GTK_CONTAINER(hbox), label);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 0, 1, 
				      tbl_len-1, tbl_len);
	    gtk_widget_show(label);
	    gtk_widget_show(hbox);

	    hbox = gtk_hbox_new(FALSE, 2);
	    button = color_patch_button(i);
	    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 1, 2, 
				      tbl_len-1, tbl_len);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(gnuplot_color_selector), 
			     GINT_TO_POINTER(i));
	    gtk_widget_show_all(button);
	    gtk_widget_show(hbox);
	}
    }
}

/* ........................................................... */

static void gpt_tab_output (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl;
    int i, tbl_len;
    GList *termlist = NULL;

    gchar *termtypes[] = {
	"postscript",
	"postscript color",
	"fig",
	"latex",
	"png",
	"plot commands",
	NULL
    };  

    for (i=0; termtypes[i] != NULL; i++) {
	termlist = g_list_append(termlist, termtypes[i]);
    }
   
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);
    
    label = gtk_label_new (_("Output to file"));
    gtk_widget_show(label);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    tbl_len++;
    label = gtk_label_new(_("output format"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(label);

    termcombo = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), termcombo, 1, 2, 
			      tbl_len-1, tbl_len);
    gtk_combo_set_popdown_strings(GTK_COMBO(termcombo), termlist);   
    gtk_widget_show(termcombo);

    /* button to generate output to file */
    filesavebutton = gtk_button_new_with_label (_("Save to file..."));
    GTK_WIDGET_SET_FLAGS(filesavebutton, GTK_CAN_DEFAULT);
    tbl_len++;
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      filesavebutton, 1, 2, tbl_len-1, tbl_len);
    g_signal_connect (G_OBJECT(filesavebutton), "clicked", 
		      G_CALLBACK(apply_gpt_changes), 
		      spec);
    gtk_widget_grab_default(filesavebutton);
    gtk_widget_show(filesavebutton);    
}

static void linetitle_callback (GtkWidget *w, GPT_SPEC *spec)
{
    set_keyspec_sensitivity(spec);
}

static void gpt_tab_lines (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl;
    int i, tbl_len, tbl_num, tbl_col;
    char label_text[32];
    GList *plot_types = NULL;
    GList *yaxis_loc = NULL;

    if (spec->flags & GPTSPEC_TS) {
	plot_types = g_list_append(plot_types, "lines");
	plot_types = g_list_append(plot_types, "points");
    } else {
	plot_types = g_list_append(plot_types, "points");
	plot_types = g_list_append(plot_types, "lines");
    }

    plot_types = g_list_append(plot_types, "linespoints"); 
    plot_types = g_list_append(plot_types, "impulses");
    plot_types = g_list_append(plot_types, "dots");

    yaxis_loc = g_list_append(yaxis_loc, "left");
    yaxis_loc = g_list_append(yaxis_loc, "right");

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);

    label = gtk_label_new(_("Lines"));
    gtk_widget_show(label);

    if (spec->n_lines > 4) {
	GtkWidget *scroller;
	
	scroller = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				       GTK_POLICY_AUTOMATIC, 
				       GTK_POLICY_AUTOMATIC);
	gtk_widget_show(scroller);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), 
					      vbox);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), scroller, label); 
    } else {
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label); 
    }  

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);

    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    for (i=0; i<spec->n_lines; i++) {
	/* identifier and key or legend text */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	sprintf(label_text, _("line %d: "), i + 1);
	label = gtk_label_new(label_text);
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 0, 1, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	label = gtk_label_new(_("legend"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	linetitle[i] = gtk_entry_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  linetitle[i], 2, 3, tbl_len-1, tbl_len);

	strip_lr(spec->lines[i].title);
	gp_string_to_entry(linetitle[i], spec->lines[i].title);

	g_signal_connect(G_OBJECT(linetitle[i]), "changed", 
			 G_CALLBACK(linetitle_callback), 
			 spec);
	g_signal_connect(G_OBJECT(linetitle[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_widget_show(linetitle[i]);

	/* line type or style */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("type"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	stylecombo[i] = gtk_combo_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  stylecombo[i], 2, 3, tbl_len-1, tbl_len);
	gtk_combo_set_popdown_strings(GTK_COMBO(stylecombo[i]), plot_types); 
	gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(stylecombo[i])->entry), 
			   spec->lines[i].style);  
	gtk_widget_show(stylecombo[i]);	

	/* scale factor for data */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("scale"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	linescale[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(linescale[i]), 6);
	gtk_entry_set_text(GTK_ENTRY(linescale[i]), spec->lines[i].scale);
#ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(linescale[i]), "activate", 
			   GTK_SIGNAL_FUNC(apply_gpt_changes), 
			   spec);
#else
	gtk_entry_set_width_chars(GTK_ENTRY(linescale[i]), 6);
	g_signal_connect(G_OBJECT(linescale[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
#endif
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  linescale[i], 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show(linescale[i]);

	/* use left or right y axis? */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("y axis"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	yaxiscombo[i] = gtk_combo_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  yaxiscombo[i], 2, 3, tbl_len-1, tbl_len);
	gtk_combo_set_popdown_strings(GTK_COMBO(yaxiscombo[i]), yaxis_loc); 
	gtk_entry_set_text (GTK_ENTRY(GTK_COMBO(yaxiscombo[i])->entry), 
			    (spec->lines[i].yaxis == 1)? "left" : "right");  
	gtk_widget_show(yaxiscombo[i]);	
    }
}

static void label_pos_to_entry (double *pos, GtkWidget *w)
{
    if (!na(pos[0]) && !na(pos[1])) {
	gchar *s = g_strdup_printf("%g %g", pos[0], pos[1]);

	gtk_entry_set_text(GTK_ENTRY(w), s);
	g_free(s);
    } else {
	gtk_entry_set_text(GTK_ENTRY(w), "");
    }
}

static void gpt_tab_labels (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl, *menu;
    int i, j, tbl_len, tbl_num, tbl_col;
    char label_text[32];
    png_plot *plot = (png_plot *) spec->ptr;

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);

    label = gtk_label_new(_("Labels"));

    gtk_widget_show(label);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	GtkWidget *hbox, *button, *image;
#ifdef OLD_GTK
	GdkPixmap *pixmap;
	GdkBitmap *mask;	
#else
	GdkPixbuf *icon;
#endif

	/* label text */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	sprintf(label_text, _("label %d: "), i + 1);
	label = gtk_label_new(label_text);
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 0, 1, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	label = gtk_label_new(_("text"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	labeltext[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(labeltext[i]), PLOT_LABEL_TEXT_LEN);
	gp_string_to_entry(labeltext[i], spec->labels[i].text);
#ifdef OLD_GTK
	gtk_signal_connect (GTK_OBJECT(labeltext[i]), "activate", 
			    GTK_SIGNAL_FUNC(apply_gpt_changes), 
			    spec);
#else
	gtk_entry_set_width_chars(GTK_ENTRY(labeltext[i]), PLOT_LABEL_TEXT_LEN);
	g_signal_connect (G_OBJECT(labeltext[i]), "activate", 
			  G_CALLBACK(apply_gpt_changes), 
			  spec);
#endif
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  labeltext[i], 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show(labeltext[i]);

	/* label placement */
	tbl_len++;

	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("position (X Y)"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	/* holder for entry and button */
	hbox = gtk_hbox_new(FALSE, 5);

	/* entry for coordinates */
	labelpos[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(labelpos[i]), PLOT_LABEL_POS_LEN);
	label_pos_to_entry(spec->labels[i].pos, labelpos[i]);
#ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(labelpos[i]), "activate", 
			   GTK_SIGNAL_FUNC(apply_gpt_changes), 
			   spec);
#else
	gtk_entry_set_width_chars(GTK_ENTRY(labelpos[i]), PLOT_LABEL_POS_LEN);
	g_signal_connect(G_OBJECT(labelpos[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
#endif
	gtk_container_add(GTK_CONTAINER(hbox), labelpos[i]);
	gtk_widget_show(labelpos[i]);

	if (plot_is_mouseable(plot)) {
	    /* button to invoke mouse-assisted placement */
	    button = gtk_button_new();
	    g_object_set_data(G_OBJECT(button), "labelpos_entry", labelpos[i]);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(plot_label_position_click), spec);
#ifdef OLD_GTK
	    pixmap = gdk_pixmap_create_from_xpm_d(mdata->w->window,
						  &mask, NULL, mini_mouse_xpm);
	    image = gtk_pixmap_new(pixmap, mask);
#else
	    icon = gdk_pixbuf_new_from_xpm_data((const char **) mini_mouse_xpm);
	    image = gtk_image_new_from_pixbuf(icon);
	    gtk_widget_set_size_request(button, 32, 26);
#endif
	    gtk_container_add (GTK_CONTAINER(button), image);
	    gtk_container_add(GTK_CONTAINER(hbox), button);
	    gtk_widget_show_all(button);
	}

	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  hbox, 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show(hbox);

	/* label justification */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("justification"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	labeljust[i] = gtk_option_menu_new();
	menu = gtk_menu_new();
	for (j=0; j<3; j++) {
	    GtkWidget *item;

	    item = gtk_menu_item_new_with_label(gp_justification_string(j));
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
#ifdef OLD_GTK
	    gtk_object_set_data(GTK_OBJECT(item), "option", GINT_TO_POINTER(j));
#endif
	}
	gtk_option_menu_set_menu(GTK_OPTION_MENU(labeljust[i]), menu);	
	gtk_option_menu_set_history(GTK_OPTION_MENU(labeljust[i]),
				    spec->labels[i].just);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  labeljust[i], 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show_all(labeljust[i]);	
    }
}

/* ........................................................... */

static void gpt_tab_XY (GtkWidget *notebook, GPT_SPEC *spec, gint axis) 
{
    GtkWidget *vbox, *manual, *tbl;
    GtkWidget *label = NULL;
    int i, tbl_len;
   
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);

    if (axis == 0) {
	label = gtk_label_new(_("X-axis"));
    } else if (axis == 1) {
	label = gtk_label_new(_("Y-axis"));
    } else if (axis == 2) {
	label = gtk_label_new(_("Y2-axis"));
    } else {
	return;
    }

    gtk_widget_show(label);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].tab == 1 + axis) {
	    GtkWidget *title_entry;

	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
            
	    label = gtk_label_new(_(gpt_titles[i].description));
	    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      label, 0, 1, tbl_len-1, tbl_len);
	    gtk_widget_show(label);

	    title_entry = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      title_entry, 1, 2, 
				      tbl_len-1, tbl_len);
	    gp_string_to_entry(title_entry, spec->titles[i]);

	    g_signal_connect(G_OBJECT(title_entry), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     spec);

	    gtk_widget_show(title_entry);
	    gpt_titles[i].widget = title_entry;
	}
    } 

    if (spec->code != PLOT_REGULAR) return;

    /* axis range: auto versus manual buttons */
    axis_range[axis].ID = axis;
    tbl_len += 3;
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);

    label = gtk_label_new("");
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      tbl_len-3, tbl_len-2);
    gtk_widget_show(label);

    axis_range[axis].isauto = 
	gtk_radio_button_new_with_label(NULL, _("auto axis range"));

    g_signal_connect(G_OBJECT(axis_range[axis].isauto), "clicked",
		     G_CALLBACK(flip_manual_range), 
		     GINT_TO_POINTER(axis_range[axis].ID));

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON
				 (axis_range[axis].isauto), TRUE);
    gtk_table_attach_defaults(GTK_TABLE(tbl), axis_range[axis].isauto, 
			      0, 1, tbl_len-2, tbl_len-1);
    gtk_widget_show(axis_range[axis].isauto);

    manual = 
	gtk_radio_button_new_with_label(gtk_radio_button_get_group 
					(GTK_RADIO_BUTTON 
					 (axis_range[axis].isauto)),
					_("manual range:")); 
    g_signal_connect(G_OBJECT(manual), "clicked",
		     G_CALLBACK(flip_manual_range), 
		     GINT_TO_POINTER(axis_range[axis].ID));

    gtk_table_attach_defaults(GTK_TABLE(tbl), manual, 0, 1, 
			      tbl_len-1, tbl_len);
    gtk_widget_show(manual);

    /* axis range min. entry */
    tbl_len++;
    label = gtk_label_new(_("minimum"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(label);
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    axis_range[axis].min = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), axis_range[axis].min, 
			      1, 2, tbl_len-1, tbl_len);
    gtk_entry_set_text(GTK_ENTRY(axis_range[axis].min), "");

    g_signal_connect(G_OBJECT(axis_range[axis].min), "activate", 
		     G_CALLBACK(apply_gpt_changes), 
		     spec);

    gtk_widget_show(axis_range[axis].min);

    /* axis range max. entry */
    tbl_len++;
    label = gtk_label_new(_("maximum"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      tbl_len-1, tbl_len);
    gtk_widget_show(label);
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    axis_range[axis].max = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), axis_range[axis].max, 
			      1, 2, tbl_len-1, tbl_len);
    gtk_entry_set_text(GTK_ENTRY(axis_range[axis].max), "");

    g_signal_connect(G_OBJECT(axis_range[axis].max), "activate", 
		     G_CALLBACK(apply_gpt_changes), 
		     spec);

    gtk_widget_show(axis_range[axis].max);

    if (na(spec->range[axis][0])) {
	flip_manual_range(NULL, GINT_TO_POINTER(axis_range[axis].ID));
    } else {
	double_to_gp_entry(spec->range[axis][0], axis_range[axis].min);
	double_to_gp_entry(spec->range[axis][1], axis_range[axis].max);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(axis_range[axis].isauto), 
				     FALSE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(manual), 
				     TRUE);
    }
}

/* ........................................................... */

int show_gnuplot_dialog (GPT_SPEC *spec) 
{
    GtkWidget *button, *notebook;
    int i;

    if (gpt_control != NULL) {
	errbox(_("You can only have one plot controller open\n"
		 "at any given time"));
	return 1;
    }

    for (i=0; i<MAX_AXES; i++) {
	axis_range[i].isauto = NULL;
    }

    for (i=0; i<NTITLES; i++) {
	gpt_titles[i].widget = NULL;
    }

    gpt_control = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(gpt_control), _("gretl plot controls"));
    gtk_container_set_border_width 
        (GTK_CONTAINER(GTK_DIALOG(gpt_control)->vbox), 10);
    gtk_container_set_border_width 
        (GTK_CONTAINER(GTK_DIALOG(gpt_control)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(gpt_control)->vbox), 2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(gpt_control), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(gpt_control), "destroy",
		     G_CALLBACK(close_plot_controller), 
		     (gpointer *) spec);
   
    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->vbox), 
		       notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    gpt_tab_main(notebook, spec);
    gpt_tab_XY(notebook, spec, 0);
    gpt_tab_XY(notebook, spec, 1);

    if (spec->flags & GPTSPEC_Y2AXIS) {
	gpt_tab_XY(notebook, spec, 2);
    }

    if (!no_edit_lines(spec)) {
	gpt_tab_lines(notebook, spec);
    }    

    gpt_tab_labels(notebook, spec); 
    gpt_tab_output(notebook, spec);

    /* "Apply" button */
    button = standard_button(GTK_STOCK_APPLY);
    GTK_WIDGET_SET_FLAGS (button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_gpt_changes), spec);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* "OK" button (apply and close) */
    button = standard_button(GTK_STOCK_OK);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_gpt_changes), spec);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(delete_widget), gpt_control);
    gtk_widget_show(button);

    /* Close button (do not apply changes */
    button = standard_button(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(delete_widget), gpt_control);
    gtk_widget_show(button);

    button = standard_button(GTK_STOCK_HELP);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(context_help), 
		     GINT_TO_POINTER(GR_PLOT));
    gtk_widget_show(button);

    set_keyspec_sensitivity(spec);

    gtk_widget_show(gpt_control);

    return 0;
}

void raise_gpt_control_window (void)
{
    if (gpt_control != NULL) {
	gdk_window_raise(gpt_control->window);
    }
}

void destroy_gpt_control_window (void)
{
    gtk_widget_destroy(gpt_control);
}


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

/* gpt_dialog.c for gretl -- GUI gnuplot controller dialog */

#include "gretl.h"
#include "plotspec.h"
#include "gpt_control.h"
#include "graphics.h"
#include "session.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "calculator.h"
#include "gpt_dialog.h"

#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION >= 2
# define USE_GTK_FONT_CHOOSER 1
#else
# define USE_GTK_FONT_CHOOSER 0
#endif

#include "../pixmaps/mouse.xpm"
#include "gppoints.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

struct gpt_titles_t {
    char *desc;        /* How the field will show up in the options dialog */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
};

struct gpt_range_t {
    gint ID;
    GtkWidget *isauto;
    GtkWidget *min;
    GtkWidget *max;
    GtkWidget *lbase;
};

#define MAX_AXES 3
#define NTITLES  4

/* apparatus for handling the case where plot components are
   added or removed via the gui dialog, but these changes are
   discarded (the user doesn't click "OK" or "Apply")
*/

#define old_lines_init(e) (e->old_n_lines = -1, e->old_lines = NULL)
#define old_labels_init(e) (e->old_n_labels = -1, e->old_labels = NULL)
#define old_arrows_init(e) (e->old_n_arrows = -1, e->old_arrows = NULL)

#define lines_not_synced(e) (e->old_n_lines < 0)
#define labels_not_synced(e) (e->old_n_labels < 0)
#define arrows_not_synced(e) (e->old_n_arrows < 0)

#define restore_lines(e) (e->old_n_lines >= 0)
#define restore_labels(e) (e->old_n_labels >= 0)
#define restore_arrows(e) (e->old_n_arrows >= 0)

typedef struct plot_editor_ plot_editor;

struct plot_editor_ {
    GtkWidget *dialog;
    GtkWidget *notebook;
    GPT_SPEC *spec;
    gchar *user_barsfile;
    gint active_bars;

    GtkWidget **linetitle;
    GtkWidget **lineformula;
    GtkWidget **stylecombo;
    GtkWidget **yaxiscombo;
    GtkWidget **linescale;
    GtkWidget **linewidth;
    GtkWidget **pointsize;
    GtkWidget **colorsel;

    GtkWidget *fitformula;
    GtkWidget *fitlegend;

    GtkWidget **labeltext;
    GtkWidget **labeljust;
    GtkWidget **labelpos;
    GtkWidget **arrowpos;

    GtkWidget *keycombo;
    GtkWidget *fitcombo;
    GtkWidget *border_check;
    GtkWidget *grid_check;
    GtkWidget *grid_combo;
    GtkWidget *y2_check;
    GtkWidget *bars_check;
    GtkWidget *fontcheck;
    GtkWidget *barscombo;

    int gui_nlines;
    int gui_nlabels;
    int gui_narrows;

    struct gpt_range_t axis_range[MAX_AXES];
    struct gpt_titles_t gpt_titles[NTITLES];

    int old_n_lines;
    int old_n_labels;
    int old_n_arrows;

    GPT_LINE *old_lines;
    GPT_LABEL *old_labels;
    GPT_ARROW *old_arrows;
};

#define PLOT_POSITION_LEN 32

const gchar *fittype_strings[] = {
    N_("none"),
    N_("linear: y = a + b*x"),
    N_("quadratic: y = a + b*x + c*x^2"),
    N_("cubic: y = a + b*x + c*x^2 + d*x^3"),
    N_("inverse: y = a + b*(1/x)"),
    N_("loess (locally weighted fit)"),
    N_("semilog: log y = a + b*x"),
    N_("linear-log: y = a + b*log(x)"),
    NULL
};

enum {
    GUI_LINE,
    GUI_LABEL,
    GUI_ARROW
};

static void gpt_tab_lines (plot_editor *ed, GPT_SPEC *spec, int ins);
static void gpt_tab_labels (plot_editor *ed, GPT_SPEC *spec, int ins);
static void gpt_tab_arrows (plot_editor *ed, GPT_SPEC *spec, int ins);
static int add_line_widget (plot_editor *ed);
static int add_label_widget (plot_editor *ed);
static int add_arrow_widget (plot_editor *ed);
static void plot_editor_set_fontname (plot_editor *ed, const char *name);
static int line_get_point_type (GPT_LINE *line, int i);

/* graph color selection apparatus */

#define XPMROWS 19
#define XPMCOLS 17

#define scale_round(v) ((v) * 255.0 / 65535.0)

static GtkWidget *get_image_for_color (const gretlRGB *color)
{
    static char **xpm = NULL;
    GdkPixbuf *icon;
    GtkWidget *image;
    char colstr[8] = {0};
    int i;

    if (color == NULL) {
	return NULL;
    }

    if (xpm == NULL) {
	xpm = strings_array_new_with_length(XPMROWS, XPMCOLS);
	if (xpm == NULL) {
	    return NULL;
	}

	/* common set-up */
	strcpy(xpm[0], "16 16 2 1");
	strcpy(xpm[1], "X      c #000000");
	strcpy(xpm[2], ".      c #000000");
	strcpy(xpm[3], "................");

	for (i=4; i<XPMROWS-1; i++) {
	    strcpy(xpm[i], ".XXXXXXXXXXXXXX.");
	}

	strcpy(xpm[XPMROWS-1], "................");
    }

    /* write in the specific color we want */
    print_rgb_hash(colstr, color);
    for (i=0; i<6; i++) {
	xpm[1][10+i] = colstr[i+1];
    }    

    icon = gdk_pixbuf_new_from_xpm_data((const char **) xpm);
    image = gtk_image_new_from_pixbuf(icon);
    g_object_unref(icon);
    
    return image;
}

static void color_select_callback (GtkWidget *button, GtkWidget *w)
{
    GtkWidget *csel;
    GtkWidget *color_button, *image;
    GdkColor gcolor;
    gpointer data;
    gretlRGB rgb;
    gint i;

    color_button = g_object_get_data(G_OBJECT(w), "color_button");
    csel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(w));

    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(csel), &gcolor);

    rgb.r = (unsigned char) (scale_round(gcolor.red));
    rgb.g = (unsigned char) (scale_round(gcolor.green));
    rgb.b = (unsigned char) (scale_round(gcolor.blue));

    i = widget_get_int(w, "colnum");
    data = g_object_get_data(G_OBJECT(color_button), "plotspec");
    
    if (data != NULL) {
	gretlRGB *prgb = malloc(sizeof *prgb);

	if (prgb != NULL) {
	    prgb->r = rgb.r;
	    prgb->g = rgb.g;
	    prgb->b = rgb.b;
	    g_object_set_data_full(G_OBJECT(color_button), "rgb",
				   prgb, free);
	}
    } else {
	set_graph_palette(i, rgb);
	update_persistent_graph_colors();
    } 

    /* update the "image" widget */
    image = g_object_get_data(G_OBJECT(color_button), "image");
    gtk_widget_destroy(image);
    image = get_image_for_color(&rgb);
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(color_button), image);
    g_object_set_data(G_OBJECT(color_button), "image", image);
  
    gtk_widget_destroy(w);
}

static void color_cancel (GtkWidget *button, GtkWidget *w)
{
    gtk_widget_destroy(w);
}

/* reset color patch button after selecting the option to
   restore the default plot colors */

static void color_patch_button_reset (GtkWidget *button, int cnum)
{
    GtkWidget *image;

    image = g_object_get_data(G_OBJECT(button), "image");
    gtk_widget_destroy(image);
    image = get_image_for_color(get_graph_color(cnum));
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(button), image);
    g_object_set_data(G_OBJECT(button), "image", image);

    if (cnum == BOXCOLOR || cnum == BOXCOLOR - 1) {
	update_persistent_graph_colors();
    }
}

static void graph_color_selector (GtkWidget *w, gpointer p)
{
    GPT_SPEC *spec;
    GtkWidget *cdlg, *csel;
    GtkWidget *button;
    gint i = GPOINTER_TO_INT(p);
    char colstr[8];
    GdkColor gcolor;

    spec = g_object_get_data(G_OBJECT(w), "plotspec");

    if (spec != NULL && spec->lines[i].rgb[0] != '\0') {
	strcpy(colstr, spec->lines[i].rgb);
    } else {
	const gretlRGB *rgb = get_graph_color(i);

	if (rgb == NULL) {
	    fprintf(stderr, "graph_color_selector: got NULL rgb\n");
	    return;
	}
	print_rgb_hash(colstr, rgb);
    }

    gdk_color_parse(colstr, &gcolor);

    cdlg = gtk_color_selection_dialog_new(_("gretl: graph color selection"));

    widget_set_int(cdlg, "colnum", i);
    g_object_set_data(G_OBJECT(cdlg), "color_button", w);

    csel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(cdlg));
    gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(csel), &gcolor);

    g_object_get(G_OBJECT(cdlg), "ok-button", &button, NULL);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_select_callback), cdlg);

    g_object_get(G_OBJECT(cdlg), "cancel-button", &button, NULL);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_cancel), cdlg);

    gtk_widget_show(cdlg);
    gtk_window_set_modal(GTK_WINDOW(cdlg), TRUE);
}

static GtkWidget *color_patch_button (int i)
{
    GtkWidget *image, *button;

    image = get_image_for_color(get_graph_color(i));

    if (image == NULL) {
	button = gtk_button_new_with_label(_("Select color"));
    } else {
	button = gtk_button_new();
	gtk_container_add(GTK_CONTAINER(button), image);
	g_object_set_data(G_OBJECT(button), "image", image);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(graph_color_selector), 
			 GINT_TO_POINTER(i));
    }	

    return button;
}

static int style_from_line_number (GPT_SPEC *spec, int i,
				   int *dotted)
{
    int j, lt = spec->lines[i].type;

    if (lt > 0) {
	j = lt - 1;
    } else if (lt == LT_AUTO) {
	j = i;
    } else if (lt == 0) {
	if (dotted != NULL) {
	    *dotted = 1;
	}
	j = i;
    } else {
	j = LT_AUTO;
    }

    if (j >= spec->n_lines) {
	/* don't go out of bounds */
	j = -1;
    }

    return j;
}

static GtkWidget *line_color_button (GPT_SPEC *spec, int i)
{
    GtkWidget *image = NULL;
    GtkWidget *button = NULL;
    int j = i, dotted = 0;

    if (spec->lines[i].type == 0) {
	return NULL;
    }

    if (spec->lines[i].rgb[0] != '\0') {
	/* we have a specific color for this line */
	gretlRGB color;

	gretl_rgb_get(&color, spec->lines[i].rgb); 
	image = get_image_for_color(&color);
    } else if (spec->lines[i].style == GP_STYLE_FILLEDCURVE) {
	image = get_image_for_color(get_graph_color(SHADECOLOR));
    } else {
	j = style_from_line_number(spec, i, &dotted);
	if (j < 0) {
	    return NULL;
	}
	image = get_image_for_color(get_graph_color(j));
    }

    if (image != NULL) {
	button = gtk_button_new();
	gtk_container_add(GTK_CONTAINER(button), image);
	g_object_set_data(G_OBJECT(button), "image", image);
	g_object_set_data(G_OBJECT(button), "plotspec", spec);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(graph_color_selector), 
			 GINT_TO_POINTER(j));
    }

    if (dotted && button != NULL) {
	gtk_widget_set_sensitive(button, FALSE);
    }

    return button;
}

static void apply_line_color (GtkWidget *cb, GPT_SPEC *spec,
			      int i)
{
    gretlRGB *rgb = NULL;
 
    if (cb != NULL && gtk_widget_is_sensitive(cb)) {
	rgb = g_object_get_data(G_OBJECT(cb), "rgb");
    }

    if (rgb != NULL && i >= 0 && i < spec->n_lines) {
	print_rgb_hash(spec->lines[i].rgb, rgb);
    }

    /* should we also update the "palette"? */
}

/* end graph color selection apparatus */

static GdkPixbuf *get_pixbuf_for_line (int dots)
{
    static char **xpm = NULL;
    int i;

    if (xpm == NULL) {
	xpm = strings_array_new_with_length(14, 31);
	if (xpm == NULL) {
	    return NULL;
	}

	/* common set-up */
	strcpy(xpm[0], "30 11 2 1");
	strcpy(xpm[1], "  c None");
	strcpy(xpm[2], ". c #000000");

	for (i=3; i<14; i++) {
	    strcpy(xpm[i], "                              ");
	}
    }

    /* write in the specific pattern we want */
    if (dots) {
	strcpy(xpm[8], ".  .  .  .  .  .  .  .  .  .  ");
    } else {
	strcpy(xpm[8], "............................  ");
    }    

    return gdk_pixbuf_new_from_xpm_data((const char **) xpm);
}

static void flip_manual_range (GtkWidget *widget, plot_editor *ed)
{
    gint i = widget_get_int(widget, "axis");
    gboolean s = button_is_active(ed->axis_range[i].isauto);

    gtk_widget_set_sensitive(ed->axis_range[i].min, !s);
    gtk_widget_set_sensitive(ed->axis_range[i].max, !s);
}

/* Take text from a gtkentry and write to gnuplot spec string */

static void entry_to_gp_string (GtkWidget *w, char *targ, size_t n)
{
    const gchar *s;

    *targ = '\0';

    g_return_if_fail(GTK_IS_ENTRY(w));

    s = gtk_entry_get_text(GTK_ENTRY(w));
    if (s != NULL && *s != '\0') {
	strncat(targ, s, n - 1);
    }
}

static void combo_to_gp_style (GtkWidget *w, int *sty)
{
    gchar *s;

    g_return_if_fail(GTK_IS_COMBO_BOX(w));

    s = combo_box_get_active_text(w);

    if (s != NULL && *s != '\0') {
	*sty = gp_style_index_from_display_name(s);
    }

    g_free(s);
}

static int ftype_from_selected (GtkWidget *w)
{
    FitType f = PLOT_FIT_NONE;

    if (GTK_IS_COMBO_BOX(w)) {
	f = gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	if (f == PLOT_FIT_LOGLIN && widget_get_int(w, "no-semilog")) {
	    /* the second-last fit option will have been
	       suppressed, so we're off by one */
	    f = PLOT_FIT_LINLOG;
	}
    }

    return f;
}

static int set_fit_type_from_combo (GtkWidget *box,
				    GPT_SPEC *spec)
{
    int oldfit = widget_get_int(box, "oldfit");
    FitType f;
    int err = 0;

    f = ftype_from_selected(box);
    
    if (f == oldfit) {
	/* no change */
	return 0;
    }

    if (f == PLOT_FIT_OLS || f == PLOT_FIT_QUADRATIC || 
	f == PLOT_FIT_CUBIC || f == PLOT_FIT_INVERSE || 
	f == PLOT_FIT_LOESS || f == PLOT_FIT_LOGLIN ||
	f == PLOT_FIT_LINLOG) {
	err = plotspec_add_fit(spec, f);
	if (err) {
	    gui_errmsg(err);
	    f = PLOT_FIT_NONE;
	} else {
	    spec->flags &= ~GPT_FIT_HIDDEN;
	}
    }
    
    if (f == PLOT_FIT_NONE) {
	if (spec->n_lines >= 2) {
	    spec->flags |= GPT_FIT_HIDDEN;
	}
	spec->fit = f;
    }

    widget_set_int(box, "oldfit", f);

    return err;
}

/* In this callback we're reacting to the user's selection of a
   "fit type" for a plot (that is, to the "changed" signal from
   the drop-down fit-type selector). We adjust other elements in
   the plot dialog (the default title, the representations under
   the "Lines" tab) to match this selection. 

   Note, however, that no changes are being "committed" at this
   point -- that will happen only if/when the user clicks "Apply"
   or "OK" in the plot editor dialog action area.
*/

static gboolean fit_type_changed (GtkComboBox *box, plot_editor *ed)
{
    GPT_SPEC *spec = ed->spec;
    const char *s1 = spec->yvarname;
    const char *s2 = spec->xvarname;
    FitType f;

    f = ftype_from_selected(GTK_WIDGET(box));

    if ((spec->flags & GPT_TS) && f != PLOT_FIT_NONE && ed->keycombo != NULL) {
	if (!gtk_widget_is_sensitive(ed->keycombo)) {
	    gtk_widget_set_sensitive(ed->keycombo, TRUE);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(ed->keycombo), 0);
	}
    }	

    if (*s1 != '\0' && *s2 != '\0') {
	/* revise the default plot title */
	gchar *title = NULL;
	
	if (f == PLOT_FIT_OLS) {
	    title = g_strdup_printf(_("%s versus %s (with least squares fit)"),
				    s1, s2);
	} else if (f == PLOT_FIT_QUADRATIC) {
	    title = g_strdup_printf(_("%s versus %s (with quadratic fit)"),
				    s1, s2);
	} else if (f == PLOT_FIT_CUBIC) {
	    title = g_strdup_printf(_("%s versus %s (with cubic fit)"),
				    s1, s2);
	} else if (f == PLOT_FIT_INVERSE) {
	    title = g_strdup_printf(_("%s versus %s (with inverse fit)"),
				    s1, s2);
	} else if (f == PLOT_FIT_LOESS) {
	    title = g_strdup_printf(_("%s versus %s (with loess fit)"),
				    s1, s2);
	} else if (f == PLOT_FIT_LOGLIN) {
	    title = g_strdup_printf(_("%s versus %s (with semilog fit)"),
				    s1, s2);
	} else if (f == PLOT_FIT_LINLOG) {
	    title = g_strdup_printf(_("%s versus %s (with linear-log fit)"),
				    s1, s2);
	} else {
	    title = g_strdup("");
	}

	if (title != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(ed->gpt_titles[0].widget), title);
	    g_free(title);
	}
    }

    /* re-jig the "Lines" tab entries for the revised fitted line */

    if (ed->fitformula != NULL && ed->fitlegend != NULL) {
	set_fit_type_from_combo(GTK_WIDGET(box), spec);
	if (f == PLOT_FIT_LOESS || f == PLOT_FIT_NONE) {
	    gtk_entry_set_text(GTK_ENTRY(ed->fitformula), "");
	} else {
	    gtk_entry_set_text(GTK_ENTRY(ed->fitformula), spec->lines[1].formula);
	}
	gtk_entry_set_text(GTK_ENTRY(ed->fitlegend), spec->lines[1].title);
    }
    
    return FALSE;
}

static void dot_callback (GtkComboBox *box, GPT_SPEC *spec)
{
    int i = widget_get_int(GTK_WIDGET(box), "linenum");
    int solid = (gtk_combo_box_get_active(box) == 0);
    GtkWidget *w;

    /* flip between colored solid line (the first option)
       and a dotted black line (second option) 
    */
    spec->lines[i].type = solid ? LT_AUTO : 0;

    /* the color selector is effective only for solid lines */
    w = g_object_get_data(G_OBJECT(box), "colorsel");
    if (w != NULL) {
	gtk_widget_set_sensitive(w, solid);
	w = g_object_get_data(G_OBJECT(box), "color-label");
	gtk_widget_set_sensitive(w, solid);
    }
}

static float spinner_get_float (GtkWidget *b)
{
    return (float) gtk_spin_button_get_value(GTK_SPIN_BUTTON(b));
}

/* take a double (which might be NA) and format it for
   a gtkentry widget */

static void double_to_gp_entry (double x, GtkWidget *w)
{
    if (w != NULL && GTK_IS_ENTRY(w)) {
	gchar *numstr;

	if (na(x)) {
	    numstr = g_strdup("NA");
	} else if (x == floor(x)) {
	    numstr = g_strdup_printf("%.1f", x);
	} else {
	    numstr = g_strdup_printf("%.10g", x);
	}
	gtk_entry_set_text(GTK_ENTRY(w), numstr);
	g_free(numstr);
    }
}

/* read a double from a gtkentry, with error checking */

static void entry_to_gp_double (GtkWidget *w, double *val)
{
    double x = NADBL;

    if (w != NULL && GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));

	if (s != NULL) {
	    if (!strcmp(s, "NA")) {
		*val = NADBL;
		return;
	    } else if (*s != '\0') {
		gchar *tmp = g_strdup(s);

		gretl_charsub(tmp, ',', '.');
		gretl_push_c_numeric_locale();	    
		if (check_atof(tmp)) {
		    errbox(gretl_errmsg_get());
		} else {
		    x = atof(tmp);
		}
		gretl_pop_c_numeric_locale();
		g_free(tmp);
	    }
	}
    }

    if (!na(x)) {
	*val = x;
    }    
}

#define gp_string_to_entry(w,s) do { \
	gtk_entry_set_text(GTK_ENTRY(w),s); } while (0)

static int get_label_pos_from_entry (GtkWidget *w, double *pos)
{
    int err = 0;

    if (GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
	int chk;
    
	chk = sscanf(s, "%lf %lf", &pos[0], &pos[1]);
	if (chk != 2) {
	    errbox(_("Invalid position, must be X Y"));
	    gtk_editable_select_region(GTK_EDITABLE(w), 0, strlen(s));
	    pos[0] = pos[1] = NADBL;
	    err = 1;
	} 
    } else {
	err = 1;
    }

    return err;
}

static int get_arrow_spec_from_entries (plot_editor *ed, int i, 
					GPT_ARROW *arrow)
{
    int j, err = 0;

    for (j=0; j<2 && !err; j++) {
	GtkWidget *entry, *b1 = NULL, *b2 = NULL;
	double *x, *y;

	if (j == 0) {
	    entry = ed->arrowpos[i];
	    b1 = g_object_get_data(G_OBJECT(entry), "arrow_check");
	    b2 = g_object_get_data(G_OBJECT(entry), "dots_check");
	    x = &arrow->x0;
	    y = &arrow->y0;
	} else {
	    entry = g_object_get_data(G_OBJECT(ed->arrowpos[i]), "pos_2");
	    x = &arrow->x1;
	    y = &arrow->y1;
	}

	if (GTK_IS_ENTRY(entry)) {
	    const gchar *s = gtk_entry_get_text(GTK_ENTRY(entry));
	    int n;
    
	    n = sscanf(s, "%lf %lf", x, y);
	    if (n != 2) {
		errbox(_("Invalid position, must be X Y"));
		gtk_editable_select_region(GTK_EDITABLE(entry), 0, strlen(s));
		*x = *y = NADBL;
		err = 1;
	    } 
	} else {
	    err = 1;
	}

	if (b1 != NULL) {
	    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b1))) {
		arrow->flags |= GP_ARROW_HEAD;
	    } else {
		arrow->flags &= ~GP_ARROW_HEAD;
	    }
	}
	if (b2 != NULL) {
	    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b2))) {
		arrow->flags |= GP_ARROW_DOTS;
	    } else {
		arrow->flags &= ~GP_ARROW_DOTS;
	    }
	}
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

static int 
set_logscale_from_entry (GPT_SPEC *spec, int i, GtkWidget *entry)
{
    double base;
    int err = 0;

    base = entry_get_numeric_value(entry, C_POS_DBL);
    if (na(base)) {
	err = 1;
    } else if (!na(base) && base < 1.1) {
	err = 1;
	errbox("bad base");
    } else {
	spec->logbase[i] = base;
	if (i == 0) {
	    *spec->xfmt = '\0';
	    *spec->xtics = '\0';
	} 
    }

    return err;
}

static void maybe_set_point_type (GPT_LINE *line, GtkWidget *w, int i)
{
    GtkWidget *ptsel = g_object_get_data(G_OBJECT(w), "pointsel");

    if (ptsel != NULL && gtk_widget_is_sensitive(ptsel)) {
	int pt = gtk_combo_box_get_active(GTK_COMBO_BOX(ptsel));
	int pt0 = line_get_point_type(line, i);

	if (pt != pt0) {
	    int ptdef = (line->type == LT_AUTO)? i : line->type - 1;

	    if (pt == ptdef) {
		line->ptype = 0;
	    } else {
		line->ptype = pt + 1;
	    }
	}
    }
}

static void plot_editor_sync (plot_editor *ed)
{
    ed->spec->flags &= ~(GPT_FILL_SWITCH | GPT_ERR_SWITCH);

    free(ed->old_lines);
    old_lines_init(ed);

    free(ed->old_labels);
    old_labels_init(ed);

    free(ed->old_arrows);
    old_arrows_init(ed);
}

static char *default_bars_filename (void)
{
    static char barsname[FILENAME_MAX];

    if (*barsname == '\0') {
	sprintf(barsname, "%sdata%cplotbars%cnber.txt",
		gretl_home(), SLASH, SLASH);
    }

    return barsname;
}

static char user_barsfile[FILENAME_MAX];

/* called on destruction of a plot dialog, if the user has 
   selected a user-defined "plotbars" file, so that is will 
   be remembered for use with the next plot
*/

static void remember_user_bars_file (const char *fname)
{
    *user_barsfile = '\0';
    strncat(user_barsfile, fname, FILENAME_MAX-1);
}

/* if there's a saved filename for a user-defined
   "plotbars" file, append it to the list of options
*/

static gboolean maybe_append_user_bars_file (GtkComboBox *combo,
					     plot_editor *ed)
{
    if (ed->user_barsfile == NULL && *user_barsfile != '\0') {
	ed->user_barsfile = g_strdup(user_barsfile);
    }

    if (ed->user_barsfile != NULL) {
	const char *p = strrchr(ed->user_barsfile, SLASH);
	const char *show = (p != NULL)? (p+1) : ed->user_barsfile;

	combo_box_append_text(combo, show);
	return TRUE;
    } else {
	return FALSE;
    }
}

static void try_adding_plotbars (GPT_SPEC *spec, plot_editor *ed)
{
    png_plot *plot = (png_plot *) spec->ptr;
    double xmin = -1, xmax = -1;
    double ymin = -1, ymax = -1;
    int err;

    err = plot_get_coordinates(plot, &xmin, &xmax, &ymin, &ymax);

    if (!err) {
	int active = gtk_combo_box_get_active(GTK_COMBO_BOX(ed->barscombo));
	const char *fname = NULL;

	if (active == 1 && ed->user_barsfile != NULL) {
	    fname = ed->user_barsfile;
	} else if (active == 0) {
	    fname = default_bars_filename();
	}

	if (fname != NULL) {
	    if (ed->axis_range[1].isauto != NULL &&
		!button_is_active(ed->axis_range[1].isauto)) {
		/* we have a manual y-axis setting that may override the
		   current (ymin, ymax) as read from the plot */
		ymin = spec->range[1][0];
		ymax = spec->range[1][1];
	    } else {
		/* freeze the current (ymin, ymax) so that the addition
		   of bars doesn't disturb the range */
		spec->range[1][0] = ymin;
		spec->range[1][1] = ymax;
	    }

	    err = plotspec_add_bars_info(spec, xmin, xmax,
					 ymin, ymax, fname);
	}
    }

    if (err) {
	gui_errmsg(err);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ed->bars_check),
				     FALSE);
    }
}

/* Respond to "OK" or "Apply", or to hitting the Enter key in
   some parts of the plot editor dialog */

static void apply_gpt_changes (GtkWidget *w, plot_editor *ed) 
{
    GPT_SPEC *spec = ed->spec;
    GPT_LINE *line;
    int suppress_y2 = 0;
    gchar *s;
    int i, k, err = 0;

    for (i=0; i<NTITLES; i++) {
	if (ed->gpt_titles[i].widget != NULL) {
	    entry_to_gp_string(ed->gpt_titles[i].widget, spec->titles[i], 
			       sizeof spec->titles[0]);
	}
    }

    if (ed->keycombo != NULL) {
        s = combo_box_get_active_text(ed->keycombo);
	spec->keyspec = gp_keypos_from_display_name(s);
	g_free(s);
    }

    spec->flags &= ~GPT_Y2AXIS;

    if (ed->y2_check != NULL && button_is_active(ed->y2_check)) {
	suppress_y2 = 1;
    } 

    for (i=0; i<ed->gui_nlines; i++) {
	line = &spec->lines[i];
	line->yaxis = 1;
	if (!suppress_y2 && ed->yaxiscombo[i] != NULL) {
	    s = combo_box_get_active_text(ed->yaxiscombo[i]);
	    if (!strcmp(s, _("right"))) {
		line->yaxis = 2;	
	    }
	    if (line->yaxis == 2) {
		spec->flags |= GPT_Y2AXIS;
	    }
	    g_free(s);
	}
    }

    if (spec->code == PLOT_REGULAR || 
	spec->code == PLOT_CURVE ||
	spec->code == PLOT_MANY_TS) {
	k = (spec->flags & GPT_Y2AXIS)? 3 : 2;
	for (i=0; i<k; i++) {
	    if (ed->axis_range[i].isauto != NULL) {
		if (button_is_active(ed->axis_range[i].isauto)) {
		    spec->range[i][0] = NADBL;
		    spec->range[i][1] = NADBL;
		} else {
		    entry_to_gp_double(ed->axis_range[i].min, &spec->range[i][0]);
		    entry_to_gp_double(ed->axis_range[i].max, &spec->range[i][1]);
		    err = validate_range(spec->range[i]);
		}
	    }
	    if (spec->code == PLOT_REGULAR && ed->axis_range[i].lbase != NULL) {
		if (gtk_widget_is_sensitive(ed->axis_range[i].lbase)) {
		    err = set_logscale_from_entry(spec, i, ed->axis_range[i].lbase);
		} else {
		    spec->logbase[i] = 0.0;
		}
	    }
	}
    }

    for (i=0; i<ed->gui_nlines && !err; i++) {
	GtkWidget *combo = ed->stylecombo[i];

	line = &spec->lines[i];
	if (combo != NULL && gtk_widget_is_sensitive(combo)) {
	    int oldalt = 0;

	    if (line->style == GP_STYLE_FILLEDCURVE ||
		line->style == GP_STYLE_ERRORBARS) {
		oldalt = line->style;
	    }
	    combo_to_gp_style(combo, &line->style);
	    if (oldalt == GP_STYLE_FILLEDCURVE &&
		line->style == GP_STYLE_ERRORBARS) {
		spec->flags |= GPT_ERR_SWITCH;
		spec->flags &= ~GPT_FILL_SWITCH;
	    } else if (oldalt == GP_STYLE_ERRORBARS &&
		       line->style == GP_STYLE_FILLEDCURVE) {
		spec->flags |= GPT_FILL_SWITCH;
		spec->flags &= ~GPT_ERR_SWITCH;
	    }
	    maybe_set_point_type(line, combo, i);
	}
	if (ed->linetitle[i] != NULL) {
	    entry_to_gp_string(ed->linetitle[i], line->title, 
			       sizeof spec->lines[0].title);
	}
	if (ed->lineformula[i] != NULL && 
	    gtk_widget_is_sensitive(ed->lineformula[i])) {
	    entry_to_gp_string(ed->lineformula[i], line->formula, 
			       sizeof spec->lines[0].formula);
	}
	if (ed->linescale[i] != NULL) {
	    entry_to_gp_double(ed->linescale[i], &line->scale);
	}
	if (ed->linewidth[i] != NULL) {
	    line->width = spinner_get_float(ed->linewidth[i]);
	}
	if (ed->colorsel[i] != NULL) {
	    apply_line_color(ed->colorsel[i], spec, i);
	}
	if (ed->pointsize[i] != NULL) {
	    line->pscale = spinner_get_float(ed->pointsize[i]);
	}	
    }

    for (i=0; i<ed->gui_nlabels && !err; i++) {
	entry_to_gp_string(ed->labeltext[i], spec->labels[i].text, 
			   sizeof spec->labels[0].text);
	if (string_is_blank(spec->labels[i].text)) {
	    continue;
	}
	err = get_label_pos_from_entry(ed->labelpos[i], spec->labels[i].pos);
	if (!err) {
	    spec->labels[i].just = 
		gtk_combo_box_get_active(GTK_COMBO_BOX(ed->labeljust[i]));
	}
    } 

    for (i=0; i<ed->gui_narrows && !err; i++) {
	err = get_arrow_spec_from_entries(ed, i, &spec->arrows[i]);
    } 

    if (!err && ed->border_check != NULL) {
	if (button_is_active(ed->border_check)) {
	    /* full border */
	    spec->border = GP_BORDER_DEFAULT;
	} else {
	    /* left and bottom only */
	    spec->border = 3;
	}
    } 

    if (!err && ed->grid_check != NULL) {
	spec->flags &= ~(GPT_GRID_Y | GPT_GRID_X);
	if (button_is_active(ed->grid_check)) {
	    int i = gtk_combo_box_get_active(GTK_COMBO_BOX(ed->grid_combo));

	    if (i == 0) {
		spec->flags |= GPT_GRID_Y;
	    } else if (i == 1) {
		spec->flags |= GPT_GRID_X;
	    } else {
		spec->flags |= (GPT_GRID_Y | GPT_GRID_X);
	    }
	} 
    }

    if (!err) {
	if (ed->fontcheck != NULL && button_is_active(ed->fontcheck)) {
	    if (ed->spec->fontstr != NULL && *ed->spec->fontstr != '\0') {
		set_gretl_png_font(ed->spec->fontstr);
		update_persistent_graph_font();
	    }
	}
    }

    if (!err && ed->fitcombo != NULL && gtk_widget_is_sensitive(ed->fitcombo)) {
	err = set_fit_type_from_combo(ed->fitcombo, spec);
	if ((spec->fit == PLOT_FIT_NONE || spec->fit == PLOT_FIT_LOESS) &&
	    ed->fitformula != NULL && spec->n_lines > 1) {
	    /* scrub irrelevant formula, if any */
	    spec->lines[1].formula[0] = '\0';
	}
    }

    if (!err && ed->bars_check != NULL) {
	if (button_is_active(ed->bars_check)) {
	    try_adding_plotbars(spec, ed);
	} else {
	    plotspec_remove_bars(spec);
	}
    }

    if (!err) {
	png_plot *plot = (png_plot *) spec->ptr;

	set_plot_has_y2_axis(plot, spec->flags & GPT_Y2AXIS);
	redisplay_edited_plot(plot);
	if (plot_is_saved(plot)) {
	    mark_session_changed();
	}
    }

    plot_editor_sync(ed);
}

static void set_keyspec_sensitivity (plot_editor *ed)
{
    gboolean state = FALSE;
    const char *p;
    int i;

    for (i=0; i<ed->gui_nlines; i++) {
	p = gtk_entry_get_text(GTK_ENTRY(ed->linetitle[i]));
	if (p != NULL && *p != 0) {
	    state = TRUE;
	    break;
	}
    }

    gtk_widget_set_sensitive(ed->keycombo, state);
}

#define TAB_MAIN_COLS 3

static void plot_editor_set_fontname (plot_editor *ed, const char *name)
{
    free(ed->spec->fontstr);
    ed->spec->fontstr = gretl_strdup(name);
}

#ifdef G_OS_WIN32

static void real_graph_font_selector (GtkButton *button, gpointer p, int type,
				      const char *default_font)
{
    char fontname[128];

    if (default_font == NULL) {
	if (type == 1) {
	    default_font = pdf_ps_saver_current_font(p);
	} else {
	    default_font = gretl_png_font();
	}
    }

    strcpy(fontname, default_font);
    win32_font_selector(fontname, APP_FONT_SELECTION);

    if (*fontname != '\0') {
	gchar *title = g_strdup_printf(_("font: %s"), fontname);

	gtk_button_set_label(button, title);
	g_free(title);
	if (type == 0) {
	    plot_editor_set_fontname(p, fontname);
	} else if (type == 1) {
	    pdf_ps_saver_set_fontname(p, fontname);
	} else if (type == 2) {
	    activate_plot_font_choice(p, fontname);
	}
    }
}

#elif USE_GTK_FONT_CHOOSER

static void graph_font_selection_ok (GtkWidget *w, GtkFontChooser *fc)
{
    gchar *fontname = gtk_font_chooser_get_font(fc);

    if (fontname != NULL && *fontname != '\0') {
	gpointer p = g_object_get_data(G_OBJECT(fc), "parent");
	gint type = widget_get_int(fc, "parent-type");

	if (type < 2) {
	    GtkWidget *b = g_object_get_data(G_OBJECT(fc), "launcher");
	    gchar *title = g_strdup_printf(_("font: %s"), fontname);

	    gtk_button_set_label(GTK_BUTTON(b), title);
	    g_free(title);
	}

	if (type == 0) {
	    plot_editor_set_fontname(p, fontname);
	} else if (type == 1) {
	    pdf_ps_saver_set_fontname(p, fontname);
	} else if (type == 2) {
	    activate_plot_font_choice(p, fontname);
	}
    }

    g_free(fontname);
    gtk_widget_destroy(GTK_WIDGET(fc));
}

static void real_graph_font_selector (GtkButton *button, gpointer p, int type,
				      const char *default_font)
{
    static GtkWidget *fontsel = NULL;
    GtkWidget *b;

    if (fontsel != NULL) {
	gtk_window_present(GTK_WINDOW(fontsel));
        return;
    }

    if (default_font == NULL) {
	if (type == 1) {
	    default_font = pdf_ps_saver_current_font(p);
	} else {
	    default_font = gretl_png_font();
	}
    }

    fontsel = gtk_font_chooser_dialog_new(_("Font for graphs"), NULL);
    gtk_font_chooser_set_font(GTK_FONT_CHOOSER(fontsel), 
			      default_font); 

    /* attach data items */
    if (button != NULL) {
	g_object_set_data(G_OBJECT(fontsel), "launcher", button);
    }
    g_object_set_data(G_OBJECT(fontsel), "parent", p);
    g_object_set_data(G_OBJECT(fontsel), "parent-type", GINT_TO_POINTER(type));

    gtk_window_set_position(GTK_WINDOW(fontsel), GTK_WIN_POS_MOUSE);

    /* destruction signals */
    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &fontsel);
    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_main_quit), NULL);

    /* button signals */
    b = gtk_dialog_get_widget_for_response(GTK_DIALOG(fontsel),
					   GTK_RESPONSE_OK);
    if (b != NULL) {
	g_signal_connect(G_OBJECT(b), "clicked",
			 G_CALLBACK(graph_font_selection_ok),
			 fontsel);
    }

    b = gtk_dialog_get_widget_for_response(GTK_DIALOG(fontsel),
					   GTK_RESPONSE_CANCEL);
    if (b != NULL) {
	g_signal_connect(G_OBJECT(b), "clicked",
			 G_CALLBACK(delete_widget),
			 fontsel);
    }

    gtk_widget_show(fontsel);

    gtk_main();
}

#else /* using GtkFontSelectionDialog */

/* callback from OK button in graph font selector */

static void graph_font_selection_ok (GtkWidget *w, GtkFontSelectionDialog *fs)
{
    gchar *fontname;

    fontname = gtk_font_selection_dialog_get_font_name(fs);

    if (fontname != NULL && *fontname != '\0') {
	gpointer p = g_object_get_data(G_OBJECT(fs), "parent");
	gint type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(fs), "parent-type"));

	if (type < 2) {
	    GtkWidget *b = g_object_get_data(G_OBJECT(fs), "launcher");
	    gchar *title = g_strdup_printf(_("font: %s"), fontname);

	    gtk_button_set_label(GTK_BUTTON(b), title);
	    g_free(title);
	}

	if (type == 0) {
	    plot_editor_set_fontname(p, fontname);
	} else if (type == 1) {
	    pdf_ps_saver_set_fontname(p, fontname);
	} else if (type == 2) {
	    activate_plot_font_choice(p, fontname);
	}
    }

    g_free(fontname);
    gtk_widget_destroy(GTK_WIDGET(fs));
}

static void real_graph_font_selector (GtkButton *button, gpointer p, int type,
				      const char *default_font)
{
    static GtkWidget *fontsel = NULL;
    GtkWidget *b;

    if (fontsel != NULL) {
	gtk_window_present(GTK_WINDOW(fontsel));
        return;
    }

    if (default_font == NULL) {
	if (type == 1) {
	    default_font = pdf_ps_saver_current_font(p);
	} else {
	    default_font = gretl_png_font();
	}
    }

    fontsel = gtk_font_selection_dialog_new(_("Font for graphs"));
    gtk_font_selection_dialog_set_font_name(GTK_FONT_SELECTION_DIALOG(fontsel), 
					    default_font);
    if (button != NULL) {
	g_object_set_data(G_OBJECT(fontsel), "launcher", button);
    }
    g_object_set_data(G_OBJECT(fontsel), "parent", p);
    g_object_set_data(G_OBJECT(fontsel), "parent-type", GINT_TO_POINTER(type));

    gtk_window_set_position(GTK_WINDOW(fontsel), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &fontsel);
    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_main_quit), NULL);

    b = gtk_font_selection_dialog_get_ok_button(GTK_FONT_SELECTION_DIALOG(fontsel));
    g_signal_connect(G_OBJECT(b), "clicked", G_CALLBACK(graph_font_selection_ok), 
		     fontsel);
    b = gtk_font_selection_dialog_get_cancel_button(GTK_FONT_SELECTION_DIALOG(fontsel));
    g_signal_connect(G_OBJECT(b), "clicked", G_CALLBACK(delete_widget), 
		     fontsel);

    gtk_widget_show(fontsel);

    gtk_main();
}

#endif /* end aternative font selectors */

static void graph_font_selector (GtkButton *button, gpointer p)
{
    real_graph_font_selector(button, p, 0, NULL);
}

void pdf_font_selector (GtkButton *button, gpointer p)
{
    real_graph_font_selector(button, p, 1, NULL);
}

void plot_add_font_selector (png_plot *plot, const char *oldfont) 
{
    real_graph_font_selector(NULL, plot, 2, oldfont);
}

static void strip_lr (gchar *txt)
{
    gchar test[16];
    gchar *p;

    sprintf(test, "(%s)", _("left"));
    p = strstr(txt, test);
    if (p != NULL) {
	*p = '\0';
    } else {
	sprintf(test, "(%s)", _("right"));
	p = strstr(txt, test);
	if (p != NULL) {
	   *p = '\0';
	}
    } 
}

static void toggle_axis_selection (GtkWidget *w, plot_editor *ed)
{
    int no_y2 = button_is_active(w);
    int i;

    for (i=0; i<ed->gui_nlines; i++) {
	if (ed->yaxiscombo[i] != NULL) {
	    gtk_widget_set_sensitive(ed->yaxiscombo[i], !no_y2);
	}
    }
}

/* re-establish the default plot colors and reset the
   color selection buttons accordingly */

static void color_default_callback (GtkWidget *w, GtkWidget *book)
{
    GtkWidget *button;
    char id[32];
    int i;

    sprintf(id, "color-button%d", BOXCOLOR);
    button = g_object_get_data(G_OBJECT(book), id);

    if (button != NULL) {
	graph_palette_reset(BOXCOLOR);
	color_patch_button_reset(button, BOXCOLOR);
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    sprintf(id, "color-button%d", i);
	    button = g_object_get_data(G_OBJECT(book), id);
	    if (button != NULL) {
		if (i == 0) {
		    graph_palette_reset(i);
		}
		color_patch_button_reset(button, i);
	    }
	}
    }
}

static void table_add_row (GtkWidget *tbl, int *rows, int cols)
{
    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, cols);    
}

static void add_color_selector (int i, GtkWidget *tbl, int cols,
				int *rows, GtkWidget *notebook)
{
    static int r0;
    int collen = (N_GP_COLORS - 1) / 2;
    GtkWidget *button, *hbox;
    GtkWidget *label;
    int row, cmin, cmax;
    char str[32];

    if (i == 0) {
	/* get baseline table row for color selectors */
	r0 = *rows;
    }

    if (i < collen || i == BOXCOLOR) {
	table_add_row(tbl, rows, cols);
	cmin = 0;
	cmax = 1;
	row = *rows;
    } else {
	/* place selector to the right */
	row = r0 + 1 + i - collen;
	cmin = 1;
	cmax = 2;
    }

    hbox = gtk_hbox_new(FALSE, 2);

    if (i == BOXCOLOR) {
	strcpy(str, _("Fill color"));
    } else {
	sprintf(str, _("Color %d"), i + 1);
    }

    label = gtk_label_new(str);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_widget_show(label);

    button = color_patch_button(i);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_widget_show_all(button);
    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, cmin, cmax, 
			      row - 1, row);
    gtk_widget_show(hbox);

    sprintf(str, "color-button%d", i);
    g_object_set_data(G_OBJECT(notebook), str, button);

    if (i == BOXCOLOR || i == BOXCOLOR - 1) {
	table_add_row(tbl, rows, TAB_MAIN_COLS);
	hbox = gtk_hbox_new(FALSE, 2);
	button = gtk_button_new_with_label(_("Reset to default"));
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 10);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(color_default_callback), 
			 notebook);
	gtk_table_attach(GTK_TABLE(tbl), hbox, 0, 2, *rows - 1, *rows,
			 GTK_FILL, 0, 0, 5);
	gtk_widget_show(button);
	gtk_widget_show(hbox);
    }
}

static void combo_set_ignore_changed (GtkComboBox *box, gboolean s)
{
    if (s) {
	g_object_set_data(G_OBJECT(box), "ignore-changed", 
			  GINT_TO_POINTER(1));
    } else {
	g_object_set_data(G_OBJECT(box), "ignore-changed", 
			  NULL);
    }
}

/* callback from the file-open selector: if the user
   cancelled then @fname will be NULL
*/

void set_plotbars_filename (const char *fname, gpointer data) 
{
    plot_editor *ed = (plot_editor *) data;
    GtkComboBox *box = GTK_COMBO_BOX(ed->barscombo);

    /* we should block response to the "changed" signal 
       for the duration of this adjustment */
    combo_set_ignore_changed(box, TRUE);

    if (fname == NULL) {
	/* canceled or error */
	if (ed->user_barsfile == NULL) {
	    gtk_combo_box_set_active(box, 0);
	} 
    } else {
	/* we got a completed file selection */
	const char *p = strrchr(fname, SLASH);
	const char *show = (p != NULL)? (p+1) : fname;
	int had_userfile = ed->user_barsfile != NULL;

	if (had_userfile) {
	    g_free(ed->user_barsfile);
	    combo_box_remove(box, 2);
	    combo_box_remove(box, 1);
	} else {
	    combo_box_remove(box, 1);
	}
	ed->user_barsfile = g_strdup(fname);
	combo_box_append_text(box, show);
	combo_box_append_text(box, _("other..."));
	gtk_combo_box_set_active(box, 1);
    }

    /* unblock the "changed" signal */
    combo_set_ignore_changed(box, FALSE);
}

/* "changed" callback from the combo box listing the default 
   plotbars file and possibly a user-specified file: the point
   of this callback is to allow the user to add a new
   plotbars file.
*/

static void plot_bars_changed (GtkComboBox *box, plot_editor *ed)
{
    int i = gtk_combo_box_get_active(box);

    if (g_object_get_data(G_OBJECT(box), "ignore-changed")) {
	return;
    }

    if (i == 0) {
	ed->active_bars = 0; /* selected the default */
    } else if (i == 1 && ed->user_barsfile != NULL) {
	ed->active_bars = 1; /* selected a previously opened own file */
    } else if (i == 2 || (i == 1 && ed->user_barsfile == NULL)) {
	/* selected "other..." */
	const char *msg = N_("To add your own \"bars\" to a plot, you must supply the\n"
			     "name of a plain text file containing pairs of dates.");
	GtkWidget *dlg;
	gint ret;

	dlg = gtk_dialog_new_with_buttons("gretl",
					  GTK_WINDOW(ed->dialog),
					  GTK_DIALOG_MODAL | 
					  GTK_DIALOG_DESTROY_WITH_PARENT,
					  GTK_STOCK_OPEN, 1,
					  _("See example"), 2,
					  GTK_STOCK_CANCEL, 3,
					  NULL);

	gretl_dialog_add_message(dlg, _(msg));
#if GTK_MAJOR_VERSION < 3
	gtk_dialog_set_has_separator(GTK_DIALOG(dlg), FALSE);
#endif
	gtk_window_set_keep_above(GTK_WINDOW(dlg), TRUE);  
	ret = gtk_dialog_run(GTK_DIALOG(dlg));
	gtk_widget_destroy(dlg);

#if 1
	if (ret != 1 && ed->active_bars >= 0) {
	    /* not opening a file: revert to default selection */
	    gtk_combo_box_set_active(box, ed->active_bars);
	}
#endif

	if (ret == 1) {
	    /* open */
	    file_selector_with_parent(OPEN_BARS, FSEL_DATA_MISC, 
				      ed, ed->dialog);
	} else if (ret == 2) {
	    /* loot at the example file */
	    const char *fname = default_bars_filename();

	    view_file(fname, 0, 0, 60, 420, VIEW_FILE);
	}
    }
} 

static GtkWidget *gp_dialog_table (int rows, int cols, 
				   GtkWidget *vbox)
{
    GtkWidget *tbl;

    tbl = gtk_table_new(rows, cols, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);

    return tbl;
}

static GtkWidget *gp_dialog_vbox (void)
{
    GtkWidget *vbox;

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);

    return vbox;
}

static GtkWidget *gp_page_vbox (GtkWidget *notebook, char *str)
{
    GtkWidget *vbox;
    GtkWidget *label;

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);
    label = gtk_label_new(str);
    gtk_widget_show(label);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label); 

    return vbox;
}

static int semilog_is_ok (GPT_SPEC *spec)
{
    const double *x = spec->data;

    if (x == NULL) {
	return 0;
    } else {
	const double *y = x + spec->nobs;

	return gretl_ispositive(0, spec->nobs - 1, y, 1);
    }
}

static int log_x_ok (GPT_SPEC *spec)
{
    const double *x = spec->data;

    if (x == NULL) {
	return 0;
    } else {
	return gretl_ispositive(0, spec->nobs - 1, x, 1);
    }
}

static int plotspec_gridval (GPT_SPEC *spec)
{
    int val = 0;

    if (spec->flags & (GPT_GRID_Y | GPT_GRID_X)) {
	if (!(spec->flags & GPT_GRID_X)) {
	    val = 1;
	} else if (!(spec->flags & GPT_GRID_Y)) {
	    val = 2;
	} else {
	    val = 3;
	}
    }    

    return val;
}

static int show_bars_check (GPT_SPEC *spec)
{
    if (!(spec->flags & GPT_TS)) {
	return 0;
    }

    return (spec->pd == 1 || 
	    spec->pd == 4 || 
	    spec->pd == 12);
}

#define plot_has_tics(s) (strcmp(s->xtics, "none") || strcmp(s->ytics, "none"))

static void gpt_tab_main (plot_editor *ed, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl;
    GtkWidget *hsep, *button;
    gchar *title;
    int i, rows = 1;
    int kactive = 0;

    vbox = gp_page_vbox(ed->notebook, _("Main"));

    tbl = gp_dialog_table(rows, TAB_MAIN_COLS, vbox);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (ed->gpt_titles[i].tab == 0) {
	    GtkWidget *entry;

	    if (i > 0) {
		table_add_row(tbl, &rows, TAB_MAIN_COLS);
	    }

	    label = gtk_label_new(_(ed->gpt_titles[i].desc));
	    gtk_table_attach_defaults(GTK_TABLE (tbl), 
				      label, 0, 1, rows-1, rows);
	    gtk_widget_show(label);

	    entry = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      entry, 1, TAB_MAIN_COLS, 
				      rows-1, rows);
				      
            if (spec->titles[i] != NULL && *spec->titles[i] != '\0') {
		gp_string_to_entry(entry, spec->titles[i]);
            }		

	    g_signal_connect(G_OBJECT(entry), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     ed);

	    gtk_widget_show(entry);
	    ed->gpt_titles[i].widget = entry;
	}
    }

    /* specify position of plot key or legend */
    table_add_row(tbl, &rows, TAB_MAIN_COLS);
    label = gtk_label_new(_("key position"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, rows-1, rows);
    gtk_widget_show(label);

    ed->keycombo = gtk_combo_box_text_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      ed->keycombo, 1, TAB_MAIN_COLS, rows-1, rows);
    for (i=0; ; i++) {
	gp_key_spec *kp = get_keypos_spec(i);

	if (kp == NULL) {
	    break;
	}
	combo_box_append_text(ed->keycombo, _(kp->str));
	if (kp->id == spec->keyspec) {
	    kactive = i;
	}
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(ed->keycombo), kactive);
    gtk_widget_show(ed->keycombo);

    if (spec->fit != PLOT_FIT_NA) {
	/* give choice of fitted line type, if applicable */
	int semilog_ok = semilog_is_ok(spec);
	int linlog_ok = log_x_ok(spec);
	
	table_add_row(tbl, &rows, TAB_MAIN_COLS);

	label = gtk_label_new(_("fitted line"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 0, 1, rows-1, rows);
	gtk_widget_show(label);

	ed->fitcombo = gtk_combo_box_text_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  ed->fitcombo, 1, TAB_MAIN_COLS, rows-1, rows);

	if (spec->flags & GPT_TS) {
	    char *p, tmp[128];

	    for (i=0; fittype_strings[i] != NULL; i++) {
		if (i == PLOT_FIT_LOGLIN && !semilog_ok) {
		    widget_set_int(ed->fitcombo, "no-semilog", 1);
		    continue;
		} else {
		    strcpy(tmp, _(fittype_strings[i]));
		    p = strchr(tmp, ':');
		    if (p != NULL) {
			gretl_charsub(tmp, 'x', 't');
		    }
		    combo_box_append_text(ed->fitcombo, tmp);
		}
	    }	    
	} else {
	    for (i=0; fittype_strings[i] != NULL; i++) {
		if (i == PLOT_FIT_LOGLIN && !semilog_ok) {
		    widget_set_int(ed->fitcombo, "no-semilog", 1);
		    continue;
		} else if (i == PLOT_FIT_LINLOG && !linlog_ok) {
		    continue;
		} else {
		    combo_box_append_text(ed->fitcombo, _(fittype_strings[i]));
		}
	    }
	}

	gtk_combo_box_set_active(GTK_COMBO_BOX(ed->fitcombo), spec->fit);
	widget_set_int(ed->fitcombo, "oldfit", spec->fit);
	g_signal_connect(G_OBJECT(ed->fitcombo), "changed",
			 G_CALLBACK(fit_type_changed), ed);
	gtk_widget_show(ed->fitcombo);
    } 

    if (spec->flags & GPT_Y2AXIS) { 
	/* give option of forcing a single y-axis */
	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	ed->y2_check = gtk_check_button_new_with_label(_("Use only one y axis"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  ed->y2_check, 0, TAB_MAIN_COLS, 
				  rows-1, rows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ed->y2_check),
				     FALSE);
	g_signal_connect(G_OBJECT(ed->y2_check), "clicked", 
			 G_CALLBACK(toggle_axis_selection), ed);
	gtk_widget_show(ed->y2_check);
    } else {	
	/* give option of removing/adding top & right border */
	if (spec->border == GP_BORDER_DEFAULT || spec->border == 3) {
	    table_add_row(tbl, &rows, TAB_MAIN_COLS);
	    ed->border_check = gtk_check_button_new_with_label(_("Show full border"));
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      ed->border_check, 0, TAB_MAIN_COLS, 
				      rows-1, rows);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ed->border_check),
					 spec->border == GP_BORDER_DEFAULT);
	    gtk_widget_show(ed->border_check);
	}
    } 

    if (plot_has_tics(spec)) {
	/* add show grid options */
	const char *grid_opts[] = {
	    N_("horizontal"),
	    N_("vertical"),
	    N_("both")
	};
	int gridval = plotspec_gridval(spec);
	GtkWidget *combo;

	/* check button */
	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	ed->grid_check = gtk_check_button_new_with_label(_("Show grid"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), ed->grid_check, 
				  0, 1, rows-1, rows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ed->grid_check),
				     gridval > 0);
	gtk_widget_show(ed->grid_check);

	/* plus combo selector */
	ed->grid_combo = combo = gtk_combo_box_text_new();
	for (i=0; i<3; i++) {
	    combo_box_append_text(combo, _(grid_opts[i]));
	}
	gtk_widget_set_sensitive(combo, gridval > 0);
	if (gridval > 0) {
	    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), gridval - 1);
	} else {
	    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	}
	gtk_table_attach(GTK_TABLE(tbl), combo,
			 1, 2, rows-1, rows,
			 GTK_FILL, 0, 0, 0);
	gtk_widget_show(combo);
	sensitize_conditional_on(combo, ed->grid_check);
    }

    if (show_bars_check(spec)) {
	/* option to display NBER "recession bars" or similar */
	gboolean userbars = FALSE;
	GtkWidget *combo;

	/* check button */
	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	ed->bars_check = gtk_check_button_new_with_label(_("Show bars"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ed->bars_check),
				     spec->nbars > 0);
	gtk_table_attach_defaults(GTK_TABLE(tbl), ed->bars_check,
				  0, 1, rows-1, rows);
	gtk_widget_show(ed->bars_check);

	/* plus combo selector */
	ed->barscombo = combo = gtk_combo_box_text_new();
	combo_box_append_text(combo, _("NBER recessions"));
	userbars = maybe_append_user_bars_file(GTK_COMBO_BOX(combo), ed);
	combo_box_append_text(combo, _("other..."));
	/* if we got a user-specific plotbars filename, make it
	   the selected value, otherwise use the NBER file */
	ed->active_bars = userbars ? 1 : 0;
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), ed->active_bars);
	gtk_widget_set_sensitive(combo, spec->nbars > 0);
	g_signal_connect(G_OBJECT(combo), "changed", 
			 G_CALLBACK(plot_bars_changed),
			 ed);	
	gtk_table_attach(GTK_TABLE(tbl), combo,
			 1, 2, rows-1, rows,
			 GTK_FILL, 0, 0, 0);
	gtk_widget_show(combo);
	sensitize_conditional_on(combo, ed->bars_check);
    }

    /* setting of graph font */
    table_add_row(tbl, &rows, TAB_MAIN_COLS);	
    hsep = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
			      rows-1, rows);  
    gtk_widget_show(hsep);

    title = g_strdup_printf(_("font: %s"), gretl_png_font());
    button = gtk_button_new_with_label(title);
    table_add_row(tbl, &rows, TAB_MAIN_COLS);
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, TAB_MAIN_COLS, 
			      rows - 1, rows); 
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(graph_font_selector), 
		     ed);
    gtk_widget_show(button);
    g_free(title);

    /* set font as default button */
    table_add_row(tbl, &rows, TAB_MAIN_COLS);
    button = gtk_check_button_new_with_label(_("Set as default"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, TAB_MAIN_COLS, 
			      rows - 1, rows); 
    gtk_widget_show(button);
    ed->fontcheck = button;

    if (frequency_plot_code(spec->code)) {
	/* give option of setting fill color */
	GtkWidget *hsep = gtk_hseparator_new();

	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
				  rows - 1, rows);  
	gtk_widget_show(hsep);
	add_color_selector(BOXCOLOR, tbl, TAB_MAIN_COLS, &rows, 
			   ed->notebook);
    }
}

static void linetitle_callback (GtkWidget *w, plot_editor *ed)
{
    set_keyspec_sensitivity(ed);
}

struct new_line_info_ {
    plot_editor *editor;
    GtkWidget *formula_entry;
    GtkWidget *dlg;
};

typedef struct new_line_info_ new_line_info;

static void gpt_tab_new_line (plot_editor *ed, new_line_info *nlinfo) 
{
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *hbox;
    int nrows = 1;
    char label_text[32];

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(nlinfo->dlg));
    hbox = gtk_hbox_new(FALSE, 5);

    tbl = gtk_table_new(nrows, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(hbox), tbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);
    gtk_widget_show(tbl);
   
    /* identifier and formula text */
    gtk_table_resize(GTK_TABLE(tbl), ++nrows, 3);
    sprintf(label_text, _("line %d: "), ed->gui_nlines + 1);
    label = gtk_label_new(label_text);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, nrows-1, nrows);
    gtk_widget_show(label);

    label = gtk_label_new(_("formula"));
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 1, 2, nrows-1, nrows);
    gtk_widget_show(label);

    nlinfo->formula_entry = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(nlinfo->formula_entry), "");
    gtk_entry_set_width_chars(GTK_ENTRY(nlinfo->formula_entry), 32);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      nlinfo->formula_entry, 2, 3, nrows-1, nrows);
    gtk_entry_set_activates_default(GTK_ENTRY(nlinfo->formula_entry), TRUE);
    gtk_widget_show(nlinfo->formula_entry);
}

/* back up existing spec components in case changes are not
   applied and we need to restore the prior state */

static int make_labels_backup (plot_editor *ed)
{
    int err = 0;

    if (labels_not_synced(ed)) {
	ed->old_labels = plotspec_clone_labels(ed->spec, &err);
	if (!err) {
	    ed->old_n_labels = ed->spec->n_labels;
	}
    }

    return err;
}

static int make_arrows_backup (plot_editor *ed)
{
    int err = 0;

    if (arrows_not_synced(ed)) {
	ed->old_arrows = plotspec_clone_arrows(ed->spec, &err);
	if (!err) {
	    ed->old_n_arrows = ed->spec->n_arrows;
	}
    }

    return err;
}

static int make_lines_backup (plot_editor *ed)
{
    int err = 0;

    if (lines_not_synced(ed)) {
	ed->old_lines = plotspec_clone_lines(ed->spec, &err);
	if (!err) {
	    ed->old_n_lines = ed->spec->n_lines;
	}
    }

    return err;
}

static int allocate_label (plot_editor *ed)
{
    int err;

    err = make_labels_backup(ed);

    if (!err) {
	err = plotspec_add_label(ed->spec);
    }

    if (!err) {
	err = add_label_widget(ed);
	if (err) {
	    ed->spec->n_labels -= 1;
	}
    }

    if (err) {
	nomem();
    }

    return err;
}

static int allocate_arrow (plot_editor *ed)
{
    int err = 0;

    err = make_arrows_backup(ed);

    if (!err) {
	err = plotspec_add_arrow(ed->spec);
    }

    if (!err) {
	err = add_arrow_widget(ed);
	if (err) {
	    ed->spec->n_arrows -= 1;
	}
    }

    if (err) {
	nomem();
    }

    return err;
}

static int allocate_line (plot_editor *ed)
{
    int err = make_lines_backup(ed);

    if (!err) {
	err = plotspec_add_line(ed->spec);
    }

    if (!err) {
	err = add_line_widget(ed);
	if (err) {
	    ed->spec->n_lines -= 1;
	}
    }

    if (err) {
	nomem();
    }

    return err;
}

static void add_label_callback (GtkWidget *w, plot_editor *ed)
{
    int err = allocate_label(ed);

    if (!err) {
	/* re-fill the "labels" notebook page */
	gint pgnum = widget_get_int(ed->notebook, "labels_page");

	gtk_notebook_remove_page(GTK_NOTEBOOK(ed->notebook), pgnum);
	gpt_tab_labels(ed, ed->spec, pgnum);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    }
}

static void add_arrow_callback (GtkWidget *w, plot_editor *ed)
{
    int err = allocate_arrow(ed);
 
    if (!err) {
	/* re-fill the "arrows" notebook page */
	gint pgnum = widget_get_int(ed->notebook, "arrows_page");

	gtk_notebook_remove_page(GTK_NOTEBOOK(ed->notebook), pgnum);
	gpt_tab_arrows(ed, ed->spec, pgnum);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    }
}

static void set_gp_formula (char *targ, const char *s)
{
    int decom = (get_local_decpoint() == ',');
    int rem = GP_MAXFORMULA - strlen(s);
    int i = 0;

    while (*s) {
	if (decom && *s == ',') {
	    targ[i++] = '.';
	} else if (*s == '^' && rem > 1) {
	    targ[i++] = '*';
	    targ[i++] = '*';
	    rem--;
	} else {
	    targ[i++] = *s;
	}
	s++;
    }

    targ[i] = '\0';
}

static void real_add_line (GtkWidget *w, new_line_info *nlinfo)
{
    plot_editor *ed = nlinfo->editor;
    GPT_SPEC *spec = ed->spec;
    GPT_LINE *line;
    char tmp[GP_MAXFORMULA];
    const gchar *s;
    gint pgnum;
    int err = 0;

    s = gtk_entry_get_text(GTK_ENTRY(nlinfo->formula_entry));
    if (s == NULL || *s == '\0') {
	errbox(_("No formula was given"));
	return;
    }

    *tmp = '\0';
    strncat(tmp, s, GP_MAXFORMULA - 1);

    err = allocate_line(ed);
    if (err) {
	gtk_widget_destroy(nlinfo->dlg);
	return;
    }

    line = &spec->lines[spec->n_lines - 1];
    set_gp_formula(line->formula, tmp);

    line->style = GP_STYLE_LINES;
    line->type = spec->n_lines; /* assign next line style */
    line->scale = NADBL;        /* mark as a non-data line */
    line->flags = GP_LINE_USER;

    /* re-fill the "lines" notebook page */
    pgnum = widget_get_int(ed->notebook, "lines_page");
    gtk_notebook_remove_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    gpt_tab_lines(ed, spec, pgnum);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    if (spec->n_lines == 2) {
	/* user-defined line has taken the place of a potential fitted line */
	spec->fit = PLOT_FIT_NA;
	if (ed->fitcombo != NULL) {
	    gtk_widget_set_sensitive(ed->fitcombo, FALSE);
	}
    }

    gtk_widget_destroy(nlinfo->dlg);
}

static void add_line_callback (GtkWidget *w, plot_editor *ed)
{
    new_line_info *nlinfo;
    GtkWidget *hbox;
    GtkWidget *button;

    nlinfo = mymalloc(sizeof *nlinfo);
    if (nlinfo == NULL) {
	return;
    }

    nlinfo->editor = ed;
    nlinfo->dlg = gretl_dialog_new(NULL, ed->dialog, GRETL_DLG_BLOCK);
    gpt_tab_new_line(ed, nlinfo);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(nlinfo->dlg));

    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), nlinfo->dlg);

    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(real_add_line), nlinfo);
    gtk_widget_grab_default(button);

    context_help_button(hbox, GPT_ADDLINE);

    gtk_widget_show_all(nlinfo->dlg);

    free(nlinfo);
}

static void remove_line (GtkWidget *w, plot_editor *ed)
{
    if (make_lines_backup(ed) != 0) {
	nomem();
    } else {
	int pgnum = widget_get_int(ed->notebook, "lines_page");

	plotspec_delete_line(ed->spec, widget_get_int(w, "linenum"));
	ed->gui_nlines -= 1;

	/* refresh the associated notebook page */
	gtk_notebook_remove_page(GTK_NOTEBOOK(ed->notebook), pgnum);
	gpt_tab_lines(ed, ed->spec, pgnum);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    }
}

static void remove_label (GtkWidget *w, plot_editor *ed)
{
    if (make_labels_backup(ed) != 0) {
	nomem();
    } else {
	int pgnum = widget_get_int(ed->notebook, "labels_page");

	plotspec_delete_label(ed->spec, widget_get_int(w, "labelnum"));
	ed->gui_nlabels -= 1;

	/* refresh the associated notebook page */
	gtk_notebook_remove_page(GTK_NOTEBOOK(ed->notebook), pgnum);
	gpt_tab_labels(ed, ed->spec, pgnum);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    }
}

static void remove_arrow (GtkWidget *w, plot_editor *ed)
{
    if (make_arrows_backup(ed) != 0) {
	nomem();
    } else {
	int pgnum = widget_get_int(ed->notebook, "arrows_page");

	plotspec_delete_arrow(ed->spec, widget_get_int(w, "arrownum"));
	ed->gui_narrows -= 1;

	/* refresh the associated notebook page */
	gtk_notebook_remove_page(GTK_NOTEBOOK(ed->notebook), pgnum);
	gpt_tab_arrows(ed, ed->spec, pgnum);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(ed->notebook), pgnum);
    }
}

static void item_remove_button (GtkWidget *tbl, int row, 
				plot_editor *ed, int i,
				int j)
{
    GtkWidget *button; 
    const gchar *key[] = {
	"linenum",
	"labelnum",
	"arrownum"
    };
    GCallback cb[] = {
	G_CALLBACK(remove_line),
	G_CALLBACK(remove_label),
	G_CALLBACK(remove_arrow)
    };

    button = gtk_button_new_with_label(_("Remove"));
    widget_set_int(button, key[j], i);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(cb[j]), ed);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      button, 0, 1, row - 1, row);
    gtk_widget_show(button);
}

static void print_line_label (GtkWidget *tbl, int row, GPT_SPEC *spec, 
			      int i)
{
    char label_text[32];
    GtkWidget *label; 

    if (spec->code == PLOT_BOXPLOTS) {
	if (i == 0) {
	    sprintf(label_text, "%s: ", _("quartiles"));
	} else if (i == 1) {
	    sprintf(label_text, "%s: ", _("median"));
	} else if (i == 2 && spec->n_lines == 3) {
	    sprintf(label_text, "%s: ", _("mean"));
	} else if (i == 2) {
	    sprintf(label_text, "%s: ", _("lower bound"));
	} else if (i == 3) {
	    sprintf(label_text, "%s: ", _("upper bound"));
	}
    } else {
	sprintf(label_text, _("line %d: "), i+1);
    }

    label = gtk_label_new(label_text);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, row - 1, row);
    gtk_widget_show(label);
}

static void print_label_label (GtkWidget *tbl, int row, GPT_SPEC *spec, 
			       int i)
{
    char label_text[32];
    GtkWidget *label; 

    sprintf(label_text, _("label %d: "), i + 1);
    label = gtk_label_new(label_text);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, row - 1, row);
    gtk_widget_show(label);
}

static GtkWidget *print_field_label (GtkWidget *tbl, int row,
				     const gchar *text)
{
    GtkWidget *label = gtk_label_new(text);

    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 1, 2, row - 1, row);
    gtk_widget_show(label);
    
    return label;
}

static GtkWidget *dash_types_combo (void)
{
    GtkWidget *dotsel;
    GtkCellRenderer *cell;
    GtkListStore *store;
    GtkTreeIter iter;
    GdkPixbuf *pbuf;
    int i;

    store = gtk_list_store_new(1, GDK_TYPE_PIXBUF);

    for (i=0; i<2; i++) {
	gtk_list_store_append(store, &iter);
	pbuf = get_pixbuf_for_line(i);
	gtk_list_store_set(store, &iter, 0, pbuf, -1);
	g_object_unref(pbuf);
    }

    dotsel = gtk_combo_box_new_with_model(GTK_TREE_MODEL(store));
    cell = gtk_cell_renderer_pixbuf_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(dotsel), cell, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(dotsel), cell,
				   "pixbuf", 0, NULL);
    
    return dotsel;
}

static GtkWidget *point_types_combo (void)
{
    GtkWidget *ptsel;
    GtkCellRenderer *cell;
    GtkListStore *store;
    GtkTreeIter iter;
    GdkPixbuf *pbuf;
    int i;

    store = gtk_list_store_new(1, GDK_TYPE_PIXBUF);

    for (i=0; i<13; i++) {
	gtk_list_store_append(store, &iter);
	pbuf = gdk_pixbuf_from_pixdata(gppoints[i], FALSE, NULL);
	gtk_list_store_set(store, &iter, 0, pbuf, -1);
	g_object_unref(pbuf);
    }

    ptsel = gtk_combo_box_new_with_model(GTK_TREE_MODEL(store));
    cell = gtk_cell_renderer_pixbuf_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(ptsel), cell, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(ptsel), cell,
				   "pixbuf", 0, NULL);
    
    return ptsel;
}

static int line_get_point_type (GPT_LINE *line, int i)
{
    if (line->ptype > 0) {
	/* a specific point-style has been selected: convert
	   to zero-based */
	return line->ptype - 1;
    } else if (line->type == LT_AUTO) {
	/* line type is set by placement of line in plot */
	return i;
    } else {
	/* a specific line-type has been selected: give the 
	   associated point-style */
	return line->type - 1;
    }
}

#define has_point(s) (s == GP_STYLE_POINTS || s == GP_STYLE_LINESPOINTS)

static void flip_pointsel (GtkWidget *box, GtkWidget *targ)
{
    GtkWidget *label = g_object_get_data(G_OBJECT(targ), "label");
    GtkWidget *psize = g_object_get_data(G_OBJECT(targ), "psize");
    GtkWidget *pslbl = g_object_get_data(G_OBJECT(targ), "pslbl");
    gchar *s = combo_box_get_active_text(box);
    int hp = has_point(gp_style_index_from_display_name(s));

    gtk_widget_set_sensitive(targ, hp);
    gtk_widget_set_sensitive(label, hp);
    gtk_widget_set_sensitive(psize, hp);
    gtk_widget_set_sensitive(pslbl, hp);
}

static GtkWidget *scroller_page (GtkWidget *vbox)
{
    GtkWidget *scroller;
	
    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_NEVER, 
				   GTK_POLICY_AUTOMATIC);
    gtk_widget_show(scroller);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), 
					  vbox);
    return scroller;
}

static int show_axis_chooser (GPT_SPEC *spec)
{
    int i, s = 0;

    if (spec->code == PLOT_REGULAR && spec->n_lines > 1) {
	s = 1;
    } else {
	for (i=0; i<spec->n_lines; i++) {
	    if (spec->lines[i].yaxis == 2) {
		s = 1;
		break;
	    }
	}
    }

    return s;
}

static int gp_style_index (int sty, GList *list)
{
    gp_style_spec *spec;
    GList *mylist = list;
    int i = 0;

    while (mylist != NULL) {
	spec = (gp_style_spec *) mylist->data;
	if (sty == spec->id) {
	    return i;
	}
	i++;
	mylist = mylist->next;
    }

    return 0;
}

static GtkWidget *gpt_hboxit (GtkWidget *w)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 0);

    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 0);
    return hbox;
}

static void gpt_linetab_add2 (GtkWidget *label,
			      GtkWidget *control,
			      GtkWidget *tbl,
			      int col, int row)
{
    GtkWidget *hb;
	    
    gtk_widget_show(label);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach(GTK_TABLE(tbl), label,
		     col-2, col-1, row-1, row,
		     0, 0, 10, 0);
    hb = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hb), control, 0, 0, 0);
    gtk_widget_show_all(hb);
    gtk_table_attach_defaults(GTK_TABLE(tbl), hb,
			      col-1, col, row-1, row);
}

static void set_combo_box_strings_from_stylist (GtkWidget *box, 
						GList *list)
{
    gp_style_spec *spec;
    GList *mylist = list;

    while (mylist != NULL) {
	spec = (gp_style_spec *) mylist->data;
	combo_box_append_text(box, _(spec->trname));
	mylist = mylist->next;
    }
}

#define IRF_plot(c) (c == PLOT_IRFBOOT || c == PLOT_MULTI_IRF)

static GList *add_style_spec (GList *list, int t)
{
    return g_list_append(list, get_style_spec(t));
}

#define line_is_formula(l) (l->formula[0] != '\0' || (l->flags & GP_LINE_USER))

static void gpt_tab_lines (plot_editor *ed, GPT_SPEC *spec, int ins)
{
    GtkWidget *notebook = ed->notebook;
    GtkWidget *color_label = NULL;
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *hbox, *sep;
    GtkWidget *page;
    GList *stylist = NULL;
    int axis_chooser;
    int i, nrows, ncols = 5;
    int pgnum = -1;

    axis_chooser = show_axis_chooser(spec);

    if (frequency_plot_code(spec->code)) {
	stylist = add_style_spec(stylist, GP_STYLE_BOXES);
    }

    if (spec->flags & GPT_TS) {
	stylist = add_style_spec(stylist, GP_STYLE_LINES);
	stylist = add_style_spec(stylist, GP_STYLE_POINTS);
    } else {
	stylist = add_style_spec(stylist, GP_STYLE_POINTS);
	stylist = add_style_spec(stylist, GP_STYLE_LINES);
    }

    stylist = add_style_spec(stylist, GP_STYLE_LINESPOINTS); 
    stylist = add_style_spec(stylist, GP_STYLE_IMPULSES);
    stylist = add_style_spec(stylist, GP_STYLE_DOTS);
    stylist = add_style_spec(stylist, GP_STYLE_STEPS);

    if (!frequency_plot_code(spec->code)) {
	stylist = add_style_spec(stylist, GP_STYLE_BOXES);
    }

    vbox = gp_dialog_vbox();
    gtk_widget_show(vbox);

    label = gtk_label_new(_("Lines"));
    gtk_widget_show(label);

    if (ed->gui_nlines > 4) {
	page = scroller_page(vbox);
    } else {
	page = vbox;
    }

    if (ins > 0) {
	pgnum = gtk_notebook_insert_page(GTK_NOTEBOOK(notebook), page, label, ins); 
    } else {
	pgnum = gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
    }  

    widget_set_int(notebook, "lines_page", pgnum);

    nrows = 1;
    tbl = gp_dialog_table(nrows, ncols, vbox);
    gtk_widget_show(tbl);
   
    for (i=0; i<ed->gui_nlines; i++) {
	GPT_LINE *line = &spec->lines[i];
	GtkWidget *ptsel = NULL;
	int label_done = 0;

	hbox = NULL;

	if (line_is_formula(line)) {
	    gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	    print_line_label(tbl, nrows, spec, i);
	    label_done = 1;
	    print_field_label(tbl, nrows, _("formula"));
	    ed->lineformula[i] = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), ed->lineformula[i],
				      2, ncols, nrows-1, nrows);
	    strip_lr(line->formula);
	    gp_string_to_entry(ed->lineformula[i], line->formula);
	    if (i == 1 && (spec->flags & GPT_AUTO_FIT)) {
		/* fitted formula: not GUI-editable */
		gtk_widget_set_sensitive(ed->lineformula[i], FALSE);
		ed->fitformula = ed->lineformula[i];
	    } else {
		g_signal_connect(G_OBJECT(ed->lineformula[i]), "activate", 
				 G_CALLBACK(apply_gpt_changes), 
				 ed);
	    }
	    gtk_widget_show(ed->lineformula[i]);
	}

	/* identifier and key or legend text */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	if (!label_done) {
	    print_line_label(tbl, nrows, spec, i);
	}
	if (line->flags & GP_LINE_USER) {
	    item_remove_button(tbl, nrows, ed, i, GUI_LINE);
	}
	print_field_label(tbl, nrows, _("legend"));
	ed->linetitle[i] = gtk_entry_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), ed->linetitle[i],
				  2, ncols, nrows-1, nrows);
	strip_lr(line->title);
	gp_string_to_entry(ed->linetitle[i], line->title);
	g_signal_connect(G_OBJECT(ed->linetitle[i]), "changed", 
			 G_CALLBACK(linetitle_callback), ed);
	g_signal_connect(G_OBJECT(ed->linetitle[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), ed);
	gtk_widget_show(ed->linetitle[i]);
	if (i == 1 && (spec->flags & GPT_AUTO_FIT)) {
	    ed->fitlegend = ed->linetitle[i];
	}

	if (line_is_formula(line) || 
	    line->style == GP_STYLE_CANDLESTICKS) {
	    goto line_width_adj;
	}

	/* line type (lines, points, etc.) */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	label = gtk_label_new(_("type"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 1, 2,
				  nrows-1, nrows);
	gtk_widget_show(label);

	ed->stylecombo[i] = gtk_combo_box_text_new();
	hbox = gpt_hboxit(ed->stylecombo[i]);

	/* the errorbars and filledcurves styles are not exchangeable
	   with other styles
	*/
	if (line->style == GP_STYLE_ERRORBARS || 
	    line->style == GP_STYLE_FILLEDCURVE) {
	    GList *altsty = NULL;

	    altsty = add_style_spec(altsty, GP_STYLE_ERRORBARS); 
	    altsty = add_style_spec(altsty, GP_STYLE_FILLEDCURVE);
	    set_combo_box_strings_from_stylist(ed->stylecombo[i], altsty);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(ed->stylecombo[i]),
				     gp_style_index(line->style, altsty));
	    g_list_free(altsty);
	} else {
	    int lt = gp_style_index(line->style, stylist);

	    set_combo_box_strings_from_stylist(ed->stylecombo[i], stylist);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(ed->stylecombo[i]), lt);
	    ptsel = point_types_combo();
	    g_object_set_data(G_OBJECT(ed->stylecombo[i]), "pointsel", ptsel);
	    g_signal_connect(G_OBJECT(ed->stylecombo[i]), "changed",
			     G_CALLBACK(flip_pointsel), ptsel);
	} 

	if (IRF_plot(spec->code) && (line->style == GP_STYLE_FILLEDCURVE ||
				     line->style == GP_STYLE_ERRORBARS)) {
	    /* styles not exchangeable in this context */
	    gtk_widget_set_sensitive(ed->stylecombo[i], FALSE);
	}

	gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, ncols,
				  nrows-1, nrows);
	gtk_widget_show_all(hbox);	

    line_width_adj:

	/* line-width adjustment */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	print_field_label(tbl, nrows, _("line width"));
	hbox = gtk_hbox_new(FALSE, 5);
	ed->linewidth[i] = gtk_spin_button_new_with_range(0.5, 6.0, 0.5);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(ed->linewidth[i]),
				  line->width);
	gtk_box_pack_start(GTK_BOX(hbox), ed->linewidth[i], FALSE, FALSE, 0);
	gtk_widget_show(ed->linewidth[i]);
	gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, ncols-2, 
				  nrows-1, nrows);
	if (spec->lines[i].style == GP_STYLE_FILLEDCURVE) {
	    gtk_widget_set_sensitive(ed->linewidth[i], FALSE);
	} else {
	    g_signal_connect(G_OBJECT(ed->linewidth[i]), "activate", 
			     G_CALLBACK(apply_gpt_changes), ed);
	}

	ed->colorsel[i] = NULL;
	color_label = NULL;

	/* line color adjustment */
	if (i < 6 && !frequency_plot_code(spec->code)) {
	    ed->colorsel[i] = line_color_button(spec, i);
	}

	if (ed->colorsel[i] != NULL) {
	    color_label = gtk_label_new(_("color"));
	    gpt_linetab_add2(color_label, ed->colorsel[i],
			     tbl, ncols, nrows);
	}

	if (line->flags & GP_LINE_USER) {
	    /* offer dotted option */
	    GtkWidget *dotcombo = dash_types_combo();
	    int dotted = (line->type == 0);

	    gtk_combo_box_set_active(GTK_COMBO_BOX(dotcombo), 
				     dotted ? 1 : 0);
	    gtk_widget_show(dotcombo);
	    widget_set_int(dotcombo, "linenum", i);
	    if (ed->colorsel[i] != NULL) {
		g_object_set_data(G_OBJECT(dotcombo), "colorsel", ed->colorsel[i]);
		g_object_set_data(G_OBJECT(dotcombo), "color-label", color_label);
		g_signal_connect(G_OBJECT(dotcombo), "changed", 
				 G_CALLBACK(dot_callback), spec);
	    }
	    gtk_box_pack_start(GTK_BOX(hbox), dotcombo, FALSE, FALSE, 5);
	}

	if (hbox != NULL) {
	    gtk_widget_show(hbox);
	    hbox = NULL;
	}

	if (ptsel != NULL) {
	    /* point type and size adjustment */
	    int pt = line_get_point_type(line, i);
	    int hp = has_point(line->style);

	    gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	    label = print_field_label(tbl, nrows, _("point"));
	    gtk_widget_set_sensitive(label, hp);
	    g_object_set_data(G_OBJECT(ptsel), "label", label);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(ptsel), pt);
	    hbox = gpt_hboxit(ptsel);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, ncols-2, 
				      nrows-1, nrows);
	    gtk_widget_show_all(hbox);
	    
	    label = gtk_label_new( _("size"));
	    g_object_set_data(G_OBJECT(ptsel), "pslbl", label);
	    gtk_widget_set_sensitive(label, hp);
	    ed->pointsize[i] = gtk_spin_button_new_with_range(0.5, 6.0, 0.5);
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(ed->pointsize[i]),
				      line->pscale);
	    gpt_linetab_add2(label, ed->pointsize[i], tbl, ncols, nrows);
	    gtk_widget_set_sensitive(ptsel, hp);
	    gtk_widget_set_sensitive(ed->pointsize[i], hp);
	    g_object_set_data(G_OBJECT(ptsel), "psize", ed->pointsize[i]);
	}

	if (axis_chooser) {
	    /* scale factor for data? */
	    GtkWidget *hbox;
	    
	    gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	    print_field_label(tbl, nrows, _("data scale"));
	    ed->linescale[i] = gtk_entry_new();
	    gtk_entry_set_max_length(GTK_ENTRY(ed->linescale[i]), 6);
	    double_to_gp_entry(line->scale, ed->linescale[i]);
	    gtk_entry_set_width_chars(GTK_ENTRY(ed->linescale[i]), 6);
	    if (line_is_formula(line)) {
		gtk_widget_set_sensitive(ed->linescale[i], FALSE);
	    } else {
		g_signal_connect(G_OBJECT(ed->linescale[i]), "activate", 
				 G_CALLBACK(apply_gpt_changes), 
				 ed);
	    }
	    hbox = gpt_hboxit(ed->linescale[i]);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, ncols-2,
				      nrows-1, nrows);
	    gtk_widget_show_all(hbox);

	    /* use left or right y axis? */
	    label = gtk_label_new( _("y axis"));
	    ed->yaxiscombo[i] = gtk_combo_box_text_new();
	    combo_box_append_text(ed->yaxiscombo[i], _("left"));
	    combo_box_append_text(ed->yaxiscombo[i], _("right"));
	    gtk_combo_box_set_active(GTK_COMBO_BOX(ed->yaxiscombo[i]), 
				     (line->yaxis == 1)? 0 : 1);
	    gpt_linetab_add2(label, ed->yaxiscombo[i], tbl, ncols, nrows);
	}	

	/* separator */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	sep = gtk_hseparator_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), sep, 0, ncols, 
				  nrows-1, nrows);
	gtk_widget_show(sep);
    } /* end iteration over lines in plot */

    if ((spec->code == PLOT_REGULAR || spec->code == PLOT_CURVE)
	&& spec->n_lines < 8) {
	/* button for adding a line (formula) */
	GtkWidget *add_button;

	gtk_table_resize(GTK_TABLE(tbl), ++nrows, ncols);
	add_button = gtk_button_new_with_label(_("Add line..."));
	g_signal_connect(G_OBJECT(add_button), "clicked", 
			 G_CALLBACK(add_line_callback), 
			 ed);
	gtk_widget_show(add_button);
	gtk_table_attach_defaults(GTK_TABLE(tbl), add_button, 
				  0, 1, nrows-1, nrows);
    }

    g_list_free(stylist);
}

static void label_pos_to_entry (double *pos, GtkWidget *w)
{
    if (!na(pos[0]) && !na(pos[1])) {
	gchar *s = g_strdup_printf("%.10g %.10g", pos[0], pos[1]);

	gtk_entry_set_text(GTK_ENTRY(w), s);
	g_free(s);
    } else {
	gtk_entry_set_text(GTK_ENTRY(w), "");
    }
}

static void gpt_tab_labels (plot_editor *ed, GPT_SPEC *spec, int ins) 
{
    GtkWidget *notebook = ed->notebook;
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *page, *sep;
    int i, j, nrows;
    png_plot *plot = (png_plot *) spec->ptr;
    int pgnum = -1;

    vbox = gp_dialog_vbox();
    gtk_widget_show(vbox);

    label = gtk_label_new(_("Labels"));
    gtk_widget_show(label);

    if (ed->gui_nlabels > 4) {
	page = scroller_page(vbox);
    } else {
	page = vbox;
    }

    if (ins > 0) {
	pgnum = gtk_notebook_insert_page(GTK_NOTEBOOK(notebook), page, label, ins); 
    } else {
	pgnum = gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
    }  

    widget_set_int(notebook, "labels_page", pgnum);
 
    nrows = 1;
    tbl = gp_dialog_table(nrows, 3, vbox);
    gtk_widget_show(tbl);
   
    for (i=0; i<ed->gui_nlabels; i++) {
	GtkWidget *hbox, *button, *image;
	GdkPixbuf *icon;

	/* label text */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, 3);
	print_label_label(tbl, nrows, spec, i);
	print_field_label(tbl, nrows, _("text"));
	ed->labeltext[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(ed->labeltext[i]), PLOT_LABEL_TEXT_LEN);
	gp_string_to_entry(ed->labeltext[i], spec->labels[i].text);
	g_signal_connect(G_OBJECT(ed->labeltext[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 ed);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  ed->labeltext[i], 2, 3, nrows-1, nrows);
	gtk_widget_show(ed->labeltext[i]);

	/* label placement */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, 3);

	item_remove_button(tbl, nrows, ed, i, GUI_LABEL);

	print_field_label(tbl, nrows, _("position (X Y)"));

	/* holder for entry and button */
	hbox = gtk_hbox_new(FALSE, 5);

	/* entry for coordinates */
	ed->labelpos[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(ed->labelpos[i]), PLOT_POSITION_LEN);
	label_pos_to_entry(spec->labels[i].pos, ed->labelpos[i]);
	g_signal_connect(G_OBJECT(ed->labelpos[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 ed);
	gtk_container_add(GTK_CONTAINER(hbox), ed->labelpos[i]);
	gtk_widget_show(ed->labelpos[i]);

	if (plot_is_mouseable(plot)) {
	    /* button to invoke mouse-assisted placement */
	    button = gtk_button_new();
	    g_object_set_data(G_OBJECT(button), "pos_entry", ed->labelpos[i]);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(plot_position_click), spec->ptr);
	    icon = gdk_pixbuf_new_from_xpm_data((const char **) mini_mouse_xpm);
	    image = gtk_image_new_from_pixbuf(icon);
	    g_object_unref(icon);
	    gtk_widget_set_size_request(button, 32, 26);
	    gtk_container_add(GTK_CONTAINER(button), image);
	    gtk_container_add(GTK_CONTAINER(hbox), button);
	    gtk_widget_show_all(button);
	}

	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  hbox, 2, 3, nrows-1, nrows);
	gtk_widget_show(hbox);

	/* label justification */
	gtk_table_resize(GTK_TABLE(tbl), ++nrows, 3);

	print_field_label(tbl, nrows, _("justification"));

	ed->labeljust[i] = gtk_combo_box_text_new();
	for (j=0; j<3; j++) {
	    combo_box_append_text(ed->labeljust[i],
				  _(gp_justification_string(j)));
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(ed->labeljust[i]), 
				 spec->labels[i].just);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  ed->labeljust[i], 2, 3, nrows-1, nrows);
	gtk_widget_show_all(ed->labeljust[i]);

	if (i < ed->gui_nlabels - 1 || spec->n_labels < 8) {
	    /* separator */
	    gtk_table_resize(GTK_TABLE(tbl), ++nrows, 3);
	    sep = gtk_hseparator_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), sep, 0, 3, 
				      nrows-1, nrows);
	    gtk_widget_show(sep);
	}
    }

    if (spec->n_labels < 8) {
	/* button for adding a label */
	GtkWidget *button;

	gtk_table_resize(GTK_TABLE(tbl), ++nrows, 3);
	button = gtk_button_new_with_label(_("Add..."));
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(add_label_callback), 
			 ed);
	gtk_widget_show(button);
	gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, 1, 
				  nrows-1, nrows);
    }
}

static void arrow_pos_to_entry (GPT_ARROW *arrow, GtkWidget *w, int j)
{
    double x = (j == 0)? arrow->x0 : arrow->x1;
    double y = (j == 0)? arrow->y0 : arrow->y1;

    if (!na(x) && !na(y)) {
	gchar *s = g_strdup_printf("%.10g %.10g", x, y);

	gtk_entry_set_text(GTK_ENTRY(w), s);
	g_free(s);
    } else {
	gtk_entry_set_text(GTK_ENTRY(w), "");
    }
}

#define ARROW_MAX 6

static void gpt_tab_arrows (plot_editor *ed, GPT_SPEC *spec, int ins) 
{
    GtkWidget *notebook = ed->notebook;
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *page, *sep;
    int i, j, r, nrows;
    png_plot *plot = (png_plot *) spec->ptr;
    int nsep = 0, pgnum = -1;

    vbox = gp_dialog_vbox();
    gtk_widget_show(vbox);

    label = gtk_label_new(_("Arrows"));
    gtk_widget_show(label);
    page = vbox;

    if (ins > 0) {
	pgnum = gtk_notebook_insert_page(GTK_NOTEBOOK(notebook), page, label, ins); 
    } else {
	pgnum = gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
    }  

    widget_set_int(notebook, "arrows_page", pgnum);

    if (ed->gui_narrows > 1) {
	nsep = ed->gui_narrows - 1;
    }
    nrows = 1 + 4 * ed->gui_narrows + nsep;
    tbl = gp_dialog_table(nrows, 3, vbox);
    r = 1;
   
    for (i=0; i<ed->gui_narrows; i++) {
	GtkWidget *button, *image, *chk;
	GdkPixbuf *icon;

	item_remove_button(tbl, r, ed, i, GUI_ARROW);

	/* entries and buttons for (start, stop) coordinates */
	for (j=0; j<2; j++) {
	    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
	    GtkWidget *apos;

	    if (j == 0) {
		label = gtk_label_new(_("from (X Y)"));
	    } else {
		label = gtk_label_new(_("to (X Y)"));
	    }
	    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 1, 2, r-1, r);

	    apos = gtk_entry_new();
	    gtk_entry_set_max_length(GTK_ENTRY(apos), PLOT_POSITION_LEN);
	    arrow_pos_to_entry(&spec->arrows[i], apos, j);
	    g_signal_connect(G_OBJECT(apos), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     ed);
	    gtk_container_add(GTK_CONTAINER(hbox), apos);

	    if (j == 0) {
		ed->arrowpos[i] = apos;
	    } else {
		g_object_set_data(G_OBJECT(ed->arrowpos[i]), "pos_2", apos);
	    }

	    if (plot_is_mouseable(plot)) {
		/* button to invoke mouse-assisted placement */
		button = gtk_button_new();
		g_object_set_data(G_OBJECT(button), "pos_entry", apos);
		g_signal_connect(G_OBJECT(button), "clicked",
				 G_CALLBACK(plot_position_click), spec->ptr);
		icon = gdk_pixbuf_new_from_xpm_data((const char **) mini_mouse_xpm);
		image = gtk_image_new_from_pixbuf(icon);
		g_object_unref(icon);
		gtk_widget_set_size_request(button, 32, 26);
		gtk_container_add(GTK_CONTAINER(button), image);
		gtk_container_add(GTK_CONTAINER(hbox), button);
		gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, 3, r-1, r);
		r++;
	    }
	}

	/* arrow head selector */
	chk = gtk_check_button_new_with_label(_("arrow has head"));
	g_object_set_data(G_OBJECT(ed->arrowpos[i]), "arrow_check", chk);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), 
				     (spec->arrows[i].flags & GP_ARROW_HEAD)?
				     TRUE : FALSE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), chk, 2, 3, r-1, r);
	r++;

	/* dotted line selector */
	chk = gtk_check_button_new_with_label(_("line is dotted"));
	g_object_set_data(G_OBJECT(ed->arrowpos[i]), "dots_check", chk);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), 
				     (spec->arrows[i].flags & GP_ARROW_DOTS)?
				     TRUE : FALSE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), chk, 2, 3, r-1, r);
	r++;	

	if (i < ed->gui_narrows - 1) {
	    /* separator */
	    sep = gtk_hseparator_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), sep, 0, 3, r-1, r);
	    r++;
	}
    }

    if (spec->n_arrows < ARROW_MAX) {
	/* button for adding an arrow */
	GtkWidget *button = gtk_button_new_with_label(_("Add..."));

	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(add_arrow_callback), ed);
	gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, 1, r-1, r);
    }

    gtk_widget_show_all(tbl);
}

static void gpt_tab_XY (plot_editor *ed, GPT_SPEC *spec, gint axis) 
{
    png_plot *plot = (png_plot *) spec->ptr;
    GtkWidget *notebook = ed->notebook;
    GtkWidget *b1, *b2, *entry, *vbox, *tbl;
    GtkWidget *label = NULL;
    char *labelstr = NULL;
    int i, nrows;
   
    if (axis == 0) {
	labelstr = _("X-axis");
    } else if (axis == 1) {
	labelstr = _("Y-axis");
    } else if (axis == 2) {
	labelstr = _("Y2-axis");
    } else {
	return;
    }

    vbox = gp_page_vbox(notebook, labelstr);

    nrows = 1;
    tbl = gp_dialog_table(nrows, 2, vbox);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (ed->gpt_titles[i].tab == 1 + axis) {
	    gtk_table_resize(GTK_TABLE(tbl), ++nrows, 2);
            
	    label = gtk_label_new(_(ed->gpt_titles[i].desc));
	    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
				      nrows-1, nrows);
	    gtk_widget_show(label);

	    entry = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, 
				      nrows-1, nrows);
	    gp_string_to_entry(entry, spec->titles[i]);
	    g_signal_connect(G_OBJECT(entry), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     ed);
	    gtk_widget_show(entry);
	    ed->gpt_titles[i].widget = entry;
	}
    } 

    if (spec->code != PLOT_REGULAR && 
	spec->code != PLOT_CURVE &&
	spec->code != PLOT_MANY_TS) {
	return;
    }

    if (axis == 0 && (spec->flags & GPT_TIMEFMT)) {
	return;
    }

    ed->axis_range[axis].ID = axis;

    /* axis range: "auto" button */
    nrows += 3;
    gtk_table_resize(GTK_TABLE(tbl), nrows, 2);

    label = gtk_label_new("");
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      nrows-3, nrows-2);
    gtk_widget_show(label);
    b1 = gtk_radio_button_new_with_label(NULL, _("auto axis range"));
    widget_set_int(b1, "axis", axis);
    g_signal_connect(G_OBJECT(b1), "clicked",
		     G_CALLBACK(flip_manual_range), ed);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 1, 
			      nrows-2, nrows-1);
    gtk_widget_show(b1);
    ed->axis_range[axis].isauto = b1;

    /* axis range: manual range button */
    b2 = gtk_radio_button_new_with_label(gtk_radio_button_get_group 
					 (GTK_RADIO_BUTTON(b1)),
					 _("manual range:")); 
    widget_set_int(b2, "axis", axis);
    g_signal_connect(G_OBJECT(b2), "clicked",
		     G_CALLBACK(flip_manual_range), ed);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 1, 
			      nrows-1, nrows);
    gtk_widget_show(b2);

    /* axis range min. entry */
    nrows++;
    label = gtk_label_new(_("minimum"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, nrows-1, nrows);
    gtk_widget_show(label);
    gtk_table_resize(GTK_TABLE(tbl), nrows, 2);
    entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, 
			      nrows-1, nrows);
    gtk_entry_set_text(GTK_ENTRY(entry), "");
    g_signal_connect(G_OBJECT(entry), "activate", 
		     G_CALLBACK(apply_gpt_changes), 
		     ed);
    gtk_widget_show(entry);
    ed->axis_range[axis].min = entry;

    /* axis range max. entry */
    nrows++;
    label = gtk_label_new(_("maximum"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      nrows-1, nrows);
    gtk_widget_show(label);
    gtk_table_resize(GTK_TABLE(tbl), nrows, 2);
    entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, 
			      nrows-1, nrows);
    gtk_entry_set_text(GTK_ENTRY(entry), "");
    g_signal_connect(G_OBJECT(entry), "activate", 
		     G_CALLBACK(apply_gpt_changes), 
		     ed);
    gtk_widget_show(entry);
    ed->axis_range[axis].max = entry;

    if (na(spec->range[axis][0])) {
	gtk_widget_set_sensitive(ed->axis_range[axis].min, FALSE);
	gtk_widget_set_sensitive(ed->axis_range[axis].max, FALSE);
    } else {
	double_to_gp_entry(spec->range[axis][0], ed->axis_range[axis].min);
	double_to_gp_entry(spec->range[axis][1], ed->axis_range[axis].max);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ed->axis_range[axis].isauto), 
				     FALSE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    }

    if (spec->code != PLOT_REGULAR) {
	return;
    }

    /* axis scale: linear vs log? */

    if ((axis == 0 && plot_get_xmin(plot) >= 0.0) ||
	(axis == 1 && plot_get_ymin(plot) >= 0.0)) {
	GtkWidget *combo;
	GList *strs = NULL;

	strs = g_list_append(strs, "e");
	strs = g_list_append(strs, "2");
	strs = g_list_append(strs, "10");

	nrows += 3;
	gtk_table_resize(GTK_TABLE(tbl), nrows, 2);

	label = gtk_label_new("");
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
				  nrows-3, nrows-2);
	gtk_widget_show(label);
	b1 = gtk_radio_button_new_with_label(NULL, _("linear scale"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 1, 
				  nrows-2, nrows-1);
	gtk_widget_show(b1);

	b2 = gtk_radio_button_new_with_label(gtk_radio_button_get_group 
					     (GTK_RADIO_BUTTON(b1)),
					     _("logarithmic scale, base:")); 
	gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 1, 
				  nrows-1, nrows);
	gtk_widget_show(b2);    

	combo = combo_box_text_new_with_entry();
	set_combo_box_strings_from_list(combo, strs);
	g_list_free(strs);
	gtk_table_attach_defaults(GTK_TABLE(tbl), combo, 1, 2, 
				  nrows-1, nrows);
	entry = gtk_bin_get_child(GTK_BIN(combo));
	g_signal_connect(G_OBJECT(entry), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 ed);
	gtk_widget_show(combo);
	ed->axis_range[axis].lbase = entry;

	if (spec->logbase[axis] > 1.1) {
	    gchar *txt = g_strdup_printf("%g", spec->logbase[axis]);

	    gtk_entry_set_text(GTK_ENTRY(entry), txt);
	    g_free(txt);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
	} else {
	    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	    gtk_widget_set_sensitive(entry, FALSE);
	}

	sensitize_conditional_on(entry, b2);
    }
}

static void gpt_tab_palette (GtkWidget *notebook) 
{
    GtkWidget *vbox, *label, *tbl;
    int i, rows = 1;

    vbox = gp_page_vbox(notebook, _("Palette"));
    tbl = gp_dialog_table(1, 2, vbox);

    label = gtk_label_new(_("These colors will be used unless overridden\n"
			    "by graph-specific choices\n"));
    gtk_table_attach(GTK_TABLE(tbl), label, 0, 2, 0, 1, 0, 0, 0, 0);
    gtk_widget_show(label);

    for (i=0; i<BOXCOLOR; i++) {
	add_color_selector(i, tbl, 2, &rows, notebook);
    }

    gtk_widget_show(tbl);
}

static void plot_editor_destroy (plot_editor *ed)
{
    if (restore_lines(ed)) {
	/* undo any unapplied changes to lines */
	free(ed->spec->lines);
	ed->spec->lines = ed->old_lines;
	ed->spec->n_lines = ed->old_n_lines;
    }

    if (restore_labels(ed)) {
	/* undo any unapplied changes to labels */
	free(ed->spec->labels);
	ed->spec->labels = ed->old_labels;
	ed->spec->n_labels = ed->old_n_labels;
    }  

    if (restore_arrows(ed)) {
	/* undo any unapplied changes to arrows */
	free(ed->spec->arrows);
	ed->spec->arrows = ed->old_arrows;
	ed->spec->n_arrows = ed->old_n_arrows;
    } 

    if (ed->user_barsfile != NULL) {
	remember_user_bars_file(ed->user_barsfile);
	g_free(ed->user_barsfile);
    }

    free(ed->linetitle);
    free(ed->lineformula);
    free(ed->stylecombo);
    free(ed->yaxiscombo);
    free(ed->linescale);
    free(ed->linewidth);
    free(ed->colorsel);
    free(ed->pointsize);

    free(ed->labeltext);
    free(ed->labelpos);
    free(ed->labeljust);

    free(ed->arrowpos);

    free(ed);
}

static GtkWidget **widget_array_new (int n, int *err)
{
    GtkWidget **pw = malloc(n * sizeof *pw);
    int i;

    if (pw != NULL) {
	for (i=0; i<n; i++) {
	    pw[i] = NULL;
	}
    } else {
	*err = E_ALLOC;
    }

    return pw;
}

static GtkWidget **widget_array_expand (GtkWidget ***ppw, int n, int *err)
{
    GtkWidget **tmp;
    int i;

    tmp = realloc(*ppw, n * sizeof *tmp);

    if (tmp != NULL) {
	for (i=0; i<n; i++) {
	    tmp[i] = NULL;
	}
	*ppw = tmp;
    } else {
	*err = E_ALLOC;
    }

    return *ppw;
}

static int add_line_widget (plot_editor *ed)
{
    int n = ed->gui_nlines + 1;
    int err = 0;

    ed->linetitle   = widget_array_expand(&ed->linetitle, n, &err);
    ed->lineformula = widget_array_expand(&ed->lineformula, n, &err);
    ed->stylecombo  = widget_array_expand(&ed->stylecombo, n, &err);
    ed->yaxiscombo  = widget_array_expand(&ed->yaxiscombo, n, &err);
    ed->linescale   = widget_array_expand(&ed->linescale, n, &err);
    ed->linewidth   = widget_array_expand(&ed->linewidth, n, &err);
    ed->colorsel    = widget_array_expand(&ed->colorsel, n, &err);
    ed->pointsize   = widget_array_expand(&ed->pointsize, n, &err);
    
    if (!err) {
	ed->gui_nlines = n;
    }

    return err;
}

static int allocate_line_widgets (plot_editor *ed, int n)
{
    int err = 0;

    if (n > 0) {
	ed->linetitle   = widget_array_new(n, &err);
	ed->lineformula = widget_array_new(n, &err);
	ed->stylecombo  = widget_array_new(n, &err);
	ed->yaxiscombo  = widget_array_new(n, &err);
	ed->linescale   = widget_array_new(n, &err);
	ed->linewidth   = widget_array_new(n, &err);
	ed->colorsel    = widget_array_new(n, &err);
	ed->pointsize   = widget_array_new(n, &err);
	if (!err) {
	    ed->gui_nlines = n;
	}
    }

    return err;
}

static int add_label_widget (plot_editor *ed)
{
    int n = ed->gui_nlabels + 1;
    int err = 0;

    ed->labeltext = widget_array_expand(&ed->labeltext, n, &err);
    ed->labelpos  = widget_array_expand(&ed->labelpos, n, &err);
    ed->labeljust = widget_array_expand(&ed->labeljust, n, &err);
    
    if (!err) {
	ed->gui_nlabels = n;
    }

    return err;
}

static int allocate_label_widgets (plot_editor *ed, int n)
{
    int err = 0;

    if (n > 0) {
	ed->labeltext = widget_array_new(n, &err);
	ed->labelpos  = widget_array_new(n, &err);
	ed->labeljust = widget_array_new(n, &err);
	if (!err) {
	    ed->gui_nlabels = n;
	}
    }

    return err;
}

static int add_arrow_widget (plot_editor *ed)
{
    int n = ed->gui_narrows + 1;
    int err = 0;

    ed->arrowpos = widget_array_expand(&ed->arrowpos, n, &err);
    
    if (!err) {
	ed->gui_narrows = n;
    }

    return err;
}

static int allocate_arrow_widgets (plot_editor *ed, int n)
{
    int err = 0;

    if (n > 0) {
	ed->arrowpos = widget_array_new(n, &err);
	if (!err) {
	    ed->gui_narrows = n;
	}
    }

    return err;
}

static plot_editor *plot_editor_new (GPT_SPEC *spec)
{
    plot_editor *ed;
    int i;

    ed = malloc(sizeof *ed);
    if (ed == NULL) {
	return NULL;
    }

    ed->spec = spec;
    ed->user_barsfile = NULL;
    ed->active_bars = -1;

    ed->dialog = NULL;

    ed->linetitle = NULL;
    ed->lineformula = NULL;
    ed->stylecombo = NULL;
    ed->yaxiscombo = NULL;
    ed->linescale = NULL;
    ed->linewidth = NULL;
    ed->colorsel = NULL;
    ed->pointsize = NULL;

    ed->labeltext = NULL;
    ed->labelpos = NULL;
    ed->labeljust = NULL;

    ed->arrowpos = NULL;

    /* also set to NULL the widgets which may or may not
       actually be used in constructing the plot editor 
       dialog 
    */
    ed->fitformula = NULL;
    ed->fitlegend = NULL;
    ed->keycombo = NULL;
    ed->fitcombo = NULL;
    ed->border_check = NULL;
    ed->grid_check = NULL;
    ed->grid_combo = NULL;
    ed->y2_check = NULL;
    ed->bars_check = NULL;
    ed->fontcheck = NULL;
    ed->barscombo = NULL;

    ed->gui_nlines = ed->gui_nlabels = ed->gui_narrows = 0;

    if (spec->code != PLOT_STACKED_BAR) {
	if (allocate_line_widgets(ed, spec->n_lines)) {
	    plot_editor_destroy(ed);
	    return NULL;
	}
    }   

    if (allocate_label_widgets(ed, spec->n_labels)) {
	plot_editor_destroy(ed);
	return NULL;
    }   

    if (allocate_arrow_widgets(ed, spec->n_arrows)) {
	plot_editor_destroy(ed);
	return NULL;
    }     

    for (i=0; i<MAX_AXES; i++) {
	ed->axis_range[i].isauto = NULL;
	ed->axis_range[i].lbase = NULL;
    }

    for (i=0; i<NTITLES; i++) {
	ed->gpt_titles[i].desc = (i == 0)? N_("Title of plot"):
	    N_("Title for axis");
	ed->gpt_titles[i].tab = i;
	ed->gpt_titles[i].widget = NULL;
    }

    old_lines_init(ed);
    old_labels_init(ed);
    old_arrows_init(ed);

    return ed;
} 

GtkWidget *plot_add_editor (png_plot *plot) 
{
    GPT_SPEC *spec = plot_get_spec(plot);
    GtkWidget *plotshell = plot_get_shell(plot);
    plot_editor *editor;
    GtkWidget *dialog, *vbox, *hbox;
    GtkWidget *button, *notebook;

    editor = plot_editor_new(spec);
    if (editor == NULL) {
	nomem();
	return NULL;
    }

    dialog = gretl_dialog_new(_("gretl plot controls"), plotshell, 
			      GRETL_DLG_RESIZE);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_set_spacing(GTK_BOX(vbox), 2);
#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(dialog), FALSE);
#endif
    
    if (plot != NULL) {
	gtk_window_set_transient_for(GTK_WINDOW(dialog), 
				     GTK_WINDOW(plotshell));
    }

    editor->dialog = dialog;

    g_signal_connect_swapped(G_OBJECT(dialog), "destroy",
			     G_CALLBACK(plot_editor_destroy), 
			     editor);
   
    editor->notebook = notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    gpt_tab_main(editor, spec);
    gpt_tab_XY(editor, spec, 0);
    gpt_tab_XY(editor, spec, 1);

    if (spec->flags & GPT_Y2AXIS) {
	gpt_tab_XY(editor, spec, 2);
    }

    if (spec->lines != NULL && spec->code != PLOT_STACKED_BAR) {
	gpt_tab_lines(editor, spec, 0);
    }

    gpt_tab_labels(editor, spec, 0);
    gpt_tab_arrows(editor, spec, 0);

    if (!frequency_plot_code(spec->code)) {
	gpt_tab_palette(notebook);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Apply" button */
    button = apply_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_gpt_changes), editor);
    gtk_widget_grab_default(button);

    /* "OK" button (apply and close) */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_gpt_changes), editor);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy), 
			     dialog);

    /* Close button (do not apply changes) */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy), 
			     dialog);

    /* Help button */
    context_help_button(hbox, GR_PLOT);

    set_keyspec_sensitivity(editor);

    gtk_widget_show_all(hbox);
    gtk_widget_show(dialog);

    return dialog;
}

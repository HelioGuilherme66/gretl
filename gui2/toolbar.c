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

/* toolbar.c: main-window toolbar, viewer window toolbars, etc. */

#include "gretl.h"
#include "console.h"
#include "session.h"
#include "datafiles.h"
#include "selector.h"
#include "textbuf.h"
#include "textutil.h"
#include "series_view.h"
#include "model_table.h"
#include "cmdstack.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "winstack.h"
#include "tabwin.h"
#include "fncall.h"
#include "fnsave.h"
#include "toolbar.h"

#include "uservar.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

/* for viewer window toolbars */
#include "../pixmaps/mini.tex.xpm"
#include "../pixmaps/mail_16.xpm"
#include "../pixmaps/mini.tsplot.xpm"
#include "../pixmaps/mini.boxplot.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.pin.xpm"
#include "../pixmaps/mini.alpha.xpm"
#include "../pixmaps/mini.en.xpm"
#include "../pixmaps/mini.split_h.xpm"
#include "../pixmaps/mini.split_v.xpm"
#include "../pixmaps/mini.join_h.xpm"
#include "../pixmaps/mini.join_v.xpm"
#include "../pixmaps/mini.winlist.xpm"
#include "../pixmaps/mini.bundle.xpm"

/* for main-window toolbar */
#include "../pixmaps/mini.calc.xpm"
#include "../pixmaps/mini.sh.xpm"
#include "../pixmaps/mini.session.xpm"
#include "../pixmaps/mini.plot.xpm"
#include "../pixmaps/mini.model.xpm"
#include "../pixmaps/mini.func.xpm"
#include "../pixmaps/mini.db.xpm"

/* for window-finder */
#include "../pixmaps/mini.gretl.xpm"
#include "../pixmaps/mini.table.xpm"
#include "../pixmaps/mini.page.xpm"
#include "../pixmaps/mini.tools.xpm"

/* for plot bar */
#include "../pixmaps/upsize.xpm"
#include "../pixmaps/downsize.xpm"

enum {
    SAVE_ITEM = 1,
    SAVE_AS_ITEM,
    EDIT_ITEM,
    PLOT_ITEM,
    EXEC_ITEM,
    COPY_ITEM,
    TEX_ITEM,
    ADD_DATA_ITEM,
    ADD_MATRIX_ITEM,
    MAIL_ITEM,
    HELP_ITEM,
    CMD_HELP_ITEM,
    GP_HELP_ITEM,
    X12A_HELP_ITEM,
    SORT_ITEM,
    SORT_BY_ITEM,
    FORMAT_ITEM,
    INDEX_ITEM,
    EDIT_SCRIPT_ITEM,
    STICKIFY_ITEM,
    ALPHA_ITEM,
    REFRESH_ITEM,
    OPEN_ITEM,
    SPLIT_H_ITEM,
    SPLIT_V_ITEM,
    EDITOR_ITEM,
    NOTES_ITEM,
    NEW_ITEM,
    CLOSE_ITEM,
    BUNDLE_ITEM,
    FIND_ITEM,
    COPY_SCRIPT_ITEM,
    BUILD_ITEM
} viewbar_flags;

struct stock_maker {
    char **xpm;
    const char *str;
};

void gretl_stock_icons_init (void)
{
    struct stock_maker stocks[] = {
	{ mini_tex_xpm, GRETL_STOCK_TEX },
	{ mail_16_xpm, GRETL_STOCK_MAIL },
	{ mini_tsplot_xpm, GRETL_STOCK_TS },
	{ mini_boxplot_xpm, GRETL_STOCK_BOX },
	{ mini_pdf_xpm, GRETL_STOCK_PDF },
	{ mini_manual_xpm, GRETL_STOCK_BOOK },
	{ mini_calc_xpm, GRETL_STOCK_CALC },
	{ mini_sh_xpm, GRETL_STOCK_CONSOLE },
	{ mini_session_xpm, GRETL_STOCK_ICONS },
	{ mini_plot_xpm, GRETL_STOCK_SCATTER },
	{ mini_model_xpm, GRETL_STOCK_MODEL },
	{ mini_func_xpm, GRETL_STOCK_FUNC },
	{ mini_pin_xpm, GRETL_STOCK_PIN },
	{ mini_alpha_xpm, GRETL_STOCK_ALPHA },
	{ mini_en_xpm, GRETL_STOCK_EN },
	{ mini_split_h_xpm, GRETL_STOCK_SPLIT_H },
	{ mini_split_v_xpm, GRETL_STOCK_SPLIT_V },
	{ mini_join_h_xpm, GRETL_STOCK_JOIN_H },
	{ mini_join_v_xpm, GRETL_STOCK_JOIN_V },
	{ mini_winlist_xpm, GRETL_STOCK_WINLIST },
	{ mini_bundle_xpm, GRETL_STOCK_BUNDLE },
	{ mini_db_xpm, GRETL_STOCK_DB},
	{ mini_gretl_xpm, GRETL_STOCK_GRETL},
	{ mini_table_xpm, GRETL_STOCK_TABLE},
	{ mini_page_xpm, GRETL_STOCK_PAGE},
	{ mini_tools_xpm, GRETL_STOCK_TOOLS},
	{ upsize_xpm, GRETL_STOCK_BIGGER},
	{ downsize_xpm, GRETL_STOCK_SMALLER}
    };
    static GtkIconFactory *gretl_factory;
    int n = G_N_ELEMENTS(stocks);

    if (gretl_factory == NULL) {
	GtkIconSet *iset;
	GdkPixbuf *pbuf;
	int i;

	gretl_factory = gtk_icon_factory_new();

	for (i=0; i<n; i++) {
	    pbuf = gdk_pixbuf_new_from_xpm_data((const char **) stocks[i].xpm);
	    iset = gtk_icon_set_new_from_pixbuf(pbuf);
	    g_object_unref(pbuf);
	    gtk_icon_factory_add(gretl_factory, stocks[i].str, iset);
	    gtk_icon_set_unref(iset);
	}

	gtk_icon_factory_add_default(gretl_factory);
    }
}

/* callbacks for viewer window toolbar */

static void copy_to_editor (GtkWidget *w, windata_t *vwin)
{
    gchar *buf = textview_get_text(vwin->text);

    if (vwin->role == VIEW_LOG) {
	/* allow for the possibility that the buffer is empty */
	gchar *s = buf;
	int n = 0;

	while (s != NULL && n < 3) {
	    s = strchr(s, '\n');
	    if (s != NULL) {
		s++;
		n++;
	    }
	}

	if (s != NULL) {
	    gchar *modbuf;

	    modbuf = g_strdup_printf("# logged commands\n%s", s);
	    do_new_script(EDIT_SCRIPT, modbuf);
	    g_free(modbuf);
	} else {
	    do_new_script(EDIT_SCRIPT, buf);
	}
    } else {
	do_new_script(EDIT_SCRIPT, buf);
    }
    
    g_free(buf);
}

static void save_as_callback (GtkWidget *w, windata_t *vwin)
{
    GtkWidget *vmain = vwin_toplevel(vwin);
    guint u = 0;

    if (g_object_get_data(G_OBJECT(vmain), "text_out")) {
	const char *opts[] = {
	    N_("Save to file"),
	    N_("Save to session as icon")
	};
	int resp;

	resp = radio_dialog(_("gretl: save text"), _("Save text"), 
			    opts, 2, 0, 0, vmain);
	if (resp < 0) {
	    return;
	} else if (resp == 1) {
	    save_output_as_text_icon(vwin);
	    return;
	} else {
	    u = SAVE_OUTPUT;
	}
    } else if (vwin->role == EDIT_SCRIPT) {
	u = SAVE_SCRIPT;
    } else if (vwin->role == EDIT_GP) {
	u = SAVE_GP_CMDS;
    } else if (vwin->role == EDIT_R) {
	u = SAVE_R_CMDS;
    } else if (vwin->role == EDIT_OX) {
	u = SAVE_OX_CMDS;
    } else if (vwin->role == EDIT_OCTAVE) {
	u = SAVE_OCTAVE_CMDS;
    } else if (vwin->role == EDIT_PYTHON) {
        u = SAVE_PYTHON_CMDS;
    } else if (vwin->role == EDIT_JULIA) {
	u = SAVE_JULIA_CODE;
    } else if (vwin->role == EDIT_STATA) {
	u = SAVE_STATA_CMDS;
    } else if (vwin->role == EDIT_SPEC) {
	u = SAVE_SPEC_FILE;
    } else if (vwin->role == VIEW_FILE) {
	u = SAVE_TEXT;
    } else {
	dummy_call();
	return;
    }

    file_save(vwin, u);
}

static void mail_script_callback (GtkWidget *w, windata_t *vwin)
{
    if (viewer_char_count(vwin) == 0) {
	infobox(_("Nothing to send"));
	return;
    }

    if (query_save_text(NULL, NULL, vwin)) {
	return;
    }
    
    send_file(vwin->fname);
}

/* callback for the "Open" icon in a script editing window,
   which enables the user to switch to a different script,
   or to open another tab if the editor is tab-enabled
*/

static void file_open_callback (GtkWidget *w, windata_t *vwin)
{
    int tabbed = window_is_tab(vwin);
    
    if (!tabbed) {
	/* Don't proceed unconditionally if there's unsaved
	   text in a single-script window; this doesn't
	   apply if we're just opening another tab.
	*/
	if (query_save_text(NULL, NULL, vwin)) {
	    return;
	}
    }

    file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);

    if (!tabbed && (vwin->flags & VWIN_CONTENT_CHANGED)) {
	mark_vwin_content_saved(vwin);
	vwin_set_filename(vwin, tryfile);
    }
}

static void open_pkg_sample (GtkWidget *w, windata_t *vwin)
{
    if (viewer_char_count(vwin) > 0) {
	int resp;

	resp = yes_no_dialog(NULL, _("Really replace content with\n"
				     "a selected file?"),
			     vwin->main);
	if (resp != GRETL_YES) {
	    return;
	}
    }

    file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);
}

static void toolbar_new_callback (GtkWidget *w, windata_t *vwin)
{
    do_new_script(vwin->role, NULL);
}

static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
#ifdef G_OS_WIN32
    /* gtksourceview printing is screwed on Windows */
    window_print(NULL, vwin);
#else
    if (textview_use_highlighting(vwin->role)) {
	int resp = yes_no_cancel_dialog(NULL,
					_("Print with syntax highlighting?"),
					vwin_toplevel(vwin));

	if (resp == GRETL_YES) {
	    sourceview_print(vwin);
	} else if (resp == GRETL_NO) {
	    window_print(NULL, vwin);
	}
    } else {
	window_print(NULL, vwin);
    }
#endif
}

static void window_help (GtkWidget *w, windata_t *vwin)
{
    show_gui_help(vwin->role);
}

static void multi_save_as_callback (GtkWidget *w, windata_t *vwin)
{
    copy_format_dialog(vwin, W_SAVE);
}

static void script_index (GtkWidget *w, windata_t *vwin)
{
    display_files(PS_FILES, NULL);
}

static void toolbar_refresh (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_SERIES) {
	series_view_refresh(w, vwin);
    }
}

static void set_matrix_name (GtkWidget *widget, dialog_t *dlg)
{
    char *vname = (char *) edit_dialog_get_data(dlg);
    GtkWidget *parent =  edit_dialog_get_window(dlg);
    const gchar *s = edit_dialog_get_text(dlg);

    if (s == NULL || gui_validate_varname(s, GRETL_TYPE_MATRIX, parent)) {
	edit_dialog_reset(dlg);
    } else {
	strcpy(vname, s);
	edit_dialog_close(dlg);
    }
}

static void add_matrix_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == XTAB) {
	char mname[VNAMELEN];
	gretl_matrix *m;
	int err, cancel = 0;

	blocking_edit_dialog(0, _("gretl: save matrix"), 
			     _("Enter a name"), NULL,
			     set_matrix_name, mname, 
			     VARCLICK_NONE, 
			     vwin_toplevel(vwin),
			     &cancel);
	if (!cancel) {
	    m = xtab_to_matrix(vwin->data);
	    if (m == NULL) {
		nomem();
	    } else {
		err = user_var_add_or_replace(mname,
					      GRETL_TYPE_MATRIX,
					      m);
		if (err) {
		    gretl_matrix_free(m);
		    gui_errmsg(err);
		} else {
		    infobox_printf(_("Saved matrix as %s"), mname);
		}
	    }
	}
    }
}

static void add_data_callback (GtkWidget *w, windata_t *vwin)
{
    int oldv = dataset->v;

    if (vwin->role == PCA) {
	add_pca_data(vwin);
    } else if (vwin->role == LEVERAGE) {
	add_leverage_data(vwin);
    } else if (vwin->role == MAHAL) {
	add_mahalanobis_data(vwin);
    } else if (vwin->role == FCAST) {
	add_fcast_data(vwin);
    } else if (vwin->role == LOESS || vwin->role == NADARWAT) {
	add_nonparam_data(vwin);
    }

    if (dataset->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }	
}

static void real_coeffint_set_alpha (GtkWidget *w, GtkWidget *dialog)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(dialog), "vwin");
    double *x = g_object_get_data(G_OBJECT(dialog), "xptr");
    CoeffIntervals *cf = vwin->data;
    GtkTextBuffer *buf;
    const char *newtext;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    reset_coeff_intervals(cf, 1.0 - *x);
    text_print_model_confints(cf, prn);
    newtext = gretl_print_get_buffer(prn);
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_set_text(buf, "", -1);
    textview_set_text(vwin->text, newtext);
    gretl_print_destroy(prn); 

    gtk_widget_destroy(dialog);
}

static void alpha_button_callback (GtkToggleButton *b, double *x)
{
    if (gtk_toggle_button_get_active(b)) {
	int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(b), "i"));

	if (i == 0) {
	    *x = 0.90;
	} else if (i == 1) {
	    *x = 0.95;
	} else if (i == 2) {
	    *x = 0.99;
	}
    } 
}

static void reformat_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_MODELTABLE) {
	format_model_table(vwin);
    } else {
	series_view_format_dialog(vwin);
    }
}

static void split_pane_callback (GtkWidget *w, windata_t *vwin)
{
    GtkWidget *hb = g_object_get_data(G_OBJECT(w), "hpane");
    GtkWidget *vb = g_object_get_data(G_OBJECT(w), "vpane");
    int vertical = 0;

    if (hb != NULL) {
	vb = w;
	vertical = 1;
    } else {
	hb = w;
    }

    /* Note: by "vertical" here we mean that the split runs vertically,
       dividing the pane into left- and right-hand sections; otherwise
       the split runs horizontally.
    */

    if (g_object_get_data(G_OBJECT(vwin->vbox), "sw") != NULL) {
	/* currently in single-view mode: so split */
	viewer_split_pane(vwin, vertical);
	gtk_widget_set_sensitive(vb, vertical);
	gtk_widget_set_sensitive(hb, !vertical);
	if (vertical) {
	    gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(vb), 
					 GRETL_STOCK_JOIN_V);
	} else {
	    gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(hb), 
					 GRETL_STOCK_JOIN_H);
	}
    } else {
	GtkWidget *paned;

	paned = g_object_get_data(G_OBJECT(vwin->vbox), "paned");

	if (paned != NULL) {
	    /* currently in split-view mode: so rejoin */
	    vertical = GTK_IS_HPANED(paned);
	    viewer_close_pane(vwin);
	    gtk_widget_set_sensitive(hb, TRUE);
	    gtk_widget_set_sensitive(vb, TRUE);
	    if (vertical) {
		gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(vb), 
					     GRETL_STOCK_SPLIT_V);
	    } else {
		gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(hb), 
					     GRETL_STOCK_SPLIT_H);
	    }
	}
    }
}

static void toggle_alpha_spin (GtkToggleButton *b, GtkWidget *w)
{
    gtk_widget_set_sensitive(w, gtk_toggle_button_get_active(b));
}

static void coeffint_set_alpha (GtkWidget *w, windata_t *vwin)
{
    CoeffIntervals *cf = vwin->data;
    GtkWidget *dialog, *tmp, *hbox;
    GtkWidget *vbox, *b, *hb2;
    GSList *group = NULL;
    GtkAdjustment *adj;
    gchar txt[16];
    double x = 1.0 - cf->alpha;
    gboolean defset = FALSE;
    int i;

    if (maybe_raise_dialog()) {
	return;
    }

    dialog = gretl_dialog_new(_("gretl: coefficient confidence intervals"), 
			      vwin_toplevel(vwin), GRETL_DLG_BLOCK);

    hbox = gtk_hbox_new(FALSE, 5);
    hb2 = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Confidence level"));
    gtk_box_pack_start(GTK_BOX(hb2), tmp, FALSE, FALSE, 0);

    tmp = gtk_label_new("1 - α :");
    gtk_box_pack_start(GTK_BOX(hb2), tmp, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), hb2, TRUE, TRUE, 10);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    vbox = gtk_vbox_new(FALSE, 5);

    /* radio button for 90%, 95%, 99% confidence */

    for (i=0; i<3; i++) {
	double a = (i == 0)? 0.90 : (i == 1)? 0.95 : 0.99;

	sprintf(txt, "%.2f", a);
	b = gtk_radio_button_new_with_label(group, txt);
	gtk_box_pack_start(GTK_BOX(vbox), b, FALSE, FALSE, 0);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
	if (a == x) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), TRUE);
	    defset = TRUE;
	} else {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), FALSE);
	}
	g_object_set_data(G_OBJECT(b), "i", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(alpha_button_callback), 
			 &x);
    }

    /* radio button for "other" confidence level, plus spinner */

    hb2 = gtk_hbox_new(FALSE, 0);
    b = gtk_radio_button_new_with_label(group, _("Other"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), !defset);
    gtk_box_pack_start(GTK_BOX(hb2), b, FALSE, FALSE, 0);
    adj = (GtkAdjustment *) gtk_adjustment_new(x, 0.60, 0.99, 0.01, 0, 0);
    tmp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.01, 2);
    g_signal_connect(G_OBJECT(tmp), "value-changed", 
		     G_CALLBACK(set_double_from_spinner), &x);
    gtk_widget_set_sensitive(tmp, !defset);
    g_signal_connect(G_OBJECT(b), "toggled", 
		     G_CALLBACK(toggle_alpha_spin),
		     tmp);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hb2), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hb2, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 10);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_delete_button(hbox, dialog);

    g_object_set_data(G_OBJECT(dialog), "vwin", vwin);
    g_object_set_data(G_OBJECT(dialog), "xptr", &x);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(real_coeffint_set_alpha), dialog);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);
}

static void stickiness_callback (GtkWidget *w, windata_t *vwin)
{
    output_policy_dialog(vwin, vwin, 1);
}

static void toolbar_plot_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_SERIES) {
	series_view_graph(w, vwin);
    } else if (vwin->role == VIEW_BUNDLE) {
	exec_bundle_plot_function(vwin->data, NULL);
    } else {
	do_nonparam_plot(vwin);
    }
}

static void editor_prefs_callback (GtkWidget *w, windata_t *vwin)
{
    preferences_dialog(TAB_EDITOR, NULL, vwin_toplevel(vwin));
}

static void build_pkg_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin_content_changed(vwin)) {
	int resp;
	
	resp = yes_no_cancel_dialog("gretl", _("Save changes?"),
				    vwin->main);
	if (resp == GRETL_CANCEL) {
	    return;
	}
	if (resp == GRETL_YES) {
	    vwin_save_callback(NULL, vwin);
	}
    }

    build_package_from_spec_file(vwin);
}

static int bundle_plot_ok (windata_t *vwin)
{
    gretl_bundle *b = vwin->data;
    gchar *pf = get_bundle_plot_function(b);
    int ret = 0;

    if (pf != NULL) {
	ret = 1;
	g_free(pf);
    }

    return ret;
}

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
    text_set_cursor(vwin->text, GDK_QUESTION_ARROW);
    set_window_help_active(vwin);
}

static int edit_script_popup_item (GretlToolItem *item)
{
    return !strcmp(item->icon, GTK_STOCK_COPY) ||
	!strcmp(item->icon, GTK_STOCK_PASTE) ||
	!strcmp(item->icon, GTK_STOCK_FIND) ||
	!strcmp(item->icon, GTK_STOCK_UNDO) ||
	!strcmp(item->icon, GTK_STOCK_FIND_AND_REPLACE);
}

static void set_plot_icon (GretlToolItem *item, int role)
{
    if (role == LOESS || role == NADARWAT || role == VIEW_BUNDLE) {
	item->icon = GRETL_STOCK_SCATTER;
    } else if (dataset_is_time_series(dataset)) {
	item->icon = GRETL_STOCK_TS;
    } else {
	item->icon = GRETL_STOCK_BOX;
    }
}

static void vwin_cut_callback (GtkWidget *w, windata_t *vwin)
{
    gtk_text_buffer_cut_clipboard(gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)),
				  gtk_clipboard_get(GDK_NONE),
				  TRUE);
}

static GretlToolItem viewbar_items[] = {
    { N_("New window"), GTK_STOCK_NEW, G_CALLBACK(toolbar_new_callback), NEW_ITEM },
    { N_("Open..."), GTK_STOCK_OPEN, G_CALLBACK(file_open_callback), OPEN_ITEM },
    { N_("Save"), GTK_STOCK_SAVE, G_CALLBACK(vwin_save_callback), SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, G_CALLBACK(save_as_callback), SAVE_AS_ITEM },
    { N_("Open in script editor"), GTK_STOCK_EDIT, G_CALLBACK(copy_to_editor), COPY_SCRIPT_ITEM },
    { N_("Save bundle content..."), GRETL_STOCK_BUNDLE, NULL, BUNDLE_ITEM },
    { N_("Print..."), GTK_STOCK_PRINT, G_CALLBACK(window_print_callback), 0 },
    { N_("Show/hide"), GRETL_STOCK_PIN, G_CALLBACK(session_notes_callback), NOTES_ITEM },
    { N_("Run"), GTK_STOCK_EXECUTE, G_CALLBACK(do_run_script), EXEC_ITEM },
    { N_("Build package"), GRETL_STOCK_TOOLS, G_CALLBACK(build_pkg_callback), BUILD_ITEM },
    { N_("Cut"), GTK_STOCK_CUT, G_CALLBACK(vwin_cut_callback), EDIT_ITEM }, 
    { N_("Copy"), GTK_STOCK_COPY, G_CALLBACK(vwin_copy_callback), COPY_ITEM }, 
    { N_("Paste"), GTK_STOCK_PASTE, G_CALLBACK(text_paste), EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, G_CALLBACK(text_find), FIND_ITEM },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, G_CALLBACK(text_replace), EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, G_CALLBACK(text_undo), EDIT_ITEM },
    { N_("Redo"), GTK_STOCK_REDO, G_CALLBACK(text_redo), EDIT_ITEM },
    { N_("Sort"), GTK_STOCK_SORT_ASCENDING, G_CALLBACK(series_view_toggle_sort), SORT_ITEM },    
    { N_("Sort by..."), GTK_STOCK_SORT_ASCENDING, G_CALLBACK(multi_series_view_sort_by), SORT_BY_ITEM },
    { N_("Preferences..."), GTK_STOCK_PREFERENCES, G_CALLBACK(editor_prefs_callback), EDIT_SCRIPT_ITEM },
    { N_("Send To..."), GRETL_STOCK_MAIL, G_CALLBACK(mail_script_callback), MAIL_ITEM },
    { N_("Scripts index"), GTK_STOCK_INDEX, G_CALLBACK(script_index), INDEX_ITEM },
    { N_("Confidence level..."), GRETL_STOCK_ALPHA, G_CALLBACK(coeffint_set_alpha), ALPHA_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, G_CALLBACK(window_tex_callback), TEX_ITEM },
    { N_("Graph"), GRETL_STOCK_TS, G_CALLBACK(toolbar_plot_callback), PLOT_ITEM },
    { N_("Reformat..."), GTK_STOCK_CONVERT, G_CALLBACK(reformat_callback), FORMAT_ITEM },
    { N_("Edit values..."), GTK_STOCK_EDIT, G_CALLBACK(series_view_edit), EDITOR_ITEM },
    { N_("Refresh"), GTK_STOCK_REFRESH, G_CALLBACK(toolbar_refresh), REFRESH_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, G_CALLBACK(add_data_callback), ADD_DATA_ITEM },
    { N_("Add as matrix..."), GTK_STOCK_ADD, G_CALLBACK(add_matrix_callback), ADD_MATRIX_ITEM },
    { N_("Stickiness..."), GRETL_STOCK_PIN, G_CALLBACK(stickiness_callback), STICKIFY_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT_H, G_CALLBACK(split_pane_callback), SPLIT_H_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT_V, G_CALLBACK(split_pane_callback), SPLIT_V_ITEM },
    { N_("Help on command"), GTK_STOCK_HELP, G_CALLBACK(activate_script_help), CMD_HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(window_help), HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(display_gnuplot_help), GP_HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(display_x12a_help), X12A_HELP_ITEM }
};

static int n_viewbar_items = G_N_ELEMENTS(viewbar_items);

#define exec_ok(r) (vwin_editing_script(r) || \
		    r == VIEW_SCRIPT || \
		    r == VIEW_PKG_SAMPLE || \
	            r == EDIT_PKG_SAMPLE)

#define open_ok(r) (vwin_editing_script(r))

#define new_ok(r) (vwin_editing_script(r))

#define edit_ok(r) (vwin_editing_script(r) || \
		    vwin_editing_buffer(r) || \
                    r == EDIT_PKG_CODE || \
		    r == EDIT_PKG_SAMPLE || \
		    r == EDIT_PKG_HELP || \
		    r == EDIT_PKG_GHLP)

#define save_as_ok(r) (r != EDIT_HEADER && \
	               r != EDIT_NOTES && \
	               r != EDIT_PKG_CODE && \
		       r != EDIT_PKG_SAMPLE && \
		       r != EDIT_PKG_HELP && \
		       r != EDIT_PKG_GHLP && \
		       r != CONSOLE && \
		       r != VIEW_BUNDLE)

#define help_ok(r) (r == LEVERAGE || \
		    r == COINT2 || \
		    r == HURST || \
		    r == RMPLOT || \
		    r == MAHAL)

#define cmd_help_ok(r) (r == EDIT_SCRIPT || \
	                r == EDIT_PKG_CODE || \
                        r == EDIT_PKG_SAMPLE || \
			r == VIEW_PKG_SAMPLE || \
			r == VIEW_PKG_CODE || \
			r == VIEW_SCRIPT || \
			r == VIEW_LOG)

/* for a non-editable script: can offer option to copy
   content into an editor window */
#define copy_script_ok(r) (r == VIEW_PKG_SAMPLE || \
			   r == VIEW_PKG_CODE || \
			   r == VIEW_SCRIPT || \
			   r == VIEW_LOG)

#define sort_ok(r) (r == VIEW_SERIES)

#define plot_ok(r) (r == VIEW_SERIES || \
		    r == LOESS || \
		    r == NADARWAT)

#define add_data_ok(r) (r == PCA || r == LEVERAGE || \
                        r == MAHAL || r == FCAST || \
			r == LOESS || r == NADARWAT)

#define split_h_ok(r) (r == SCRIPT_OUT || r == FNCALL_OUT || \
		       r == VIEW_LOG || r == VIEW_PKG_CODE || \
		       vwin_editing_script(r))

#define split_v_ok(r) (r == SCRIPT_OUT || r == FNCALL_OUT)

/* Screen out unwanted menu items depending on the context; also
   adjust the callbacks associated with some items based on
   context.
*/

static GCallback tool_item_get_callback (GretlToolItem *item, windata_t *vwin, 
					 int latex_ok, int sortby_ok,
					 int format_ok, int save_ok)
{
    GCallback func = item->func;
    int f = item->flag;
    int r = vwin->role;

    if (r == EDIT_SPEC) {
	/* This is a "special" that should maybe be regularized:
	   a bit like editing a script, but different...
	*/
	if (f == NEW_ITEM || f == OPEN_ITEM || f == EXEC_ITEM) {
	    return NULL;
	}
    } else if (f == BUILD_ITEM) {
	return NULL;
    }

    if (use_toolbar_search_box(r) && f == FIND_ITEM) {
	/* using an "inline" search box: skip the
	   "Find" button */
	return NULL;
    }

    if (copy_script_ok(r)) {
	if (f == SAVE_AS_ITEM) {
	    return NULL;
	}
    } else if (f == COPY_SCRIPT_ITEM) {
	return NULL;
    }

    if (r == EDIT_PKG_SAMPLE && f == OPEN_ITEM) {
	return G_CALLBACK(open_pkg_sample);
    } else if (!edit_ok(r) && f == EDIT_ITEM) {
	return NULL;
    } else if (!open_ok(r) && f == OPEN_ITEM) {
	return NULL;
    } else if (!new_ok(r) && f == NEW_ITEM) {
	return NULL;
    } else if (!exec_ok(r) && f == EXEC_ITEM) {
	return NULL;
    } else if (!cmd_help_ok(r) && f == CMD_HELP_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && f == MAIL_ITEM) {
	return NULL;
    } else if (!help_ok(r) && f == HELP_ITEM) {
	return NULL;
    } else if ((!latex_ok || !multiple_formats_ok(vwin)) && f == TEX_ITEM) {
	return NULL;
    } else if (!add_data_ok(r) && f == ADD_DATA_ITEM) {
	return NULL;
    } else if (r != XTAB && f == ADD_MATRIX_ITEM) {
	return NULL;
    } else if (!sort_ok(r) && f == SORT_ITEM) {
	return NULL;
    } else if (!sortby_ok && f == SORT_BY_ITEM) {
	return NULL;
    } else if (!plot_ok(r) && f == PLOT_ITEM) {
	if (r == VIEW_BUNDLE && bundle_plot_ok(vwin)) {
	    ; /* alright then */
	} else {
	    return NULL;
	}
    } else if (!split_h_ok(r) && f == SPLIT_H_ITEM) {
	return NULL;
    } else if (!split_v_ok(r) && f == SPLIT_V_ITEM) {
	return NULL;
    } else if (!format_ok && f == FORMAT_ITEM) {
	return NULL;
    } else if (r != VIEW_SERIES && f == EDITOR_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && r != EDIT_PKG_CODE && 
	       r != EDIT_PKG_SAMPLE && f == EDIT_SCRIPT_ITEM) {
	return NULL;
    } else if (r != VIEW_SCRIPT && f == INDEX_ITEM) {
	return NULL;
    } else if (r != SCRIPT_OUT && f == STICKIFY_ITEM) {
	return NULL;
    } else if (r != COEFFINT && f == ALPHA_ITEM) {
	return NULL;
    } else if (r != VIEW_SERIES && f == REFRESH_ITEM) {
	return NULL;
    } else if (r != EDIT_GP && f == GP_HELP_ITEM) {
	return NULL;
    } else if (r != EDIT_X12A && f == X12A_HELP_ITEM) {
	return NULL;
    } else if (f == SAVE_ITEM && !save_ok) {
	return NULL;
    } else if (r != EDIT_NOTES && f == NOTES_ITEM) {
	return NULL;
    } else if (r != VIEW_BUNDLE && f == BUNDLE_ITEM) {
	return NULL;
    } else if (f == SAVE_AS_ITEM) {
	if (!save_as_ok(r)) {
	    return NULL;
	} else if (multiple_formats_ok(vwin) || (r == PRINT && vwin->data != NULL)) {
	    func = G_CALLBACK(multi_save_as_callback);
	}
    }

    return func;
}

static GtkWidget *tool_item_get_menu (GretlToolItem *item, windata_t *vwin)
{
    GtkWidget *menu = NULL;

    if (vwin->role == VIEW_BUNDLE) {
	if (item->flag == BUNDLE_ITEM) {
	    menu = make_bundle_content_menu(vwin);
	} else if (item->flag == PLOT_ITEM) {
	    menu = make_bundle_plot_menu(vwin);
	} else if (item->flag == SAVE_ITEM) {
	    menu = make_bundle_save_menu(vwin);
	    if (menu != NULL) {
		item->tip = N_("Save...");
	    }
	}
    }

    if (menu != NULL) {
	/* don't leak: record pointer to menu so it can
	   be destroyed when the window is closed */
	vwin_record_toolbar_popup(vwin, menu);
    }

    return menu;
}

static void gretl_toolbar_flat (GtkWidget *w)
{
    static int style_done;

    gtk_widget_set_name(w, "gretl_toolbar");

    if (!style_done) {
	gtk_rc_parse_string("style \"gretl-tb-style\"\n{\n"
			    "  GtkToolbar::shadow-type = GTK_SHADOW_NONE\n"
			    "}\n"
			    "widget \"*.gretl_toolbar\" style \"gretl-tb-style\"");
	style_done = 1;
    }
}

GtkWidget *gretl_toolbar_new (GtkWidget *sibling)
{
    GtkWidget *tb = gtk_toolbar_new();

    gtk_toolbar_set_icon_size(GTK_TOOLBAR(tb), GTK_ICON_SIZE_MENU);
    gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(tb), FALSE);

    if (sibling == NULL) {
	/* if we're not alongside a menu bar ("sibling"),
	   show the toolbar without a shadow 
	*/
	gretl_toolbar_flat(tb);
    }       

    return tb;
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_widget_set_tooltip_text(w, str);
}

static GtkToolItem *gretl_menu_button (const char *icon,
				       const char *tip,
				       GtkWidget **pw)
{
    GtkWidget *img, *button = gtk_button_new();
    GtkToolItem *item = gtk_tool_item_new();

    gtk_widget_set_tooltip_text(GTK_WIDGET(item), _(tip));
    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    img = gtk_image_new_from_stock(icon, GTK_ICON_SIZE_MENU);
    gtk_container_add(GTK_CONTAINER(button), img);
    gtk_container_add(GTK_CONTAINER(item), button);
    *pw = button;

    return item;
}

static void gretl_tool_item_set_tip (GtkWidget *item,
				     GretlToolItem *tool)
{
    const char *accel = NULL;

    if (tool->flag == EXEC_ITEM) {
	accel = "Ctrl-R";
    } else if (tool->flag == COPY_ITEM) {
	accel = "Ctrl-C";
    } else if (tool->flag == SAVE_ITEM) {
	accel = "Ctrl-S";
    } else if (tool->flag == FIND_ITEM) {
	accel = "Ctrl-F";
    }

    if (accel != NULL) {
	gchar *s = g_strdup_printf("%s (%s)", _(tool->tip), accel);

	gtk_widget_set_tooltip_text(item, s);
	g_free(s);
    } else {
	gtk_widget_set_tooltip_text(item, _(tool->tip));
    }    
}

GtkWidget *gretl_toolbar_insert (GtkWidget *tbar,
				 GretlToolItem *tool,
				 GCallback func,
				 gpointer data,
				 gint pos)
{
    GtkToolItem *item;

    item = gtk_tool_button_new_from_stock(tool->icon);
    gretl_tool_item_set_tip(GTK_WIDGET(item), tool);
    g_signal_connect(G_OBJECT(item), "clicked", func, data);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, pos);

    return GTK_WIDGET(item);
}

static void button_menu_pos (GtkMenu *menu,
			     gint *x,
			     gint *y,
			     gboolean *push_in,
			     gpointer data)
{
    GtkWidget *button = data;
    gint wx, wy, tx, ty;

    gdk_window_get_origin(gtk_widget_get_window(button), &wx, &wy);
    gtk_widget_translate_coordinates(button, gtk_widget_get_toplevel(button), 
				     0, 0, &tx, &ty);
    *x = wx + tx;
    *y = wy + ty + 26;
    *push_in = TRUE;
}

static void tool_item_popup (GtkWidget *button, GdkEvent *event, 
			     GtkWidget *menu)
{
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, 
		   button_menu_pos, button,
		   event->button.button, event->button.time);
}

static GtkWidget *vwin_toolbar_insert (GretlToolItem *tool,
				       GCallback func,
				       GtkWidget *menu,
				       windata_t *vwin,
				       gint pos)
{
    GtkToolItem *item;

    if (menu != NULL) {
	/* make and insert a button that pops down a menu */
	GtkWidget *button;
	
	item = gretl_menu_button(tool->icon, tool->tip, &button);
	g_signal_connect(G_OBJECT(button), "button-press-event", 
			 G_CALLBACK(tool_item_popup), menu);
    } else {
	/* make and insert a regular callback button */
	item = gtk_tool_button_new_from_stock(tool->icon);
	g_signal_connect(G_OBJECT(item), "clicked", func, vwin);
	if (tool->flag == NEW_ITEM && window_is_tab(vwin)) {
	    gtk_widget_set_tooltip_text(GTK_WIDGET(item), _("New tab"));
	} else {
	    gretl_tool_item_set_tip(GTK_WIDGET(item), tool);
	}
    }

    gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), item, pos);

    return GTK_WIDGET(item);
}

static void viewbar_add_items (windata_t *vwin, ViewbarFlags flags)
{
    int sortby_ok = has_sortable_data(vwin);
    int format_ok = can_format_data(vwin);
    int latex_ok = latex_is_ok();
    int save_ok = (flags & VIEWBAR_EDITABLE);
    GtkWidget *hpane = NULL, *vpane = NULL;
    GtkWidget *button;
    GtkWidget *menu;
    GretlToolItem *item;
    GCallback func;
    int i;
 
    for (i=0; i<n_viewbar_items; i++) {
	func = NULL;
	menu = NULL;

	item = &viewbar_items[i];

	/* Is there anything to hook up, in context? We
	   try first for a menu to attach to the toolbar
	   button; failing that we test for a "direct"
	   callback function.
	*/
	menu = tool_item_get_menu(item, vwin);
	if (menu == NULL && item->func != NULL) {
	    func = tool_item_get_callback(item, vwin, latex_ok, sortby_ok,
					  format_ok, save_ok);
	}
	if (func == NULL && menu == NULL) {
	    continue;
	}

	if (item->flag == PLOT_ITEM) {
	    set_plot_icon(item, vwin->role);
	}

	button = vwin_toolbar_insert(item, func, menu, vwin, -1);

	if (func == (GCallback) split_pane_callback) {
	    if (hpane == NULL) {
		hpane = button;
	    } else {
		vpane = button;
	    }
	}

	if (item->flag == SAVE_ITEM) { 
	    if (vwin->role != CONSOLE && vwin->role != VIEW_BUNDLE) {
		/* nothing to save just yet */
		g_object_set_data(G_OBJECT(vwin->mbar), "save_button", button); 
		gtk_widget_set_sensitive(button, FALSE);
	    }
	} else if (item->flag == SAVE_AS_ITEM) {
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_as_button", button);
	    if (strstr(vwin->fname, "script_tmp")) {
		gtk_widget_set_sensitive(button, FALSE);
	    }
	}
    }

    if (hpane != NULL) {
	g_object_set_data(G_OBJECT(hpane), "vpane", vpane);
    }

    if (vpane != NULL) {
	g_object_set_data(G_OBJECT(vpane), "hpane", hpane);
    }
}

void vwin_add_viewbar (windata_t *vwin, ViewbarFlags flags)
{
    if ((flags & VIEWBAR_HAS_TEXT) || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->main), "text_out", 
			  GINT_TO_POINTER(1));
    }

    vwin->mbar = gretl_toolbar_new(NULL);
    viewbar_add_items(vwin, flags);
    vwin_pack_toolbar(vwin);
}

GtkWidget *build_text_popup (windata_t *vwin)
{
    GtkWidget *pmenu = gtk_menu_new();
    GretlToolItem *item;
    GCallback func;
    GtkWidget *w;
    int i;

    for (i=0; i<n_viewbar_items; i++) {
	item = &viewbar_items[i];
	if (item->flag == SPLIT_H_ITEM || item->flag == SPLIT_V_ITEM) {
	    continue;
	}
	if (vwin->role == EDIT_SCRIPT) {
	    /* the script editor popup may have some special stuff
	       added: don't clutter it up */
	    if (edit_script_popup_item(item)) {
		func = item->func;
	    } else {
		func = NULL;
	    }
	} else {
	    func = tool_item_get_callback(item, vwin, 0, 0, 0, 0);
	}
	if (func != G_CALLBACK(NULL)) {
	    if (func == G_CALLBACK(text_paste)) {
		GtkClipboard *cb = gtk_clipboard_get(GDK_NONE);

		if (!gtk_clipboard_wait_is_text_available(cb)) {
		    continue;
		}
	    } else if (func == G_CALLBACK(text_undo) && !text_can_undo(vwin)) {
		continue;
	    }
	    w = gtk_menu_item_new_with_label(_(item->tip));
	    g_signal_connect(G_OBJECT(w), "activate", func, vwin);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), w);
	}
    }

    if (vwin->role != EDIT_SCRIPT) {
	if (window_is_undockable(vwin)) {
	    add_undock_popup_item(pmenu, vwin);
	} else if (window_is_dockable(vwin)) {
	    add_dock_popup_item(pmenu, vwin);
	}
    }

    return pmenu;
}

/* callbacks for main-window toolbar icons */

static void tbar_calc (void)
{
#ifdef G_OS_WIN32
    create_child_process(calculator);
#else
    gretl_fork("calculator", NULL, NULL);
#endif 
}

static void tbar_open_data (void)
{
    display_files(TEXTBOOK_DATA, NULL);
}

static void tbar_command_ref (void)
{
    plain_text_cmdref(NULL);
}

static void tbar_xy_graph (void)
{
    if (data_status) {
	if (dataset->v == 2) {
	    do_graph_var(mdata->active_var);
	} else if (mdata_selection_count() == 2) {
	    plot_from_selection(GR_XY);
	} else {
	    selection_dialog(GR_XY, _("gretl: define graph"), 
			     do_graph_from_selector);
	}
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_model (void)
{
    if (data_status) {
	selection_dialog(OLS, _("gretl: specify model"), do_model);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_new_script (void)
{
    do_new_script(EDIT_SCRIPT, NULL);
}

static void tbar_show_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, mdata);
}

/* end toolbar icon callbacks */

static GretlToolItem mainbar_items[] = {
    { N_("launch calculator"),  GRETL_STOCK_CALC,    G_CALLBACK(tbar_calc), 0 },
    { N_("new script"),         GTK_STOCK_EDIT,      G_CALLBACK(tbar_new_script), 0 },
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, G_CALLBACK(gretl_console), 0 },
    { N_("session icon view"),  GRETL_STOCK_ICONS,   G_CALLBACK(view_session), 0 },
    { N_("function packages"),  GRETL_STOCK_FUNC,    G_CALLBACK(tbar_show_funcs), 0 },
    { N_("command reference"),  GTK_STOCK_HELP,      G_CALLBACK(tbar_command_ref), 0 },
    { N_("X-Y graph"),          GRETL_STOCK_SCATTER, G_CALLBACK(tbar_xy_graph), 0 },
    { N_("OLS model"),          GRETL_STOCK_MODEL,   G_CALLBACK(tbar_model), 0 },
    { N_("gretl database"),     GRETL_STOCK_DB,      G_CALLBACK(show_native_dbs), 0 },
    { N_("open dataset"),       GTK_STOCK_OPEN,      G_CALLBACK(tbar_open_data), 0 },
};

void add_mainwin_toolbar (GtkWidget *vbox)
{
    GretlToolItem *item;
    GtkWidget *hbox;
    int i, n = G_N_ELEMENTS(mainbar_items);

    mdata->mbar = gretl_toolbar_new(NULL);

    for (i=0; i<n; i++) {
	item = &mainbar_items[i];
	gretl_toolbar_insert(mdata->mbar, item, item->func, mdata, -1);
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), mdata->mbar, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

/* Add a temporary menubar for use in a script output
   window, while we're waiting for the output. If the
   output window is being reused this is a bit more
   complicated; we have to "hide" the regular menubar
   before inserting the temporary one.
 */

void vwin_add_tmpbar (windata_t *vwin)
{
    GretlToolItem item = {
	N_("Stop"), 
	GTK_STOCK_STOP, 
	G_CALLBACK(do_stop_script), 
	0
    };
    GtkWidget *hbox, *tmp;

    hbox = g_object_get_data(G_OBJECT(vwin->main), "top-hbox");

    if (hbox != NULL) {
	/* We're replacing a "real" menubar temporarily: ref. the
	   widgets in @hbox before removing them so we can put
	   them back later.
	*/
	GtkWidget *winlist = g_object_get_data(G_OBJECT(hbox), "winlist");

	g_object_ref(G_OBJECT(vwin->mbar));
	gtk_container_remove(GTK_CONTAINER(hbox), vwin->mbar);
	if (winlist != NULL) {
	    g_object_ref(G_OBJECT(winlist));
	    gtk_container_remove(GTK_CONTAINER(hbox), winlist);
	}
    } else {
	/* starting from scratch */
	hbox = gtk_hbox_new(FALSE, 0);
	g_object_set_data(G_OBJECT(vwin->main), "top-hbox", hbox);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
    }

    tmp = gretl_toolbar_new(NULL);
    gretl_toolbar_insert(tmp, &item, item.func, NULL, 0);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    start_wait_for_output(vwin, hbox);
    gtk_widget_show_all(hbox);
}

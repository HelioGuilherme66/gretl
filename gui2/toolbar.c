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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* toolbar.c: the gretl toolbar */

#include "gretl.h"
#include "console.h"
#include "session.h"

#ifdef G_OS_WIN32
# include "webget.h"
#endif

/* pixmaps for gretl toolbar */
#include "../pixmaps/mini.calc.xpm"
#include "../pixmaps/mini.edit.xpm"
#include "../pixmaps/mini.sh.xpm"
#include "../pixmaps/mini.session.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.netscape.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.plot.xpm"
#include "../pixmaps/mini.model.xpm"
#include "../pixmaps/mini.ofolder.xpm"

static GtkWidget *toolbar_box;

/* callbacks for gretl toolbar icons */

static void show_calc (void)
{
#ifdef G_OS_WIN32
    create_child_process(calculator, NULL);
#else
    gretl_fork(calculator, NULL);
#endif 
}

static void open_textbook_data (void)
{
    display_files(NULL, TEXTBOOK_DATA, NULL);
}

#ifndef G_OS_WIN32

static void netscape_open (const char *url)
{
# ifdef USE_GNOME
#  ifndef OLD_GTK
    gnome_url_show(url, NULL); 
#  else
    gnome_url_show(url); 
#  endif  
# else
    int err;
    char ns_cmd[128];

    sprintf(ns_cmd, "netscape -remote \"openURLNewWindow(%s)\"", url);
    err = gretl_spawn(ns_cmd);
    if (err) gretl_fork("netscape", url);
# endif /* USE_GNOME */
}

#endif /* ! G_OS_WIN32 */

static void gretl_website (void)
{
#ifdef G_OS_WIN32
    if (goto_url("http://gretl.sourceforge.net/")) {
	errbox("Failed to open URL");
    }
#else
    netscape_open("http://gretl.sourceforge.net/");
#endif
}

static void gretl_pdf (void)
{
    char manurl[64];

    sprintf(manurl, "http://gretl.sourceforge.net/%s", _("manual.pdf"));

#ifdef G_OS_WIN32
    if (goto_url(manurl)) {
	errbox(_("Failed to open URL"));
    }
#else
    netscape_open(manurl);
#endif
}

static void xy_graph (void)
{
    if (data_status) {
	if (mdata_selection_count() == 2) {
	    plot_from_selection(NULL, GR_XY, NULL);
	} else {
	    selector_callback(NULL, GR_XY, NULL);
	}
    } else {
	errbox(_("Please open a data file first"));
    }
}

static void ols_model (void)
{
    if (data_status) {
	model_callback(NULL, OLS, NULL);
    } else {
	errbox(_("Please open a data file first"));
    }
}

static void go_session (void)
{
    if (data_status) {
	view_session();
    } else {
	errbox(_("Please open a data file first"));
    }
}

static void new_script_callback (void)
{
    do_new_script(NULL, 0, NULL);
}

/* end toolbar icon callbacks */

#ifndef OLD_GTK

static GtkWidget *image_button_new (GdkPixbuf *pix, void (*toolfunc)())
{
    GtkWidget *image = gtk_image_new_from_pixbuf(pix);
    GtkWidget *button = gtk_button_new();

    gtk_widget_set_size_request(button, 26, 24); /* 26, 24 */

    gtk_container_add (GTK_CONTAINER(button), image);
    g_signal_connect (G_OBJECT(button), "clicked",
                      G_CALLBACK(toolfunc), NULL);

    return button;
}

static void gretl_toolbar_append_tool (GtkWidget *tbar,
				       GtkWidget *w,
				       const char *toolstr)
{
    gtk_box_pack_start(GTK_BOX(tbar), w, FALSE, FALSE, 0);
    gretl_tooltips_add(w, toolstr);
}

#endif /* !OLD_GTK */

static void make_toolbar (GtkWidget *w, GtkWidget *box)
{
    GtkWidget *button;
    GtkWidget *toolbar;
#ifdef OLD_GTK
    GtkWidget *iconw;
    GdkPixmap *icon;
    GdkBitmap *mask;
    GdkColormap *cmap;
#else
    GtkWidget *hbox;
    GdkPixbuf *icon;
#endif
    int i;
    const char *toolstrings[] = {
	N_("launch calculator"), 
	N_("new script"), 
	N_("open gretl console"),
	N_("session icon view"),
	N_("gretl website"), 
	N_("gretl manual (PDF)"),
	N_("show help"), 
	N_("X-Y graph"), 
	N_("OLS model"),
	N_("open dataset"),
	NULL
    };
    gchar **toolxpm = NULL;
    void (*toolfunc)() = NULL;
    const char *toolstr;

#ifdef OLD_GTK
    cmap = gdk_colormap_get_system();
    toolbar_box = gtk_handle_box_new();
    gtk_handle_box_set_shadow_type(GTK_HANDLE_BOX(toolbar_box), 
				   GTK_SHADOW_NONE);
    gtk_box_pack_start(GTK_BOX(box), toolbar_box, FALSE, FALSE, 0);

    toolbar = gtk_toolbar_new(GTK_ORIENTATION_HORIZONTAL,
			      GTK_TOOLBAR_ICONS);
    gtk_container_set_border_width(GTK_CONTAINER(toolbar), 0);
    gtk_toolbar_set_space_size(GTK_TOOLBAR(toolbar), 0);
    gtk_container_add(GTK_CONTAINER(toolbar_box), toolbar);

    colorize_tooltips(GTK_TOOLBAR(toolbar)->tooltips);
#else
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);

    toolbar_box = gtk_handle_box_new();
    gtk_box_pack_start(GTK_BOX(hbox), toolbar_box, FALSE, FALSE, 0);

    toolbar = gtk_hbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(toolbar_box), toolbar);
#endif

    for (i=0; toolstrings[i] != NULL; i++) {
	switch (i) {
	case 0:
	    toolxpm = mini_calc_xpm;
	    toolfunc = show_calc;
	    break;
	case 1:
	    toolxpm = mini_edit_xpm;
	    toolfunc = new_script_callback;
	    break;
	case 2:
	    toolxpm = mini_sh_xpm;
	    toolfunc = show_gretl_console;
	    break;
	case 3:
	    toolxpm = mini_session_xpm;
	    toolfunc = go_session;
	    break;
	case 4:
	    toolxpm = mini_netscape_xpm;
	    toolfunc = gretl_website;
	    break;  
	case 5:
	    toolxpm = mini_pdf_xpm;
	    toolfunc = gretl_pdf;
	    break;    
	case 6:
	    toolxpm = mini_manual_xpm;
	    toolfunc = do_gui_help;
	    break;
	case 7:
	    toolxpm = mini_plot_xpm;
	    toolfunc = xy_graph;
	    break;
	case 8:
	    toolxpm = mini_model_xpm;
	    toolfunc = ols_model;
	    break;
	case 9:
	    toolxpm = mini_ofolder_xpm;
	    toolfunc = open_textbook_data;
	    break;
	default:
	    break;
	}

	toolstr = _(toolstrings[i]);
#ifdef OLD_GTK
	icon = gdk_pixmap_colormap_create_from_xpm_d(NULL, cmap, &mask, 
						     NULL, toolxpm);
	iconw = gtk_pixmap_new(icon, mask);
	button = gtk_toolbar_append_item(GTK_TOOLBAR(toolbar),
					 NULL, toolstr, NULL,
					 iconw, toolfunc, NULL);
#else
	icon = gdk_pixbuf_new_from_xpm_data((const char **) toolxpm);
	button = image_button_new(icon, toolfunc);
	gretl_toolbar_append_tool(toolbar, button, toolstr);
	gdk_pixbuf_unref(icon);
#endif
    }

#ifdef OLD_GTK
    gtk_widget_show(toolbar);
    gtk_widget_show(toolbar_box);
#else
    gtk_widget_show_all(hbox);
#endif
}

/* public interface */

void show_or_hide_toolbar (int want_toolbar)
{
    if (want_toolbar && toolbar_box == NULL) {
	GtkWidget *vbox = g_object_get_data(G_OBJECT(mdata->w), "vbox");

	make_toolbar(mdata->w, vbox);
    } else if (!want_toolbar && toolbar_box != NULL) {
	gtk_widget_destroy(toolbar_box);
	toolbar_box = NULL;
    }
}


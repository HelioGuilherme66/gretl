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

/* treeutils.c for gretl */

#include "gretl.h"
#include "treeutils.h"
#include "datafiles.h"

/* these live in dialogs.c */
extern GtkWidget *active_edit_id; 
extern GtkWidget *active_edit_name;
extern GtkWidget *active_edit_text;

static gint list_alpha_compare (GtkTreeModel *model, 
				GtkTreeIter *a, GtkTreeIter *b,
				gpointer p)
{
    gchar *vname_a, *vname_b;
    gint ret;

    gtk_tree_model_get(model, a, 1, &vname_a, -1);
    gtk_tree_model_get(model, b, 1, &vname_b, -1);

    if (strcmp(vname_a, "const") == 0) {
	ret = 0;
    } else if (strcmp(vname_b, "const") == 0) {
	ret = 1;
    } else {
	ret = strcmp(vname_a, vname_b);
    }

    g_free(vname_a);
    g_free(vname_b);
    
    return ret;
}

static gint list_id_compare (GtkTreeModel *model, 
			     GtkTreeIter *a, GtkTreeIter *b,
			     gpointer p)
{
    gchar *vnum_a, *vnum_b;
    gint ret;

    gtk_tree_model_get(model, a, 0, &vnum_a, -1);
    gtk_tree_model_get(model, b, 0, &vnum_b, -1);

    ret = strcmp(vnum_a, vnum_b);

    g_free(vnum_a);
    g_free(vnum_b);
    
    return ret;
}

static gboolean no_select_row_zero (GtkTreeSelection *selection,
				    GtkTreeModel *model,
				    GtkTreePath *path,
				    gboolean path_currently_selected,
				    gpointer data)
{
    if (tree_path_get_row_number(path) == 0) {
	return FALSE;
    }

    return TRUE;
}

static void my_gtk_entry_append_text (GtkEntry *entry, gchar *add)
{
    const gchar *old = gtk_entry_get_text(entry);
    gchar *new = g_strdup_printf("%s%s", old, add);
    
    gtk_entry_set_text(entry, new);
    gtk_editable_set_position(GTK_EDITABLE(entry), -1);
    g_free(new);
}

static void update_dialogs_from_varclick (int active_var)
{
    const gchar *edttext;

    if (active_edit_id != NULL) {
	gchar addvar[9];

	edttext = gtk_entry_get_text(GTK_ENTRY(active_edit_id));
	if (*edttext != '\0') {
	    sprintf(addvar, " %d", active_var);
	} else {
	    sprintf(addvar, "%d", active_var);
	}
	my_gtk_entry_append_text(GTK_ENTRY(active_edit_id), addvar);
    } else if (active_edit_name != NULL) {
	edttext = gtk_entry_get_text(GTK_ENTRY(active_edit_name));
	my_gtk_entry_append_text(GTK_ENTRY(active_edit_name), 
				 datainfo->varname[active_var]);
	my_gtk_entry_append_text(GTK_ENTRY(active_edit_name), " ");
    } else if (active_edit_text != NULL) {
	GtkTextBuffer *tbuf;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(active_edit_text));
	if (tbuf != NULL) {
	    gtk_text_buffer_insert_at_cursor(tbuf, 
					     datainfo->varname[active_var], 
					     -1);
	}
    }  
}

gboolean main_varclick (GtkWidget *widget, GdkEventButton *event,
			windata_t *win)
{
    GtkTreeView *view = GTK_TREE_VIEW(win->listbox);
    GtkTreePath *path;
    gint row = 0;

    if (gtk_tree_view_get_path_at_pos(view, event->x, event->y, &path, 
				      NULL, NULL, NULL)) {
	row = tree_path_get_row_number(path);

	if (row != 0) {
	    gchar *varnum;

	    g_object_set_data(G_OBJECT(win->listbox), "active_row",
			      GINT_TO_POINTER(row));
	    tree_view_get_string(view, row, 0, &varnum);
	    win->active_var = atoi(varnum);
	    g_free(varnum);
	    update_dialogs_from_varclick(win->active_var);
	}
	gtk_tree_path_free(path);
    } else {
	/* clicked below the lines representing variables */
	return FALSE;
    }

    return (row == 0);
}

static void
bool_col_toggled (GtkCellRendererToggle *cell, gchar *path_str, windata_t *vwin)
{
    GtkTreeView *treeview = GTK_TREE_VIEW(vwin->listbox);
    GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreePath *path = gtk_tree_path_new_from_string(path_str);
    GtkTreeIter iter;
    gboolean val;
    gint col;

    col = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(model), "boolcol"));
    gtk_tree_model_get_iter(model, &iter, path);
    gtk_tree_model_get(model, &iter, col, &val, -1);

    if (val) {
	return;
    }

    if (vwin->role == FUNC_FILES) {
	browser_load_func(NULL, vwin);
    }

    gtk_list_store_set(GTK_LIST_STORE(model), &iter, col, TRUE, -1);
    gtk_tree_path_free(path);
}

void vwin_add_list_box (windata_t *vwin, GtkBox *box, 
			int ncols, gboolean hidden_col,
			GType *types, const char **titles) 
{
    GtkListStore *store; 
    GtkWidget *view, *scroller;
    GtkCellRenderer *renderer;
    GtkCellRenderer *bool_renderer = NULL;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    int i, viscols = ncols;

    if (hidden_col) viscols--;

    store = gtk_list_store_newv(ncols, types);

    vwin->listbox = view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 0, NULL);

    for (i=0; i<viscols; i++) {
	if (types[i] == G_TYPE_BOOLEAN) {
	    bool_renderer = gtk_cell_renderer_toggle_new();
	    g_object_set_data(G_OBJECT(store), "boolcol", GINT_TO_POINTER(i));
	    g_signal_connect(bool_renderer, "toggled",
			     G_CALLBACK(bool_col_toggled), vwin);
	    column = gtk_tree_view_column_new_with_attributes(_(titles[i]),
							      bool_renderer,
							      "active", i,
							      NULL);
#if 1
	    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column),
					    GTK_TREE_VIEW_COLUMN_FIXED);
	    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
#endif
	    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	} else {
	    column = gtk_tree_view_column_new_with_attributes(_(titles[i]),
							      renderer,
							      "text", i, 
							      NULL);
	    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	    if (vwin != mdata) {
		g_object_set(G_OBJECT(column), "resizable", TRUE, NULL);
	    }
	}	
    }

    if (hidden_col) {
	column = gtk_tree_view_column_new_with_attributes(NULL,
							  renderer,
							  "text", i, 
							  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
	gtk_tree_view_column_set_visible(column, FALSE);
    }

    /* set the selection properties on the tree view */
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));

    if (vwin == mdata) { 
	/* gretl main window */
	gtk_tree_selection_set_mode(select, GTK_SELECTION_EXTENDED);
	gtk_tree_selection_set_select_function(select, 
					       (GtkTreeSelectionFunc)
					       no_select_row_zero,
					       NULL, NULL);

	gtk_widget_set_events(view, GDK_POINTER_MOTION_MASK 
			      | GDK_POINTER_MOTION_HINT_MASK);

        g_signal_connect(G_OBJECT(view), "motion_notify_event",
			 G_CALLBACK(listbox_drag), NULL);

    } else {
	gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
	g_signal_connect(G_OBJECT(select), "changed",
			 G_CALLBACK(listbox_select_row),
			 vwin);
    }

    g_signal_connect(G_OBJECT(view), "key_press_event",
		     G_CALLBACK(catch_listbox_key),
		     vwin);

    g_signal_connect(G_OBJECT(view), "button_press_event",
		     G_CALLBACK(listbox_double_click),
		     vwin);

    /* set sort properties on the tree model */
    gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), 0, 
				    (GtkTreeIterCompareFunc) 
				    list_id_compare,
				    NULL, NULL);

    gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), 1, 
				    (GtkTreeIterCompareFunc) 
				    list_alpha_compare,
				    NULL, NULL);

    scroller = gtk_scrolled_window_new(NULL, NULL);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroller),
					GTK_SHADOW_IN);    

    gtk_container_add(GTK_CONTAINER(scroller), view);

    gtk_box_pack_start(box, scroller, TRUE, TRUE, TRUE);

    gtk_widget_show(view);
    gtk_widget_show(scroller);
}

void tree_view_get_string (GtkTreeView *view, int row, int col, gchar **val)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *path;

    model = gtk_tree_view_get_model(view);
    path = g_strdup_printf("%d", row);
    gtk_tree_model_get_iter_from_string(model, &iter, path);
    gtk_tree_model_get(model, &iter, col, val, -1);
    g_free(path);
}

int tree_path_get_row_number (GtkTreePath *path)
{
    return gtk_tree_path_get_indices(path)[0];
}



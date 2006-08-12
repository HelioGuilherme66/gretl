#ifndef TREEUTILS_H
#define TREEUTILS_H

void vwin_add_list_box (windata_t *vwin, GtkBox *box, 
			int ncols, gboolean hidden_col,
			GType *types, const char **titles,
			int tree);

void tree_view_get_string (GtkTreeView *view, int row, int col, gchar **val);

void tree_view_get_int (GtkTreeView *view, int row, int col, int *val);

int tree_path_get_row_number (GtkTreePath *path); 

gboolean main_varclick (GtkWidget *widget, GdkEventButton *event,
			windata_t *win);

int *main_window_selection_as_list (void); 

#endif /* TREEUTILS_H */

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

#ifndef FNCALL_H
#define FNCALL_H

int open_function_package (const char *pkgname,
			   const char *fname,
			   windata_t *vwin);

void function_call_cleanup (void);

void destroy_gui_package_info (void);

void gui_define_list (void);

void fncall_register_genr (int addv, gpointer p);

gchar *get_bundle_plot_function (gretl_bundle *b);

int exec_bundle_plot_function (gretl_bundle *b, const char *funname);

int try_exec_bundle_print_function (gretl_bundle *b, PRN *prn);

void get_fncall_param_info (GtkWidget *dlg, int *series_ok,
			    char **pname);

void maybe_add_packages_to_model_menus (windata_t *vwin);

void maybe_add_packages_to_menus (windata_t *vwin);

int package_is_available_for_menu (const gchar *pkgname,
				   const char *fname);

int gui_function_pkg_query_register (const char *fname, 
				     GtkWidget *parent);

void gui_function_pkg_unregister (const gchar *pkgname);

int gui_function_pkg_revise_status (const gchar *pkgname,
				    const gchar *fname,
				    const gchar *label,
				    const gchar *mpath,
				    gboolean uses_subdir);

int n_registered_packages (void);

int n_user_handled_packages (void);

void get_registered_pkg_info (int i, char **name, char **path,
			      char **label, int *modelwin);

int query_addons (void);

int download_addon (const char *pkgname, char **local_path);

char *installed_addon_status_string (const char *path,
				     const char *svstr);

#endif /* FNCALL_H */

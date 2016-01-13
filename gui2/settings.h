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

#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef G_OS_WIN32
int read_win32_config (int debug);
#else
int gretl_config_init (void);
#endif

#ifdef HAVE_TRAMO
int get_tramo_ok (void);
#endif

#ifdef HAVE_X12A
int get_x12a_ok (void);
#endif

#if defined(MAC_NATIVE) && defined(PKGBUILD)
void set_up_mac_look (void);
#endif

#ifdef G_OS_WIN32
int using_wimp (void);
void set_wimp_preferred (int s);
void set_up_windows_look (void);
#endif

void set_gretl_startdir (void);

int using_hc_by_default (void);

int get_manpref (void);

int autoicon_on (void);

int use_tabbed_editor (void);

int use_tabbed_model_viewer (void);

int session_prompt_on (void);

void set_session_prompt (int val);

int get_keep_folder (void);

void set_script_output_policy (int p, windata_t *vwin);

int get_script_output_policy (void);

int get_thread_warn (void);

void set_thread_warn (int s);

int write_rc (void);

void dump_rc (void);

void force_english_help (void);

int preferences_dialog (int page, const char *varname, GtkWidget *parent);

void font_selector (GtkAction *action);

void set_fixed_font (const char *fontname);

void update_persistent_graph_colors (void);

void update_persistent_graph_font (void);

void set_app_font (const char *fontname);

const char *get_app_fontname (void);

const char *get_fixed_fontname (void);

void get_default_dir_for_action (char *s, int action);

void working_dir_dialog (void);

int gui_set_working_dir (char *dirname);

void set_working_dir_callback (GtkWidget *w, char *path);

void set_path_callback (char *setvar, char *setting);

void set_datapage (const char *str);

void set_scriptpage (const char *str);

const char *get_datapage (void);

const char *get_scriptpage (void);

const char *get_default_hc_string (int ci);

int check_for_prog (const char *prog);

void get_model_table_prefs (int *colheads,
			    int *use_tstats,
			    int *do_pvals,
			    int *do_asts,
			    int *figs,
			    char *fmt);

void set_model_table_prefs (int colheads,
			    int use_tstats,
			    int do_pvals,
			    int do_asts,
			    int figs,
			    char fmt);

void set_author_mail (const char *s);

const char *get_author_mail (void);

const char *get_sourceview_style (void);

#endif /* SETTINGS_H */

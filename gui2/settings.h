#ifndef SETTINGS_H
#define SETTINGS_H

#define COLOR_MAX 3

#ifdef G_OS_WIN32
void read_rc (void);
#endif

int using_olddat (void);

int using_hc_by_default (void);

int get_manpref (void);

void set_rcfile (void);

void write_rc (void);

void dump_rc (void);

void force_english_help (void);

int options_dialog (int page);

void options_dialog_callback (gpointer p, guint u, GtkWidget *w);

void font_selector (gpointer data, guint which, GtkWidget *widget);

void set_fixed_font (void);

#ifndef USE_GNOME

void set_app_font (const char *fontname);

const char *get_app_fontname (void);

#endif

void gnuplot_color_selector (GtkWidget *w, gpointer p);

GtkWidget *color_patch_button (int colnum);

void get_default_dir (char *s, int action);

void filesel_set_path_callback (const char *setting, char *strvar);

void set_datapage (const char *str);

void set_scriptpage (const char *str);

const char *get_datapage (void);

const char *get_scriptpage (void);

#ifdef HAVE_TRAMO
void set_tramo_ok (int set);
#endif

#ifdef HAVE_X12A
void set_x12a_ok (int set);
#endif

#ifndef G_OS_WIN32
void first_time_set_user_dir (void);
#endif

int check_for_prog (const char *prog);

#endif /* SETTINGS_H */

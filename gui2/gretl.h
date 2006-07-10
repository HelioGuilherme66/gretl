/*
 *   Copyright (c) by Allin Cottrell
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

/* gretl.h, main header file for gretl gui */

#ifndef GRETL_H
#define GRETL_H

#ifdef _WIN32
# include "winconfig.h"
#else
# include "config.h"
# ifdef USE_GNOME
#   include <gnome.h>
# endif
#endif

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include "libgretl.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#if GTK_MAJOR_VERSION < 2
# define OLD_GTK
# include <gtkextra/gtkiconfilesel.h>
#else
# define GNULL (gconstpointer) NULL
#endif

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
#endif

#include "gretltypes.h"
#include "gui_utils.h"
#include "callbacks.h"
#include "dialogs.h"
#include "library.h"
#include "settings.h"
#include "helpfiles.h"

#if GTK_MAJOR_VERSION >= 2
# include "gui_recode.h"
#endif

#ifdef ENABLE_NLS
# include "locale.h"
#endif

#define GRETL_BUFSIZE 8192
#define MAXSTR FILENAME_MAX

/* basic global program vars */
extern double **Z;
extern DATAINFO *datainfo;
extern PATHS paths; 
extern MODEL **models;

/* global counters */
extern int plot_count;
extern int orig_vars;

/* global state variables */
extern int data_status;
extern char *storelist;
extern float gui_scale;

/* global filenames */
extern char cmdfile[MAXLEN];
extern char scriptfile[MAXLEN];
extern char sessionfile[MAXLEN];
extern char tryfile[MAXLEN];

/* global option-related vars */
extern int expert;
extern int updater;
extern int winsize;
extern int main_x;
extern int main_y;
extern int mainwin_width;
extern int mainwin_height;

#if defined(G_OS_WIN32)
extern int wimp;
#else
# ifndef OSX_BUILD
extern char viewps[MAXSTR];
extern char viewpdf[MAXSTR];
extern char Browser[MAXSTR];
# endif
#endif

extern char calculator[MAXSTR];
extern char latex[MAXSTR];
extern char viewdvi[MAXSTR];

#if defined(HAVE_AUDIO) && !defined(G_OS_WIN32)
extern char midiplayer[MAXSTR];
#endif

/* global GUI equipment */
extern windata_t *mdata;
extern GtkTargetEntry gretl_drag_targets[];
#ifdef OLD_GTK
extern GdkFont *fixed_font;
#else
extern PangoFontDescription *fixed_font;
#endif

#include "gretl_enums.h"

/* functions follow */

#ifndef WIN32
int gretl_fork (const char *prog, const char *arg);
void set_wm_icon (GtkWidget *w, gpointer data);
#endif

int mdata_selection_count (void);
int mdata_active_var (void);
void populate_varlist (void);
void clear_varlist (GtkWidget *widget);
int gui_restore_sample (void);

/* functions defined in files other than gretl.c */
void about_dialog (gpointer data);
void stats_calculator (gpointer p, guint code, GtkWidget *w);

#endif /* GRETL_H */

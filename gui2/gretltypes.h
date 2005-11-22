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

#ifndef GRETLTYPES_H
#define GRETLTYPES_H

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
#endif

#define GRETL_STOCK_TEX    "gretl-tex"
#define GRETL_STOCK_MAIL   "gretl-mail"
#define GRETL_STOCK_XY     "gretl-plot"
#define GRETL_STOCK_TS     "gretl-tsplot"
#define GRETL_STOCK_BOX    "gretl-boxplot"

enum windata_flags {
    VWIN_HELP_ACTIVE = 1 << 0,
    VWIN_BUSY        = 1 << 1
};

typedef struct _windata_t windata_t;

struct _windata_t {
    GtkWidget *dialog;
    GtkWidget *vbox;
    GtkWidget *listbox; 
    GtkWidget *mbar;
    GtkWidget *w;
    GtkWidget *status;
    GtkWidget *popup;
    GtkItemFactory *ifac; 
    windata_t *gretl_parent;
    windata_t **gretl_children;
    gpointer data;
    int active_var; 
    int role;
    int n_model_tests;
    int n_gretl_children;
    unsigned char flags;
    char fname[MAXLEN];
#ifdef USE_GTKSOURCEVIEW
    GtkSourceBuffer *sbuf;
#endif
};

#define window_is_busy(w)    (w->flags & VWIN_BUSY)
#define set_window_busy(w)   (w->flags |= VWIN_BUSY)
#define unset_window_busy(w) (w->flags &= ~VWIN_BUSY)

#define window_help_is_active(w)    (w->flags & VWIN_HELP_ACTIVE)
#define set_window_help_active(w)   (w->flags |= VWIN_HELP_ACTIVE)
#define unset_window_help_active(w) (w->flags &= ~VWIN_HELP_ACTIVE)

#endif /* GRETLTYPES_H */

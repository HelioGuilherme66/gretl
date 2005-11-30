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

#ifndef HELPFILES_H
#define HELPFILES_H

void helpfile_init (void);

void context_help (GtkWidget *widget, gpointer data);

void plain_text_cmdref (gpointer p, guint cmdnum, GtkWidget *w);

gint edit_script_help (GtkWidget *widget, GdkEventButton *b,
		       windata_t *vwin);

void display_pdf_help (gpointer p, guint uguide, GtkWidget *w);

void datafile_find (GtkWidget *widget, gpointer data);

void menu_find (gpointer data, guint dbfind, GtkWidget *widget);

void text_find_callback (GtkWidget *w, gpointer data);

char *quoted_help_string (const char *s);

void gretl_tooltips_init (void);

void gretl_tooltips_add (GtkWidget *w, const gchar *str);

#ifdef OLD_GTK
void colorize_tooltips (GtkTooltips *tip);
#endif

#endif /* HELPFILES_H */

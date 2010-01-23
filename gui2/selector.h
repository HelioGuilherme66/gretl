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

#ifndef SELECTOR_H
#define SELECTOR_H

typedef struct iterinfo_t iterinfo;

struct iterinfo_t {
    int ci;
    int maxiters;
    double tol;
};

void clear_selector (void);

selector *selection_dialog (const char *title, int (*callback)(), 
			    guint cmdcode);

selector *simple_selection (const char *title, int (*callback)(), 
			    guint cmdcode, gpointer p);

void modelspec_dialog (int ci);

void selector_set_varnum (int v);

void selector_from_model (void *ptr, int ci);

char *main_window_selection_as_string (void);

void data_save_selection_wrapper (int file_code);

void add_remove_functions_dialog (char **pubnames, int npub,
				  char **privnames, int npriv,
				  void *p1, void *p2);

int selector_code (const selector *sr);

const char *selector_list (const selector *sr);

int selector_list_hasconst (const selector *sr);

gpointer selector_get_data (const selector *sr);

gretlopt selector_get_opts (const selector *sr);

int selector_get_depvar_number (const selector *sr);

int selector_get_VAR_order (const selector *sr);

const char *selector_entry_text (const selector *sr);

int selector_error (const selector *sr);

void maybe_clear_selector (const int *dlist);

GtkWidget *selector_get_window (const selector *sr);

#endif /* SELECTOR_H */

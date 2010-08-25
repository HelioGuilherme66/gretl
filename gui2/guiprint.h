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

/*  guiprint.h for gretl */ 

#ifndef GUIPRINT_H
#define GUIPRINT_H

#ifdef NATIVE_PRINTING
void print_window_content (char *fullbuf, char *selbuf, 
			   const char *fname);
# ifdef G_OS_WIN32
int win32_print_graph (char *emfname);
# else
void gtk_print_graph (const char *fname);
# endif
#endif

void special_print_summary (const Summary *summ,
			    const DATAINFO *pdinfo,
			    PRN *prn);

void special_print_vmatrix (const VMatrix *corr,
			    const DATAINFO *pdinfo, 
			    PRN *prn);

void special_print_fit_resid (const FITRESID *fr, 
			      const DATAINFO *pdinfo, 
			      PRN *prn);

void special_print_forecast (const FITRESID *fr, 
			     const DATAINFO *pdinfo, 
			     PRN *prn);

void special_print_confints (const CoeffIntervals *cf, 
			     PRN *prn);

int text_print_equation (const MODEL *pmod, const DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn);

int text_equation_ok (const MODEL *pmod);

int csv_to_clipboard (void);

int csv_selected_to_clipboard (void);

int copy_vars_formatted (windata_t *vwin, int fmt, int action);

int scalars_to_clipboard_as_csv (void);

int matrix_to_clipboard_as_csv (const gretl_matrix *m);

int font_has_minus (PangoFontDescription *desc);

#endif /* GUIPRINT_H */

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

#ifndef GRETL_MULTIPLOT_H_
#define GRETL_MULTIPLOT_H_

int gretl_multiplot_collecting (void);

int gretl_multiplot_start (gretlopt opt);

int gretl_multiplot_add_plot (gchar *buf);

int gretl_multiplot_finalize (gretlopt opt);

int gretl_multiplot_revise (gretlopt opt);

int gretl_multiplot_from_array (gretlopt opt);

int check_multiplot_options (gretlopt opt);

#endif /* GRETL_MULTIPLOT_H_ */

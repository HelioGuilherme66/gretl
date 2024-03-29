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

#ifndef ESTIM_PRIVATE_H
#define ESTIM_PRIVATE_H

double dwstat (int order, MODEL *pmod, const DATASET *dset);

double rhohat (int order, int t1, int t2, const double *uhat);

int check_for_effective_const (MODEL *pmod, const DATASET *dset);

void maybe_shift_ldepvar (MODEL *pmod, DATASET *dset);

MODEL ivreg_via_gmm (const int *list, DATASET *dset, gretlopt opt);

int get_x13as_maxpd (void);

#endif /* ESTIM_PRIVATE_H */

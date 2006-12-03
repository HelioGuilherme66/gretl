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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef MISSING_PRIVATE_H
#define MISSING_PRIVATE_H

#define masked(m,t) (m != NULL && m[t] == '1')
#define has_missing_obs(m)  ((m)->missmask != NULL)

void set_reference_missmask (const MODEL *pmod);

int apply_reference_missmask (MODEL *pmod);

int reference_missmask_present (void);

int undo_daily_repack (MODEL *pmod, double **Z, 
		       const DATAINFO *pdinfo);

int repack_missing_daily_obs (MODEL *pmod, double **Z, 
			      const DATAINFO *pdinfo);

int adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		 int n, const double **Z, int *misst);

#endif /* MISSING_PRIVATE_H */

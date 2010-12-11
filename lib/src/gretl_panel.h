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

#ifndef GRETL_PANEL_H_
#define GRETL_PANEL_H_

int panel_diagnostics (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn);

MODEL real_panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn);

MODEL panel_wls_by_unit (const int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn);

int panel_autocorr_test (MODEL *pmod, int order, 
			 double **Z, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn);

int groupwise_hetero_test (MODEL *pmod, DATAINFO *pdinfo,
			   gretlopt opt, PRN *prn);

int panel_tsls_robust_vcv (MODEL *pmod, const double **Z, 
			   const DATAINFO *pdinfo);

int set_panel_structure_from_vars (int uv, int tv, 
				   double **Z, 
				   DATAINFO *pdinfo);

int set_panel_structure_from_line (const char *line, 
				   double **Z, 
				   DATAINFO *pdinfo);

int switch_panel_orientation (double **Z, DATAINFO *pdinfo);

int balanced_panel (const DATAINFO *pdinfo);

int *panel_list_omit (const MODEL *orig, const int *drop, int *err);

int *panel_list_add (const MODEL *orig, const int *add, int *err);

int panel_variance_info (const double *x, const DATAINFO *pdinfo,
			 double xbar, double *psw, double *psb);

int plausible_panel_time_var (const double **Z, const DATAINFO *pdinfo);

int panel_isconst (int t1, int t2, int pd, const double *x,
		   int bygroup);

#endif /* GRETL_PANEL_H_ */

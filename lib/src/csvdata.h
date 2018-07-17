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

#ifndef CSVDATA_H
#define CSVDATA_H

typedef enum {
    AGGR_NONE,
    AGGR_COUNT,
    AGGR_AVG,
    AGGR_SUM,
    AGGR_MIN,
    AGGR_MAX,
    AGGR_SEQ,
    AGGR_MIDAS
} AggrType;

#define NON_NUMERIC 1.0e99

void import_na_init (void);

int import_na_string (const char *s);

int import_obs_label (const char *s);

int test_markers_for_dates (DATASET *dset, 
			    int *reversed, 
			    char *skipstr, 
			    PRN *prn);

int non_numeric_check (DATASET *dset, int **plist,
		       gretl_string_table **pst,
		       PRN *prn);

void reverse_data (DATASET *dset, PRN *prn);

int gretl_join_data (const char *fname,
		     const char **vnames,
		     int nvars,
		     DATASET *dset, 
		     const int *ikeyvars,
		     const char *okey,
		     const char *filter,
		     const char *srcname,
		     AggrType agg,
		     int seqval,
		     const char *auxname,
		     const char *tconvstr,
		     const char *tconvfmt,
		     int midas_pd,
		     gretlopt opt,
		     PRN *prn);

int probe_csv (const char *fname, char ***varnames,
	       int *nvars, gretlopt *opt);

void normalize_join_colname (char *targ, const char *src,
			     int underscore, int k);

int csv_open_needs_matrix (gretlopt opt);

#endif /* CSVDATA_H */

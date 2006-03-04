/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2005 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef DATAIO_H
#define DATAIO_H

#include <stdio.h>
#include <string.h>

typedef enum {
    GRETL_DATA_FLOAT = 1, /* single-precision binary data */
    GRETL_DATA_DOUBLE,    /* double-precision binary data */
    GRETL_DATA_OCTAVE,    /* data in Gnu Octave format */
    GRETL_DATA_CSV,       /* data in Comma Separated Values format */
    GRETL_DATA_R,         /* data in Gnu R format */
    GRETL_DATA_R_TS,      /* data in Gnu R format (time series) */
    GRETL_DATA_GZIPPED,   /* gzipped data */
    GRETL_DATA_TRAD,      /* traditional (ESL-style) data */
    GRETL_DATA_DAT,       /* data in PcGive format */
    GRETL_DATA_DB         /* gretl native database format */
} GretlDataFormat;

typedef enum {
    GRETL_NATIVE_DATA,    /* gretl native format data file */
    GRETL_XML_DATA,       /* gretl xml format data file */
    GRETL_CSV_DATA,       /* comma-separated data file */
    GRETL_BOX_DATA,       /* BOX1 format data file */
    GRETL_OCTAVE,         /* GNU octave ascii data file */
    GRETL_GNUMERIC,       /* gnumeric workbook data */
    GRETL_EXCEL,          /* MS Excel spreadsheet data */
    GRETL_WF1,            /* Eviews workfile data */
    GRETL_DTA,            /* Stata .dta data */
    GRETL_SCRIPT,         /* file containing gretl commands */
    GRETL_NATIVE_DB,      /* gretl database */
    GRETL_RATS_DB,        /* RATS 4.0 database */
    GRETL_UNRECOGNIZED    /* none of the above */
} GretlFileType;

typedef enum {
    CLEAR_FULL,           /* fully clear the dataset */
    CLEAR_SUBSAMPLE       /* dataset is sub-sampled: clear partially */
} DataClearCode;

typedef enum {
    DATA_NONE,    /* no dataset is currently open */
    DATA_CLEAR,   /* dataset is open: dataset info should be cleared */
    DATA_APPEND   /* dataset is open: attempt to append new data */
} DataOpenCode;

typedef enum {
    VARNAME_RESERVED = 1, /* vername is a gretl reserved name */
    VARNAME_FIRSTCHAR,    /* first character is not alphabetical */
    VARNAME_BADCHAR       /* illegal character in second or subsequent place */
} GretlVarnameError;

#define WORKSHEET_IMPORT(f) (f == GRETL_GNUMERIC || f == GRETL_EXCEL || \
                             f == GRETL_WF1 || f == GRETL_DTA)

#define free_datainfo(p) do { if (p != NULL) { clear_datainfo(p, 0); free(p); } \
                            } while (0);

#define DBNA  -999.0 /* missing value code for gretl databases */

#define GRETL_SCALAR_DIGITS 12

/* functions follow */

int dateton (const char *date, const DATAINFO *pdinfo);

char *ntodate (char *datestr, int nt, const DATAINFO *pdinfo);

char *ntodate_full (char *datestr, int t, const DATAINFO *pdinfo);

int get_info (const char *hdrfile, PRN *prn);

int get_precision (const double *x, int n, int placemax);

double get_date_x (int pd, const char *obs);

int write_data (const char *fname, const int *list, 
		const double **Z, const DATAINFO *pdinfo, 
	        gretlopt opt, PATHS *ppaths);

int data_report (const DATAINFO *pdinfo, PATHS *ppaths, PRN *prn);

int is_gzipped (const char *fname);

void gz_switch_ext (char *targ, char *src, char *ext);

int merge_data (double ***pZ, DATAINFO *pdinfo,
		double **addZ, DATAINFO *addinfo,
		PRN *prn);

int gretl_get_data (double ***pZ, DATAINFO **ppdinfo, 
		    char *datfile, PATHS *ppaths, 
		    DataOpenCode code, PRN *prn);

int open_nulldata (double ***pZ, DATAINFO *pdinfo, 
		   int data_status, int length,
		   PRN *prn);

int import_csv (double ***pZ, DATAINFO **ppdinfo, 
                const char *fname, PRN *prn);

int import_octave (double ***pZ, DATAINFO **ppdinfo, 
		   const char *fname, PRN *prn);

int import_box (double ***pZ, DATAINFO **ppdinfo, 
		const char *fname, PRN *prn);

int import_other (double ***pZ, DATAINFO **ppdinfo, 
		  int ftype, const char *fname, PRN *prn);

int add_obs_markers_from_file (DATAINFO *pdinfo, const char *fname);

GretlFileType detect_filetype (char *fname, PATHS *ppaths, PRN *prn);

int check_varname (const char *varname);

int check_atof (const char *numstr);

int transpose_data (double ***pZ, DATAINFO *pdinfo);

#endif /* DATAIO_H */

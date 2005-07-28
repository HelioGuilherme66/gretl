/*
 *  Copyright (c) by Allin Cottrell 2002-2004
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

/* import data from Eviews workfiles */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgretl.h"

static int get_data (FILE *fp, long pos, double **Z, int i, int n)
{
    double x;
    int t, nobs = 0;
    int got;
    int err = 0;

    fseek(fp, pos, SEEK_SET);
    got = fread(&nobs, sizeof nobs, 1, fp);
    if (got <= 0) {
	return 1;
    }

    /* should we be able to handle an offset here? */
    if (nobs != n) {
	fputs("problem: this does not match the specification "
	      " for the dataset\n", stderr);
    }

    fseek(fp, pos + 22, SEEK_SET);
    for (t=0; t<nobs; t++) {
	got = fread(&x, sizeof x, 1, fp);
	if (got <= 0) {
	    err = 1;
	    break;
	}
	if (x == 1e-37) {
	    Z[i][t] = NADBL;
	} else {
	    Z[i][t] = x;
	}
    }

    return err;
}

static int read_wf1_variables (FILE *fp, long pos, double ***pZ,
			       DATAINFO *dinfo, PRN *prn)
{
    int nv = dinfo->v + 1; /* RESID */
    char vname[32];
    short code;
    long u;
    int i, j = 1;
    int err = 0;

    for (i=0; i<nv && !err; i++, pos += 70) {
	/* read the 'code' for the 'object' (should be 44 for a regular
	   variable?) */
	fseek(fp, pos + 60, SEEK_SET);
	fread(&code, sizeof code, 1, fp);
	if (code == 43) {
	    /* constant: skip */
	    continue;
	} else if (code != 44) {
	    pprintf(prn, "byte %ld: unknown object code %d\n", 
		    pos + 60, (int) code);
	    continue;
	}

	/* grab the variable name */
	fseek(fp, pos + 20, SEEK_SET);
	fscanf(fp, "%31s", vname);
	if (!strcmp(vname, "C") || !strcmp(vname, "RESID")) {
	    continue;
	}
	fprintf(stderr, "Got variable '%s'\n", vname);
	dinfo->varname[j][0] = 0;
	strncat(dinfo->varname[j], vname, 8);

	/* get stream position for the data */
	fseek(fp, pos + 12, SEEK_SET);
	fread(&u, sizeof u, 1, fp);
	if (u > 0) {
	    /* follow up at the pos given above, if non-zero */
	    err = get_data(fp, u, *pZ, j++, dinfo->n);
	} else {
	    fputs("Couldn't find the data: skipping this variable\n", stderr);
	}
    }

    if (j < dinfo->v) {
	dataset_drop_last_variables(dinfo->v - j, pZ, dinfo);
    }

    return err;
}

static int parse_wf1_header (FILE *fp, DATAINFO *dinfo)
{
    int nvars = 0, nobs = 0, startyr = 0;
    short pd = 0, startper = 0;
    int err = 0;

    fseek(fp, 114, SEEK_SET);
    fread(&nvars, sizeof nvars, 1, fp);

    fseek(fp, 124, SEEK_SET);
    fread(&pd, sizeof pd, 1, fp);

    fseek(fp, 126, SEEK_SET);
    fread(&startper, sizeof startper, 1, fp);

    fseek(fp, 128, SEEK_SET);
    fread(&startyr, sizeof startyr, 1, fp);

    fseek(fp, 140, SEEK_SET);
    fread(&nobs, sizeof nobs, 1, fp);

    if (nvars <= 2 || nobs <= 0 || startyr <= 0 ||
	pd <= 0 || startper < 0) {
	err = E_DATA;
    } else {
	dinfo->v = nvars - 2; /* skip C and RESID */
	dinfo->n = nobs;
	dinfo->pd = pd;
    }

    fprintf(stderr, "header info:\n"
	    " number of variables = %d\n"
	    " number of observations = %d\n"
	    " data frequency = %d\n"
	    " starting year or major = %d\n"
	    " starting sub-period or minor = %d\n",
	    dinfo->v, dinfo->n, dinfo->pd,
	    startyr, startper);

    if (!err) {
	if (startper > 0) {
	    sprintf(dinfo->stobs, "%d:%d", startyr, startper);
	} else {
	    sprintf(dinfo->stobs, "%d", startyr);
	}

	if (dinfo->pd > 1 || startyr > 10) {
	    dinfo->structure = TIME_SERIES;
	}

	dinfo->sd0 = get_date_x(dinfo->pd, dinfo->stobs);
    }

    return err;
}

static int check_file_type (FILE *fp)
{
    char test[22] = {0};
    int err = 0;

    fread(test, 1, 21, fp);

    if (strcmp(test, "New MicroTSP Workfile")) {
	err = 1;
    }

    return err;
}

int wf1_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn)
{
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    int err = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    if (check_file_type(fp)) {
	fclose(fp);
	pputs(prn, "This file does not seem to be an Eviews workfile");
	return E_DATA;
    }

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    err = parse_wf1_header(fp, newinfo);

    if (!err) {
	err = start_new_Z(&newZ, newinfo, 0);
    }

    if (!err) {
	/* is the position 172 (always) right? */
	err = read_wf1_variables(fp, 172L, &newZ, newinfo, prn);
    }

    if (!err) {
	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}	

	if (*pZ == NULL) {
	    *pZ = newZ;
	    *pdinfo = *newinfo;
	} else {
	    err = merge_data(pZ, pdinfo, newZ, newinfo, prn);
	}
    }

    fclose(fp);

    return err;
}  

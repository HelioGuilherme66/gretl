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

/* subsample.c for gretl */

#include "libgretl.h"
#include "internal.h"

/* .......................................................... */

int attach_subsample_to_model (MODEL *pmod, double **fullZ, 
			       const DATAINFO *fullinfo)
     /* if the data set is currently subsampled, record the
	subsample info in pmod->subdum */
{
    int i, t, n = fullinfo->n;

    /* no subsample currently in force */
    if (fullZ == NULL) return 0;

    pmod->subdum = malloc(n * sizeof(double));
    if (pmod->subdum == NULL) return E_ALLOC;

    i = varindex(fullinfo, "subdum");
    if (i == fullinfo->v) { /* safety measure: should be impossible */
	fprintf(stderr, "mystery failure in attach_subsample_to_model\n");
	return 1;   
    } 

    for (t=0; t<n; t++)
	pmod->subdum[t] = (*fullZ)[i*n + t];
    
    return 0;
}

/* .......................................................... */

static int subsampled (const double *Z, const DATAINFO *pdinfo, 
		       const int subnum)
     /* Is the data set currently "sub-sampled" via selection of 
	cases?  Use this func _only_ if subnum tests < pdinfo->v */
{
    int t, n = pdinfo->n;

    for (t=0; t<n; t++)
	if (floatneq(Z[n*subnum + t], 0.0)) return 1;
    return 0;
}

/* .......................................................... */

int model_sample_issue (const MODEL *pmod, MODELSPEC *spec, 
			const double *Z, const DATAINFO *pdinfo)
     /* check a model (or modelspec) against the data info to see if 
	it may have been estimated on a different (subsampled) data 
	set from the current one */
{
    int i, n = pdinfo->n;
    double *subdum;
    extern int _identical (const double *x, const double *y, const int n);

    if (pmod == NULL && spec == NULL) return 0;

    /* if no sub-sampling has been done, we're OK */
    if ((i = varindex(pdinfo, "subdum")) == pdinfo->v) return 0;

    if (pmod != NULL) subdum = pmod->subdum;
    else subdum = spec->subdum;

    /* case: model has no sub-sampling info recorded */
    if (subdum == NULL) {
	/* if data set is not currently sub-sampled, we're OK */
	if (!subsampled(Z, pdinfo, i)) return 0;
	/* data set is subsampled, model is not: problem */
	else {
	    fprintf(stderr, "dataset is subsampled, model is not\n");
	    return 1;
	}
    }

    /* case: model has sub-sampling info recorded */
    if (!subsampled(Z, pdinfo, i)) {
	/* data set not sub-sampled: problem */
	fprintf(stderr, "model is subsampled, dataset is not\n");
	return 1;
    } else { /* do the subsamples (model and current data set) agree? */
	if (_identical(&Z[n*i], subdum, n))
	    return 0;
	else {
	    fprintf(stderr, "model and dataset subsamples not the same\n");
	    return 1;
	}
    }

    /* can't be reached */
    return 1;
}

/* .......................................................... */

int allocate_case_markers (char ***S, int n)
{
    int t;

    *S = malloc(n * sizeof(char *));
    if (*S == NULL) {
	return E_ALLOC;
    }
    for (t=0; t<n; t++) {
	(*S)[t] = malloc(9);
	if ((*S)[t] == NULL) {
	    free(*S);
	    return E_ALLOC;
	}
    }
    return 0;
}

/* .......................................................... */

static void prep_subdinfo (DATAINFO *dinfo, int markers, int n)
{
    dinfo->sd0 = 1.;
    dinfo->pd = 1;
    dinfo->time_series = 0;
    dinfo->extra = 0;
    if (markers) dinfo->markers = 1;
    else dinfo->markers = 0;
    strcpy(dinfo->stobs, "1");
    sprintf(dinfo->endobs, "%d", n);
}

#ifdef NOT_YET
/* .......................................................... */

int case_sample_direct (double **oldZ, double **newZ,
			DATAINFO *olddinfo, DATAINFO *newdinfo,
			double *dummy)
     /* subsample directly from a full-length dummy var */
{
    int i, t, st, sn, n = olddinfo->n;
    char **S = NULL;

    sn = 0;
    /* how many cases in sub-sample? */
    for (t=0; t<n; t++) {
	if (dummy[t] == 1.0) sn++;
    }
    /* FIXME handle pathological cases for sn? */

    newdinfo->n = sn;
    newdinfo->v = olddinfo->v;
    if (start_new_Z(newZ, newdinfo, 1)) return E_ALLOC;

    newdinfo->varname = olddinfo->varname;
    newdinfo->label = olddinfo->label;

    if (olddinfo->markers && allocate_case_markers(&S, sn)) {
	free(*newZ);
	return E_ALLOC;
    }

    /* copy across data and case markers, if any */
    st = 0;
    for (t=0; t<n; t++) {
	if (dummy[t] == 1.) {
	    for (i=1; i<olddinfo->v; i++) 
		(*newZ)[i * sn + st] = (*oldZ)[i * n + t];
	    if (olddinfo->markers) 
		strcpy(S[st], olddinfo->S[t]);
	    st++;
	}
    }

    prep_subdinfo(newdinfo, olddinfo->markers, sn);
    if (olddinfo->markers) newdinfo->S = S;

    return 0;
}
#endif

/* .......................................................... */

int set_sample_dummy (const char *line, 
		      double **oldZ, double **newZ,
		      DATAINFO *oldinfo, DATAINFO *newinfo,
		      char *msg, const int opt)
     /* sub-sample the data set, based on the criterion of skipping
	all observations with missing data values, or using as a
	mask a specified dummy variable, or masking with a specified
	boolean condition */
{
    double xx, *dum = NULL;
    char **S, dumv[9];
    int missobs = 0, subnum, dumnum;
    int i, t, st, sn, n = oldinfo->n;

    *msg = '\0';
    dumv[0] = '\0';
    if (opt == OPT_O && 
	(line == NULL || sscanf(line, "%*s %s", dumv) <= 0))
	missobs = 1; 

    if (missobs) {
	/* construct missing obs dummy on the fly */
	dum = malloc(n * sizeof *dum);
	if (dum == NULL) return E_ALLOC;
	sn = 0;
	for (t=0; t<n; t++) {
	    dum[t] = 1.0;
	    for (i=1; i<oldinfo->v; i++) {
		if (na((*oldZ)[i*n + t])) {
		    dum[t] = 0.;
		    break;
		}
	    }
	    if (floateq(dum[t], 1.0)) sn++;
	}
    } else if (opt == OPT_O) {  
	/* the name of a dummy variable was passed in */
	dumnum = varindex(oldinfo, dumv);
	if (dumnum == oldinfo->v) {
	    sprintf(msg, "Variable '%s' not defined", dumv);
	    return 1;
	} 
	sn = isdummy(dumnum, oldinfo->t1, oldinfo->t2, *oldZ, n);
    } else if (opt == OPT_R) {
	/* construct dummy from boolean expression */	
	GENERATE genr;
	char formula[MAXLEN];

	/* + 4 below to omit the word "smpl" */
	sprintf(formula, "genr subdum=%s", line + 4);
	genr = genr_func(oldZ, oldinfo, formula, 0, NULL, 1);
	if (genr.errcode) {
	    strcpy(msg, genr.errmsg);
	    return 1;
	}
	if (add_new_var(oldinfo, oldZ, &genr)) {
	    strcpy(msg, "Failed to add sub-sampling dummy variable");
	    return 1;
	}
	subnum = varindex(oldinfo, "subdum");
	dumnum = subnum;
	sn = isdummy(subnum, oldinfo->t1, oldinfo->t2, *oldZ, n);
    } else {
	/* impossible */
	strcpy(msg, "Sub-sample command failed mysteriously");
	return 1;
    }

    /* does this policy lead to an empty sample, or no change
       in the sample, perchance? */
    if (sn == 0) {
	if (opt == OPT_O && !missobs)
	    sprintf(msg, "'%s' is not a dummy variable", dumv);
	else if (missobs)
	    strcpy(msg, "No observations would be left!");
	else { /* case of boolean expression */
	    if ((*oldZ)[subnum * n + oldinfo->t1] == 0)
		strcpy(msg, "No observations would be left!");
	    else
		strcpy(msg, "No observations were dropped!");
	}
	return 1;
    }
    if (sn == n) {
	strcpy(msg, "No observations were dropped!");
	return 1;
    }

    /* create or reuse "hidden" dummy to record sub-sample */
    subnum = varindex(oldinfo, "subdum");
    if (subnum == oldinfo->v) {
	if (grow_Z(1, oldZ, oldinfo)) return E_ALLOC;
	strcpy(oldinfo->varname[subnum], "subdum");
	strcpy(oldinfo->label[subnum], "automatic sub-sampling dummy");
    }
    for (t=0; t<n; t++) {
	if (missobs) 
	    (*oldZ)[subnum * n + t] = dum[t];
	else if (opt == OPT_O)
	    /* ?possibility of missing values here? */
	    (*oldZ)[subnum * n + t] = (*oldZ)[dumnum * n + t];
    }

    newinfo->n = sn;
    newinfo->v = oldinfo->v;
    if (start_new_Z(newZ, newinfo, 1)) {
	if (dum != NULL) free(dum);
	return E_ALLOC;
    }

    /* link varnames and descriptions (not dependent on series length) */
    newinfo->varname = oldinfo->varname;
    newinfo->label = oldinfo->label;

    /* case markers */
    if (oldinfo->markers && allocate_case_markers(&S, sn)) {
	free(*newZ);
	free(dum);
	return E_ALLOC;
    }

    /* copy across data and case markers, if any */
    st = 0;
    for (t=0; t<n; t++) {
	xx = (missobs)? dum[t] : (*oldZ)[dumnum * n + t];
	if (xx == 1.) {
	    for (i=1; i<oldinfo->v; i++) 
		(*newZ)[i * sn + st] = (*oldZ)[i * n + t];
	    if (oldinfo->markers) 
		strcpy(S[st], oldinfo->S[t]);
	    st++;
	}
    }

    prep_subdinfo(newinfo, oldinfo->markers, sn);
    if (oldinfo->markers) newinfo->S = S;

    free(dum);

    return 0;
}

/* .......................................................... */

int set_sample (const char *line, DATAINFO *pdinfo, char *msg)
{
    int nf, new_t1, new_t2;
    char cmd[5], newstart[8], newstop[8];

    nf = count_fields(line);

    if (nf == 1) return 0;
	
    if (nf == 2) {
	if (sscanf(line, "%s %s", cmd, newstart) != 2) {
	    sprintf(msg, "error reading smpl line");
	    return 1;
	} else {
	    new_t1 = dateton(newstart, pdinfo->pd, pdinfo->stobs, msg);
	    if (new_t1 < 0 || strlen(msg)) return 1;
	    if (new_t1 > pdinfo->n) {
		sprintf(msg, "error in new starting obs");
		return 1;
	    }
	    pdinfo->t1 = new_t1;
	    return 0;
	}
    }
    if (sscanf(line, "%s %s %s", cmd, newstart, newstop) != 3) {
	sprintf(msg, "error reading smpl line");
	return 1;
    }
    if (strcmp(newstart, ";")) {
	new_t1 = dateton(newstart, pdinfo->pd, pdinfo->stobs, msg);
	if (new_t1 < 0 || strlen(msg)) {
	    return 1;
	}
	pdinfo->t1 = new_t1;
    }
    if (strcmp(newstop, ";")) {
	new_t2 = dateton(newstop, pdinfo->pd, pdinfo->stobs, msg);
	if (strlen(msg)) return 1;
	if (new_t2 >= pdinfo->n) {
	    sprintf(msg, "error in new ending obs");
	    return 1;
	}
	pdinfo->t2 = new_t2;
    }
    return 0;
}

/* ........................................................... */

static int datamerge (double **fullZ, DATAINFO *fullinfo,
		      double **subZ, DATAINFO *subinfo)
{
    int i, t, dumn, subt;
    int newvars = subinfo->v - fullinfo->v;
    int subn = subinfo->n, n = fullinfo->n;
    double *newZ = NULL;

    if (newvars <= 0) return 0;

    dumn = varindex(subinfo, "subdum");
    if (dumn == subinfo->v) return E_NOMERGE;

    newZ = realloc(*fullZ, subinfo->v * n * sizeof **fullZ);  
    if (newZ == NULL) return E_ALLOC;
    else *fullZ = newZ;

    subt = 0;
    for (t=0; t<n; t++) {
	if ((*fullZ)[dumn * n + t] == 1.0) {
	    for (i=fullinfo->v; i<subinfo->v; i++) {
		(*fullZ)[i*n + t] = (*subZ)[i*subn + subt];
	    }
	    subt++;
	} else {
	    for (i=fullinfo->v; i<subinfo->v; i++) { 
		(*fullZ)[i*n + t] = NADBL;
	    }
	}
    }

    fullinfo->v = subinfo->v;
    fullinfo->varname = subinfo->varname;
    fullinfo->label = subinfo->label;

    return 0;
}

/* ........................................................... */

int restore_full_sample (double **subZ, double **fullZ, double **Z,
			 DATAINFO **subinfo, DATAINFO **fullinfo,
			 DATAINFO **datainfo, char *msg)
{
    int i, t, n, err = 0;

    if (*subZ == NULL) {
        (*datainfo)->t1 = 0;
        (*datainfo)->t2 = (*datainfo)->n - 1;
        return 0;
    }

    n = (*fullinfo)->n;

    /* in case any new vars added, try to merge them in */
    err = datamerge(fullZ, *fullinfo, Z, *subinfo);
    if (err == E_ALLOC)
        sprintf(msg, "Out of memory expanding data set\n");
    if (err == E_NOMERGE)
        sprintf(msg, "Missing sub-sample information; can't merge data\n");

    /* zero out the "subdum" dummy variable */
    i = varindex(*fullinfo, "subdum");
    if (i < (*fullinfo)->v)
        for (t=0; t<n; t++) (*fullZ)[i*n + t] = 0.;

    /* reorganize pointers for data set */
    *subZ = *Z;
    *Z = *fullZ;
    free(*subZ);
    *subZ = NULL;
    *fullZ = NULL;
    /* and data info struct */
    *subinfo = *datainfo;
    *datainfo = *fullinfo;
    clear_datainfo(*subinfo, 1);
    /* FIXME should free subinfo and/or fullinfo? Maybe not! */
    free(*subinfo);
    *subinfo = NULL;
    *fullinfo = NULL;
    
    return 0;
}

/* ........................................................... */

int count_missing_values (double **pZ, DATAINFO *pdinfo, print_t *prn)
{
    int i, v, t, n = pdinfo->n;
    int missval = 0, missobs = 0, oldmiss = 0, tmiss;
    int year, yearmiss = 0, totvals = 0, yearbak;

    v = varindex(pdinfo, "year");
    if (v == pdinfo->v) v = varindex(pdinfo, "YEAR");
    if (v == pdinfo->v) v = 0;
    else yearbak = (int) (*pZ)[v*n + pdinfo->t1];

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	tmiss = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (hidden_var(i, pdinfo)) continue;
	    if (na((*pZ)[i*n + t])) missval++;
	    totvals++;
	}
	if ((tmiss = missval - oldmiss)) missobs++;
	if (v) {
	    year = (int) (*pZ)[v*n + t];
	    if (year != yearbak) {
		pprintf(prn, "%d: %4d missing data values\n", 
			yearbak, yearmiss);
		yearmiss = tmiss;
		yearbak = year;
	    } else
		yearmiss += tmiss;
	}
	oldmiss = missval;
    }
    if (v) 
	pprintf(prn, "%d: %4d missing data values\n", 
		year, yearmiss);
    
    pprintf(prn, "\nNumber of observations (rows) with missing data "
	    "values = %d (%d%%)\n", missobs, 
	    (int) (100.0 * missobs / (pdinfo->t2 - pdinfo->t1 + 1)));
    pprintf(prn, "Total number of missing data values = %d (%d%% "
	    "of total data values)\n", missval, 
	    (int) (100.0 * missval / totvals));
    return missval;
}

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

/* subsample.c for gretl */

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_panel.h"
#include "cmd_private.h"
#include "dbread.h"
#include "uservar.h"
#include "gretl_xml.h"

#define SUBDEBUG 0
#define FULLDEBUG 0

typedef enum {
    SUBSAMPLE_NONE,
    SUBSAMPLE_DROP_MISSING,
    SUBSAMPLE_DROP_EMPTY,
    SUBSAMPLE_DROP_WKENDS,
    SUBSAMPLE_USE_DUMMY,
    SUBSAMPLE_BOOLEAN,
    SUBSAMPLE_RANDOM,
    SUBSAMPLE_UNKNOWN
} SubsampleMode;

/*
  The purpose of the static pointers below: When the user subsamples
  the current dataset in a non-trivial way -- i.e., by selecting cases
  rather than just moving the starting or ending points of the data
  range -- we create a new sub-dataset, and we need to keep the full
  dataset around so that it can be restored later.  The pointer
  @fullset is used to record the address of the full dataset.

  In addition, @peerset keeps track of the location of the DATASET
  struct associated with the backed-up full dataset; by means of this,
  we can know when to free the full dataset and when not to (for
  instance, if we're freeing an auxiliary dataset).
*/

static DATASET *fullset;
static DATASET *peerset;

#define SUBMASK_SENTINEL 127

static int smpl_get_int (const char *s, DATASET *dset, int *err);
static int test_set_sample (const char *s, DATASET *dset,
			    int *t1, int *t2);

/* accessors for full dataset, when sub-sampled */

DATASET *fetch_full_dataset (void)
{
    return fullset;
}

static int get_submask_length (const char *s)
{
    int n = 1;

    if (s == NULL || s == RESAMPLED) {
	n = 0;
    } else {
	while (*s != SUBMASK_SENTINEL) {
	    n++;
	    s++;
	}
    }

    return n;
}

void free_subsample_mask (char *s)
{
    if (s != RESAMPLED && s != NULL) {
	free(s);
    }
}

char *copy_subsample_mask (const char *src, int *err)
{
    char *ret = NULL;

    if (src == RESAMPLED) {
	ret = RESAMPLED;
    } else if (src != NULL) {
	int n = get_submask_length(src);

	ret = malloc(n * sizeof *ret);
	if (ret != NULL) {
	    memcpy(ret, src, n);
	} else {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

char *copy_datainfo_submask (const DATASET *dset, int *err)
{
    char *mask = NULL;

    if (complex_subsampled()) {
	mask = copy_subsample_mask(dset->submask, err);
    }

    return mask;
}

int write_model_submask (const MODEL *pmod, FILE *fp)
{
    int ret = 0;

    if (pmod->submask == RESAMPLED) {
	fputs("<submask length=\"0\"></submask>\n", fp);
	ret = 1;
    } else if (pmod->submask != NULL) {
	int i, n = get_submask_length(pmod->submask);

	fprintf(fp, "<submask length=\"%d\">", n);
	for (i=0; i<n; i++) {
	    fprintf(fp, "%d ", (int) pmod->submask[i]);
	}
	fputs("</submask>\n", fp);
	ret = 1;
    }

    return ret;
}

int write_datainfo_submask (const DATASET *dset, FILE *fp)
{
    int ret = 0;

    if (dset->submask == RESAMPLED) {
	unsigned int seed = get_resampling_seed();

	fprintf(fp, "<resample seed=\"%u\" n=\"%d\"/>\n", seed, dset->n);
	ret = 1;
    } else if (complex_subsampled()) {
	int i, n = get_submask_length(dset->submask);

	fprintf(fp, "<submask length=\"%d\">", n);
	for (i=0; i<n; i++) {
	    fprintf(fp, "%d ", (int) dset->submask[i]);
	}
	fputs("</submask>\n", fp);

	if (dset->restriction != NULL) {
	    gretl_xml_put_tagged_string("restriction", dset->restriction, fp);
	}

	ret = 1;
    }

    return ret;
}

int submask_cmp (const char *m1, const char *m2)
{
    if (m1 == NULL && m2 == NULL) {
	return 0;
    } else if (m1 == NULL || m2 == NULL) {
	return 1;
    }
	
    if (m1 == RESAMPLED || m2 == RESAMPLED) {
	return m1 != RESAMPLED || m2 != RESAMPLED;
    }

    while (*m1 != SUBMASK_SENTINEL && *m2 != SUBMASK_SENTINEL) {
	if (*m1 != *m2) {
	    return 1;
	}
	m1++;
	m2++;
    }

    return 0;
}

/* all values apart from the sentinel are initialized to zero; once
   the mask is used, 1s will indicate included observations and
   0s will indicate excluded observations
*/

static char *make_submask (int n)
{
    char *mask = calloc(n + 1, 1);
    
    if (mask != NULL) {
	mask[n] = SUBMASK_SENTINEL;
    }

    return mask;
}

#if 0

static void debug_print_submask (char *mask, char *msg)
{
    if (mask == NULL) {
	fprintf(stderr, "%s: NULL\n", msg);
    } else if (mask == RESAMPLED) {
	fprintf(stderr, "%s: RESAMPLED\n", msg);
    } else {
	char *s = mask;

	fprintf(stderr, "%s: ", msg);
	
	while (*s != SUBMASK_SENTINEL) {
	    if (*s == 0) {
		fputc('0', stderr);
	    } else if (*s == 1) {
		fputc('1', stderr);
	    } else {
		fputc('?', stderr);
	    }
	    s++;
	}
	fputc('\n', stderr);
    }
}

#endif

void set_dataset_resampled (DATASET *dset, unsigned int seed)
{
    dset->submask = RESAMPLED;
    dset->rseed = seed;
}

int dataset_is_resampled (const DATASET *dset)
{
    return (dset != NULL && dset->submask == RESAMPLED);
}

void maybe_free_full_dataset (const DATASET *dset)
{
    if (dset == peerset) {
	if (fullset != NULL) {
#if SUBDEBUG
	    fprintf(stderr, "maybe_free_full_dataset: freeing fullset at %p (Z at %p)\n",
		    (void *) fullset, (void *) fullset->Z);
#endif
	    if (fullset->Z != NULL) {
		free_Z(fullset);
	    }
	    clear_datainfo(fullset, CLEAR_SUBSAMPLE);
	    free(fullset);
	    fullset = NULL;
	}
	peerset = NULL;
    }
}

/* we do this on "restore full sample" */

static void relink_to_full_dataset (DATASET *dset)
{
#if SUBDEBUG
    fprintf(stderr, "relink_to_full_dataset: fullset = %p (freeing and nulling)\n",
	    (void *) fullset);
    fprintf(stderr, "fullset: v = %d, n = %d\n", fullset->v, fullset->n);
#endif

    *dset = *fullset;
    free(fullset);
    fullset = NULL;
    peerset = NULL;
}

/* sync malloced elements of the fullset struct that might
   have been moved via realloc
*/

static void sync_datainfo_members (const DATASET *dset)
{
    if (fullset->v > dset->v) {
	int i;

#if FULLDEBUG
	fprintf(stderr, "*** sync_datainfo_members: fullset->v = %d but dset->v = %d\n",
		fullset->v, dset->v);
	fprintf(stderr, " deleting the last %d element(s) of fullZ\n", 
		fullset->v - dset->v);
#endif
	for (i=dset->v; i<fullset->v; i++) {
	    free(fullset->Z[i]);
	    fullset->Z[i] = NULL;
	}
	fullset->v = dset->v;
    }

    fullset->varname = dset->varname;
    fullset->varinfo = dset->varinfo;
    fullset->descrip = dset->descrip;
}

/* attach_subsample_to_model:
 * @pmod: model to which subsample should be attached.
 * @dset: pointer to current dataset info.
 *
 * If the dataset is currently subsampled, record the subsample
 * information with the model so that it can be retrieved later.
 * 
 * Returns: 0 if the recording is not needed, or on success; non-zero
 * error code failure.
 */

int attach_subsample_to_model (MODEL *pmod, const DATASET *dset)
{
    int err = 0;

#if SUBDEBUG
    fprintf(stderr, "attach_subsample_to_model: fullset = %p\n",
	    (void *) fullset);
#endif

    if (fullset != NULL) {
	/* sync, in case anything has moved */
	sync_datainfo_members(dset);

	if (pmod->submask != NULL) {
	    free_subsample_mask(pmod->submask);
	}

	pmod->submask = copy_subsample_mask(dset->submask, &err);
    }

    return err;
}

#define MDEBUG 0

/* This is called from objstack.c for each saved model, to
   determine whether a proposed permanent shrinkage of the
   dataset will invalidate the model.
*/

int subsample_check_model (MODEL *pmod, char *mask)
{
    if (submask_cmp(pmod->submask, mask)) {
	return E_DATA;
    } else {
	return 0;
    }
}

/* Called from objstack.c for each saved model, when the
   user calls for the imposition of a permanent subsample
   restriction on the dataset. We've already checked that
   the "new" subsample matches that on which model was
   estimated, so all we have to do now is sync by removing
   the model's subsample mask.
*/

int remove_model_subsample_info (MODEL *pmod)
{
    free_subsample_mask(pmod->submask);
    pmod->submask = NULL;
    
    return 0;
}

/* If series have been added to a resampled dataset, we can't
   bring these back to the "full" dataset, which may have a
   longer or shorter series length, and from which there is
   no definite mapping by row. So we just delete them. In
   this function we destroy their varnames and varinfo
   structures; the numerical arrays get deleted later.
*/

static int resample_sync_dataset (DATASET *dset)
{
    if (dset->v > fullset->v) {
	char **varname;
	VARINFO **varinfo;
	int i, nv = fullset->v;

	for (i=fullset->v; i<dset->v; i++) {
	    free(dset->varname[i]);
	    free(dset->varinfo[i]);
	}

	varname = realloc(dset->varname, nv * sizeof *varname);
	if (varname == NULL) {
	    return E_ALLOC;
	}
	dset->varname = varname;

	varinfo = realloc(dset->varinfo, nv * sizeof *varinfo);
	if (varinfo == NULL) {
	    return E_ALLOC;
	}
	dset->varinfo = varinfo;
    }

    /* sync */
    fullset->varname = dset->varname;
    fullset->varinfo = dset->varinfo;
    fullset->descrip = dset->descrip;

    return 0;
}

/* Apparatus for updating full dataset when restoring full sample
   after sub-sampling.  
*/

static void
update_full_data_values (const DATASET *dset)
{
    int i, s, t;

#if SUBDEBUG
    fprintf(stderr, "update_full_data_values: fullset->Z=%p, dset->Z=%p, dset=%p\n",
	    (void *) fullset->Z, (void *) dset->Z, (void *) dset);
#endif

    for (i=1; i<fullset->v && i<dset->v; i++) {
	s = 0;
	for (t=0; t<fullset->n; t++) {
	    if (dset->submask[t] == 1) {
		fullset->Z[i][t] = dset->Z[i][s++];
	    } else if (dset->submask[t] == 'p') {
		/* skip panel padding (?) */
		s++;
	    }
	}
    }
}

static int update_case_markers (const DATASET *dset)
{
    int err = 0;

    if (fullset->markers == DAILY_DATE_STRINGS) {
	; /* don't mess with them */
    } else if (dated_daily_data(fullset)) {
	; /* again, don't mess! */
    } else if (dset->markers && !fullset->markers) {
	dataset_allocate_obs_markers(fullset);
	if (fullset->S == NULL) {
	    err = 1;
	} else {
	    int t, subt = 0;

	    for (t=0; t<fullset->n; t++) {
		if (dset->submask[t]) {
		    strcpy(fullset->S[t], dset->S[subt++]);
		} else {
		    sprintf(fullset->S[t], "%d", t + 1);
		}
	    }
	}
    }	
	
    return err;
}

static int add_new_vars_to_full (DATASET *dset)
{
    int V1 = dset->v;
    int V0 = fullset->v;
    int N = fullset->n;
    double **newZ = NULL;
    int i, t, s;
    int err = 0;

    if (V1 <= V0) {
	return 0;
    }

    if (dset->submask == NULL) {
	return E_NOMERGE;
    }

#if SUBDEBUG
    fprintf(stderr, "add_new_vars_to_full:\n");
    fprintf(stderr, " V1 = dset->v = %d; V0 = fullset->v = %d\n",
	    V1, V0);
    fprintf(stderr, " dset->Z = %p, fullset->Z = %p\n", (void *) dset->Z, 
	    (void *) fullset->Z);
#endif

    /* allocate expanded data array */
    newZ = realloc(fullset->Z, V1 * sizeof *fullset->Z);

    if (newZ == NULL) {
	return E_ALLOC;
    } 

    fullset->Z = newZ;

    for (i=V0; i<dset->v && !err; i++) {
#if FULLDEBUG
	fprintf(stderr, "adding to full: var %d (%s, level %d)\n",
		i, dset->varname[i], series_get_stack_level(dset, i));
#endif
	fullset->Z[i] = malloc(N * sizeof **fullset->Z);
	if (fullset->Z[i] == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "After expansion, fullset->Z = %p (%d vars)\n", 
	    (void *) fullset->Z, dset->v);
#endif

    for (i=V0; i<dset->v; i++) {
	s = 0;
	for (t=0; t<N; t++) {
	    fullset->Z[i][t] = (dset->submask[t])? 
		dset->Z[i][s++] : NADBL;
	}
    }

    fullset->v = V1;

    return 0;
}

static int sync_data_to_full (DATASET *dset)
{
    int err;

    /* update values for pre-existing series, which may have been
       modified via "genr" etc. */
    update_full_data_values(dset);

    /* if case markers were added when subsampled, carry them back */
    update_case_markers(dset);

    /* delete any newly added hidden vars */
    err = dataset_destroy_hidden_variables(dset, fullset->v);

    /* in case any new vars were added when subsampled, try to merge
       them into the full dataset */
    if (!err) {
	err = add_new_vars_to_full(dset);
    }

    return err;
}

/* Here we make a mask representing the "complete" sample
   restriction currently in force, for use in cumulating
   restrictions. "Complete" means that we take into account
   both an existing subsampling restriction, if any, and
   possible setting of the t1 and t2 members of @dset to
   exclude certain observation ranges.

   The mask returned, if non-NULL, will be the length of
   the full dataset. It will be NULL (without error) if
   the current dataset is neither subsampled nor range-
   restricted.
*/

static char *make_current_sample_mask (DATASET *dset, int *err)
{
    char *currmask = NULL;
    int range_set;
    int s, t;

    range_set = (dset->t1 > 0 || dset->t2 < dset->n - 1);

    if (dset->submask == NULL) {
	/* no pre-existing mask, so not subsampled, but
	   we should restrict currmask to observations
	   included in the current (t1, t2) range, if
	   applicable
	*/
	if (range_set) {
	    currmask = make_submask(dset->n);
	    if (currmask == NULL) {
		*err = E_ALLOC;
	    } else {
		for (t=dset->t1; t<=dset->t2; t++) {
		    currmask[t] = 1;
		}
	    }
	}
    } else {
	/* there's a pre-existing mask: in addition we
	   mask out observations outside of the current
	   (t1, t2) range, if applicable
	*/
	currmask = copy_subsample_mask(dset->submask, err);
	if (!*err && range_set) {
	    s = -1;
	    for (t=0; t<fullset->n; t++) {
		if (dset->submask[t]) s++;
		if (s < dset->t1 || s > dset->t2) {
		    currmask[t] = 0;
		} 
	    }
	}
    }

    return currmask;
}

/* Deal with the case where sampling has been done simply by
   moving the sample-range endpoints
*/

static int restore_full_easy (DATASET *dset, ExecState *state)
{
    int t1min, t2max;

    if (state == NULL) {
	/* not inside a function */
	t1min = 0;
	t2max = dset->n - 1;
    } else {
	/* don't go outside the bounds set on entry to
	   a function */
	sample_range_get_extrema(dset, &t1min, &t2max);
    }

    if (dset->t1 != t1min || dset->t2 != t2max) {
	dset->t1 = t1min;
	dset->t2 = t2max;
#if SUBDEBUG
	fprintf(stderr, "restore_full_sample: just set t1=%d and t2=%d\n",
		t1min, t2max);
#endif
    }

    return 0;
}

/* restore_full_sample: 
 * @dset: dataset struct.
 * @state: structure representing program execution state,
 * or %NULL. 
 *
 * Restore the full data range, undoing any sub-sampling that was
 * previously performed (either by shifting the endpoints of the
 * sample range or by selecting cases on some criterion or other).
 * However, if @state is not %NULL, and if it contains a submask,
 * the "full" range is relative to that mask.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restore_full_sample (DATASET *dset, ExecState *state)
{
    int err = 0;

    if (dset == NULL) {
	return E_NODATA;
    } else if (!complex_subsampled()) {
	return restore_full_easy(dset, state);
    }

    if (dset != peerset) {
	fprintf(stderr, "restore_full_sample: dset is not peerset!\n");
	return E_DATA;
    }

#if FULLDEBUG || SUBDEBUG
    fprintf(stderr, "\nrestore_full_sample: dset=%p, state=%p, fullset=%p\n", 
	    (void *) dset, (void *) state, (void *) fullset);
#endif

    /* Beyond this point we are doing a non-trivial restoration
       of a stored "full" dataset which has previously been
       subsampled, e.g., by some boolean criterion.
    */

    if (dset->submask == RESAMPLED) {
	err = resample_sync_dataset(dset);
    } else {
	if (dset->padmask != NULL) {
	    fprintf(stderr, "restore_full_sample: first undo panel padding\n");
	    err = undo_panel_padding(dset);
	}
	if (!err) {
	    sync_datainfo_members(dset);
	    err = sync_data_to_full(dset);
	}
    }

    if (err == E_NOMERGE) {
        gretl_errmsg_set(_("Missing sub-sample information; can't merge data\n"));
    }

    if (err) {
	return err;
    }

    /* destroy sub-sampled data array */
#if SUBDEBUG
    fprintf(stderr, "restore_full_sample: freeing sub-sampled Z at %p (v = %d, n = %d)\n"
	    " and clearing dset at %p\n", (void *) dset->Z, dset->v, dset->n,
	    (void *) dset);
#endif
    free_Z(dset);
    clear_datainfo(dset, CLEAR_SUBSAMPLE);

    relink_to_full_dataset(dset);

    if (state != NULL) {
	/* in this case restoring the "full" sample really means, relative 
	   to the original state 
	*/
	if (state->submask != NULL) {
	    err = restrict_sample_from_mask(state->submask, dset, OPT_NONE);
	} else {
	    int t1min, t2max;

	    sample_range_get_extrema(dset, &t1min, &t2max);
	    if (dset->t1 < t1min) {
		dset->t1 = t1min;
	    }
	    if (dset->t2 > t2max) {
		dset->t2 = t2max;
	    }
	}
    }

    return err;
}

static int overlay_masks (char *targ, const char *src, int n)
{
    int i, sn = 0;
    
    for (i=0; i<n; i++) {
	if (targ[i] == 1 && src[i] == 1) {
	    targ[i] = 1;
	    sn++;
	} else {
	    targ[i] = 0;
	}
    }
	
    return sn;
}

static int make_weekday_mask (const DATASET *dset, char *mask)
{
    if (!dated_daily_data(dset) || 
	dset->markers != 0 || dset->pd < 6) {
	return E_PDWRONG;
    } else {
	char datestr[OBSLEN];
	int t, wd;

	for (t=0; t<dset->n; t++) {
	    ntodate(datestr, t, dset);
	    wd = weekday_from_date(datestr);
	    mask[t] = (wd >= 1 && wd <= 5);
	}
	
	return 0;
    }
}

/* write into @mask: 0 for observations at which _all_ variables
   (or all variables in @list, if @list is non-NULL) have missing 
   values, 1 for all other observations (that have at least one 
   non-NA). (Refinement: if no list is given, we ignore the generic
   variables "time" and "index".)
*/

static int 
make_empty_mask (const int *list, const DATASET *dset, char *mask)
{
    int i, vi, t, vt;

    if (list != NULL && list[0] > 0) {
	/* check specified list of variables */
	for (t=0; t<dset->n; t++) {
	    mask[t] = 0;
	    for (i=1; i<=list[0]; i++) {
		vi = list[i];
		if (!na(dset->Z[vi][t])) {
		    mask[t] = 1;
		    break;
		}
	    }
	}
    } else {
	/* check (almost) all variables */
	vt = current_series_index(dset, "time");
	vi = current_series_index(dset, "index");

	for (t=0; t<dset->n; t++) {
	    mask[t] = 0;
	    for (i=1; i<dset->v; i++) {
		if (!series_is_hidden(dset, i) &&
		    i != vt && i != vi &&
		    !na(dset->Z[i][t])) {
		    mask[t] = 1;
		    break;
		}
	    }
	}
    }

    return 0;
}

/* write into @mask: 0 for observations at which _any_ variable
   (or any variable in @list, if @list is non-NULL) has a missing 
   value, 1 for all other observations.
*/

static int 
make_missing_mask (const int *list, const DATASET *dset, char *mask)
{
    int i, vi, t;

    if (list != NULL && list[0] > 0) {
	/* check specified list of variables */
	for (t=0; t<dset->n; t++) {
	    mask[t] = 1;
	    for (i=1; i<=list[0]; i++) {
		vi = list[i];
		if (na(dset->Z[vi][t])) {
		    mask[t] = 0;
		    break;
		}
	    }
	}
    } else {	
	/* check all variables */
	for (t=0; t<dset->n; t++) {
	    mask[t] = 1;
	    for (i=1; i<dset->v; i++) {
		if (!series_is_hidden(dset, i) && na(dset->Z[i][t])) {
		    mask[t] = 0;
		    break;
		}
	    }
	}
    }

    return 0;
}

static int copy_dummy_to_mask (char *mask, const double *x, int n)
{
    int t, err = 0;

    for (t=0; t<n && !err; t++) {
	if (x[t] == 1.0) {
	    mask[t] = 1;
	} else if (!na(x[t]) && x[t] != 0.0) { /* NA? */
	    err = 1;
	}
    }

    return err;
}

static int mask_from_temp_dummy (const char *s, DATASET *dset, 
				 char *mask, PRN *prn)
{
    char formula[MAXLINE];
    double *x;
    int err = 0;

    *formula = '\0';
    strncat(formula, s, MAXLINE - 1);

    x = generate_series(formula, dset, prn, &err);

    if (!err) {
	err = copy_dummy_to_mask(mask, x, dset->n);
	if (err) {
	    gretl_errmsg_sprintf(_("'%s' is not a dummy variable"), "mask");
	}
    }

    free(x);

    return err;
}

static int mask_from_dummy (const char *s, const DATASET *dset,
			    char *mask)
{
    char dname[VNAMELEN] = {0};
    int dnum, err = 0;

    gretl_scan_varname(s, dname);
    dnum = series_index(dset, dname);

    if (dnum == dset->v) {
	gretl_errmsg_sprintf(_("Variable '%s' not defined"), dname);
	err = 1;
    } else {
	err = copy_dummy_to_mask(mask, dset->Z[dnum], dset->n);
	if (err) {
	    gretl_errmsg_sprintf(_("'%s' is not a dummy variable"), dname);
	}
    }

    return err;
}

/* how many observations are selected by the given 
   subsample mask? */

static int 
count_selected_cases (const char *mask, const DATASET *dset)
{
    int i, n = 0;

    for (i=0; i<dset->n; i++) {
	if (mask[i]) {
	    n++;
	}
    }

    return n;
}

/* panel: how many distinct cross-sectional units are included 
   in the masked subset of observations? */

static int 
count_panel_units (const char *mask, const DATASET *dset)
{
    int u, ubak = -1;
    int i, n = 0;

    for (i=0; i<dset->n; i++) {
	if (mask[i]) {
	    u = i / dset->pd;
	    if (u != ubak) {
		n++;
		ubak = u;
	    }
	}
    }

#if SUBDEBUG
    fprintf(stderr, "count_panel_units: got n = %d\n", n);
#endif

    return n;
}

/* construct mask for taking random sub-sample from dataset:
   we're selecting 'subn' cases without replacement, and
   there may or may not be an existing mask in place that
   has to be respected.
*/

static int make_random_mask (const char *s, const char *oldmask, 
			     DATASET *dset, char *mask)
{
    unsigned u;
    int targ, avail, oldn = dset->n;
    int i, subn, cases, rejn;
    int err = 0;

    /* how many cases are requested? */
    subn = smpl_get_int(s, dset, &err);

    if (subn <= 0 || subn >= dset->n) {
	err = 1;
    } else if (oldmask != NULL) {
	/* dataset is already sub-sampled */
	oldn = 0;
	for (i=0; i<dset->n; i++) {
	    if (oldmask[i]) {
		oldn++;
	    }
	}
	if (subn >= oldn) {
	    err = 1;
	}
    }	

    if (err && err != E_PARSE) {
	gretl_errmsg_sprintf(_("Invalid number of cases %d"), subn);
	return err;
    }	

    /* Which is smaller: the number of cases to be selected or the
       complement, the number to be discarded?  For the sake of
       efficiency we'll go for the smaller value.
    */
    rejn = oldn - subn;
    if (rejn < subn) {
	/* select rejn observations to discard */
	targ = rejn;
	avail = 1;
    } else {
	/* select subn observations to include */
	targ = subn;
	avail = 0;
    }

    for (i=0; i<dset->n; i++) {
	if (oldmask != NULL && oldmask[i] == 0) {
	    /* obs is already excluded, hence not selectable */
	    mask[i] = -1;
	} else {
	    mask[i] = avail;
	}
    }

    cases = 0;

    while (cases < targ) {
	u = gretl_rand_int_max(dset->n);
	if (mask[u] == avail) {
	    /* obs is available, not yet selected */
	    mask[u] = !avail;
	    cases++;
	}
    }

    if (oldmask != NULL) {
	/* undo the 'already excluded' coding above */
	for (i=0; i<dset->n; i++) {
	    if (mask[i] == -1) {
		mask[i] = 0;
	    }
	}
    }

    return err;
}

int backup_full_dataset (DATASET *dset)
{
#if SUBDEBUG
    int newfull = 0;
#endif

    if (fullset == NULL) {
	fullset = malloc(sizeof *fullset);
	if (fullset == NULL) {
	    return E_ALLOC;
	}
#if SUBDEBUG
	newfull = 1;
#endif
    }

    if (dset != NULL) {
	*fullset = *dset;
	peerset = dset;
    } 

#if SUBDEBUG
    fprintf(stderr, "backup_full_dataset: fullset = %p (%s)\n",
	    (void *) fullset, newfull ? "new" : "old");
#endif

    return 0;
}

static void destroy_full_dataset (DATASET *dset)
{
    if (fullset != NULL) {
	if (fullset->varname == dset->varname) {
	    fullset->varname = NULL;
	}
	if (fullset->varinfo == dset->varinfo) {
	    fullset->varinfo = NULL;
	}
	if (fullset->descrip == dset->descrip) {
	    fullset->descrip = NULL;
	}
	destroy_dataset(fullset);
	fullset = NULL;
    }

    dset->varname = NULL;
    dset->varinfo = NULL;
    dset->descrip = NULL;

    free_Z(dset);
    clear_datainfo(dset, CLEAR_SUBSAMPLE); 

    peerset = NULL;
}

int complex_subsampled (void)
{
    return (fullset != NULL && fullset->Z != NULL);
}

int get_full_length_n (void)
{
    return (fullset != NULL) ? fullset->n : 0;
}

int dataset_is_complex_subsampled (const DATASET *dset)
{
    return (fullset != NULL && fullset->Z != NULL &&
	    dset == peerset);
}

/* When sub-sampling on some boolean criterion, check to see if we can
   meet the criterion by simply adjusting the endpoints of the sample
   range: life will be simpler if that is so.
*/

static int mask_contiguous (const char *mask,
			    const DATASET *dset,
			    int *pt1, int *pt2)
{
    int t, t1 = 0, t2 = dset->n - 1;
    int contig = 1;

    for (t=0; t<dset->n; t++) {
	if (mask[t] == 0) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=dset->n - 1; t>=0; t--) {
	if (mask[t] == 0) {
	    t2--;
	} else {
	    break;
	}
    }

    for (t=t1; t<=t2; t++) {
	if (mask[t] == 0) {
	    /* there's a hole inside the range */
	    contig = 0;
	    break;
	}
    }

    if (contig && dataset_is_panel(dset)) {
	int n = t2 - t1 + 1;

	/* sample must leave a whole number of panel units; moreover,
	   to retain "panelness" this number must be greater than 1 
	*/
	if (t1 % dset->pd != 0 || n % dset->pd != 0) {
	    contig = 0;
	} else if (n == dset->pd) {
	    contig = 0;
	}
    }

    if (contig) {
	*pt1 = t1;
	*pt2 = t2;
    }

    return contig;
}

static void 
copy_data_to_subsample (DATASET *subset, const DATASET *dset,
			int maxv, const char *mask)
{
    int i, t, s;

#if SUBDEBUG
    fprintf(stderr, "copy_data_to_subsample: subset = %p, dset = %p\n",
	    (void *) subset, (void *) dset);
#endif

    /* copy data values */
    for (i=1; i<maxv; i++) {
	s = 0;
	for (t=0; t<dset->n; t++) {
	    if (mask == NULL) {
		subset->Z[i][s++] = dset->Z[i][t];
	    } else if (mask[t] == 1) {
		subset->Z[i][s++] = dset->Z[i][t];
	    } else if (mask[t] == 'p') {
		/* panel padding */
		subset->Z[i][s++] = NADBL;
	    }
	}
    }

    /* copy observation markers, if any */
    if (dset->markers && subset->markers) {
	s = 0;
	for (t=0; t<dset->n; t++) {
	    if (mask == NULL || mask[t] == 1 || mask[t] == 'p') {
		strcpy(subset->S[s++], dset->S[t]);
	    }
	}
    }

    /* copy panel time info? */
    if (dataset_is_panel(subset)) {
	if (subset->pd == dset->pd) {
	    subset->panel_pd = dset->panel_pd;
	    subset->panel_sd0 = dset->panel_sd0;
	}
    }

    if (subset->stobs[0] == '\0' || subset->endobs[0] == '\0') {
	/* impose simple default if not already handled */
	strcpy(subset->stobs, "1");
	sprintf(subset->endobs, "%d", subset->n);
    }
}

int get_restriction_mode (gretlopt opt)
{
    int mode = SUBSAMPLE_UNKNOWN;

    if (opt & (OPT_M | OPT_C)) {
	mode = SUBSAMPLE_DROP_MISSING;
    } else if (opt & OPT_R) {
	mode = SUBSAMPLE_BOOLEAN;
    } else if (opt & OPT_N) {
	mode = SUBSAMPLE_RANDOM;
    } else if (opt & OPT_O) {
	mode = SUBSAMPLE_USE_DUMMY;
    } else if (opt & OPT_A) {
	mode = SUBSAMPLE_DROP_EMPTY;
    } else if (opt & OPT_W) {
	mode = SUBSAMPLE_DROP_WKENDS;
    }

    return mode;
}

/* Copy list of panel periods from Ti to T0; return 1 if members of
   the list Ti are consecutive, 0 if they are not. 
*/

static int copy_periods_list (int *T0, const int *Ti)
{
    int j, ret = 1;

    for (j=0; j<=Ti[0] && ret; j++) {
	if (j > 1 && Ti[j] != Ti[j-1] + 1) {
	    ret = 0;
	} else {
	    T0[j] = Ti[j];
	}
    }

    return ret;
}

/* When sub-sampling panel data on some boolean criterion: see if the
   exclusion of certain rows leaves a balanced panel.  Note that the
   requirement is not simply that each unit has the same number of
   temporal observations -- they must have the _same_ observations,
   and in addition the observations must be temporally contiguous.  
*/

static int mask_leaves_balanced_panel (char *mask, const DATASET *dset)
{
    int T = dset->pd;
    int *T0, *Ti;
    int i, u, ubak = -1;
    int ret = 1;

    T0 = gretl_list_new(T);
    Ti = gretl_list_new(T);

    if (T0 == NULL || Ti == NULL) {
	free(T0);
	free(Ti);
	return 0;
    }

    T0[0] = Ti[0] = 0;

    for (i=0; i<dset->n && ret; i++) {
	if (mask[i]) {
	    u = i / T;
	    if (u != ubak) {
		if (Ti[0] > 0) {
		    if (T0[0] == 0) {
			/* we haven't made the reference list, T0, yet */
			ret = copy_periods_list(T0, Ti);
		    } else if (gretl_list_cmp(Ti, T0)) {
			/* the current list differs from the reference one */
			ret = 0;
		    }
		}
		Ti[0] = 1;
		Ti[1] = i % T;
	    } else {
		Ti[0] += 1;
		Ti[Ti[0]] = i % T;
	    }
	    ubak = u;
	}
    }

    free(T0);
    free(Ti);

    return ret;
}

static int 
make_panel_submask (char *mask, const DATASET *dset, int *err)
{
    int T = dset->pd;
    int N = dset->n / T;
    char *umask, *pmask;
    int i, np = 0;

    umask = calloc(N + T, 1);
    if (umask == NULL) {
	*err = E_ALLOC;
	return 0;
    }

    pmask = umask + N;
 
    for (i=0; i<dset->n; i++) {
	if (mask[i]) {
	    umask[i / T] = 1;
	    pmask[i % T] = 1;
	}
    }

    for (i=0; i<dset->n; i++) {
	if (!mask[i]) {
	    if (umask[i / T] && pmask[i % T]) {
		mask[i] = 'p'; /* mark as padding row */
		np++;
	    }
	}
    }

#if SUBDEBUG
    fprintf(stderr, "make_panel_submask: number of padding rows = %d\n",
	    np);
#endif

    free(umask);

    return np;
}

static int add_daily_date_strings (char *selected, 
				   const DATASET *dset,
				   DATASET *subset)
{
    int err, t, i = 0;
    
    err = dataset_allocate_obs_markers(subset);

    if (!err) {
	subset->markers = DAILY_DATE_STRINGS;
	for (t=0; t<dset->n; t++) {
	    if (selected[t]) {
		ntodate(subset->S[i++], t, dset);
	    }
	}
    }

    return err;
}

/* See if sub-sampling dated daily data leaves a useable daily dataset
   (as in dropping weekends and/or holidays).
*/

static int try_for_daily_subset (char *selected, 
				 const DATASET *dset,
				 DATASET *subset,
				 gretlopt *optp)
{
    char datestr[OBSLEN];
    long ed, ed0 = 0, edbak = 0;
    int t, wd, ngaps = 0;
    int delta, mon_delta;
    int delta_max = 0;
    int started = 0;
    int any_sat = 0;
    int newpd = 5;
    int err = 0;

    if (dset->pd > 5) {
	/* first pass: look at weekend status of included obs */
	for (t=0; t<dset->n; t++) {
	    if (selected[t]) {
		ntodate(datestr, t, dset);
		wd = weekday_from_date(datestr);
		if (wd == 0) {
		    /* got a Sunday: result must be 7-day */
		    newpd = 7;
		    break;
		} else if (wd == 6) {
		    any_sat = 1;
		}
	    }
	}
	if (newpd == 5 && any_sat) {
	    /* no Sundays, but got a Saurday: result
	       must be 6-day data */
	    newpd = 6;
	}
    }

#if SUBDEBUG
    fprintf(stderr, "try_for_daily_subset: newpd = %d\n", newpd);
#endif

    /* determine expected "delta days" for Mondays */
    mon_delta = (newpd == 7) ? 1 : (newpd == 6)? 2 : 3;

    /* second pass: count the calendar gaps */
    for (t=0; t<dset->n; t++) {
	if (selected[t]) {
	    ntodate(datestr, t, dset);
	    wd = weekday_from_date(datestr);
	    ed = get_epoch_day(datestr);
	    if (started) {
		delta = ed - edbak;
		ngaps += (wd == 1)? (delta > mon_delta) : (delta > 1);
		if (delta > delta_max) {
		    delta_max = delta;
		}
	    } else {
		ed0 = ed;
		started = 1;
	    }
	    edbak = ed;
	}
    }

    /* Now, we should probably restrict the proportion of
       missing days for this treatment: for exclusion of
       non-trading days 10 percent should be generous. 
       But we'll also require that the maximum "daily
       delta" be less than 10.
    */

    if (delta_max < 10) {
	double keepfrac = subset->n / (double) dset->n;
	double keepmin = 0.9 * newpd / (double) dset->pd;

#if SUBDEBUG
	fprintf(stderr, " daily: delta_max=%d, keepfrac=%.4f\n"
		" keepmin=%.4f, ngaps=%d, newpd=%d\n", delta_max, keepfrac,
		keepmin, ngaps, newpd);
#endif

	if (keepfrac >= keepmin) {
	    if (dset->S == NULL && ngaps > 0) {
		err = add_daily_date_strings(selected, dset, subset);
		/* flag not to destroy subset date strings */
		*optp |= OPT_P;
	    }
	    if (!err) {
		subset->structure = TIME_SERIES;
		subset->pd = newpd;
		subset->sd0 = (double) ed0;
		subset->t2 = subset->n - 1;
		ntodate(subset->stobs, 0, subset);
		ntodate(subset->endobs, subset->t2, subset);
	    }
	}
    }

    return err;
}

/* does the given policy lead to an empty sample, or no change
   in the sample, perchance? 
*/

static int check_subsample_n (int n, DATASET *dset,
			      char **pmask, PRN *prn)
{
    if (n == 0) {
	gretl_errmsg_set(_("No observations would be left!"));
	return E_DATA;
    } else if (n == dset->n) {
	/* not really an error, just a no-op */
	if (gretl_messages_on()) {
	    pputs(prn, _("No observations were dropped!"));
	    pputc(prn, '\n');
	}
	free(*pmask);
	*pmask = NULL;
    }

    return 0;
}

#define needs_string_arg(m) (m == SUBSAMPLE_RANDOM || \
	                     m == SUBSAMPLE_USE_DUMMY || \
                             m == SUBSAMPLE_BOOLEAN)

static int 
make_restriction_mask (int mode, const char *s,
		       const int *list, DATASET *dset,
		       const char *oldmask, char **pmask,
		       PRN *prn)
{
    char *mask = NULL;
    int sn = 0, err = 0;

    if (needs_string_arg(mode) && (s == NULL || *s == '\0')) {
	return E_ARGS;
    }

    mask = make_submask(dset->n);
    if (mask == NULL) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "make_restriction_mask: oldmask = %p\n", (void *) oldmask);
#endif

    /* construct subsample mask in one of several possible ways */

    if (mode == SUBSAMPLE_DROP_MISSING) { 
	err = make_missing_mask(list, dset, mask);
    } else if (mode == SUBSAMPLE_DROP_EMPTY) {
	err = make_empty_mask(list, dset, mask);
    } else if (mode == SUBSAMPLE_DROP_WKENDS) {
	err = make_weekday_mask(dset, mask);
    } else if (mode == SUBSAMPLE_RANDOM) {
	err = make_random_mask(s, oldmask, dset, mask);
    } else if (mode == SUBSAMPLE_USE_DUMMY) {
	err = mask_from_dummy(s, dset, mask);
    } else if (mode == SUBSAMPLE_BOOLEAN) {
	err = mask_from_temp_dummy(s, dset, mask, prn);
    } else {
	gretl_errmsg_set(_("Sub-sample command failed mysteriously"));
	err = 1;
    }

    if (err) {
	/* exit now on unrecoverable error */
	free(mask);
	return err;
    }

    /* cumulate sample restrictions, if appropriate */
    if (oldmask != NULL && mode != SUBSAMPLE_RANDOM) {
	sn = overlay_masks(mask, oldmask, dset->n);
    } else {
	sn = count_selected_cases(mask, dset);
    }

    err = check_subsample_n(sn, dset, &mask, prn);

    if (err) {
	free(mask);
    } else {
	*pmask = mask;
    }

    return err;
}

/* Here the mask we're after is "mixed" in the sense that
   we are starting from a t1, t2 setting but possibly
   compounding this with a prior "subsampling".
*/

static int 
make_mixed_mask (int t1, int t2, DATASET *dset,
		 const char *oldmask, char **pmask,
		 PRN *prn)
{
    char *mask = NULL;
    int t, sn = 0, err = 0;

    mask = make_submask(dset->n);
    if (mask == NULL) {
	return E_ALLOC;
    }

    for (t=t1; t<=t2; t++) {
	mask[t] = 1;
    }

    if (oldmask != NULL) {
	sn = overlay_masks(mask, oldmask, dset->n);
    } else {
	sn = t2 - t1 + 1;
    }

    err = check_subsample_n(sn, dset, &mask, prn);

    if (err) {
	free(mask);
    } else {
	*pmask = mask;
    }

    return err;
}

static void finalize_panel_subset (DATASET *subset,
				   DATASET *dset,
				   int npad)
{
    int pdp = subset->pd;
    int den = 10.0;

    while ((pdp = pdp / 10)) {
	den *= 10;
    }

    subset->sd0 = 1.0 + 1.0 / den;
    ntodate(subset->stobs, 0, subset); 
    ntodate(subset->endobs, subset->n - 1, subset);

    if (dset->pangrps != NULL && npad == 0) {
	/* carry over panel group names from full dataset */
	subset->pangrps = gretl_strdup(dset->pangrps);
    }
}    

/* This is also used elsewhere: 

   in gui2/session.c, when re-establishing a previously sub-sampled data 
   state on reopening a saved session

   in gretl_func.c, on exit from a user function when the dataset was 
   sub-sampled on entry to the function (and we need to re-establish the
   original sub-sample)
*/

int 
restrict_sample_from_mask (char *mask, DATASET *dset, gretlopt opt)
{
    DATASET *subset;
    gretlopt zopt = OPT_R;
    int err = 0;

    if (dset->auxiliary) {
	fprintf(stderr, "restrict_sample_from_mask: attempting to restrict "
		"an auxiliary data\n");
	return E_DATA;
    }

    if (mask == RESAMPLED) {
	fprintf(stderr, "restrict_sample_from_mask: got RESAMPLED!\n");
	return E_DATA;
    }

    if ((opt & OPT_B) && !dataset_is_panel(dset)) {
	gretl_errmsg_set(_("--balanced option is invalid: the (full) "
			   "dataset is not a panel"));
	return E_BADOPT;
    }

    subset = datainfo_new();
    if (subset == NULL) {
	return E_ALLOC;
    }

    subset->n = count_selected_cases(mask, dset);
    subset->v = dset->v;

#if SUBDEBUG
    fprintf(stderr, "restrict_sample_from_mask:\n new subset=%p, "
	    "%d obs in subsample vs %d in full dataset\n",
	    (void *) subset, subset->n, dset->n);
#endif

    if (dataset_is_panel(dset)) {
	/* are we able to reconstitute a panel? */
	int nunits = count_panel_units(mask, dset);
	int ok = 0, npad = 0;

	if (nunits > 1 && subset->n > nunits) {
	    /* there's some possibility of doing so */
	    if (opt & OPT_B) {
		/* add padding rows? only if this was requested */
		npad = make_panel_submask(mask, dset, &err);
		if (err) {
		    free(subset);
		    return err;
		}
		ok = 1;
	    } else {		
		ok = mask_leaves_balanced_panel(mask, dset);
	    } 
	    if (ok) {
		subset->structure = STACKED_TIME_SERIES;
		subset->n += npad;
		subset->pd = subset->n / nunits;
		finalize_panel_subset(subset, dset, npad);
	    }
	} else if (nunits == 1 && subset->n == dset->pd) {
	    /* time series for single panel unit */
	    subset->structure = SPECIAL_TIME_SERIES;
	}
    } else if (dated_daily_data(dset)) {
	/* see if we can preserve daily time series */
	err = try_for_daily_subset(mask, dset, subset, &zopt);
    }

    if (!err) {
	/* set up the sub-sampled dataset */
	err = start_new_Z(subset, zopt);
    }

    if (err) { 
	free(subset);
	return err;
    }

#if SUBDEBUG
    fprintf(stderr, "started new Z for subset (v=%d, n=%d, Z=%p)\n", 
	    subset->v, subset->n, (void *) subset->Z);
#endif

    /* link (don't copy) varnames and descriptions, since these are
       not dependent on the series length */
    subset->varname = dset->varname;
    subset->varinfo = dset->varinfo;
    subset->descrip = dset->descrip;

    /* set up case markers? */
    if (dset->markers) {
	err = dataset_allocate_obs_markers(subset);
	if (err) {
	    free_Z(subset);
	    free(subset);
	    return E_ALLOC;
	}
    }

    /* copy across data (and case markers, if any) */
    copy_data_to_subsample(subset, dset, dset->v, mask);

    if (opt & OPT_T) {
	/* --permanent */
	check_models_for_subsample(mask, NULL);
	destroy_full_dataset(dset);
    } else {
	err = backup_full_dataset(dset);
	subset->submask = copy_subsample_mask(mask, &err);
    }

    /* switch pointers */
    *dset = *subset;
    free(subset);

    return err;
}

static char *expand_mask (char *tmpmask, const char *oldmask,
			  int *err)
{
    char *newmask = make_submask(fullset->n);

    if (newmask == NULL) {
	*err = E_ALLOC;
    } else {
	/* map @tmpmask onto the full data range */
	int t, i = 0;
	
	for (t=0; t<fullset->n; t++) {
	    if (oldmask[t]) {
		newmask[t] = tmpmask[i++];
	    }
	}
    }

    return newmask;
}

/* Below: we do this "precompute" thing if the dataset is already
   subsampled and the user wants to compound the restriction with a
   new one of the form "obs=x" or "obs!=x".  The "obs" references may
   get out of whack if we restore the full dataset first, as we
   usually do.  For example, say the spec is "obs!=50" to exclude
   observation 50: presumably the user means to exclude the 50th
   observation in the current, subsampled dataset, which may not be
   the same as the 50th observation in the full dataset.
*/

static char *precompute_mask (const char *s, const char *oldmask,
			      DATASET *dset, PRN *prn, int *err)
{
    char *tmpmask = make_submask(dset->n);
    char *newmask = NULL;

#if SUBDEBUG
    fprintf(stderr, "restrict_sample: precomputing new mask\n");
#endif

    if (tmpmask == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	/* fill out mask relative to current, restricted dataset */
	*err = mask_from_temp_dummy(s, dset, tmpmask, prn);
    }

    if (!*err) {
	if (fullset != NULL) {
	    newmask = expand_mask(tmpmask, oldmask, err);
	} else {
	    newmask = tmpmask;
	    tmpmask = NULL;
	}
    }

    free(tmpmask);

    return newmask;
}

/* Intended for time series data: trim any missing values
   at the start and end of the current sample range, then
   check the remaining range for missing values and flag
   an error if any are found. If @rt1 and @rt2 are non-
   NULL we record the resulting sample range in them,
   otherwise we assign in to the t1 and t2 members of @dset.
*/

static int set_contiguous_sample (const int *list,
				  DATASET *dset,
				  int *rt1, int *rt2)
{
    int ct1 = dset->t1;
    int ct2 = dset->t2;
    int err = 0;

    if (list != NULL && list[0] > 0) {
	err = list_adjust_sample(list, &ct1, &ct2, dset, NULL);
    } else {
	int *biglist = NULL;
	int nvars = 0;

	biglist = full_var_list(dset, &nvars);
	if (nvars == 0) {
	    ; /* no-op */
	} else if (biglist == NULL) {
	    err = E_ALLOC;
	} else {
	    err = list_adjust_sample(biglist, &ct1, &ct2, 
				     dset, NULL);
	    free(biglist);
	}
    }

    if (!err) {
	if (rt1 != NULL && rt2 != NULL) {
	    *rt1 = ct1;
	    *rt2 = ct2;
	} else {
	    dset->t1 = ct1;
	    dset->t2 = ct2;
	}
    }

    return err;
}

static void destroy_restriction_string (DATASET *dset)
{
    if (dset->restriction != NULL) {
	free(dset->restriction);
	dset->restriction = NULL;
    }
}

/* Make a string representing the sample restriction. This is
   for reporting purposes. Note that in some cases the
   incoming @restr string may be NULL, for instance if the
   no-missing option is chosen via the gretl GUI.
*/

static int make_restriction_string (DATASET *dset, char *old, 
				    const char *restr, int mode)
{
    char *s = NULL;
    int n = 0, err = 0;

    destroy_restriction_string(dset);

    if (old != NULL) {
	n = strlen(old);
    }

    if (mode == SUBSAMPLE_RANDOM) {
	n += strlen("random");
    } else if (mode == SUBSAMPLE_DROP_MISSING) {
	n += strlen("no-missing");
    } else if (restr != NULL) {
	n += strlen(restr);
    }

    if (n > 0) {
	s = malloc(n + 5);
	if (s == NULL) {
	    err = E_ALLOC;
	} 
    }

    if (s != NULL) {
	*s = '\0';
	if (old != NULL) {
	    strcpy(s, old);
	    strcat(s, " && ");
	}
	if (mode == SUBSAMPLE_RANDOM) {
	    strcat(s, "random");
	} else if (mode == SUBSAMPLE_DROP_MISSING) {
	    strcat(s, "no-missing");
	} else if (restr != NULL) {
	    strcat(s, restr);
	}
	dset->restriction = s;
    }

    return err;
}

static int full_sample (const DATASET *dset)
{
    if (complex_subsampled()) {
	return 0;
    } else if (dset->t1 > 0 || dset->t2 < dset->n - 1) {
	return 0;
    } else {
	return 1;
    }
}

static int handle_ts_restrict (char *mask, DATASET *dset,
			       gretlopt opt, int t1)
{
    char stobs[OBSLEN];
    int pd = dset->pd;
    double sd0;
    int err;

    ntodate(stobs, t1, dset);
    sd0 = get_date_x(dset->pd, stobs);

    err = restrict_sample_from_mask(mask, dset, opt);

    if (!err) {
	/* re-establish time-series characteristics */
	dset->structure = TIME_SERIES;
	dset->pd = pd;
	dset->sd0 = sd0;
	strcpy(dset->stobs, stobs);
	ntodate(dset->endobs, dset->n - 1, dset);
    }

    return err;
}

/* under what conditions will be bother testing for a
   contiguous subsample mask? */

static int try_for_contig (gretlopt opt, DATASET *dset)
{
    if (opt & (OPT_T | OPT_N)) {
	/* permanent restriction, or random subsample: no */
	return 0;
    }

    if (opt & OPT_C) {
	/* we're already forcing contiguity */
	return 0;
    }

    if ((opt & OPT_Z) && !dataset_is_time_series(dset)) {
	/* forcing a resize and not time-series: no */
	return 0;
    }

    return 1;
}

static int handle_resize_option (gretlopt *opt,
				 const char *param,
				 DATASET *dset,
				 int *rt1, int *rt2)
{
    int err = 0;

    if (*opt & (OPT_O | OPT_M | OPT_A | OPT_N | OPT_C)) {
	/* besides --resize we got another option
	   that implies a restriction: so there's
	   nothing to be done here
	*/
	return 0;
    }

    /* add implicit "restrict" flag */
    *opt |= OPT_R;

    /* and inspect @param */
    if (param != NULL) {
	if (strchr(param, '=') ||
	    strchr(param, '<') ||
	    strchr(param, '>')) {
	    ; /* must be a restriction spec */
	} else {
	    /* may be obs1 [ obs2 ] or similar? */
	    err = test_set_sample(param, dset, rt1, rt2);
	}
    }

    return err;
}

static int check_restrict_options (gretlopt opt, const char *param)
{
    /* We'll accept the redundant combination of the options --dummy
       (OPT_O) and --restrict (OPT_R), but other than that the options
       --dummy, 
       --restrict, 
       --no-missing (OPT_M), 
       --no-all-missing (OPT_A) and
       --random (OPT_N) 
       are all mutually incompatible.
    */
    if (incompatible_options(opt, OPT_O | OPT_M | OPT_A | OPT_N) ||
	incompatible_options(opt, OPT_R | OPT_M | OPT_A | OPT_N)) {
	return E_BADOPT;
    }

    /* The --contiguous option (OPT_C) is compatible only with 
       --no-missing (implied) and possibly --resize (OPT_Z) */
    if (opt & OPT_C) {
	if ((opt & ~(OPT_C | OPT_M | OPT_Z)) != OPT_NONE) {
	    return E_BADOPT;
	}
    }

    if (opt & (OPT_O | OPT_R | OPT_N)) {
	/* parameter is required */
	if (param == NULL || *param == '\0') {
	    return E_ARGS;
	}
    }

    if (opt & OPT_U) {
	/* --current */
	if (!(opt & (OPT_T | OPT_Z))) {
	    /* requires --permanent or --resize */
	    return E_ARGS;
	}
    }
	
    return 0;
}

static int check_permanent_option (gretlopt opt,
				   DATASET *dset,
				   int *n_models)
{
    if (gretl_function_depth() > 0) {
	/* can't permanently shrink the dataset within a function */
	gretl_errmsg_set(_("The dataset cannot be modified at present"));
	return E_DATA;
    } else if (!(opt & (OPT_O | OPT_M | OPT_A | OPT_N | OPT_R | OPT_U))) {
	/* we need some kind of restriction flag */
	return E_ARGS;
    } else if (gretl_in_gui_mode() && !(opt & OPT_F)) {
	/* GUI program, without internal "force" option */
	*n_models = n_stacked_models();
	if (*n_models > 0 && !full_sample(dset)) {
	    /* too difficult to recover gracefully on error */
	    gretl_errmsg_set(_("The full dataset must be restored before "
			       "imposing a permanent sample restriction"));
	    return E_DATA;
	}
    }

    return 0;
}

/* "Precomputing" a mask means computing a mask based on a
   current subsampled dataset (before restoring the full dataset,
   which is a part of the subsampling process).
*/

#if 1 /* modified 2015-12-06 */

/* Do this if we're cumulating a boolean restriction on top of
   an existing restriction.
*/

static int do_precompute (int mode, char *oldmask, const char *param)
{
    return oldmask != NULL && mode == SUBSAMPLE_BOOLEAN;
}

#else

/* Do this only if we're cumulating (rather than replacing) sample
   restrictions _and_ the present resampling parameter makes reference
   to the built-in "obs" series. Given the last point, the user might
   reasonably expect the "obs" values to refer to the currently
   subsampled dataset rather than the full dataset.
*/

static int do_precompute (int mode, char *oldmask, const char *param)
{
    int ret = 0;
    
    if (fullset == NULL) {
	; /* not subsampled: no */
    } else if (oldmask == NULL) {
	; /* not cumulating restrictions: no */
    } else if (mode != SUBSAMPLE_BOOLEAN) {
	; /* can't be keying off "obs" series: no */
    } else {
	/* look for "obs" as a "word" in itself */
	const char *s = param;

	while (!ret && strstr(s, "obs") != NULL) {
	    if (gretl_namechar_spn(s) == 3 &&
		(s == param || !isalpha(*(s-1)))) {
		ret = 1;
	    }
	    s += 3;
	}
    }

    return ret;
}

#endif /* new versus old */

/* restrict_sample: 
 * @param: restriction string (or %NULL).
 * @list: list of variables in case of OPT_M (or %NULL).     
 * @dset: dataset struct.
 * @state: structure representing program state (or %NULL).
 * @opt: option flags.
 * @prn: printing apparatus.
 * @n_dropped: location to receive count of dropped
 * observations, or NULL.
 *
 * Sub-sample the data set, based on the criterion of skipping all
 * observations with missing data values (OPT_M); or using as mask a
 * specified dummy variable (OPT_O); or masking with a specified
 * boolean condition (OPT_R); or selecting at random (OPT_N).
 *
 * In case OPT_M or OPT_C a @list of variables may be supplied; in 
 * cases OPT_O, OPT_R and OPT_N, @param must contain specifics.
 *
 * In case OPT_P is included, the restriction will rePlace any
 * existing sample restriction, otherwise the resulting restriction
 * will be the logical product of the new restriction and any
 * existing restriction.
 *
 * In case the original dataset was a panel and OPT_B was given,
 * we'll pad with missing values if necessary, to try to reconstitute 
 * a balanced panel.
 *
 * In case OPT_T is included, the sample restriction will be
 * permanent (the full datasset is destroyed), otherwise the
 * restriction can be undone via the command "smpl full".
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restrict_sample (const char *param, const int *list,
		     DATASET *dset, ExecState *state, 
		     gretlopt opt, PRN *prn,
		     int *n_dropped)
{
    char *oldrestr = NULL;
    char *oldmask = NULL;
    char *mask = NULL;
    int free_oldmask = 0;
    int force_resize = 0;
    int permanent = 0;
    int n_models = 0;
    int rt1 = -1, rt2 = -1;
    int n_orig = 0;
    int mode, err = 0;
    
    if (dset == NULL || dset->Z == NULL) {
	return E_NODATA;
    }

    n_orig = dset->n;

    /* general check on incoming options */
    err = check_restrict_options(opt, param);

    /* now take a closer look at some finicky options */

    if (!err && (opt & OPT_Z)) {
	err = handle_resize_option(&opt, param, dset, &rt1, &rt2);
	if (!err) {
	    force_resize = 1;
	}
    }

    if (!err && (opt & OPT_T)) {
	err = check_permanent_option(opt, dset, &n_models);
	if (!err) {
	    permanent = 1;
	}
    }

    if (err) {
	return err;
    }    

    gretl_error_clear();

#if FULLDEBUG || SUBDEBUG
    fprintf(stderr, "\nrestrict_sample: param='%s'\n", param);
    fprintf(stderr, " dset=%p, state=%p, fullset=%p\n", (void *) dset, 
	    (void *) state, (void *) fullset);
#endif

    if (opt & OPT_C) {
	/* --contiguous */
	if (force_resize) {
	    err = set_contiguous_sample(list, dset, &rt1, &rt2);
	    if (err) {
		return err;
	    }
	} else {
	    /* we should be done at this point */
	    return set_contiguous_sample(list, dset, NULL, NULL);
	}
    }	

    mode = get_restriction_mode(opt);
    if (mode == SUBSAMPLE_UNKNOWN) {
	gretl_errmsg_set("Unrecognized sample command");
	return 1;
    } 

    if (!(opt & OPT_P)) {
	/* not replacing but cumulating any existing restrictions */
	oldmask = make_current_sample_mask(dset, &err);
	if (err) {
	    return err;
	}
	free_oldmask = 1;
	if (dset->restriction != NULL) {
	    oldrestr = gretl_strdup(dset->restriction);
	}
    } else if (state != NULL && state->submask != NULL) {
	/* subsampling within a function: this necessarily
	   cumulates with the incoming sample restriction 
	   recorded in state->submask
	*/
	oldmask = state->submask;
    }

    if (do_precompute(mode, oldmask, param)) {
	/* we come in here only if cumulating restrictions */
	mask = precompute_mask(param, oldmask, dset, prn, &err);
    }

    if (!err && !full_sample(dset)) {
	/* restore the full data range, for housekeeping purposes */
	err = restore_full_sample(dset, NULL);
    }

    if (err) {
	/* clean up and get out */
	free(mask);
	if (free_oldmask) {
	    free(oldmask);
	}
	return err;
    }

    if (mask == NULL) {
	/* not already handled by "precompute" above */
	if (rt1 >= 0 && rt2 >= 0) {
	    err = make_mixed_mask(rt1, rt2, dset, oldmask,
				  &mask, prn);
	} else {
	    err = make_restriction_mask(mode, param, list, dset, 
					oldmask, &mask, prn);
	}
    }

    if (!err && n_models > 0) {
	err = check_models_for_subsample(mask, n_dropped);
    }    

    if (!err && mask != NULL) {
	int contiguous = 0;
	int t1 = 0, t2 = 0;

	if (try_for_contig(opt, dset)) {
	    contiguous = mask_contiguous(mask, dset, &t1, &t2);
	} else if (dataset_is_time_series(dset) && (opt & OPT_C)) {
	    contiguous = 1;
	    t1 = rt1;
	    t2 = rt2;
	}

#if SUBDEBUG
	fprintf(stderr, "restrict sample: contiguous range? %s\n",
		contiguous ? "yes" : "no");
#endif

	if (contiguous) {
	    if (force_resize && dataset_is_time_series(dset)) {
		/* apply the restriction, but then re-establish the
		   time-series character of the dataset
		*/
		err = handle_ts_restrict(mask, dset, opt, t1);
	    } else {
		/* just move the sample range pointers, avoiding
		   the overhead of creating a parallel dataset
		*/
		dset->t1 = t1;
		dset->t2 = t2;
	    }
	} else {
	    err = restrict_sample_from_mask(mask, dset, opt);
	}
    }

    free(mask);

    if (free_oldmask) {
	free(oldmask);
    }

    if (!err) {
	if (permanent) {
	    destroy_restriction_string(dset);
	} else {
	    make_restriction_string(dset, oldrestr, param, mode);
	}
	if (n_dropped != NULL) {
	    *n_dropped = n_orig - sample_size(dset);
	}	
    }

    free(oldrestr);

#if SUBDEBUG
    fprintf(stderr, "restrict sample: dset: t1=%d, t2=%d, n=%d\n",
	    dset->t1, dset->t2, dset->n);
#endif    

    return err;
}

/* perma_sample: 
 * @dset: dataset struct.
 * @opt: option flags: must be OPT_T | OPT_U.
 * @prn: printing apparatus.
 * @n_dropped: location to receive count of dropped models,
 * or NULL.
 *
 * Make the current sub-sampling of the dataset permanent.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int perma_sample (DATASET *dset, gretlopt opt, PRN *prn,
		  int *n_dropped)
{
    gretlopt testopt = OPT_T | OPT_U;
    
    if (dset->submask == NULL) {
	pputs(prn, "smpl: nothing to be done\n");
	return 0;
    } else if (dset->submask == RESAMPLED) {
	pputs(prn, "smpl: dataset is resampled\n");
	return E_DATA;
    } else if (opt != testopt) {
	return E_BADOPT;
    }

    if (n_dropped != NULL) {
	int err;

	err = check_models_for_subsample(dset->submask, n_dropped);
	if (err) {
	    return err;
	}
    } else {
	check_models_for_subsample(dset->submask, NULL);
    }

    free(dset->submask);
    dset->submask = NULL;
    free(dset->restriction);
    dset->restriction = NULL;

    if (fullset->varname == dset->varname) {
	fullset->varname = NULL;
    }
    if (fullset->varinfo == dset->varinfo) {
	fullset->varinfo = NULL;
    }
    if (fullset->descrip == dset->descrip) {
	fullset->descrip = NULL;
    }
    
    destroy_dataset(fullset);
    fullset = peerset = NULL;

    return 0;
}

enum {
    SMPL_T1,
    SMPL_T2
};

static int panel_round (const DATASET *dset, int t, int code)
{
    if (code == SMPL_T1) {
	while ((t + 1) % dset->pd != 1) {
	    t++;
	}
    } else {
	while ((t + 1) % dset->pd != 0) {
	    t--;
	}
    }

    return t;
}

static int smpl_get_int (const char *s, DATASET *dset, int *err)
{
    int k = -1;

    if (integer_string(s)) {
	k = atoi(s);
    } else {
	double x;

	if (gretl_is_scalar(s)) {
	    x = gretl_scalar_get_value(s, NULL);
	} else {
	    x = generate_scalar(s, dset, err);
	}

	if (!xna(x) && x < INT_MAX) {
	    k = (int) x;
	}
    }

    return k;
}

static int get_sample_limit (const char *s, DATASET *dset, int code)
{
    int ret = -1;
    int err = 0;

#if SUBDEBUG
    fprintf(stderr, "get_sample_limit: s = '%s'\n", s);
#endif

    if (*s == '-' || *s == '+') {
	/* increment/decrement form */
	int incr = smpl_get_int(s + 1, dset, &err);

	if (!err) {
	    if (*s == '-') {
		incr = -incr;
	    }
	    if (dataset_is_panel(dset)) {
		incr *= dset->pd;
	    }
	    if (code == SMPL_T1) {
		ret = dset->t1 + incr;
	    } else {
		ret = dset->t2 + incr;
	    }
	}
    } else {
	/* absolute form */
	ret = get_t_from_obs_string(s, dset);
	if (ret < 0) {
	    gretl_error_clear();
	    ret = smpl_get_int(s, dset, &err);
	    /* convert to base 0 */
	    if (!err) ret--;
	}
	if (ret >= 0 && dataset_is_panel(dset)) {
	    ret = panel_round(dset, ret, code);
	}
    }

    return ret;
}

/* Catch the case where we're in a function and attempting to
   move t1 or t2 out of the range established on entry to the
   function: we don't want to carry forward any error message
   that implies t was out of the full data range.
*/

static void maybe_clear_range_error (int t, DATASET *dset)
{
    if (t >= 0 && t < dset->n) {
	gretl_error_clear();
    }
}

static int real_set_sample (const char *start,
			    const char *stop,
			    DATASET *dset,
			    int *t1, int *t2)
{
    int nf, new_t1 = dset->t1, new_t2 = dset->t2;
    int tmin = 0, tmax = 0;
    int testing = 0;

    if (dset == NULL) {
	return E_NODATA;
    }

    if (t1 != NULL && t2 != NULL) {
	testing = 1;
    }

    gretl_error_clear();

    nf = (start != NULL) + (stop != NULL);

#if SUBDEBUG
    fprintf(stderr, "set_sample: start='%s', stop='%s', dset=%p\n", 
	    start, stop, (void *) dset);
    if (dset != NULL) {
	fprintf(stderr, "dset->v = %d, dset->n = %d, pd = %d\n",
		dset->v, dset->n, dset->pd);
    }
#endif

    if (nf == 2 && dset->n == 0 && !testing) {
	/* database special */
	return db_set_sample(start, stop, dset);
    }

    if (nf == 0) {
	/* no-op, just print the current sample */
	return 0;
    }

    sample_range_get_extrema(dset, &tmin, &tmax);

#if SUBDEBUG
    fprintf(stderr, "sample extrema: lo = %d, hi = %d\n", tmin, tmax);
#endif
	
    if (nf == 1) {
	/* implicitly just setting the starting observation? */
	if (start == NULL) {
	    return E_ARGS;
	}
	new_t1 = get_sample_limit(start, dset, SMPL_T1);
	if (new_t1 < tmin || new_t1 > tmax) {
	    maybe_clear_range_error(new_t1, dset);
	    gretl_errmsg_set(_("error in new starting obs"));
	    return 1;
	}
	if (testing) {
	    *t1 = new_t1;
	} else {
	    dset->t1 = new_t1;
	}
	return 0;
    }

    /* now we're looking at the 2 fields case */

    if (strcmp(start, ";")) {
	new_t1 = get_sample_limit(start, dset, SMPL_T1);
	if (new_t1 < tmin || new_t1 > tmax) {
	    maybe_clear_range_error(new_t1, dset);
	    gretl_errmsg_set(_("error in new starting obs"));
	    return 1;
	}	
    }

    if (strcmp(stop, ";")) {
	new_t2 = get_sample_limit(stop, dset, SMPL_T2);
	if (new_t2 < tmin || new_t2 > tmax) {
	    maybe_clear_range_error(new_t2, dset);
	    gretl_errmsg_set(_("error in new ending obs"));
	    return 1;
	}
    }

    if (new_t1 < tmin || new_t1 > new_t2) {
	gretl_error_clear();
	gretl_errmsg_set(_("Invalid null sample"));
	return 1;
    }

    if (testing) {
	*t1 = new_t1;
	*t2 = new_t2;
    } else {
	dset->t1 = new_t1;
	dset->t2 = new_t2;
    }

    return 0;
}

int set_sample (const char *start, const char *stop, DATASET *dset)
{
    return real_set_sample(start, stop, dset, NULL, NULL);
}

static int test_set_sample (const char *s, DATASET *dset,
			    int *t1, int *t2)
{
    char start[OBSLEN], stop[OBSLEN];
    int err = 0;

    if (sscanf(s, "%15s %15s", start, stop) == 2) {
	err = real_set_sample(start, stop, dset, t1, t2);
    } else if (sscanf(s, "%15s", start) == 1) {
	err = real_set_sample(start, NULL, dset, t1, t2);
    }

    return err;
}

/**
 * count_missing_values:
 * @dset: dataset struct.
 * @opt: use %OPT_V for verbose operation, %OPT_A to 
 * examine all data.
 * @prn: printing struct.
 * @err: location to receive error code.
 *
 * Prints a count of missing values (if any) in the current
 * dataset over the currently defined sample range (or the 
 * entire data range if %OPT_A is given). If %OPT_V is given 
 * this includes a count of missing values at each observation; 
 * otherwise it just includes global and per-variable counts.
 *
 * Returns: 0 if no missing values are found (or on error),
 * otherwise the total number of missing values.
 */

int count_missing_values (const DATASET *dset, gretlopt opt, 
			  PRN *prn, int *err)
{
    int missval = 0, missobs = 0, totvals = 0, oldmiss = 0;
    int T, t1, t2;
    int *missvec;
    double missfrac;
    int i, t, tmiss;

    if (opt & OPT_A) {
	t1 = 0;
	t2 = dset->n - 1;
    } else {
	t1 = dset->t1;
	t2 = dset->t2;
    }

    T = t2 - t1 + 1;

    missvec = malloc(dset->v * sizeof missvec);

    if (missvec == NULL) {
	*err = E_ALLOC;
	return 0;
    }

    for (i=0; i<dset->v; i++) {
	missvec[i] = 0;
    }

    for (t=t1; t<=t2; t++) {
	tmiss = 0;
	for (i=1; i<dset->v; i++) {
	    if (series_is_hidden(dset, i)) {
		continue;
	    }
	    if (na(dset->Z[i][t])) {
		if (missvec[i] == 0) {
		    missvec[0] += 1;
		}
		missvec[i] += 1;
		missval++;
	    }
	    totvals++;
	}

	tmiss = missval - oldmiss;

	if (tmiss) {
	    missobs++;

	    if (opt & OPT_V) {
		/* verbose print by observation */
		if (dset->markers) { 
		    pprintf(prn, "%8s %4d %s\n", dset->S[t], tmiss,
			    _("missing values"));
		} else {
		    char tmp[OBSLEN];

		    ntodate(tmp, t, dset);
		    pprintf(prn, "%8s %4d %s\n", tmp, tmiss,
			    _("missing values"));
		}
	    }
	}
	oldmiss = missval;
    }

    missfrac = 100.0 * (double) missobs / T;

    pprintf(prn, _("\nNumber of observations (rows) with missing data "
		   "values = %d (%.2f%%)\n"), missobs, missfrac);

    pprintf(prn, _("Total number of missing data values = %d (%.2f%% "
	    "of total data values)\n"), missval, 
	    (100.0 * (double) missval / totvals));

    if (missvec[0] > 0) {
	pputc(prn, '\n');
	for (i=1; i<dset->v; i++) {
	    if (missvec[i] > 0) {
		missfrac = 100.0 * (double) missvec[i] / T;
		pprintf(prn, "%8s: %d %s (%.2f%%); %d %s (%.2f%%)\n", 
			dset->varname[i], missvec[i], _("missing values"),
			missfrac, T - missvec[i], _("valid values"),
			100.0 - missfrac);
	    }
	}
    }

    free(missvec);

    return missval;
}

static void copy_series_info (DATASET *dest, const DATASET *src,
			      int maxv)
{
    int i;

    for (i=1; i<maxv; i++) {
	strcpy(dest->varname[i], src->varname[i]);
	if (src->varinfo != NULL) {
	    copy_varinfo(dest->varinfo[i], src->varinfo[i]);
	}
    }
}

/* Below: For a model that was estimated using a data sample which
   differs from the current sample, reconstitute the dataset
   appropriate to the model.  If the model used a subsample, we use
   the special dummy variable (submask) saved with the model to select
   observations from the full dataset.  Another possibility is that
   the model used the full dataset, which has by now been subsampled.

   We attach this dataset to the model, so that it can be used for,
   e.g., residual or fitted versus actual plots.  The attached data
   set will be destroyed when the model itself is destroyed.

   By default we create a dataset containing series up to the
   highest numbered series associated with the model (regressand
   and regressors).  If OPT_F is given we include all series; if
   OPT_G is given we include only the constant.
*/

int add_dataset_to_model (MODEL *pmod, const DATASET *dset,
			  gretlopt opt)
{
    const DATASET *srcset;
    char *mask = NULL;
    int maxv, sn = 0;

    if (pmod->dataset != NULL) {
	/* FIXME? */
	destroy_dataset(pmod->dataset);
	pmod->dataset = NULL;
    }

#if SUBDEBUG
    fprintf(stderr, "add_dataset_to_model: fullset=%p, pmod->submask=%p\n",
	    (void *) fullset, (void *) pmod->submask);
#endif

    if (fullset != NULL) {
	sync_datainfo_members(dset);
	srcset = fullset;
    } else {
	srcset = dset;
    }

    if (pmod->submask == NULL) {
	/* no subsample info: pmod was estimated on the full dataset,
	   so we'll reconstruct the full dataset */
	sn = srcset->n;
    } else {
	/* pmod was estimated on a subsample, which has to 
	   be reconstructed */
	int t;

	mask = calloc(srcset->n, 1);
	if (mask == NULL) {
	    return E_ALLOC;
	}
	for (t=0; t<srcset->n; t++) {
	    if (pmod->submask[t] > 0) {
		mask[t] = 1;
		sn++;
	    }
	}
	if (sn == 0) {
	    free(mask);
	    return 1;
	}
    }

    if (opt & OPT_F) {
	maxv = srcset->v;
    } else if (opt & OPT_G) {
	maxv = 1;
    } else {
	maxv = highest_numbered_var_in_model(pmod, dset) + 1;
    }

    /* allocate auxiliary dataset */
    pmod->dataset = create_auxiliary_dataset(maxv, sn, 0);
    if (pmod->dataset == NULL) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "pmod->dataset allocated at %p\n", 
	    (void *) pmod->dataset);
#endif

    /* copy across info on series */
    copy_series_info(pmod->dataset, srcset, maxv);

    /* copy across data */
    copy_data_to_subsample(pmod->dataset, srcset, maxv, mask);

    /* dataset characteristics such as pd: if we're rebuilding the
       full dataset copy these across; but if we're reconstructing a
       subsampled dataset these features from the full dataset will be
       wrong, and we stay with the simple defaults
    */
    if (pmod->submask == NULL) {
	pmod->dataset->pd = srcset->pd;
	pmod->dataset->sd0 = srcset->sd0;
	strcpy(pmod->dataset->stobs, srcset->stobs);
	strcpy(pmod->dataset->endobs, srcset->endobs);
	pmod->dataset->structure = srcset->structure;
    }

    free(mask);

#if SUBDEBUG
    fputs("add_subsampled_dataset_to_model: success\n", stderr);
#endif

    return 0;
}

static int submasks_match (const DATASET *dset, const MODEL *pmod)
{
    char *s1 = dset->submask;
    char *s2 = pmod->submask;
    
    if (s1 == RESAMPLED && s2 == RESAMPLED) {
	return dset->n == pmod->full_n &&
	    dset->rseed == pmod->smpl.rseed;
    } else if (s1 == RESAMPLED || s2 == RESAMPLED) {
	return 0;
    } else {
	return submask_cmp(s1, s2) == 0;
    }
}

/* check the subsample mask from a model against datainfo to 
   see if it may have been estimated on a different
   (subsampled) data set from the current one
*/

int model_sample_problem (const MODEL *pmod, const DATASET *dset)
{
    int ret = 1;

    if (pmod->submask == NULL) {
	/* the model has no sub-sampling info recorded */
	if (dset->submask == NULL) {
	    /* data set is not sub-sampled either, OK */
	    ret = 0;
	} else {
	    fputs(I_("dataset is subsampled, model is not\n"), stderr);
	    gretl_errmsg_set(_("dataset is subsampled, model is not\n"));
	    ret = 1;
	}
    } else {
	/* the model does have sub-sampling info recorded */
	if (dset->submask == NULL) {
	    fputs(I_("model is subsampled, dataset is not\n"), stderr);
	    gretl_errmsg_set(_("model is subsampled, dataset is not\n"));
	    ret = 1;
	} else if (submasks_match(dset, pmod)) {
	    /* the subsamples (model and current data set) agree, OK */
	    ret = 0;
	} else {
	    /* the subsamples differ */
	    gretl_errmsg_set(_("model and dataset subsamples not the same\n"));
	    ret = 1;
	}
    }

    return ret;
}

static void dataset_type_string (char *str, const DATASET *dset)
{
    if (dataset_is_time_series(dset)) {
	strcpy(str, _("time series"));
    } else if (dataset_is_panel(dset)) {
        strcpy(str, _("panel"));
    } else {
        strcpy(str, _("undated"));
    }
}

static void pd_string (char *str, const DATASET *dset)
{
    if (custom_time_series(dset)) {
	strcpy(str, _("special"));
    } else {
	switch (dset->pd) {
	case 1:
	    strcpy(str, _("annual")); break;
	case 4:
	    strcpy(str, _("quarterly")); break;
	case 12:
	    strcpy(str, _("monthly")); break;
	case 24:
	    strcpy(str, _("hourly")); break;
	case 52:
	    strcpy(str, _("weekly")); break;
	case 5:
	    strcpy(str, _("daily (5 days)")); break;
	case 6:
	    strcpy(str, _("daily (6 days)")); break;
	case 7:
	    strcpy(str, _("daily (7 days)")); break;
	case 10:
	    strcpy(str, _("decennial")); break;
	default:
	    strcpy(str, _("unknown")); break;
	}
    }
}

void print_sample_obs (const DATASET *dset, PRN *prn)
{
    char d1[OBSLEN], d2[OBSLEN];

    ntodate(d1, dset->t1, dset);
    ntodate(d2, dset->t2, dset);

    pprintf(prn, "%s: %s - %s", _("Current sample"), d1, d2);
    pprintf(prn, " (n = %d)\n", dset->t2 - dset->t1 + 1);
}

void print_sample_status (const DATASET *dset, PRN *prn)
{
    char tmp[128];

    if (complex_subsampled()) {
	pprintf(prn, "%s\n\n", _("Full dataset"));
	dataset_type_string(tmp, fullset);
	pprintf(prn, "%s: %s\n", _("Type"), tmp);
	if (dataset_is_time_series(fullset)) {
	    pd_string(tmp, fullset);
	    pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
	} else if (dataset_is_panel(fullset)) {
	    int nu = fullset->n / fullset->pd;

	    pprintf(prn, "%s: %d\n", _("Number of cross-sectional units"), nu);
	    pprintf(prn, "%s: %d\n", _("Number of time periods"), fullset->pd);
	}
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Range"), 
		fullset->stobs, fullset->endobs, fullset->n);

	pprintf(prn, "\n%s\n", _("Subsampled data"));
	if (dset->restriction != NULL) {
	    pprintf(prn, "(%s: %s)\n\n", _("restriction"), dset->restriction);
	} else {
	    pputc(prn, '\n');
	}
    }	

    dataset_type_string(tmp, dset);
    pprintf(prn, "%s: %s\n", _("Type"), tmp);
    if (dataset_is_time_series(dset)) {
	pd_string(tmp, dset);
	pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
    } else if (dataset_is_panel(dset)) {
	int nu = dset->n / dset->pd;

	pprintf(prn, "%s: %d\n", _("Number of cross-sectional units"), nu);
	pprintf(prn, "%s: %d\n", _("Number of time periods"), dset->pd);
    }
    if (dset->t1 == 0 && dset->t2 == dset->n - 1) {
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Range"), 
		dset->stobs, dset->endobs, dset->n);
    } else {
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Range"), 
		dset->stobs, dset->endobs, dset->n);
	print_sample_obs(dset, prn);
    } 
}

/**
 * data_report:
 * @dset: data information struct.
 * @fname: filename for current datafile.
 * @prn: gretl printing struct.
 * 
 * Write out a summary of the content of the current data set.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 * 
 */

int data_report (const DATASET *dset, const char *fname, PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], tmp[MAXLEN];
    char tstr[48];
    int i;

    ntodate(startdate, 0, dset);
    ntodate(enddate, dset->n - 1, dset);

    sprintf(tmp, _("Data file %s\nas of"), 
	    (*fname != '\0')? fname : _("(unsaved)"));

    print_time(tstr);
    pprintf(prn, "%s %s\n\n", tmp, tstr);

    if (dset->descrip != NULL && *dset->descrip != '\0') {
	pprintf(prn, "%s:\n\n", _("Description"));
	pputs(prn, dset->descrip);
	pputs(prn, "\n\n");
    }

    dataset_type_string(tmp, dset);
    pprintf(prn, "%s: %s\n", _("Type of data"), tmp);
    
    if (dataset_is_time_series(dset)) {
	pd_string(tmp, dset);
	pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
    }	

    pprintf(prn, "%s: %s - %s (n = %d)\n\n", _("Range"),
	    startdate, enddate, dset->n);

    pprintf(prn, "%s:\n\n", _("Listing of variables"));

    for (i=1; i<dset->v; i++) {
	pprintf(prn, "%*s  %s\n", VNAMELEN, dset->varname[i], 
		series_get_label(dset, i));
    }

    return 0;
}



#include "gretl_func.h"
#include "gretl_string_table.h" /* for csvdata */
#include "csvdata.h"

/* call hansl code from dbnomics.gfn to get a series bundle */

static gretl_bundle *get_dbn_series_bundle (const char *datacode,
					    int *err)
{
    gretl_bundle *b = NULL;
    fncall *fc;

    fc = get_pkg_function_call("dbnomics_get_series", "dbnomics", NULL);
    if (fc == NULL) {
	*err = E_DATA;
    } else {
	*err = push_anon_function_arg(fc, GRETL_TYPE_STRING,
				      (void *) datacode);
	if (!*err) {
	    *err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, NULL,
				       &b, NULL, NULL);
	}
	if (b != NULL) {
	    int dberr = gretl_bundle_get_int(b, "error", NULL);

	    if (dberr) {
		const char *msg =
		    gretl_bundle_get_string(b, "errmsg", NULL);

		*err = E_DATA;
		if (msg != NULL) {
		    gretl_errmsg_set(msg);
		} else {
		    gretl_errmsg_sprintf(_("%s: no data found"), datacode);
		}
		gretl_bundle_destroy(b);
		b = NULL;
	    }
	} else if (!*err) {
	    gretl_errmsg_sprintf(_("%s: no data found"), datacode);
	    *err = E_DATA;
	}
    }

    return b;
}

/* write the info from @P (periods) and @v (values) to
   CSV, then grab it back using the gretl CSV reader
   to populate @dbset
*/

static int dbn_dset_from_csv (DATASET *dbset,
			      gretl_array *P,
			      gretl_matrix *v)
{
    gchar *fname;
    FILE *fp;
    int T, err = 0;

    fname = g_strdup_printf("%sdnomics_tmp.txt", gretl_dotdir());
    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	char **S = gretl_array_get_strings(P, &T);
	int t;

	gretl_push_c_numeric_locale();
	fputs("obs dbnomics_data\n", fp);
	for (t=0; t<T; t++) {
	    if (na(v->val[t])) {
		fprintf(fp, "%s NA\n", S[t]);
	    } else {
		fprintf(fp, "%s %.12g\n", S[t], v->val[t]);
	    }
	}
	gretl_pop_c_numeric_locale();

	fclose(fp);
	err = import_csv(fname, dbset, OPT_NONE, NULL);
	gretl_remove(fname);
    }
    
    g_free(fname);

    return err;
}

/* obtain a dbnomics series bundle and process the info
   it contains into a SERIESINFO struct
*/

static int
get_dbnomics_series_info (const char *id, SERIESINFO *sinfo)
{
    gretl_bundle *b;
    DATASET dbset = {0};
    gretl_array *P;
    gretl_matrix *v;
    int T, err = 0;

    /* FIXME check for required form PROV/DSET/SERIES */

    b = get_dbn_series_bundle(id, &err);
    if (err) {
	fprintf(stderr, "get_dbn_series_bundle: err=%d\n", err);
	goto bailout;
    }

    T = gretl_bundle_get_int(b, "T", &err);
    P = gretl_bundle_get_array(b, "period", &err);
    v = gretl_bundle_get_matrix(b, "value", &err);
    if (!err && (T <= 0 || P == NULL || v == NULL)) {
	fprintf(stderr, "get_dbnomics_series_info: invalid bundle content\n");
	err = E_DATA;
	goto bailout;
    }

    /* write bundle content as CSV and use CSV reader to
       construct a one-series dataset */
    err = dbn_dset_from_csv(&dbset, P, v);

    if (!err) {
	/* transcribe info to SERIESINFO format */
	const char *s2 = gretl_bundle_get_string(b, "series_name", NULL);
	char *rawname = strrchr(id, '/') + 1;
	gchar *descrip;

	sinfo->t1 = dbset.t1;
	sinfo->t2 = dbset.t2;
	sinfo->nobs = dbset.n;
	sinfo->pd = dbset.pd;
	strcpy(sinfo->stobs, dbset.stobs);
	strcpy(sinfo->endobs, dbset.endobs);
	/* set up name and description */
	normalize_join_colname(sinfo->varname, rawname, 0, 0);
	descrip = g_strdup_printf("%s: %s", id, s2);
	series_info_set_description(sinfo, descrip);
	g_free(descrip);
	/* steal the data array */
	sinfo->data = dbset.Z[1];
	dbset.Z[1] = NULL;
    }

    free_Z(&dbset);
    clear_datainfo(&dbset, CLEAR_FULL);

 bailout:

    gretl_bundle_destroy(b);

    return err;
}

/* transfer the data stored on @sinfo into @Z */

static int get_dbnomics_data (const char *fname,
			      SERIESINFO *sinfo,
			      double **Z)
{
    memcpy(Z[1], sinfo->data, sinfo->nobs * sizeof(double));
    free(sinfo->data);
    sinfo->data = NULL;

    return 0;
}

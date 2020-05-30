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

/* Handling of ESRI shapfiles and GeoJSON files for gretl */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shapefile.h"
#include "libgretl.h"
#include "gretl_typemap.h"
#include "libset.h"

#define GEODEBUG 0

enum { DBF, SHP, GEO };

enum {
    PRJ0,      /* geoplot default: quasi-Mercator */
    WGS84,     /* null projection */
    EPSG3857,  /* "proper" Mercator */
    EPSG2163,  /* US National Atlas Equal Area */
    EPSG3035,  /* Europe Equal Area */
};

static int proj;

#define GEOHUGE 1.0e100

#ifdef WIN32

static void geoplot_reslash (char *fname)
{
    char *s = fname;

    while (*s) {
	if (*s == '\\') *s = '/';
	s++;
    }
}

#endif

static char *get_full_read_path (char *fname)
{
    if (!g_path_is_absolute(fname)) {
	gretl_addpath(fname, 0);
    }

    return fname;
}

static gchar *ensure_full_write_path (const char *fname)
{
    gchar *ret;

    if (g_path_is_absolute(fname)) {
	ret = g_strdup(fname);
    } else {
	ret = g_build_filename(gretl_workdir(), fname, NULL);
    }

#ifdef WIN32
    geoplot_reslash(ret);
#endif

    return ret;
}

static char *put_ext (char *fname, const char *ext)
{
    char *p = strrchr(fname, '.');

    *p = '\0';
    strcat(p, ext);
    return fname;
}

static int vec_length_error (const char *which, int n,
			     int n_wanted)
{
    gretl_errmsg_sprintf("%s: expected %d rows but got %d",
			 which, n_wanted, n);
    return E_DATA;
}

static int skip_object (int i, const gretl_matrix *z,
			const gretl_matrix *m,
			double *pzi)
{
    int ret = 0;

    if ((z != NULL && i >= z->rows) ||
	(m != NULL && i >= m->rows)) {
	/* safety first! */
	ret = 1;
    } else {
	if (z != NULL) {
	    /* skip on payload value is missing */
	    *pzi = z->val[i];
	    ret = na(*pzi);
	}
	if (!ret && m != NULL && m->val[i] == 1) {
	    /* or skip on entity masked out */
	    ret = 1;
	}
    }

    return ret;
}

static void record_extrema (double x, double y,
			    double *fmin, double *gmax)
{
    if (x < fmin[0]) {
	fmin[0] = x;
    }
    if (x > gmax[0]) {
	gmax[0] = x;
    }
    if (y < fmin[1]) {
	fmin[1] = y;
    }
    if (y > gmax[1]) {
	gmax[1] = y;
    }
}

static gretl_matrix *ring2matrix (gretl_array *ring)
{
    int i, n = gretl_array_get_length(ring);
    gretl_matrix *ret = gretl_matrix_alloc(n, 2);
    GretlType rtype = gretl_array_get_type(ring);
    const char *sx, *sy;
    gretl_matrix *mi;
    gretl_array *ri;

    for (i=0; i<n; i++) {
	if (rtype == GRETL_TYPE_MATRICES) {
	    mi = gretl_array_get_data(ring, i);
	    gretl_matrix_set(ret, i, 0, mi->val[0]);
	    gretl_matrix_set(ret, i, 1, mi->val[1]);
	} else if (rtype == GRETL_TYPE_STRINGS) {
	    ri = gretl_array_get_data(ring, i);
	    sx = gretl_array_get_data(ri, 0);
	    sy = gretl_array_get_data(ri, 1);
	    gretl_matrix_set(ret, i, 0, atof(sx));
	    gretl_matrix_set(ret, i, 1, atof(sy));
	} else {
	    fprintf(stderr, "ring2matrix: i=%d of %d, rtype=%d!\n",
		    i, n, rtype);
	}
    }

    return ret;
}

/* See https://source.opennews.org/articles/choosing-right-map-projection/
   and "Map Projections -- A Working Manual", U.S. Geological Survey
   Professional Paper 1395, John P. Snyder, 1987.
*/

#define d2r (M_PI / 180.0)
#define Radius 1000.0
#define s45 sin(45.0 * d2r)
#define s52 sin(52.0 * d2r)
#define c45 cos(45.0 * d2r)
#define c52 cos(52.0 * d2r)
#define lm100 (-100.0 * d2r)
#define l10 (10.0 * d2r)

/* EPSG:2163, U.S. National Atlas Equal Area,
   and EPSG:3035, Single CRS for all Europe
*/

static void lambert_azimuthal (double *px, double *py)
{
    static double sphivec[2];
    static double cphivec[2];
    static double lam0[2];
    static int filled;
    double lat = *py, lon = *px;
    double phi = lat * d2r;
    double sphi = sin(phi);
    double cphi = cos(phi);
    double lam = lon * d2r;
    double ldiff, cldiff;
    double sphi0, cphi0;
    double k;
    int i = proj == EPSG3035 ? 1 : 0;

    if (!filled) {
	sphivec[0] = s45; sphivec[1] = s52;
	cphivec[0] = c45; cphivec[1] = c52;
	lam0[0] = lm100; lam0[1] = l10;
	filled = 1;
    }

    sphi0 = sphivec[i];
    cphi0 = cphivec[i];
    ldiff = lam - lam0[i];
    cldiff = cos(ldiff);
    k = Radius * sqrt(2.0 / (1 + sphi0*sphi + cphi0*cphi*cldiff));
    *px = k * cphi * sin(ldiff);
    *py = k * (cphi0*sphi - sphi0*cphi*cldiff);
}

/* EPSG:3857 (Mercator)
   x = R * (lam - lam0)
   y = R * ln tan (pi/4 + phi/2)
*/

static void mercator (double *px, double *py)
{
    double lat = *py, lon = *px;
    double phi = lat * d2r;

    *px = Radius * lon * d2r;
    *py = Radius * log(tan(G_PI_4 + 0.5*phi));
}

static int crs_is_nonstandard (gretl_bundle *crs)
{
    gretl_bundle *props;
    const char *crsname;
    int err = 0;
    int ret = 0;

    props = gretl_bundle_get_bundle(crs, "properties", &err);
    if (!err) {
	crsname = gretl_bundle_get_string(props, "name", &err);
    }
    if (!err) {
	const char *s = strstr(crsname, "crs:");

	if (s != NULL) {
	    /* RFC 7946: anything but OGC::CRS84 is non-conforming */
	    if (strcmp(s + 4, "OGC::CRS84")) {
		fprintf(stderr, "Got non-standard crs %s\n", s+4);
		ret = 1;
	    }
	}
    }

    return ret;
}

static gretl_array *geojson_get_features (const char *fname,
					  int *non_standard,
					  int *err)
{
    gretl_bundle *(*jfunc) (const char *, const char *, int *);
    gretl_array *a = NULL;
    GError *gerr = NULL;
    gchar *JSON = NULL;
    gsize len = 0;
    gboolean ok;

    ok = g_file_get_contents(fname, &JSON, &len, &gerr);

    if (ok) {
	gretl_bundle *b = NULL;
	gretl_bundle *c = NULL;
	GretlType type = 0;

	jfunc = get_plugin_function("json_get_bundle");
	if (jfunc == NULL) {
	    *err = E_FOPEN;
	} else {
	    b = jfunc(JSON, NULL, err);
	    if (!*err) {
		a = gretl_bundle_steal_data(b, "features", &type, NULL, err);
	    }
	    if (!*err && gretl_bundle_has_key(b, "crs")) {
		c = gretl_bundle_get_data(b, "crs", &type, NULL, err);
		if (!*err && type == GRETL_TYPE_BUNDLE) {
		    *non_standard = crs_is_nonstandard(c);
		}
	    }
	}
	gretl_bundle_destroy(b);
	g_free(JSON);
    } else if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
    }

    return a;
}

static gretl_array *features_from_bundle (gretl_bundle *b,
					  int *non_standard,
					  int *err)
{
    gretl_array *a = NULL;

    a = gretl_bundle_get_array(b, "features", err);

    if (!*err && gretl_bundle_has_key(b, "crs")) {
	gretl_bundle *c = NULL;
	GretlType type = 0;

	c = gretl_bundle_get_data(b, "crs", &type, NULL, err);
	if (!*err && type == GRETL_TYPE_BUNDLE) {
	    *non_standard = crs_is_nonstandard(c);
	}
    }

    return a;
}

static gretl_matrix *make_bbox (double *gmin, double *gmax)
{
    gretl_matrix *bbox = gretl_matrix_alloc(2, 2);

    if (bbox != NULL) {
	gretl_matrix_set(bbox, 0, 0, gmin[0]);
	gretl_matrix_set(bbox, 0, 1, gmin[1]);
	gretl_matrix_set(bbox, 1, 0, gmax[0]);
	gretl_matrix_set(bbox, 1, 1, gmax[1]);
    }

    return bbox;
}

/* Writes the coordinates content of the GeoJSON file
   identified by @geoname to the gnuplot-compatible
   plain text data file @datname. Returns the bounding
   box for the set of coordinates.
*/

static gretl_matrix *geo2dat (gretl_array *features,
			      const char *datname,
			      const gretl_matrix *zvec,
			      const gretl_matrix *mask,
			      int *non_standard)
{
    gretl_array *AC, *ACj, *ACjk;
    gretl_matrix *X, *bbox = NULL;
    gretl_bundle *fi, *geom;
    double gmin[2] = {GEOHUGE, GEOHUGE};
    double gmax[2] = {-GEOHUGE, -GEOHUGE};
    FILE *fp = NULL;
    const char *gtype;
    int nf, mp, nac, ncj;
    int i, j, k, p;
    int err = 0;

    nf = gretl_array_get_length(features);
    if (zvec != NULL && zvec->rows != nf) {
	err = vec_length_error("payload", zvec->rows, nf);
    } else if (mask != NULL && mask->rows != nf) {
	err = vec_length_error("mask", mask->rows, nf);
    }

    if (!err) {
	fp = gretl_fopen(datname, "wb");
	if (fp == NULL) {
	    err = E_FOPEN;
	}
    }

    if (err) {
	return NULL;
    }

    for (i=0; i<nf; i++) {
	double x, y, z = 0;

	if (skip_object(i, zvec, mask, &z)) {
	    continue;
	}

	fi = gretl_array_get_data(features, i);
	geom = gretl_bundle_get_bundle(fi, "geometry", NULL);
	gtype = gretl_bundle_get_string(geom, "type", NULL);
	if (!strcmp(gtype, "Polygon")) {
	    mp = 0;
	} else if (!strcmp(gtype, "MultiPolygon")) {
	    mp = 1;
	} else {
	    gretl_errmsg_sprintf("can't handle geometry type '%s'", gtype);
	    err = E_DATA;
	    break;
	}
	AC = gretl_bundle_get_array(geom, "coordinates", NULL);
	nac = gretl_array_get_length(AC);
	if (mp == 0) {
	    /* got Polygon */
	    for (j=0; j<nac; j++) {
		ACj = gretl_array_get_data(AC, j);
		X = ring2matrix(ACj);
		for (k=0; k<X->rows; k++) {
		    x = gretl_matrix_get(X, k, 0);
		    y = gretl_matrix_get(X, k, 1);
		    if (proj == EPSG3857) {
			mercator(&x, &y);
		    } else if (proj >= EPSG2163) {
			lambert_azimuthal(&x, &y);
		    }
		    if (zvec != NULL) {
			fprintf(fp, "%.8g %.8g %.8g\n", x, y, z);
		    } else {
			fprintf(fp, "%#.8g %#.8g\n", x, y);
		    }
		    record_extrema(x, y, gmin, gmax);
		}
		gretl_matrix_free(X);
		if (j < nac-1) {
		    fputc('\n', fp);
		}
	    }
	} else {
	    /* got MultiPolygon */
	    for (j=0; j<nac; j++) {
		ACj = gretl_array_get_data(AC, j);
		ncj = gretl_array_get_length(ACj);
		for (k=0; k<ncj; k++) {
		    ACjk = gretl_array_get_data(ACj, k);
		    X = ring2matrix(ACjk);
		    for (p=0; p<X->rows; p++) {
			x = gretl_matrix_get(X, p, 0);
			y = gretl_matrix_get(X, p, 1);
			if (proj == EPSG3857) {
			    mercator(&x, &y);
			} else if (proj >= EPSG2163) {
			    lambert_azimuthal(&x, &y);
			}
			if (zvec != NULL) {
			    fprintf(fp, "%.8g %.8g %.8g\n", x, y, z);
			} else {
			    fprintf(fp, "%#.8g %#.8g\n", x, y);
			}
			record_extrema(x, y, gmin, gmax);
		    }
		    gretl_matrix_free(X);
		    if (k < ncj-1) {
			fputc('\n', fp);
		    }
		}
		if (j < nac-1) {
		    fputc('\n', fp);
		}
	    }
	}
	if (i < nf-1) {
	    fputs("\n\n", fp); /* end of entity block */
	}
    }

    fputc('\n', fp);
    fclose(fp);

    if (!err) {
	bbox = make_bbox(gmin, gmax);
    }

    return bbox;
}

static int json_get_char (gchar **ps, int targ)
{
    gchar *s = *ps;
    int ret = 0;

    while (isspace(*s)) s++;
    if (*s == targ) {
	ret = 1;
	s++;
	while (isspace(*s)) s++;
    }

    *ps = s;

    return ret;
}

/* fast version of geo2dat() -- see above -- which does
   not rely on full JSON parsing.
*/

static gretl_matrix *fast_geo2dat (const char *geoname,
				   const char *datname,
				   const gretl_matrix *zvec,
				   const gretl_matrix *mask)
{
    GError *gerr = NULL;
    gsize len = 0;
    const char *targ = "\"coordinates\"";
    double gmin[2] = {GEOHUGE, GEOHUGE};
    double gmax[2] = {-GEOHUGE, -GEOHUGE};
    gretl_matrix *bbox = NULL;
    FILE *fp;
    char *test;
    int i, nlbr, nrbr;
    gchar *JSON, *p, *s;
    gboolean ok;
    int err = 0;

    ok = g_file_get_contents(geoname, &JSON, &len, &gerr);
    if (!ok) {
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	}
	return NULL;
    }

    fp = gretl_fopen(datname, "wb");
    if (fp == NULL) {
	g_free(JSON);
	return NULL;
    }

    s = JSON;

    for (i=0; !err; i++) {
	double x, y, z = 0;

	p = strstr(s, targ);
	if (p == NULL) break;

	if (skip_object(i, zvec, mask, &z)) {
	    s = p + 14;
	    continue;
	}

	if (i > 0) {
	    fputc('\n', fp);
	}
	p += 13;
	if (!json_get_char(&p, ':')) {
	    gretl_errmsg_set("geo2dat: expected ':'");
	    err = E_DATA;
	    break;
	}
	nlbr = 0;
	while (json_get_char(&p, '[')) {
	    nlbr++;
	}
	if (nlbr < 3 || nlbr > 4) {
	    gretl_errmsg_sprintf("geo2dat: unexpected count of '[': %d", nlbr);
	    err = E_DATA;
	}
	while (!err) {
	    x = strtod(p, &test);
	    if (json_get_char(&test, ',')) {
		p = test;
		y = strtod(p, &test);
		if (!json_get_char(&test, ']')) {
		    gretl_errmsg_sprintf("geo2dat: found '%c', expected ']'", *test);
		    err = E_DATA;
		    break;
		}
		nrbr = 1;
		if (proj == EPSG3857) {
		    mercator(&x, &y);
		} else if (proj >= EPSG2163) {
		    lambert_azimuthal(&x, &y);
		}
		if (zvec != NULL) {
		    fprintf(fp, "%.8g %.8g %.8g\n", x, y, z);
		} else {
		    fprintf(fp, "%#.8g %#.8g\n", x, y);
		}
		record_extrema(x, y, gmin, gmax);
		p = test;
		while (json_get_char(&p, ']')) {
		    nrbr++;
		}
		if (nrbr == nlbr) {
		    /* reached the end of a feature */
		    fputc('\n', fp);
		    break;
		} else if (nrbr > 1) {
		    /* reached the end of a sub-feature */
		    fputc('\n', fp);
		}
		json_get_char(&p, ',');
		while (json_get_char(&p, '[')) {
		    ;
		}
	    } else if (json_get_char(&test, '}')) {
		break;
	    } else {
		gretl_errmsg_sprintf("geo2dat: unexpected char '%c'", *test);
		err = E_DATA;
		break;
	    }
	}
	s = p;
    }

    fclose(fp);
    g_free(JSON);

    if (!err) {
	bbox = make_bbox(gmin, gmax);
    }

    return bbox;
}

static void output_dbf_string (const char *s, FILE *fp)
{
    fputc('"', fp);
    while (*s) {
	if (*s != 0x0d && *s != 0x0a) {
	    fputc(*s, fp);
	}
	s++;
    }
    fputc('"', fp);
}

DBFHandle open_dbf (const char *dbfname,
		    int *fcount, int *rcount,
		    int *err)
{
    DBFHandle DBF = DBFOpen(dbfname, "rb");

    if (DBF == NULL) {
	gretl_errmsg_sprintf("DBFOpen(%s) failed", dbfname);
	*err = E_FOPEN;
    } else {
	*fcount = DBFGetFieldCount(DBF);
	if (*fcount == 0) {
	    gretl_errmsg_set("There are no fields in this DBF table!");
	    *err = E_DATA;
	} else {
	    *rcount = DBFGetRecordCount(DBF);
	    if (*rcount == 0) {
		gretl_errmsg_set("There are no records in this DBF table!");
		*err = E_DATA;
	    }
	}
	if (*err) {
	    DBFClose(DBF);
	}
    }

    return DBF;
}

/* Fill out "properties" bundles within array @ff with info
   from a DBF file.
*/

int dbf_get_properties (gretl_array *ff, const char *dbfname)
{
    DBFHandle DBF;
    int width, decimals;
    int n, fcount, rcount;
    DBFFieldType etype;
    char title[32];
    int i, j;
    int err = 0;

    DBF = open_dbf(dbfname, &fcount, &rcount, &err);
    if (err) {
	return E_FOPEN;
    }

    n = gretl_array_get_length(ff);
    if (rcount != n) {
	gretl_errmsg_sprintf("Number of DBF records (%d) doesn't match "
			     "number of SHP entities (%d)", rcount, n);
	DBFClose(DBF);
	return E_DATA;
    }

    for (j=0; j<rcount && !err; j++) {
	/* access feature bundle */
	gretl_bundle *bf = gretl_array_get_data(ff, j);
	/* create properties bundle */
	gretl_bundle *pj = gretl_bundle_new();

	if (pj == NULL) {
	    err = E_ALLOC;
	    break;
	}

	for (i=0; i<fcount; i++) {
            etype = DBFGetFieldInfo(DBF, i, title, &width, &decimals);
	    if (etype == FTInvalid) {
		continue;
	    }
	    if (DBFIsAttributeNULL(DBF, j, i)) {
		; /* do what? */
	    } else if (etype == FTString) {
		const char *s = DBFReadStringAttribute(DBF, j, i);

		gretl_bundle_set_string(pj, title, s);
	    } else if (etype == FTInteger) {
		int k = DBFReadIntegerAttribute(DBF, j, i);

		gretl_bundle_set_int(pj, title, k);
	    } else if (etype == FTDouble) {
		double x = DBFReadDoubleAttribute(DBF, j, i);

		gretl_bundle_set_scalar(pj, title, x);
	    }
	}

	gretl_bundle_donate_data(bf, "properties", pj,
				 GRETL_TYPE_BUNDLE, 0);
    }

    DBFClose(DBF);

    return err;
}

/* dbf2csv, adapted from dbfdump by Frank Warmerdam.
   Outputs the content of the .dbf component of a shapefile
   (metadata) as CSV.
*/

static int dbf2csv (const char *dbfname,
		    const char *csvname)
{
    DBFHandle DBF;
    FILE *fp;
    int width, decimals;
    int fcount, rcount;
    DBFFieldType etype;
    char title[32];
    int i, j;
    int err = 0;

    DBF = open_dbf(dbfname, &fcount, &rcount, &err);
    if (err) {
	return err;
    }

    fp = gretl_fopen(csvname, "wb");
    if (fp == NULL) {
	DBFClose(DBF);
	return E_FOPEN;
    }

    /* print column headings */
    for (i=0; i<fcount; i++) {
	etype = DBFGetFieldInfo(DBF, i, title, &width, &decimals);
	fputs(title, fp);
	if (i < fcount - 1) {
	    fputc(',', fp);
	}
    }
    fputc('\n', fp);

    /* Read all the records */
    for (j=0; j<rcount; j++) {
	for (i=0; i<fcount; i++) {
            etype = DBFGetFieldInfo(DBF, i, title, &width, &decimals);

	    if (DBFIsAttributeNULL(DBF, j, i)) {
		fputs("(NULL)", fp);
	    } else {
		const char *s;

		switch (etype) {
		case FTString:
		    s = DBFReadStringAttribute(DBF, j, i);
		    output_dbf_string(s, fp);
		    break;
		case FTInteger:
		    fprintf(fp, "%d", DBFReadIntegerAttribute(DBF, j, i));
		    break;
		case FTDouble:
		    fprintf(fp, "%.8g", DBFReadDoubleAttribute(DBF, j, i));
		    break;
		default:
		    break;
		}
	    }
	    if (i < fcount - 1) {
		fputc(',', fp);
	    }
	    fflush(fp);
	}

        if (DBFIsRecordDeleted(DBF, j)) {
            fputs("(DELETED)", fp);
	}
	fputc('\n', fp);
    }

    fclose(fp);
    DBFClose(DBF);

    return 0;
}

#define SHP_DEBUG 0

gretl_bundle *shp_get_bundle (const char *shpname, int *err)
{
    gretl_bundle *ret = NULL;
    gretl_array *ff = NULL;
    SHPHandle SHP;
    double gmin[4], gmax[4];
    char *dbfname;
    int n_shapetype, n_entities;
    int i, j, k, p;

    dbfname = gretl_strdup(shpname);
    put_ext(dbfname, ".dbf");
    *err = gretl_test_fopen(dbfname, "rb");
    if (*err) {
	return NULL;
    }

    SHP = SHPOpen(shpname, "rb");
    if (SHP == NULL) {
	*err = E_FOPEN;
	free(dbfname);
	return NULL;
    }

    ret = gretl_bundle_new();
    if (ret == NULL) {
	*err = E_ALLOC;
	SHPClose(SHP);
	free(dbfname);
	return NULL;
    }

    SHPGetInfo(SHP, &n_entities, &n_shapetype, gmin, gmax);
    SHPSetFastModeReadObject(SHP, TRUE);

#if SHP_DEBUG
    fprintf(stderr, "shp: n_entities = %d\n", n_entities);
#endif
    gretl_bundle_set_string(ret, "type", "FeatureCollection");

    /* array for "feature" bundles */
    ff = gretl_array_new(GRETL_TYPE_BUNDLES, n_entities, err);

    for (i=0; i<n_entities && !*err; i++) {
	gretl_bundle *bfi = NULL;
	gretl_bundle *bgi = NULL;
	gretl_array *mi = NULL;
	gretl_matrix *xy = NULL;
        SHPObject *obj;

	obj = SHPReadObject(SHP, i);

        if (obj == NULL) {
	    fprintf(stderr, "Unable to read shape %d, terminating.\n", i);
	    *err = E_DATA;
        } else if (obj->nParts > 0 && obj->PartStart[0] != 0) {
            fprintf(stderr, "PartStart[0] = %d, not zero as expected.\n",
		    obj->PartStart[0]);
	    *err = E_DATA;
        }

#if SHP_DEBUG
	fprintf(stderr, "*** entity %d: %d vertices in %d part(s) ***\n",
		i+1, obj->nVertices, obj->nParts);
#endif

	if (obj->nParts > 1) {
	    for (p=1; p < obj->nParts && !*err; p++) {
		if (obj->PartStart[p] <= obj->PartStart[p-1]) {
		    gretl_errmsg_set("SHP parts are not in order!");
		    *err = E_DATA;
		}
	    }
	}

	if (!*err) {
	    bfi = gretl_bundle_new(); /* features[i] */
	    bgi = gretl_bundle_new(); /* geometry for features[i] */
	    if (bfi == NULL || bgi == NULL) {
		*err = E_ALLOC;
	    }
	}

	if (!*err) {
	    mi = gretl_array_new(GRETL_TYPE_MATRICES, obj->nParts, err);
	    if (!*err) {
		/* Note 2020-05-10: even when an SHP entity has more
		   than one "part", this does not seem to map onto the
		   GeoJSON "MultiPolygon" type.
		*/
		gretl_bundle_set_string(bgi, "type", "Polygon");
	    }
	}

	for (p=0, j=0; p < obj->nParts && !*err; p++) {
	    int rows;

	    if (p == obj->nParts - 1) {
		rows = obj->nVertices - obj->PartStart[p];
	    } else {
		rows = obj->PartStart[p+1] - obj->PartStart[p];
	    }
	    xy = gretl_matrix_alloc(rows, 2);
	    if (xy == NULL) {
		*err = E_ALLOC;
	    } else {
		for (k=0; k<rows && !*err; k++) {
		    if (j >= obj->nVertices) {
			gretl_errmsg_set("Reading off the end of shp array!");
			*err = E_DATA;
			break;
		    }
		    gretl_matrix_set(xy, k, 0, obj->fX[j]);
		    gretl_matrix_set(xy, k, 1, obj->fY[j]);
		    j++;
		}
		/* push matrix onto array */
		gretl_array_set_data(mi, p, xy);
	    }
	}
#if SHP_DEBUG
	fprintf(stderr, " got %d x,y pairs\n", j);
#endif

	if (!*err) {
	    /* attach array of coordinate matrices under geometry */
	    gretl_bundle_donate_data(bgi, "coordinates", mi, GRETL_TYPE_ARRAY, 0);
	    /* attach geometry bundle under feature[i] */
	    gretl_bundle_donate_data(bfi, "geometry", bgi, GRETL_TYPE_BUNDLE, 0);
	    gretl_bundle_set_string(bfi, "type", "Feature");
	    /* and attach feature[i] to features array */
	    gretl_array_set_data(ff, i, bfi);
	} else {
	    gretl_bundle_destroy(bfi);
	    gretl_bundle_destroy(bgi);
	    gretl_array_destroy(mi);
	}

        SHPDestroyObject(obj);
    }

    SHPClose(SHP);

    if (!*err) {
	*err = dbf_get_properties(ff, dbfname);
    }
    free(dbfname);

    if (!*err) {
	gretl_matrix *bbox;

	gretl_bundle_donate_data(ret, "features", ff,
				 GRETL_TYPE_ARRAY, 0);
	bbox = make_bbox(gmin, gmax);
	if (bbox != NULL) {
	    gretl_bundle_donate_data(ret, "bbox", bbox,
				     GRETL_TYPE_MATRIX, 0);
	}
    } else {
	gretl_array_destroy(ff);
	gretl_bundle_destroy(ret);
	ret = NULL;
    }

    return ret;
}

/* Writes the coordinates content of the shapefile
   identified by @shpname to the gnuplot-compatible
   plain text data file @datname. Returns the bounding
   box for the set of coordinates.
*/

static gretl_matrix *shp2dat (const char *shpname,
			      const char *datname,
			      const gretl_matrix *zvec,
			      const gretl_matrix *mask)
{
    gretl_matrix *bbox = NULL;
    SHPHandle SHP;
    FILE *fp;
    int n_shapetype, n_entities, i, part;
    double gmin[4], gmax[4];
    int get_extrema = 0;
    int nskip = 0;
    int prec = 8;
    int err = 0;

    SHP = SHPOpen(shpname, "rb");
    if (SHP == NULL) {
	return NULL;
    }

    SHPGetInfo(SHP, &n_entities, &n_shapetype, gmin, gmax);

    if (zvec != NULL || mask != NULL) {
	if (zvec != NULL && zvec->rows != n_entities) {
	    err = vec_length_error("payload", zvec->rows, n_entities);
	} else if (mask != NULL && mask->rows != n_entities) {
	    err = vec_length_error("mask", mask->rows, n_entities);
	}
	if (err) {
	    SHPClose(SHP);
	    return NULL;
	}
	for (i=0; i<n_entities; i++) {
	    if ((zvec != NULL && na(zvec->val[i])) ||
		(mask != NULL && mask->val[i] == 1)) {
		nskip++;
		break;
	    }
	}
    }

    if (nskip > 0 || proj > WGS84) {
	for (i=0; i<2; i++) {
	    gmin[i] = GEOHUGE;
	    gmax[i] = -GEOHUGE;
	}
	get_extrema = 1;
    }

    fp = gretl_fopen(datname, "wb");
    if (fp == NULL) {
	SHPClose(SHP);
	return NULL;
    }

    SHPSetFastModeReadObject(SHP, TRUE);

    for (i=0; i<n_entities && !err; i++) {
        SHPObject *obj;
	double x, y, z = 0;
	int j;

	if (skip_object(i, zvec, mask, &z)) {
	    continue;
	}

	obj = SHPReadObject(SHP, i);

        if (obj == NULL) {
	    fprintf(stderr, "Unable to read shape %d, terminating.\n", i);
	    err = E_DATA;
        } else if (obj->nParts > 0 && obj->PartStart[0] != 0) {
            fprintf(stderr, "PartStart[0] = %d, not zero as expected.\n",
		    obj->PartStart[0]);
	    err = E_DATA;
        }

	if (err) {
	    if (obj != NULL) {
		SHPDestroyObject(obj);
	    }
	    break;
	}

        for (j=0, part=1; j<obj->nVertices && !err; j++) {
	    if (part < obj->nParts && obj->PartStart[part] == j) {
		part++;
		fputc('\n', fp);
            }
	    x = obj->fX[j];
	    y = obj->fY[j];
	    if (proj == EPSG3857) {
		mercator(&x, &y);
	    } else if (proj >= EPSG2163) {
		lambert_azimuthal(&x, &y);
	    }
	    if (zvec != NULL) {
		fprintf(fp, "%.*g %.*g %.*g\n", prec, x, prec, y, prec, z);
	    } else {
		fprintf(fp, "%.*g %.*g\n", prec, x, prec, y);
	    }
	    if (get_extrema) {
		record_extrema(x, y, gmin, gmax);
	    }
        }

	if (i < n_entities - 1) {
	    fputs("\n\n", fp);
	}

        SHPDestroyObject(obj);
    }

    fputc('\n', fp);
    fclose(fp);
    SHPClose(SHP);

    if (err) {
	gretl_remove(datname);
    } else {
	bbox = make_bbox(gmin, gmax);
    }

    return bbox;
}

/* Write metadata from GeoJSON file @fname to @csvname */

static int geojson_to_csv (const char *fname,
			   const char *csvname,
			   char **mapname)
{
    GError *gerr = NULL;
    gchar *JSON = NULL;
    gsize len = 0;
    gboolean ok;
    int err = 0;

#if GEODEBUG
    fprintf(stderr, "geojson_to_csv: starting\n");
#endif

    ok = g_file_get_contents(fname, &JSON, &len, &gerr);

#if GEODEBUG
    fprintf(stderr, " g_file_get_contents: ok = %d\n", ok);
#endif

    if (!ok) {
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	} else {
	    fprintf(stderr, "g_file_get_contents failed for '%s'\n", fname);
	}
	err = E_DATA;
    } else {
	gretl_bundle *(*jfunc) (const char *, const char *, int *);
	FILE *fp = NULL;
	void *ptr = NULL;
	GretlType type;
	gretl_bundle *jb, *pp, *fi;
	gretl_array *features = NULL;
	gretl_array *keys = NULL;
	const char *key;
	int nf = 0, nk = 0;
	int i, j;

	jfunc = get_plugin_function("json_get_bundle");
	if (jfunc == NULL) {
	    return E_DATA;
	}
	fp = gretl_fopen(csvname, "wb");
	if (fp == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't open %s for writing"), csvname);
	    return E_FOPEN;
	}
#if GEODEBUG
	fprintf(stderr, " calling json_get_bundle()\n");
#endif
	jb = jfunc(JSON, NULL, &err);
	if (!err) {
	    features = gretl_bundle_get_array(jb, "features", &err);
	    if (err) {
		gretl_errmsg_sprintf(_("Couldn't read '%s'"), "features");
	    }
	}
#if GEODEBUG
	fprintf(stderr, " after json_get_bundle, err = %d\n", err);
#endif
	if (!err) {
	    nf = gretl_array_get_length(features);
	    fi = gretl_array_get_element(features, 0, NULL, &err);
	}
	if (!err) {
	    pp = gretl_bundle_get_bundle(fi, "properties", &err);
	}
	if (!err) {
	    keys = gretl_bundle_get_keys(pp, &err);
	}
	if (!err) {
	    nk = gretl_array_get_length(keys);
	}
	for (j=0; j<nk && !err; j++) {
	    key = gretl_array_get_data(keys, j);
	    fprintf(fp, "%s%c", key, j < nk-1 ? ',' : '\n');
	}
	for (i=0; i<nf && !err; i++) {
	    fi = gretl_array_get_element(features, i, NULL, &err);
	    pp = gretl_bundle_get_bundle(fi, "properties", &err);
	    for (j=0; j<nk && !err; j++) {
		key = gretl_array_get_data(keys, j);
		ptr = gretl_bundle_get_data(pp, key, &type, NULL, &err);
		if (err) {
		    fprintf(stderr, "error at feature %d, propkey %d\n", i, j);
		    break;
		}
		if (type == GRETL_TYPE_STRING) {
		    fprintf(fp, "\"%s\"", (char *) ptr);
		} else if (type == GRETL_TYPE_INT) {
		    fprintf(fp, "%d", *(int *) ptr);
		} else if (type == GRETL_TYPE_DOUBLE) {
		    fprintf(fp, "%g", *(double *) ptr);
		} else {
		    fprintf(stderr, "Got property type %s\n", gretl_type_get_name(type));
		    fprintf(fp, "\"\"");
		}
		fputc(j < nk-1 ? ',' : '\n', fp);
	    }
	}
#if GEODEBUG
	fprintf(stderr, " after getting properties, err = %d\n", err);
#endif
	gretl_array_destroy(keys);
	gretl_bundle_destroy(jb);
	fclose(fp);
    }

    g_free(JSON);

    if (!err && mapname != NULL) {
	*mapname = gretl_strdup(fname);
    }

    return err;
}

/* Write metadata from shapefile @fname to @csvname */

static int shapefile_to_csv (const char *fname,
			     const char *csvname,
			     int ftype,
			     char **mapname)
{
    char *dbfname = gretl_strdup(fname);
    char *shpname = gretl_strdup(fname);
    char *shxname = gretl_strdup(fname);
    int err;

    if (ftype == SHP) {
	put_ext(dbfname, ".dbf");
    } else if (ftype == DBF) {
	put_ext(shpname, ".shp");
    }
    put_ext(shxname, ".shx");

    if (gretl_stat(dbfname, NULL) != 0) {
	gretl_errmsg_sprintf(_("Couldn't open '%s'"), dbfname);
	err = E_FOPEN;
	goto bailout;
    }
    if (gretl_stat(shpname, NULL) != 0) {
	gretl_errmsg_sprintf(_("Couldn't open '%s'"), shpname);
	err = E_FOPEN;
	goto bailout;
    }
    if (gretl_stat(shxname, NULL) != 0) {
	gretl_errmsg_sprintf(_("Couldn't open '%s'"), shxname);
	err = E_FOPEN;
	goto bailout;
    }

    err = dbf2csv(dbfname, csvname);

    if (!err && mapname != NULL) {
	*mapname = shpname;
	shpname = NULL;
    }

 bailout:

    free(dbfname);
    free(shpname);
    free(shxname);

    return err;
}

/* Driver for extracting coordinates data from either
   a shapefile or a GeoJSON file.
*/

static gretl_matrix *map2dat (const char *mapname,
			      gretl_bundle *map,
			      const char *datname,
			      const gretl_matrix *zvec,
			      const gretl_matrix *mask,
			      int *non_standard)
{
    gretl_matrix *ret = NULL;
    char infile[MAXLEN];

    if (mapname != NULL) {
	strcpy(infile, mapname);
	get_full_read_path(infile);
    } else {
	infile[0] = '\0';
    }

    if (zvec != NULL && zvec->cols > 1) {
	gretl_errmsg_set("Invalid payload");
	return NULL;
    }
    if (mask != NULL && mask->cols > 1) {
	gretl_errmsg_set("Invalid mask");
	return NULL;
    }

    gretl_push_c_numeric_locale();

    if (has_suffix(infile, ".shp")) {
	ret = shp2dat(infile, datname, zvec, mask);
    } else if (*infile && libset_get_bool(GEOJSON_FAST)) {
	/* FIXME revisit this */
	ret = fast_geo2dat(infile, datname, zvec, mask);
    } else {
	/* Regular GeoJSON procedure, either reading from
	   @infile or working from pre-loaded @map bundle.
	*/
	gretl_array *features;
	int err = 0;

	if (map != NULL) {
	    features = features_from_bundle(map, non_standard, &err);
	} else {
	    features = geojson_get_features(infile, non_standard, &err);
	}
	if (features != NULL) {
	    ret = geo2dat(features, datname, zvec, mask, non_standard);
	    if (map == NULL) {
		gretl_array_destroy(features);
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

static int real_map_to_csv (const char *fname,
			    const char *csvname,
			    char **mapname)
{
    int ftype;
    int ret;

    if (has_suffix(fname, ".dbf")) {
	ftype = DBF;
    } else if (has_suffix(fname, ".shp")) {
	ftype = SHP;
    } else {
	ftype = GEO;
    }

#if GEODEBUG
    fprintf(stderr, "real_map_to_csv...\n");
#endif

    gretl_push_c_numeric_locale();

    if (ftype == GEO) {
	ret = geojson_to_csv(fname, csvname, mapname);
    } else {
	ret = shapefile_to_csv(fname, csvname, ftype, mapname);
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

/* Get metadata from map file for importation to
   gretl dataset.
*/

int map_get_data (const char *fname, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    gchar *csvname = NULL;
    gchar *base = NULL;
    char *mapname = NULL;
    int err = 0;

    base = g_path_get_basename(fname);
    csvname = gretl_make_dotpath(base);
    put_ext(csvname, ".csv");

    err = real_map_to_csv(fname, csvname, &mapname);

    if (!err) {
	err = import_csv(csvname, dset, opt, prn);
	if (!err) {
	    dset->mapfile = mapname;
	    mapname = NULL;
	}
    }

    g_free(base);
    g_free(csvname);
    free(mapname);

    return err;
}

/* turn a numeric property into a column vector */

static gretl_matrix *payload_from_prop (gretl_bundle *b,
					const char *key,
					int *err)
{
    gretl_matrix *ret = NULL;
    gretl_array *f;
    int i;

    f = gretl_bundle_get_array(b, "features", err);
    if (f != NULL) {
	if (gretl_array_get_type(f) != GRETL_TYPE_BUNDLES) {
	    *err = E_TYPES;
	}
    }

    if (!*err) {
	gretl_bundle *fi = gretl_array_get_data(f, 0);
	int n = gretl_array_get_length(f);
	gretl_bundle *prop;

	ret = gretl_matrix_alloc(n, 1);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
	for (i=0; i<n && !*err; i++) {
	    fi = gretl_array_get_data(f, i);
	    prop = gretl_bundle_get_bundle(fi, "properties", err);
	    ret->val[i] = gretl_bundle_get_scalar(prop, key, err);
	}
	if (*err && ret != NULL) {
	    gretl_matrix_free(ret);
	    ret = NULL;
	}
    }

    return ret;
}

static void set_projection (const char *s)
{
    if (!strcmp(s, "Mercator") || !strcmp(s, "mercator") ||
	!strcmp(s, "EPSG3857")) {
	proj = EPSG3857;
    } else if (!strcmp(s, "WGS84") || !strcmp(s, "EPSG4326")) {
	proj = WGS84;
    } else if (!strcmp(s, "EPSG2163")) {
	proj = EPSG2163;
    } else if (!strcmp(s, "EPSG3035")) {
	proj = EPSG3035;
    } else {
	fprintf(stderr, "set_projection: '%s'?\n", s);
	proj = PRJ0;
    }
}

static gretl_matrix *vector_minmax (const gretl_matrix *z)
{
    gretl_matrix *ret = gretl_matrix_alloc(1, 2);
    double zi;
    int i;

    ret->val[0] = +GEOHUGE;
    ret->val[1] = -GEOHUGE;

    for (i=0; i<z->rows; i++) {
	zi = z->val[i];
	if (!na(zi)) {
	    if (zi < ret->val[0]) {
		ret->val[0] = zi;
	    }
	    if (zi > ret->val[1]) {
		ret->val[1] = zi;
	    }
	}
    }

    return ret;
}

static int is_image_filename (const char *s)
{
    if (has_suffix(s, ".pdf") ||
	has_suffix(s, ".eps") ||
	has_suffix(s, ".png") ||
	has_suffix(s, ".svg") ||
	has_suffix(s, ".emf")) {
	return 1;
    } else {
	/* should really be .gp or .plt? */
	return 0;
    }
}

int geoplot (const char *mapfile,
	     gretl_bundle *map,
	     gretl_matrix *payload,
	     gretl_bundle *opts)
{
    const gretl_matrix *mask = NULL;
    gretl_matrix *bbox = NULL;
    gchar *plotfile = NULL;
    gchar *datfile = NULL;
    const char *sval;
    int plotfile_is_image = 0;
    int non_standard = 0;
    int free_payload = 0;
    int show = 1;
    int err = 0;

    if (opts != NULL) {
	if (gretl_bundle_has_key(opts, "show")) {
	    show = gretl_bundle_get_int(opts, "show", &err);
	}
	if (gretl_bundle_has_key(opts, "plotfile")) {
	    sval = gretl_bundle_get_string(opts, "plotfile", &err);
	    if (sval != NULL) {
		/* note: respect workdir as default output path */
		plotfile = ensure_full_write_path(sval);
		if (is_image_filename(plotfile)) {
		    plotfile_is_image = 1;
		    show = 0;
		}
	    }
	}
	if (!err && map != NULL && payload == NULL &&
	    gretl_bundle_has_key(opts, "payload")) {
	    sval = gretl_bundle_get_string(opts, "payload", &err);
	    if (sval != NULL) {
		payload = payload_from_prop(map, sval, &err);
		free_payload = 1;
	    }
	}
	if (err) {
	    goto bailout;
	}
    }

    /* catch the case of no output */
    if (!show && plotfile == NULL) {
        gretl_errmsg_set("geoplot: no output was specified");
	err = E_ARGS;
	goto bailout;
    }

    /* do we have a sub-sampling mask? */
    if (gretl_bundle_has_key(opts, "mask")) {
	mask = gretl_bundle_get_matrix(opts, "mask", &err);
    }

    /* specific projection wanted? */
    sval = gretl_bundle_get_string(opts, "projection", NULL);
    if (sval != NULL) {
	set_projection(sval);
    }

    if (!err) {
	/* output filenames */
	if (plotfile != NULL && !plotfile_is_image) {
	    datfile = g_strdup_printf("%s.dat", plotfile);
	} else {
	    datfile = gretl_make_dotpath("geoplot_tmp.dat");
	}
#ifdef WIN32
	geoplot_reslash(datfile);
#endif
	/* write out the polygons data for gnuplot */
	bbox = map2dat(mapfile, map, datfile, payload, mask, &non_standard);
	if (bbox == NULL) {
	    err = E_DATA;
	}
    }

    if (!err) {
	gretl_matrix *zrange = NULL;

	if (payload != NULL) {
	    zrange = vector_minmax(payload);
	}
	if (proj > 0) {
	    non_standard = 1;
	}
	err = write_map_gp_file(plotfile, plotfile_is_image, datfile,
				bbox, zrange, opts, non_standard, show);
	gretl_matrix_free(zrange);
	gretl_matrix_free(bbox);
    }

 bailout:

    g_free(plotfile);
    g_free(datfile);
    if (free_payload) {
	gretl_matrix_free(payload);
    }

    return err;
}

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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef WIN32
# include "winconfig.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include <unistd.h> /* for getcwd */

#include "strutils.h"
#include "zipunzip.h"

static gchar *gretl_zipfile_get_topdir (const char *fname)
{
    zipinfo *zinfo;
    gchar *topdir = NULL;

    zinfo = zipfile_get_info(fname, 0, NULL);

    if (zinfo != NULL) {
	int i, n, gotit = 0;
	const gchar *s;

	for (i=0; i<zinfo->nfiles && !gotit; i++) {
	    s = zinfo->fnames[i];
	    if (s != NULL) {
		n = strlen(s);
		if (n > 13 && !strcmp(s + n - 11, "session.xml")) {
		    topdir = g_strndup(s, n - 11);
		    if (topdir != NULL) {
			n = strlen(topdir);
			if (topdir[n-1] == '/' || topdir[n-1] == '\\') {
			    topdir[n-1] = '\0';
			}
		    }
		}
	    }
	}
	zipinfo_destroy(zinfo);
    }

    return topdir;
}

int gretl_native_unzip (const char *fname,
			const char *path,
			gchar **zdirname,
			GError **gerr)
{
    int err;

    if (zdirname != NULL) {
	*zdirname = gretl_zipfile_get_topdir(fname);
	if (*zdirname == NULL) {
	    return 1;
	}
    }

    /* for verbose operation, make 4th arg ZIP_VERBOSE or ZIP_TRACE */
    err = zipfile_extract_files(fname, NULL, path, 0, gerr);

    if (*gerr != NULL && !err) {
	/* shouldn't happen */
	err = 1;
    }

    /* don't let ZIP error codes get confused with gretl ones */
    return (err != 0);
}

/*
 * @fname: full path to zipfile to be created.
 * @path: path relative to userdir for files to be picked up
 * and zipped.
 */

int gretl_native_make_zipfile (const char *fname,
			       const char *path,
			       GError **gerr)
{
    const char *array[2] = { path, NULL };
    int err;

    err = zipfile_archive_files(fname, array, 6, 
				ZIP_RECURSE_DIRS,
				gerr);

    if (*gerr != NULL && !err) {
	/* shouldn't happen */
	err = 1;
    }    

    /* don't let ZIP error codes get confused with gretl ones */
    return (err != 0);
}

int gretl_native_zip_datafile (const char *fname, const char *path,
			       int level, GError **gerr)
{
    char thisdir[FILENAME_MAX];
    int err = 0;

#if ZDEBUG
    fprintf(stderr, "gretl_native_zip_datafile\n"
	    " fname = '%s', path = '%s'\n", fname, path);
#endif

    if (getcwd(thisdir, FILENAME_MAX - 1) == NULL) {
	err = E_FOPEN; /* ?? */
    } else {
	char zipname[FILENAME_MAX];
	const char *array[3] = {
	    "data.xml", "data.bin", NULL
	};

#if ZDEBUG
	fprintf(stderr, " cwd = '%s'\n", thisdir);
#endif
	if (!g_path_is_absolute(fname)) {
	    build_path(zipname, thisdir, fname, NULL);
	} else {
	    strcpy(zipname, fname);
	}
#if ZDEBUG
	fprintf(stderr, " zipname = '%s'\n", zipname);
#endif
	gretl_chdir(path);
	err = zipfile_archive_files(zipname, array, level, 0, gerr);
	gretl_chdir(thisdir);
    }

    if (*gerr != NULL && !err) {
	/* shouldn't happen */
	err = 1;
    }    

    /* don't let ZIP error codes get confused with gretl ones */
    return (err != 0);
}




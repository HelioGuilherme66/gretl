#include <gretl/libgretl.h>
#include <gretl/gretl_func.h>
#include <gretl/gretl_zip.h>

int blacklisted (const char *fname)
{
    char *s, line[256];
    FILE *fp;
    int ret = 0;

    fp = fopen("blacklist.txt", "r");

    if (fp == NULL) {
	return 0;
    }

    while (!ret && fgets(line, sizeof line, fp)) {
	s = line;
	s += strspn(s, " ");
	if (*s != '#') {
	    if (!strncmp(s, fname, strlen(fname))) {
		ret = 1;
	    }
	}
    }

    fclose(fp);

    return ret;
}

int has_suffix (const char *fname, const char *sfx)
{
    int len = strlen(fname);
    
    if (len > 4) {
	const char *p = fname + len - 4;

	return strcmp(p, sfx) == 0;
    }

    return 0;
}

int misspelled_include (const char *targ, const char *found)
{
    char s1[128], s2[128];
    
    *s1 = '\0';
    strncat(s1, targ, 127);
    *s2 = '\0';
    strncat(s2, found, 127);

    if (strlen(s1) == strlen(s2)) {
	gretl_charsub(s1, '-', '_');
	gretl_charsub(s2, '-', '_');
	if (!strcmp(s1, s2)) {
	    return 1;
	}
    }
    
    return 0;
}

/* In printing the script for testing purposes, we'll ensure that
   it "includes" the gfn file in our testing location by printing
   that line first. Then as we print subsequent lines we'll check
   for an inclusion line: if that is found we'll bypass it, but if
   it's not found, that is an error in the package.
*/

void print_sample_script (char *s, const char *pkgname,
			  int zipfile,
			  FILE *fp, FILE *ferr)
{
    char line[2048];
    char gfntest[128];
    char incstr[128];
    int inc_found = 0;
    int inc_ok = 0;

    if (zipfile) {
	fprintf(fp, "include ./%s/%s.gfn\n", pkgname, pkgname);
    } else {
	fprintf(fp, "include ./%s.gfn\n", pkgname);
    }

    sprintf(gfntest, "%s.gfn", pkgname);

    bufgets_init(s);

    while (bufgets(line, sizeof line, s)) {
	char *p;

	if (!inc_found) {
	    p = line + strspn(line, " \t\r\n");
	    if (!strncmp(p, "include ", 8)) {
		*incstr = '\0';
		sscanf(p + 8, "%127s", incstr);
		if (!strcmp(incstr, gfntest)) {
		    inc_found = inc_ok = 1;
		} else if (strstr(incstr, gfntest) != NULL) {
		    inc_found = 1;
		} else if (misspelled_include(gfntest, incstr)) {
		    fprintf(ferr, "*** %s: found this 'include' "
			    "line:\n '%s'\n", gfntest, p);
		    inc_found = 1;
		} else {
		    fprintf(ferr, "*** %s: found this 'include' "
			    "line:\n '%s'\n", gfntest, p);
		    fputs(line, fp);
		}
	    } else {
		fputs(line, fp);
	    }
	} else {
	    fputs(line, fp);
	}
    }

    bufgets_finalize(s);

    if (!inc_found) {
	fprintf(ferr, "%s: missing \"include\"\n", gfntest);
    } else if (!inc_ok) {
	fprintf(ferr, "%s: broken \"include\"\n", gfntest);
    }
}

int main (int argc, char **argv)
{
    char *p, *fnarg;
    char gfnname[128];
    char pkgname[128];
    fnpkg *pkg = NULL;
    char *sample = NULL;
    int zipfile = 0;
    int err = 0;

    if (argc < 2) {
	fprintf(stderr, "%s: please supply a function package filename\n",
		argv[0]);
	exit(EXIT_FAILURE);
    }

    fnarg = argv[1];

    if (has_suffix(fnarg, ".zip")) {
	zipfile = 1;
    } else if (!has_suffix(fnarg, ".gfn")) {
	fprintf(stderr, "%s: expected the name of a .gfn or .zip file\n",
		argv[0]);
	exit(EXIT_FAILURE);
    }

    if (blacklisted(fnarg)) {
	fprintf(stderr, "%s is blacklisted, skipping it\n", fnarg);
	return 0;
    }

    libgretl_init();

    if (zipfile) {
	err = gretl_unzip(fnarg);
	if (err) {
	    fprintf(stderr, "%s: error unzipping %s\n", argv[0], fnarg);
	    exit(EXIT_FAILURE);
	}
    }

    *pkgname = '\0';
    strncat(pkgname, fnarg, 127);
    p = strrchr(pkgname, '.');
    *p = '\0';

    if (zipfile) {
	sprintf(gfnname, "./%s/%s.gfn", pkgname, pkgname);
    } else {
	strcpy(gfnname, "./");
	strncat(gfnname, fnarg, 125);
    }

    pkg = get_function_package_by_filename(gfnname, &err);

    if (!err && pkg != NULL) {
	err = function_package_get_properties(pkg,
					      "sample-script", 
					      &sample,
					      NULL);
    }

    if (err || sample == NULL) {
	fprintf(stderr, "*** %s: no sample script\n", fnarg);
    } else {
	char outname[FILENAME_MAX];
	FILE *fp = NULL;
	FILE *ferr = NULL;

	sprintf(outname, "%s.gfn.inp", pkgname);
	fp = fopen(outname, "w");
	ferr = fopen("pkg_include_errors.txt", "a");

	if (fp == NULL || ferr == NULL) {
	    fprintf(stderr, "%s: couldn't open files for writing\n", argv[0]);
	    exit(EXIT_FAILURE);
	}

	if (zipfile) {
	    FILE *fzip = fopen("zipdirs.txt", "a");
	    
	    if (fzip != NULL) {
		fprintf(fzip, "%s\n", pkgname);
		fclose(fzip);
	    }
	}	

	print_sample_script(sample, pkgname, zipfile, fp, ferr);
	fprintf(stderr, "wrote %s\n", outname);

	fclose(fp);
	fclose(ferr);
    }
    
    libgretl_cleanup();
    
    return 0;
}

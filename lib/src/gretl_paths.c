/*
 *  Copyright (c) Allin Cottrell
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

/* gretl_paths.c for gretl  */

#include "libgretl.h"
#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#else
# include <dirent.h>
#endif

#ifndef WIN32
# include <glib.h>
# if GLIB_CHECK_VERSION(2,0,0)
#  define GLIB2
# endif /* GLIB_CHECK_VERSION */
#endif /* ! WIN32 */

enum {
    CURRENT_DIR,
    DATA_SEARCH,
    SCRIPT_SEARCH,
    USER_SEARCH
};

/* .......................................................... */

static int add_gdt_suffix (char *fname)
{
    int added = 0;

    if (strrchr(fname, '.') == NULL) {
	strcat(fname, ".gdt");
	added = 1;
    }

    return added;
}

/* .......................................................... */

int path_append (char *file, const char *path)
{
    char temp[MAXLEN];
    int n, pathlen = strlen(file) + strlen(path) + 1;

    if (pathlen > MAXLEN) return 1;

    strcpy(temp, path);
    n = strlen(temp);

    if (temp[n - 1] != SLASH && n < MAXLEN - 1) {
	temp[n] = SLASH;
	temp[n + 1] = '\0';
    }

    strcat(temp, file);
    strcpy(file, temp);

    return 0;
}

#ifdef WIN32

static int try_open_file (char *targ, const char *finddir, 
			  WIN32_FIND_DATA *fdata, int code)
{
    FILE *fp = NULL;
    char tmp[MAXLEN];
    int n = strlen(finddir);
    int found = 0;
    
    strcpy(tmp, finddir);
    tmp[n-1] = '\0';
    strcat(tmp, fdata->cFileName);
    strcat(tmp, "\\");
    strcat(tmp, targ);

    fp = fopen(tmp, "r");
    if (fp == NULL && code == DATA_SEARCH) {
	if (add_gdt_suffix(tmp)) {
	    fp = fopen(tmp, "r");
	}
    }

    if (fp != NULL) {
	fclose(fp);
	strcpy(targ, tmp);
	found = 1;
    }	

    return found;
}

static void make_finddir (char *targ, const char *src)
{
    int n = strlen(src);

    strcpy(targ, src);

    if (targ[n-1] != '\\') {
	strcat(targ, "\\*");
    } else {
	strcat(targ, "*");
    }
}

static int got_subdir (WIN32_FIND_DATA *fdata)
{
    int ret = 0;

    if (fdata->dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
	if (strcmp(fdata->cFileName, ".") &&
	    strcmp(fdata->cFileName, "..")) {
	    ret = 1;
	}
    }

    return ret;
}

static int find_in_subdir (const char *topdir, char *fname, int code)
{
    HANDLE handle;
    WIN32_FIND_DATA fdata;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_finddir(finddir, topdir);

    handle = FindFirstFile(finddir, &fdata); 
    if (handle != INVALID_HANDLE_VALUE) {
	if (got_subdir(&fdata)) {
	    found = try_open_file(fname, finddir, &fdata, code);
	} 
	while (!found && FindNextFile(handle, &fdata)) {
	    if (got_subdir(&fdata)) {
		found = try_open_file(fname, finddir, &fdata, code);
	    }
	} 
	FindClose(handle);
    }

    return found;
}

#else /* end of win32 file-finding, on to posix */

static int try_open_file (char *targ, const char *finddir, 
			  struct dirent *dirent, int code)
{
    FILE *fp = NULL;
    char tmp[MAXLEN];
    int found = 0;
    
    strcpy(tmp, finddir);
    strcat(tmp, dirent->d_name);
    strcat(tmp, "/");
    strcat(tmp, targ);

    fp = fopen(tmp, "r");
    if (fp == NULL && code == DATA_SEARCH) {
	if (add_gdt_suffix(tmp)) {
	    fp = fopen(tmp, "r");
	}
    }

    if (fp != NULL) {
	fclose(fp);
	strcpy(targ, tmp);
	found = 1;
    }	

    return found;
}

static void make_finddir (char *targ, const char *src)
{
    int n = strlen(src);

    strcpy(targ, src);

    if (targ[n-1] != '/') {
	strcat(targ, "/");
    } 
}

static int got_subdir (const char *topdir, struct dirent *dirent)
{
    int ret = 0;

    if (strcmp(dirent->d_name, ".") && strcmp(dirent->d_name, "..")) {
	char tmp[MAXLEN];
	DIR *sub;

	strcpy(tmp, topdir);
	strcat(tmp, dirent->d_name);
	sub = opendir(tmp);
	if (sub != NULL) {
	    closedir(sub);
	    ret = 1;
	}
    }

    return ret;
}

static int find_in_subdir (const char *topdir, char *fname, int code)
{
    DIR *dir;
    struct dirent *dirent;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_finddir(finddir, topdir);

    dir = opendir(finddir);
    if (dir != NULL) {
	while (!found && (dirent = readdir(dir))) {
	    if (got_subdir(finddir, dirent)) {
		found = try_open_file(fname, finddir, dirent, code);
	    }
	}
	closedir(dir);
    }

    return found;
}

#endif /* win32 vs posix */

/* .......................................................... */

static char *search_dir (char *fname, const char *topdir, int code)
{
    FILE *test;
    char orig[MAXLEN];

    strcpy(orig, fname);

    if (path_append(fname, topdir) == 0) {
	test = fopen(fname, "r");
	if (test != NULL) {
	    fclose(test);
	    return fname;
	}
	if (code == DATA_SEARCH && add_gdt_suffix(fname)) {
	    test = fopen(fname, "r");
	    if (test != NULL) {
		fclose(test);
		return fname;
	    }
	}	    
	strcpy(fname, orig);
	if (code != CURRENT_DIR && find_in_subdir(topdir, fname, code)) {
	    return fname;
	}
    }

    return NULL;
}

/* .......................................................... */

static int path_is_absolute (const char *fname)
{
    int ret = 0;

#ifdef WIN32
    if (fname[1] == ':') ret = 1; /* drive letter? */
#endif

    /* we'll count as absolute paths specified using "." */
    if (*fname == '.' || *fname == SLASH) {
	ret = 1;
    }

    return ret;
}

/* .......................................................... */

static void make_path_absolute (char *fname, const char *orig)
{
    char thisdir[MAXLEN];
    int offset = 0;

    if (getcwd(thisdir, MAXLEN-1) != NULL) {
#ifdef WIN32		
	lower(thisdir); /* hmmm */
#endif
	if (strstr(fname, thisdir) == NULL) {
	    strcpy(fname, thisdir);
	    strcat(fname, SLASHSTR);
	    if (*orig == '.' && orig[1] == SLASH && strlen(orig) > 2) {
		offset = 2;
	    }
	    strcat(fname, orig + offset);
	}
    }
}

/**
 * addpath:
 * @fname: initially given file name.
 * @ppaths: path information struct.
 * @script: if non-zero, suppose the file is a command script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.  Usually called by getopenfile().
 *
 * Returns: the full name of the file that was found, or NULL if no
 * file could be found.
 */

char *addpath (char *fname, PATHS *ppaths, int script)
{
    char orig[MAXLEN];
    char *tmp = fname;
    FILE *test;

    strcpy(orig, fname);

    /* try opening filename as given */
    test = fopen(fname, "r");
    if (test != NULL) { 
	fclose(test); 
	if (!path_is_absolute(fname)) {
	    make_path_absolute(fname, orig);
	}
	return fname;
    } else if (path_is_absolute(fname)) {  
	/* unable to open file as given: if the path was absolute, fail */
	return NULL;
    }

    /* try looking where script was found */
    if (*ppaths->currdir) {
	if ((fname = search_dir(fname, ppaths->currdir, CURRENT_DIR))) {
	    return fname;
	}
    }

    fname = tmp;
    strcpy(fname, orig);

    if (script) {
	/* for a script, try system script dir (and subdirs) */
	if ((fname = search_dir(fname, ppaths->scriptdir, SCRIPT_SEARCH))) { 
	    return fname;
	}
    } else {
	/* for a data file, try system data dir (and subdirs) */
	if ((fname = search_dir(fname, ppaths->datadir, DATA_SEARCH))) { 
	    return fname;
	}
    } 

    /* or try looking in user's dir (and subdirs) */
    fname = tmp;
    strcpy(fname, orig);
    if ((fname = search_dir(fname, ppaths->userdir, USER_SEARCH))) { 
	return fname;
    }

    fname = tmp;
    strcpy(fname, orig);

    return NULL;
}

static int get_quoted_filename (const char *line, char *fname)
{
    char *p;
    int quote = '"';
    int ret = 0;

    p = strchr(line, quote);
    if (p == NULL) {
	quote = '\'';
	p = strchr(line, quote);
    }

    if (p != NULL) {
	char *q = strrchr(line, quote);

	if (q != NULL) {
	    size_t len = q - p;

	    if (len > 0) {
		*fname = 0;
		strncat(fname, p+1, len-1);
		ret = 1;
	    } 
	}
    }

    return ret;
}

/**
 * getopenfile:
 * @line: command line (e.g. "open foo").
 * @fname: filename to be filled out.
 * @ppaths: path information struct.
 * @setpath: if non-zero, set @ppaths->currdir based on the file
 * that is found (if any).
 * @script: if non-zero, suppose the file is a command script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.
 *
 * Returns: 0 on successful parsing of @line, 1 on error.
 */

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 int setpath, int script)
{
    char *fullname;

    /* get the initial filename off the command line */
    if (get_quoted_filename(line, fname)) return 0; 

    if (sscanf(line, "%*s %s", fname) != 1) return 1;

    /* try a basic path search on this filename */
    fullname = addpath(fname, ppaths, script);

    if (fullname != NULL && setpath) {
	int n, spos = slashpos(fname);

	if (spos) {
	    strncpy(ppaths->currdir, fname, (size_t) spos);
	    n = strlen(ppaths->currdir);
	    ppaths->currdir[n] = SLASH;
	    ppaths->currdir[n+1] = '\0';
	} else {
	    ppaths->currdir[0] = '.';
	    ppaths->currdir[1] = SLASH;
	    ppaths->currdir[2] = '\0';
	}	    
    }

    return 0;
}

/* .......................................................... */

static char *internal_path_stuff (int code, const char *path)
{
    static char gretl_lib_path[MAXLEN];

    if (code == 1) {
#ifdef WIN32
	strcpy(gretl_lib_path, path);
#else
# ifdef GLIB2 
	const char *sfx = "-gtk2/";
# else
	const char *sfx = "-gtk1/";
# endif
	char *p = strstr(path, "/share");

	if (p) {
	    size_t len = p - path;

	    *gretl_lib_path = 0;
	    strncat(gretl_lib_path, path, len);
	    strcat(gretl_lib_path, "/lib/gretl");
	    strcat(gretl_lib_path, sfx);
	} else {
	    sprintf(gretl_lib_path, "%s/lib/gretl%s", path, sfx);
	}
#endif /* WIN32 */
	return NULL;
    } 
    else if (code == 0) {
	return gretl_lib_path;
    }
    return NULL;
}

/* .......................................................... */

const char *fetch_gretl_lib_path (void)
{
    return internal_path_stuff (0, NULL);
}

void set_string_table_written (PATHS *ppaths)
{
    ppaths->status |= STRING_TABLE_WRITTEN;
}

int gretl_string_table_written (PATHS *ppaths)
{
    int ret = 0;

    if (ppaths->status & STRING_TABLE_WRITTEN) ret = 1;

    ppaths->status &= ~STRING_TABLE_WRITTEN;

    return ret;
}

int gretl_using_gui (const PATHS *ppaths)
{
    return ppaths->status & GRETL_USING_GUI;
}

static void ensure_slash (char *str)
{
    if (str[strlen(str) - 1] != SLASH) {
	strcat(str, SLASHSTR);
    }
}

/* .......................................................... */

void show_paths (PATHS *ppaths)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", ppaths->gretldir);
    printf("userdir: %s\n", ppaths->userdir);
    printf("datadir: %s\n", ppaths->datadir);
    printf("scriptdir: %s\n", ppaths->scriptdir);
    printf("gnuplot: %s\n", ppaths->gnuplot);
}

/* .......................................................... */

#ifdef WIN32

int set_paths (PATHS *ppaths, int defaults, int gui)
{
    char envstr[MAXLEN];

    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL) {
	    strcpy(ppaths->gretldir, home);
	    ensure_slash(ppaths->gretldir);
	} else {
	    strcpy(ppaths->gretldir, "c:\\userdata\\gretl\\");
	}

	sprintf(ppaths->binbase, "%sdb\\", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "f:\\"); 

	strcpy(ppaths->x12a, "c:\\userdata\\x12arima\\x12a.exe");
	strcpy(ppaths->x12adir, "c:\\userdata\\x12arima");

	if (gui) {
	    strcpy(ppaths->dbhost_ip, "152.17.150.2");
	} else {
	    ppaths->dbhost_ip[0] = '\0';
	}
	ppaths->currdir[0] = '\0';

	strcpy(ppaths->pngfont, "verdana 8");
    } else {
	ensure_slash(ppaths->gretldir);
    }

    sprintf(ppaths->datadir, "%sdata\\", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts\\", ppaths->gretldir);

    if (gui) {
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretl_hlp.txt"));
	sprintf(ppaths->cmd_helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli_hlp.txt"));
	ppaths->status = GRETL_USING_GUI;
    } else { 
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli_hlp.txt"));
	ppaths->status = 0;
    }

    ensure_slash(ppaths->userdir);

    *ppaths->plotfile = '\0';

    sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%sshare\\gtksourceview-1.0"
	    "\\language-specs", ppaths->gretldir);
    putenv(envstr);

    internal_path_stuff (1, ppaths->gretldir);

    return 0;
}

#else /* not Windows */

int set_paths (PATHS *ppaths, int defaults, int gui)
{
    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL) {
	    strcpy(ppaths->gretldir, home);
	    ensure_slash(ppaths->gretldir);
	} else {
	    strcpy(ppaths->gretldir, GRETL_PREFIX);
	    strcat(ppaths->gretldir, "/share/gretl/");
	} 

	sprintf(ppaths->binbase, "%sdb/", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "/mnt/dosc/userdata/rats/oecd/");

	if (gui) {
	    strcpy(ppaths->dbhost_ip, "152.17.150.2");
	} else {
	    ppaths->dbhost_ip[0] = '\0';
	}

	strcpy(ppaths->gnuplot, "gnuplot");
	*ppaths->pngfont = 0;

	ppaths->currdir[0] = '\0';	

	/* try to set a default userdir */
	home = getenv("HOME");
	if (home != NULL) {
	    strcpy(ppaths->userdir, home);
	    strcat(ppaths->userdir, "/gretl/");
	} else {
	    *ppaths->userdir = '\0';
	}

#ifdef HAVE_X12A
	strcpy(ppaths->x12a, "x12a");
	sprintf(ppaths->x12adir, "%sx12arima", ppaths->userdir);
#endif
    } else {
	ensure_slash(ppaths->gretldir);
    }

    sprintf(ppaths->datadir, "%sdata/", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts/", ppaths->gretldir);
    
    if (gui) {
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretl.hlp"));
	sprintf(ppaths->cmd_helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli.hlp"));
	ppaths->status = GRETL_USING_GUI;
    } else {
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli.hlp"));
	ppaths->status = 0;
    }

    ensure_slash(ppaths->userdir);

    *ppaths->plotfile = '\0';

    internal_path_stuff (1, ppaths->gretldir);

    return 0;
}

#endif /* win32 versus unix */

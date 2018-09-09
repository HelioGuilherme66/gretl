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

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_string_table.h"
#include "texprint.h"

#ifdef USE_RLIB
# include "gretl_foreign.h"
#endif
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <unistd.h>

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h> /* for 'open' */

#include <glib/gstdio.h>

#if defined(WIN32) && defined(PKGBUILD)
# define PLUGIN_SFX "plugins"
#elif defined(USE_GTK3)
# define PLUGIN_SFX "gretl-gtk3"
#else
# define PLUGIN_SFX "gretl-gtk2"
#endif

struct INTERNAL_PATHS {
    char gretldir[MAXLEN];
    char dotdir[MAXLEN];
    char workdir[MAXLEN];
    char gnuplot[MAXLEN];
    char plotfile[MAXLEN];
    char plugpath[MAXLEN];
    char binbase[MAXLEN+4];
    char x12a[MAXLEN];
    char x12adir[MAXLEN];
    char tramo[MAXLEN];
    char tramodir[MAXLEN];
    char rbinpath[MAXLEN];
    char rlibpath[MAXLEN];
    char oxlpath[MAXLEN];
    char octpath[MAXLEN];
    char statapath[MAXLEN];
    char pypath[MAXLEN];
    char jlpath[MAXLEN];
    char mpiexec[MAXLEN];
    char mpi_hosts[MAXLEN];
    char dbhost[32];
    char pngfont[128];
    unsigned char status;
};

static struct INTERNAL_PATHS paths;

static char current_dir[MAXLEN];

static int force_en_cmdref;
static int force_en_fnref;

static void set_helpfile_option (gretlopt opt)
{
    if (opt & OPT_N) {
	force_en_cmdref = 1;
	force_en_fnref = 1;
    }
}

static const char *helpfiles[] = {
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_commands_en.xml has been translated for your language,
       as in gretl_cli_cmdref.pt -- otherwise leave it untranslated
    */
    N_("gretl_cli_cmdref.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_commands_en.xml has been translated for your language,
       as in gretl_gui_cmdref.pt -- otherwise leave it untranslated
    */
    N_("gretl_gui_cmdref.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_commands_en.xml has been translated for your language,
       as in gretl_gui_help.pt -- otherwise leave it untranslated
    */
    N_("gretl_gui_help.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_functions_en.xml has been translated for your language,
       as in gretl_cli_fnref.pt -- otherwise leave it untranslated
    */
    N_("gretl_cli_fnref.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_functions_en.xml has been translated for your language,
       as in gretl_gui_fnref.pt -- otherwise leave it untranslated
    */
    N_("gretl_gui_fnref.en")
};

const char *helpfile_path (int id, int cli, int en)
{
    static char hpath[MAXLEN+18];
    int i = -1;

    *hpath = '\0';

    if ((id == GRETL_CMDREF && force_en_cmdref) ||
	(id == GRETL_FUNCREF && force_en_fnref)) {
	en = 1;
    }

    if (cli) {
	/* Command-line program */
	if (id == GRETL_CMDREF) {
	    i = 0;
	} else if (id == GRETL_FUNCREF) {
	    i = 3;
	}
    } else {
	/* GUI program */
	if (id == GRETL_CMDREF) {
	    i = 1;
	} else if (id == GRETL_GUI_HELP) {
	    i = 2;
	} else if (id == GRETL_FUNCREF) {
	    i = 4;
	}
    }

    if (i >= 0) {
	sprintf(hpath, "%s%s", paths.gretldir, en ? helpfiles[i] :
		_(helpfiles[i]));
    }

    return hpath;
}

int using_translated_helpfile (int id)
{
    int ret = 0;
    int i = 0;

    if (id == GRETL_CMDREF) {
	if (force_en_cmdref) return 0;
	i = 1;
    } else if (id == GRETL_FUNCREF) {
	if (force_en_fnref) return 0;
	i = 4;
    } else {
	return 0;
    }

    /* If we're not forcing English help, the criterion is
       that the relevant help file has a "translated" filename
       and the translation can actually be opened.
    */

    if (strcmp(helpfiles[i], _(helpfiles[i]))) {
	char test[MAXLEN];
	int err;

	sprintf(test, "%s%s", paths.gretldir, _(helpfiles[i]));
	err = gretl_test_fopen(test, "r");
	if (err) {
	    if (id == GRETL_CMDREF) {
		force_en_cmdref = 1;
	    } else {
		force_en_fnref = 1;
	    }
	} else {
	    ret = 1;
	}
    }

    return ret;
}

/* If @fname does not already have suffix @sfx, add it.
   With the qualification that if the @fname bears either of
   the standard gretl data-file suffixes, ".gdt" or ".gdtb",
   we won't stick the other one onto the end.
*/

static int maybe_add_suffix (char *fname, const char *sfx)
{
    if (has_suffix(fname, ".gdtb") && !strcmp(sfx, ".gdt")) {
	return 0;
    } else if (has_suffix(fname, ".gdt") && !strcmp(sfx, ".gdtb")) {
	return 0;
    } else if (!has_suffix(fname, sfx)) {
	strcat(fname, sfx);
	return 1;
    }

    return 0;
}

/* Convenience wrapper macro for the GLib UTF-8 validation
   function. If and only if this returns non-zero for a
   given filename can that name be passed to GLib's gstdio
   functions on MS Windows.
*/

#define valid_utf8(s) g_utf8_validate(s, -1, NULL)

/**
 * utf8_encoded:
 * @s: the string to examine.
 *
 * The primary use of this function is to determine
 * whether a filename can be passed to regular C-library
 * functions on MS Windows: if it's in UTF-8 the answer
 * is No -- unless it's in the ASCII subset of UTF-8.
 *
 * Returns: non-zero if @s validates as UTF-8 and
 * contains bytes that are not printable ASCII,
 * otherwise zero.
 */

int utf8_encoded (const char *s)
{
    int ret = 0;

    if (g_utf8_validate(s, -1, NULL)) {
	const unsigned char *p = (const unsigned char *) s;

	while (*p) {
	    if (*p < 32 || *p > 126) {
		/* not printable ASCII */
		ret = 1;
		break;
	    }
	    p++;
	}
    }

    return ret;
}

#define FDEBUG 0

/**
 * gretl_fopen:
 * @fname: name of file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for the C library's fopen(), making allowance for
 * the possibility that @fname has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE *gretl_fopen (const char *fname, const char *mode)
{
    FILE *fp = NULL;

    gretl_error_clear();

#if FDEBUG
    fprintf(stderr, "gretl_fopen: got '%s'\n", fname);
#endif

#ifdef WIN32
    if (valid_utf8(fname)) {
	fp = g_fopen(fname, mode);
    } else {
	fp = fopen(fname, mode);
    }
#else
    fp = fopen(fname, mode);
#endif

#if FDEBUG
    fprintf(stderr, "  after fopen, errno = %d\n", errno);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno(fname, errno);
    }

    return fp;
}

/**
 * gretl_mktemp:
 * @pattern: template for filename; must end with "XXXXXX".
 * @mode: e.g. "w" for text use or "wb" for binary mode.
 *
 * A wrapper for the combination of mkstemp() and fdopen(),
 * making allowance for the possibility that @pattern has to
 * be converted from UTF-8 to the locale encoding or vice versa.
 * On successful exit @pattern holds the name of the newly
 * created file.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE *gretl_mktemp (char *pattern, const char *mode)
{
    FILE *fp = NULL;
    int fd;

    gretl_error_clear();

#ifdef WIN32
    if (valid_utf8(pattern)) {
	fd = g_mkstemp(pattern);
    } else {
	fd = mkstemp(pattern);
    }
#else
    fd = mkstemp(pattern);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno(NULL, errno);
    } else if (fd != -1) {
	fp = fdopen(fd, mode);
    }

#if 0
    fprintf(stderr, "gretl_mktemp: name='%s', fd=%d, fp=%p\n",
	    pattern, fd, (void *) fp);
#endif

    return fp;
}

/**
 * gretl_test_fopen:
 * @fname: name of file to be opened.
 * @mode: mode as used with fopen().
 *
 * Attempts to open @fname in the given mode; if the opening
 * is successful the stream is then closed.
 *
 * Returns: 0 on success, -1 on filename encoding
 * failure, or the system errno on failure of fopen().
 */

int gretl_test_fopen (const char *fname, const char *mode)
{
    FILE *fp = NULL;
    int err = 0;

    gretl_error_clear();

#ifdef WIN32
    if (valid_utf8(fname)) {
	fp = g_fopen(fname, mode);
	if (fp == NULL) {
	    err = errno;
	} else {
	    fclose(fp);
	    if (*mode == 'w') {
		g_remove(fname);
	    }
	}
	return err;
    }
#endif

    fp = fopen(fname, mode);
    if (fp == NULL) {
	err = errno;
    } else {
	fclose(fp);
	if (*mode == 'w') {
	    remove(fname);
	}
    }

    return err;
}

#ifdef WIN32

static int win32_open_fchdir (int fd)
{
    static char savedir[FILENAME_MAX];
    int ret = 0;

    if (fd == 0) {
	/* call to record current directory */
	if (getcwd(savedir, FILENAME_MAX - 1) == NULL) {
	    *savedir = '\0';
	    /* return invalid value */
	    ret = -1;
	} else {
	    /* return fake positive file id! */
	    ret = 1;
	}
    } else {
	/* call to change back to prior directory */
	if (*savedir == '\0') {
	    ret = -1;
	} else {
	    ret = gretl_chdir(savedir);
	}
    }

    return ret;
}

#endif

/**
 * gretl_open:
 * @pathname: name of file to be opened.
 * @flags: flags to pass to the system open().
 * @mode: ignored unless @flags contains O_CREAT
 * or O_TMPFILE.
 *
 * A wrapper for the C library's open(), making allowance for
 * the possibility that @pathname has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: new file descriptor, or -1 on error.
 */

int gretl_open (const char *pathname, int flags, int mode)
{
    mode_t m = 0;
    int fd = -1;

    gretl_error_clear();

    if (flags & O_CREAT) {
	m = (mode_t) mode;
    }

#ifdef WIN32
    if (!strcmp(pathname, ".")) {
	return win32_open_fchdir(0);
    } else if (valid_utf8(pathname)) {
	fd = g_open(pathname, flags, m);
    } else {
	fd = open(pathname, flags, m);
    }
#else
    fd = open(pathname, flags, m);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno(pathname, errno);
    }

    return fd;
}

/* On MS Windows support cheap fakery of fchdir, otherwise
   call the real fchdir() */

int gretl_fchdir (int fd)
{
#ifdef WIN32
    return win32_open_fchdir(fd);
#else
    return fchdir(fd);
#endif
}

/**
 * gretl_stat:
 * @fname: name of file to be examined.
 * @buf: pointer to a C struct stat (or NULL is OK if
 * the caller just wants the return value from the stat
 * call).
 *
 * A wrapper for the C library's stat(), making allowance for
 * the possibility that @fname has to be converted from UTF-8
 * to the locale encoding or vice versa.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_stat (const char *fname, struct stat *buf)
{
    struct stat tmp = {0};

    gretl_error_clear();

#ifdef WIN32
    if (utf8_encoded(fname)) {
	/* A native stat() call won't work with such a filename:
	   we should either call g_stat(), which expects UTF-8
	   on Windows, or convert @fname before calling stat().
	   Unfortunately g_stat() from GLib 2.36.4 crashes on
	   (some variants of) 32-bit Windows, so it seems we need
	   to do the conversion ourselves.
	*/
	gunichar2 *wname;

	wname = g_utf8_to_utf16(fname, -1, NULL, NULL, NULL);
	if (wname != NULL) {
	    int ret = wstat(wname, buf == NULL ? &tmp : buf);

	    g_free(wname);
	    return ret;
	}
    }
#endif

    return stat(fname, buf == NULL ? &tmp : buf);
}

/**
 * gretl_file_exists:
 * @fname: name of file to be examined.
 *
 * Uses the C library's stat() function, making allowance for
 * the possibility that @fname has to be converted from UTF-8
 * to the locale encoding or vice versa.
 *
 * Returns: 1 if @fname is the name of an existing file,
 * otherwise 0.
 */

int gretl_file_exists (const char *fname)
{
    return gretl_stat(fname, NULL) == 0;
}

#ifdef WIN32

static int win32_rename (const char *oldpath,
			 const char *newpath)
{
    int u_old = valid_utf8(oldpath);
    int u_new = valid_utf8(newpath);
    int ret = -1;

    if (u_old && u_new) {
	/* OK, both names are UTF-8 */
	ret = g_rename(oldpath, newpath);
    } else if (u_old || u_new) {
	/* let's get both names into the locale charset */
	gchar *tmp = NULL;

	if (u_old) {
	    /* only the original name is UTF-8 */
	    tmp = g_win32_locale_filename_from_utf8(oldpath);
	} else {
	    /* only the new name is UTF-8 */
	    tmp = g_win32_locale_filename_from_utf8(newpath);
	}

	if (tmp == NULL) {
	    /* couldn't convert a filename */
	    ret = -1;
	} else {
	    /* rename() on Windows sets EEXIST if the target
	       is already present */
	    remove(newpath);
	    if (u_old) {
		ret = rename(tmp, newpath);
	    } else {
		ret = rename(oldpath, tmp);
	    }
	    g_free(tmp);
	}
    } else {
	/* neither name is in UTF-8 */
	remove(newpath);
	ret = rename(oldpath, newpath);
    }

    return ret;
}

#endif

/**
 * gretl_rename:
 * @oldpath: name of file to be opened.
 * @newpath: new name to give the file.
 *
 * A wrapper for the C library's rename(), making allowance for
 * the possibility that @oldpath and/or @newpath have to be
 * converted from UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_rename (const char *oldpath, const char *newpath)
{
    int err = 0;

    gretl_error_clear();

#ifdef WIN32
    err = win32_rename(oldpath, newpath);
#else
    err = rename(oldpath, newpath);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno("gretl_rename", errno);
    }

    return err;
}

/**
 * gretl_remove:
 * @path: name of file or directory to remove.
 *
 * A wrapper for remove(), making allowance for
 * the possibility that @path has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 0 on sucess, non-zero on failure.
 */

int gretl_remove (const char *path)
{
    int ret = -1;

#ifdef WIN32
    if (valid_utf8(path)) {
	ret = g_remove(path);
    } else {
	/* @path is in locale encoding */
	ret = remove(path);
	if (ret == -1) {
	    /* allow for the possibility that we're trying to remove a
	       directory on win32 -> recode to UTF-8 and use g_remove
	    */
	    GError *gerr = NULL;
	    gchar *pconv = NULL;
	    gsize sz;

	    pconv = g_locale_to_utf8(path, -1, NULL, &sz, &gerr);
	    if (pconv == NULL) {
		if (gerr != NULL) {
		    gretl_errmsg_set(gerr->message);
		    g_error_free(gerr);
		}
	    } else {
		ret = g_remove(pconv);
		g_free(pconv);
	    }
	}
    }
#else
    ret = remove(path);
#endif

    return ret;
}

/**
 * gretl_gzopen:
 * @fname: name of gzipped file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for zlib's gzopen(), making allowance for
 * the possibility that @fname has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: pointer to gzip stream, or %NULL on failure.
 */

gzFile gretl_gzopen (const char *fname, const char *mode)
{
    gzFile fz = NULL;

    gretl_error_clear();

#ifdef WIN32
    if (utf8_encoded(fname)) {
	/* have to convert to locale */
	gchar *tmp = g_win32_locale_filename_from_utf8(fname);

	if (tmp != NULL) {
	    fz = gzopen(tmp, mode);
	    g_free(tmp);
	}
    } else {
	fz = gzopen(fname, mode);
    }
#else
    fz = gzopen(fname, mode);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno("gzopen", errno);
    }

    return fz;
}

/**
 * gretl_chdir:
 * @path: name of directory.
 *
 * A wrapper for POSIX chdir(), making allowance for
 * the possibility that @path has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 0 on sucess, non-zero on failure.
 */

int gretl_chdir (const char *path)
{
    int err = 0;

    gretl_error_clear();

#ifdef WIN32
    if (1) {
	int len = strlen(path);
	char *ptmp = NULL;

	if (len > 1 && path[len-1] == '\\' && path[len-2] != ':') {
	    /* trim trailing slash for non-root dir */
	    ptmp = gretl_strndup(path, len - 1);
	    path = ptmp;
	}
	if (valid_utf8(path)) {
	    err = g_chdir(path);
	} else {
	    err = chdir(path);
	}
	free(ptmp);
    }
#else
    err = chdir(path);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno("chdir", errno);
    }

    return err;
}

/**
 * gretl_isdir:
 * @path: path to check.
 *
 * A test for whether or not @path is the name of a directory,
 * allowing for the possibility that @path has to be converted
 * from UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 1 if @path is the name of a directory, else 0.
 */

int gretl_isdir (const char *path)
{
    struct stat buf;
    int err;

    err = gretl_stat(path, &buf);

    return err ? 0 : S_ISDIR(buf.st_mode);
}

/**
 * gretl_mkdir:
 * @path: name of directory to be created.
 *
 * Calls the underlying library function to create the
 * specified directory with mode 0755.  If the directory in
 * question already exists, this does not count as an error.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_mkdir (const char *path)
{
    int err, tried = 0;

    errno = 0;

#ifdef WIN32
    if (valid_utf8(path)) {
	err = g_mkdir_with_parents(path, 0755);
    } else {
	gchar *pconv = g_locale_to_utf8(path, -1, NULL, NULL, NULL);

	if (pconv != NULL) {
	    err = g_mkdir_with_parents(path, 0755);
	    g_free(pconv);
	    tried = 1;
	}
    }
#endif

    if (!tried) {
	err = g_mkdir_with_parents(path, 0755);
    }

    if (err) {
	fprintf(stderr, "%s: %s\n", path, gretl_strerror(errno));
	err = 1;
    }

    return err;
}

#ifndef WIN32

static const char *gretl_readd (DIR *d)
{
    struct dirent *e = readdir(d);

    return (e == NULL)? NULL : e->d_name;
}

static int real_deltree (const char *path)
{
    DIR *dir;
    int err = 0;

    errno = 0;
    dir = gretl_opendir(path);

    if (dir == NULL) {
	err = 1;
    } else {
	const char *fname;

	err = chdir(path);
	while ((fname = gretl_readd(dir)) != NULL && !err) {
	    /* recursively delete dir's contents */
	    if (strcmp(fname, ".") && strcmp(fname, "..")) {
		if (gretl_isdir(fname)) {
		    err = real_deltree(fname);
		} else {
		    err = remove(fname);
		}
	    }
	}
	if (!err) {
	    closedir(dir);
	    /* delete the directory itself */
	    if (chdir("..") == 0) {
		err = gretl_remove(path);
	    }
	}
    }

    if (err) {
	gretl_errmsg_set_from_errno(path, errno);
	err = E_FOPEN;
    }

    return err;
}

#endif /* WIN32 or not */

/**
 * gretl_deltree:
 * @path: name of directory to be deleted.
 *
 * Carries out recursive deletion of the specified directory.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_deltree (const char *path)
{
#ifdef WIN32
    return win32_delete_dir(path);
#else
    char tmp[FILENAME_MAX];
    char *savedir = NULL;
    int err;

    savedir = getcwd(tmp, FILENAME_MAX - 1);

    err = real_deltree(path);

    if (savedir != NULL) {
	gretl_chdir(savedir);
    }

    return err;
#endif
}

DIR *gretl_opendir (const char *name)
{
#ifdef WIN32
    int n = strlen(name);
    int fixit = 0;

    if (n > 0 && name[n-1] == ':') {
	/* opendir doesn't work on e.g. "f:" */
	fixit = 1;
    } else if (n > 3 && name[n-1] == '\\') {
	/* and neither does it work on e.g. "c:\foo\" */
	fixit = 2;
    }

    if (fixit) {
	char tmp[MAXLEN];

	*tmp = '\0';
	strncat(tmp, name, MAXLEN - 2);
	if (fixit == 1) {
	    /* append backslash */
	    strcat(tmp, "\\");
	} else {
	    /* chop trailing backslash */
	    tmp[strlen(tmp)-1] = '\0';
	}
	return opendir(tmp);
    }
#endif
    return opendir(name);
}

/**
 * gretl_setenv:
 * @name: name of variable to be set.
 * @value: value to set.
 *
 * Cross-platform wrapper for setenv().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_setenv (const char *name, const char *value)
{
#ifdef WIN32
    gchar *estr;
    int ok;

    /* belt and braces */
    estr = g_strdup_printf("%s=%s", name, value);
    putenv(estr);

    ok = SetEnvironmentVariable(name, value);

    return !ok;
#else
    return setenv(name, value, 1);
#endif
}

/**
 * gretl_write_access:
 * @fname: name of file to test.
 *
 * Returns: 0 on success (meaning that the current user has
 * write access to @fname), non-zero on failure.
 */

int gretl_write_access (char *fname)
{
    int err;

    gretl_error_clear();

#ifdef WIN32
    if (utf8_encoded(fname)) {
	gchar *tmp = g_win32_locale_filename_from_utf8(fname);

	if (tmp != NULL) {
	    err = !win32_write_access(tmp);
	    g_free(tmp);
	} else {
	    err = 1;
	}
    } else {
	err = !win32_write_access(fname);
    }
#else
    err = access(fname, W_OK);
    if (errno != 0) {
	gretl_errmsg_set_from_errno(fname, errno);
    }
#endif

    return err;
}

/**
 * gretl_is_xml_file:
 * @fname: name of file to test.
 *
 * Returns: 1 if @fname appears to be a (possibly gzipped) XML file,
 * otherwise 0.
 */

int gretl_is_xml_file (const char *fname)
{
    gzFile fz;
    char test[6];
    int ret = 0;

    fz = gretl_gzopen(fname, "rb");
    if (fz != Z_NULL) {
	if (gzread(fz, test, 5)) {
	    test[5] = '\0';
	    if (!strcmp(test, "<?xml")) ret = 1;
	}
	gzclose(fz);
    }

    gretl_error_clear();

    return ret;
}

/**
 * gretl_path_prepend:
 * @file: target filename.
 * @path: path to prepend.
 *
 * Creates a path string by prepending @path, plus an appropriate
 * separator if needed, to @file.  The result is written back into
 * @file: this variable is assumed to have storage for at least
 * #MAXLEN characters.
 *
 * Returns: 0 on success, or 1 if the final path string would
 * exceed #MAXLEN characters (including nul-termination).
 */

int gretl_path_prepend (char *file, const char *path)
{
    int n = strlen(file) + strlen(path) + 1;
    char temp[MAXLEN];

    if (n > MAXLEN) {
	return 1;
    }

    strcpy(temp, path);
    n = strlen(temp);

    if (temp[n-1] != SLASH && n < MAXLEN - 1) {
	strcat(temp, SLASHSTR);
    }

    strcat(temp, file);
    strcpy(file, temp);

    return 0;
}

enum {
    ADD_GDT = 1 << 0,
    ADD_GFN = 1 << 1,
    SUBDIRS = 1 << 2
};

#ifdef WIN32

static int try_open_file (char *targ, const char *finddir,
			  WIN32_FIND_DATA *fdata, int flags)
{
    char tmp[MAXLEN];
    int n = strlen(finddir);
    int err, found = 0;

    strcpy(tmp, finddir);
    tmp[n-1] = '\0';
    strcat(tmp, fdata->cFileName);
    strcat(tmp, "\\");
    strcat(tmp, targ);

    err = gretl_test_fopen(tmp, "r");

    if (err && (flags & ADD_GDT)) {
	if (maybe_add_suffix(tmp, ".gdt")) {
	    err = gretl_test_fopen(tmp, "rb");
	    if (err) {
		/* try .gdtb also */
		strncat(tmp, "b", 1);
		err = gretl_test_fopen(tmp, "rb");
	    }
	}
    }

    if (!err) {
	strcpy(targ, tmp);
	found = 1;
    }

    return found;
}

static void make_findname (char *targ, const char *src)
{
    strcpy(targ, src);

    if (utf8_encoded((const unsigned char *) targ)) {
	/* FIXME : do we want this? */
	gchar *tmp;
	gsize sz;

	tmp = g_locale_from_utf8(src, -1, NULL, &sz, NULL);
	if (tmp != NULL) {
	    strcpy(targ, tmp);
	    g_free(tmp);
	}
    }

    if (targ[strlen(targ)-1] != '\\') {
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

static int find_in_subdir (const char *topdir, char *fname, int flags)
{
    HANDLE handle;
    WIN32_FIND_DATA fdata;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_findname(finddir, topdir);

    handle = FindFirstFile(finddir, &fdata);
    if (handle != INVALID_HANDLE_VALUE) {
	if (got_subdir(&fdata)) {
	    found = try_open_file(fname, finddir, &fdata, flags);
	}
	while (!found && FindNextFile(handle, &fdata)) {
	    if (got_subdir(&fdata)) {
		found = try_open_file(fname, finddir, &fdata, flags);
	    }
	}
	FindClose(handle);
    }

    return found;
}

#else /* end of win32 file-finding, on to posix */

static int try_open_file (char *targ, const char *finddir,
			  struct dirent *dirent, int flags)
{
    char tmp[MAXLEN];
    int err, found = 0;

    strcpy(tmp, finddir);
    strcat(tmp, dirent->d_name);
    strcat(tmp, "/");
    strcat(tmp, targ);

    err = gretl_test_fopen(tmp, "r");

    if (err && (flags & ADD_GDT)) {
	if (maybe_add_suffix(tmp, ".gdt")) {
	    err = gretl_test_fopen(tmp, "r");
	    if (err) {
		/* try .gdtb also */
		strcat(tmp, "b");
		err = gretl_test_fopen(tmp, "r");
	    }
	}
    }

    if (!err) {
	strcpy(targ, tmp);
	found = 1;
    }

    return found;
}

static void make_findname (char *targ, const char *src)
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

static int find_in_subdir (const char *topdir, char *fname, int flags)
{
    DIR *dir;
    struct dirent *dirent;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_findname(finddir, topdir);

    dir = opendir(finddir);
    if (dir != NULL) {
	while (!found && (dirent = readdir(dir))) {
	    if (got_subdir(finddir, dirent)) {
		found = try_open_file(fname, finddir, dirent, flags);
	    }
	}
	closedir(dir);
    }

    return found;
}

#endif /* win32 vs posix */

#define SEARCH_DEBUG 0

char *search_dir (char *fname, const char *topdir, int flags)
{
    char orig[MAXLEN];
    int err;

#if SEARCH_DEBUG
    fprintf(stderr, "search_dir: trying '%s' for '%s' (ADD_GDT=%d)\n",
	    topdir, fname, flags & ADD_GDT);
#endif

    strcpy(orig, fname);

    if (gretl_path_prepend(fname, topdir) == 0) {
	err = gretl_test_fopen(fname, "r");
#if SEARCH_DEBUG
	fprintf(stderr, " trying '%s'\n", fname);
#endif
	if (!err) {
	    return fname;
	}
	if (flags & ADD_GDT) {
	    if (maybe_add_suffix(fname, ".gdt")) {
#if SEARCH_DEBUG
		fprintf(stderr, " trying '%s'\n", fname);
#endif
		err = gretl_test_fopen(fname, "r");
		if (!err) {
		    return fname;
		}
	    } else if (maybe_add_suffix(fname, ".gdtb")) {
		err = gretl_test_fopen(fname, "r");
		if (!err) {
		    return fname;
		}
	    }
	} else if (flags & ADD_GFN) {
	    if (maybe_add_suffix(fname, ".gfn")) {
		err = gretl_test_fopen(fname, "r");
		if (!err) {
		    return fname;
		}
	    }
	}
	strcpy(fname, orig);
	if (flags & SUBDIRS) {
	    if (find_in_subdir(topdir, fname, flags)) {
		return fname;
	    }
	}
    }

    return NULL;
}

/**
 * get_plausible_search_dirs:
 * @stype: DATA_SEARCH for data file packages, DB_SEARCH for
 * gretl databases, FUNCS_SEARCH for function packages, or
 * SCRIPT_SEARCH for hansl scripts.
 * @n_dirs: location to receive the number of directories.
 *
 * Returns: an array of plausible search paths, depending on the
 * @type of search. The array should be freed when you are done
 * with it, using strings_array_free().
 */

char **get_plausible_search_dirs (SearchType stype, int *n_dirs)
{
    char **dirs = NULL;
    const char *subdir;
    char dirname[MAXLEN];
    int err = 0;

    *n_dirs = 0;

    if (stype == DATA_SEARCH) {
	subdir = "data";
    } else if (stype == DB_SEARCH) {
	subdir = "db";
    } else if (stype == FUNCS_SEARCH) {
	subdir = "functions";
    } else if (stype == SCRIPT_SEARCH) {
	subdir = "scripts";
    } else {
	fprintf(stderr, "get_plausible_search_dir: no type specified\n");
	return NULL;
    }

    /* system dir first */
    gretl_build_path(dirname, gretl_home(), subdir, NULL);
    err = strings_array_add(&dirs, n_dirs, dirname);

#ifdef OS_OSX
    if (!err) {
	/* the user's ~/Library */
	gretl_build_path(dirname, gretl_app_support_dir(), subdir, NULL);
	err = strings_array_add(&dirs, n_dirs, dirname);
    }
#else
    if (!err) {
	/* the user's dotdir */
	gretl_build_path(dirname, gretl_dotdir(), subdir, NULL);
	err = strings_array_add(&dirs, n_dirs, dirname);
    }
#endif

    if (!err) {
	/* the user's working dir */
	gretl_build_path(dirname, gretl_workdir(), subdir, NULL);
	err = strings_array_add(&dirs, n_dirs, dirname);
    }

    if (!err) {
	/* working dir, no subdir */
	strcpy(dirname, gretl_workdir());
	err = strings_array_add(&dirs, n_dirs, dirname);
    }

    if (!err) {
	/* a legacy thing: some files may have been written to
	   the "default workdir" in the past
	*/
	const char *wd = maybe_get_default_workdir();

	if (wd != NULL) {
	    gretl_build_path(dirname, wd, subdir, NULL);
	    err = strings_array_add(&dirs, n_dirs, dirname);
	    if (!err && stype != FUNCS_SEARCH) {
		strcpy(dirname, wd);
		err = strings_array_add(&dirs, n_dirs, dirname);
	    }
	}
    }

    return dirs;
}

/* it's a dirent thing */
#ifndef NAME_MAX
# define NAME_MAX 255
#endif

/* Search relevant paths for an addon, trying to ensure we
   pick up the most recent version in case there are any
   duplicates.
*/

char *gretl_addon_get_path (const char *name)
{
    char *ret = NULL;
    char path[FILENAME_MAX];
    double version, maxver = 0;
    char **dirs;
    int i, n_dirs, err;

    *path = '\0';

    dirs = get_plausible_search_dirs(FUNCS_SEARCH, &n_dirs);

    gretl_push_c_numeric_locale();

    for (i=0; i<n_dirs; i++) {
	const char *fndir = dirs[i];
	struct dirent *dirent;
	const char *dname;
	DIR *dir;
	int found = 0;

	if ((dir = gretl_opendir(fndir)) == NULL) {
	    continue;
	}

	while ((dirent = readdir(dir)) != NULL && !found) {
	    dname = dirent->d_name;
	    if (!strcmp(dname, name)) {
		gretl_build_path(path, fndir, dname, dname, NULL);
		strcat(path, ".gfn");
		err = gretl_test_fopen(path, "r");
		if (!err) {
		    version = function_package_get_version(path);
		    if (!na(version) && version > maxver) {
			maxver = version;
			free(ret);
			ret = gretl_strdup(path);
		    }
		}
		found = 1;
	    }
	}

	closedir(dir);
    }

    gretl_pop_c_numeric_locale();

    strings_array_free(dirs, n_dirs);

    return ret;
}

static int is_addon (const char *name)
{
    if (!strcmp(name, "gig") ||
	!strcmp(name, "SVAR") ||
	!strcmp(name, "HIP") ||
	!strcmp(name, "ivpanel") ||
	!strcmp(name, "dbnomics")) {
	return 1;
    } else {
	return 0;
    }
}

/**
 * gretl_function_package_get_path:
 * @name: the name of the package to find, without the .gfn extension.
 * @type: %PKG_SUBDIR for a package that lives in its own subdirectory,
 * %PKG_TOPLEV for a package file not in a subdirectory, or %PKG_ALL
 * for a package that may be at either level.
 *
 * Searches a list of directories in which we might expect to find
 * function packages, and, if the package in question is found,
 * returns a newly allocated string holding the full path to
 * the package's gfn file. Public (system) directories are
 * searched first, then directories in the user's filespace.
 *
 * Returns: allocated path on success, otherwise NULL.
 */

char *gretl_function_package_get_path (const char *name,
				       PkgType type)
{
    char *ret = NULL;
    char path[FILENAME_MAX];
    char **dirs;
    int err, found = 0;
    int i, n_dirs;

    if (type == PKG_ALL && is_addon(name)) {
	type = PKG_SUBDIR;
    }

    if (type == PKG_SUBDIR) {
	return gretl_addon_get_path(name);
    }

    *path = '\0';

    dirs = get_plausible_search_dirs(FUNCS_SEARCH, &n_dirs);

    for (i=0; i<n_dirs && !found; i++) {
	const char *fndir = dirs[i];
	struct dirent *dirent;
	const char *dname;
	char *p, test[NAME_MAX+1];
	DIR *dir;

	if ((dir = gretl_opendir(fndir)) == NULL) {
	    continue;
	}

	if (type != PKG_TOPLEV) {
	    /* look preferentially for .gfn files in their own
	       subdirectories */
	    while ((dirent = readdir(dir)) != NULL && !found) {
		dname = dirent->d_name;
		if (!strcmp(dname, name)) {
		    sprintf(path, "%s%c%s%c%s.gfn", fndir, SLASH,
			    dname, SLASH, dname);
		    err = gretl_test_fopen(path, "r");
		    if (!err) {
			found = 1;
		    } else {
			*path = '\0';
		    }
		}
	    }
	}

	if (!found && type != PKG_SUBDIR) {
	    /* look for .gfn files in the top-level functions
	       directory */
	    rewinddir(dir);
	    while ((dirent = readdir(dir)) != NULL && !found) {
		dname = dirent->d_name;
		if (has_suffix(dname, ".gfn")) {
		    strcpy(test, dname);
		    p = strrchr(test, '.');
		    *p = '\0';
		    if (!strcmp(test, name)) {
			sprintf(path, "%s%c%s", fndir, SLASH, dname);
			found = 1;
		    }
		}
	    }
	}

	closedir(dir);
    }

    strings_array_free(dirs, n_dirs);

    if (*path != '\0') {
	ret = gretl_strdup(path);
    }

    return ret;
}

/* Search for file with basename @fname in directory
   @dirname, or in any subdirectory of same up to
   depth @maxdepth. If the file is found, write its
   path into @fullname and return 1, otherwise
   return 0 (= not found).
*/

static int find_file_in_dir (const char *fname,
			     const char *dirname,
			     char *fullname,
			     int maxdepth,
			     int depth)
{
    char tmp[FILENAME_MAX];
    struct dirent *entry;
    struct stat sbuf;
    DIR *dir;
    int found = 0;

    dir = gretl_opendir(dirname);

    if (dir == NULL) {
	return 0;
    }

    /* look for top-level plain file first */
    while (!found && (entry = readdir(dir))) {
	if (!strcmp(entry->d_name, ".") ||
	    !strcmp(entry->d_name, "..")) {
	    continue;
	}
	sprintf(tmp, "%s%c%s", dirname, SLASH, entry->d_name);
	if (gretl_stat(tmp, &sbuf) < 0) {
	    continue;
	} else if ((sbuf.st_mode & S_IFREG) &&
		   !strcmp(entry->d_name, fname)) {
	    strcpy(fullname, tmp);
	    found = 1;
	}
    }

    if (!found && depth < maxdepth) {
	/* then look in subdirs */
	rewinddir(dir);
	depth++;
	while (!found && (entry = readdir(dir))) {
	    if (!strcmp(entry->d_name, ".") ||
		!strcmp(entry->d_name, "..")) {
		continue;
	    }
	    sprintf(tmp, "%s%c%s", dirname, SLASH, entry->d_name);
	    if (gretl_stat(tmp, &sbuf) < 0) {
		continue;
	    } else if (sbuf.st_mode & S_IFDIR) {
		found = find_file_in_dir(fname, tmp, fullname,
					 maxdepth, depth);
	    }
	}
    }

    closedir(dir);

    return found;
}

/**
 * get_package_data_path:
 * @fname: the basename of the file whose full path is wanted.
 * @fullname: location to which the full path should be written
 * (should be at least FILENAME_MAX bytes).
 *
 * Looks for @fname in association with the name of a function
 * package, which must have been set previously using the
 * --frompkg option with the "open" command.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int get_package_data_path (const char *fname, char *fullname)
{
    const char *pkgname;
    int err = 0;

    *fullname = '\0';
    pkgname = get_optval_string(OPEN, OPT_K);

    if (pkgname == NULL) {
	err = E_DATA;
    } else {
	const char *ppath;
	char *gfnpath;

	ppath = get_function_package_path_by_name(pkgname);

	if (ppath != NULL) {
	    gfnpath = gretl_strdup(ppath);
	} else {
	    gfnpath = gretl_addon_get_path(pkgname);
	}

	if (gfnpath == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't find package %s"),
				 pkgname);
	    err = E_DATA;
	} else {
	    char *p = strrchr(gfnpath, SLASH);
	    const char *needle;

	    if (p != NULL) {
		*p = '\0';
	    }

	    /* trim path from @fname if present */
	    needle = strrchr(fname, SLASH);
	    if (needle != NULL) {
		needle++;
	    } else {
		needle = fname;
	    }

	    if (!find_file_in_dir(needle, gfnpath, fullname, 1, 0)) {
		gretl_errmsg_sprintf(_("Couldn't find file %s for package %s"),
				     needle, pkgname);
		*fullname = '\0';
		err = E_FOPEN;
	    }
	    free(gfnpath);
	}
    }

    return err;
}

/**
 * gretl_addpath:
 * @fname: on input, the initially given file name; on output
 * a path may be prepended and/or a suffix may be appended.
 * This variable must be of size at least #MAXLEN bytes to allow
 * for possible additions.
 * @script: if non-zero, assume the file we're looking for
 * is a hansl script.
 *
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened. Usually called by get_full_filename().
 * If @fname does not already have a dot-extension we may also try adding
 * an appropriate gretl extension in case no file is found.
 *
 * Returns: the path to the file that was found (in @fname), or
 * NULL if no file could be found even allowing for prepending
 * a path and/or adding a suffix.
 */

char *gretl_addpath (char *fname, int script)
{
    char orig[MAXLEN];
    char *test;
    int err;

#if SEARCH_DEBUG
    fprintf(stderr, "gretl_addpath: fname='%s', script=%d\n",
	    fname, script);
#endif

    /* keep a backup of the original input */
    strcpy(orig, fname);

    if (g_path_is_absolute(fname)) {
	err = gretl_test_fopen(fname, "r");
	if (err && !script && maybe_add_suffix(fname, ".gdt")) {
	    err = gretl_test_fopen(fname, "r");
	    if (err) {
		strcpy(fname, orig);
	    }
	}
	return err ? NULL : fname;
    }

    /* try workdir first */
    gretl_build_path(fname, paths.workdir, orig, NULL);
    err = gretl_test_fopen(fname, "r");

    if (!err) {
	return fname;
    } else {
	const char *gpath = gretl_current_dir();
	char trydir[MAXLEN];

	strcpy(fname, orig);

	if (*gpath != '\0') {
	    /* try looking where the last-opened script was found */
	    int flags = script ? 0 : ADD_GDT;

	    test = search_dir(fname, gpath, flags);
	    if (test != NULL) {
		return fname;
	    }
	}

	strcpy(fname, orig);

	/* now try gretl installation dir */
	gpath = gretl_home();

	if (*gpath != '\0') {
	    /* try searching some standard gretl paths */
	    if (script) {
		sprintf(trydir, "%sscripts", gpath);
		test = search_dir(fname, trydir, SUBDIRS);
		if (test != NULL) {
		    return fname;
		} else {
		    strcpy(fname, orig);
		    sprintf(trydir, "%sfunctions", gpath);
		    test = search_dir(fname, trydir, ADD_GFN | SUBDIRS);
		    if (test != NULL) {
			return fname;
		    }
		}
	    } else if (has_suffix(fname, ".bin")) {
		/* database? */
		sprintf(trydir, "%sdb", gpath);
		test = search_dir(fname, trydir, 0);
		if (test != NULL) {
		    return fname;
		}
	    } else {
		/* data file */
		sprintf(trydir, "%sdata", gpath);
		test = search_dir(fname, trydir, ADD_GDT | SUBDIRS);
		if (test != NULL) {
		    return fname;
		}
	    }
	}

	strcpy(fname, orig);

	/* now try user's personal filespace */
#ifdef OS_OSX
	gpath = gretl_app_support_dir();
#else
	gpath = gretl_dotdir();
#endif

	if (*gpath != '\0') {
	    /* try looking in ~/Library or dotdir */
	    if (script) {
		sprintf(trydir, "%sscripts", gpath);
		test = search_dir(fname, trydir, SUBDIRS);
		if (test != NULL) {
		    return fname;
		} else {
		    strcpy(fname, orig);
		    sprintf(trydir, "%sfunctions", gpath);
		    test = search_dir(fname, trydir, ADD_GFN | SUBDIRS);
		    if (test != NULL) {
			return fname;
		    }
		}
	    } else if (has_suffix(fname, ".bin")) {
		/* database? */
		sprintf(trydir, "%sdb", gpath);
		test = search_dir(fname, trydir, 0);
		if (test != NULL) {
		    return fname;
		}
	    } else {
		/* data file? */
		sprintf(trydir, "%sdata", gpath);
		test = search_dir(fname, trydir, ADD_GDT | SUBDIRS);
		if (test != NULL) {
		    return fname;
		}
	    }
	}

	strcpy(fname, orig);
	gpath = gretl_workdir();

	if (*gpath != '\0') {
	    /* try looking in user's dir (and subdirs) */
	    test = search_dir(fname, gpath, SUBDIRS);
	    if (test != NULL) {
		return fname;
	    }
	}

	strcpy(fname, orig);
	gpath = maybe_get_default_workdir();

	if (gpath != NULL && *gpath != '\0') {
	    /* try looking in default workdir? */
	    test = search_dir(fname, gpath, SUBDIRS);
	    if (test != NULL) {
		return fname;
	    }
	}
    }

#ifdef WIN32
    /* try looking on the desktop? */
    if (1) {
	char *dtdir = desktop_path();

	strcpy(fname, orig);
	if (dtdir != NULL) {
	    test = search_dir(fname, dtdir, 0);
	    free(dtdir);
	}
	if (test != NULL) {
	    return fname;
	}
    }
#endif

    strcpy(fname, orig);

    gretl_error_clear();

    return NULL;
}

/* It is assumed here that @fname starts with "~/" */

char *gretl_prepend_homedir (const char *fname, int *err)
{
    char *homedir = getenv("HOME");
    char *ret = NULL;

    if (homedir != NULL) {
	ret = malloc(strlen(homedir) + strlen(fname));
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    strcpy(ret, homedir);
	    strcat(ret, fname + 1);
	}
    } else {
	*err = E_DATA;
    }

    return ret;
}

static int substitute_homedir (char *fname)
{
    char *homedir = getenv("HOME");
    int err = 0;

    if (homedir != NULL) {
	int len = strlen(fname);
	int homelen = strlen(homedir);

	if (len + homelen > MAXLEN) {
	    err = 1;
	} else {
	    char tmp[MAXLEN];

	    strcpy(tmp, homedir);
	    strcat(tmp, fname + 1);
	    strcpy(fname, tmp);
	}
    }

    return err;
}

static int get_gfn_special (char *fname)
{
    int ok = 0;

    if (!strchr(fname, '/') && !strchr(fname, '\\')) {
	char *p, pkgname[64];
	char *pkgpath;

	*pkgname = '\0';
	strncat(pkgname, fname, 63);
	p = strstr(pkgname, ".gfn");
	*p = '\0';
	pkgpath = gretl_function_package_get_path(pkgname, PKG_ALL);
	if (pkgpath != NULL) {
	    strcpy(fname, pkgpath);
	    free(pkgpath);
	    ok = 1;
	}
    }

    return ok;
}

/**
 * get_full_filename:
 * @fname: input filename.
 * @fullname: filename to be filled out: must be at least #MAXLEN bytes.
 * @opt: if OPT_S, treat as a script; if OPT_I we're responding
 * to the "include" command; if OPT_W, pass @fname through as is.
 *
 * Includes elementary path-searching: try adding various paths to the
 * given @fname, if appropriate, and see if it can be opened. For
 * internal gretl use.
 *
 * Returns: 0 on success, non-zero on error.
 */

int get_full_filename (const char *fname, char *fullname, gretlopt opt)
{
    int script = (opt & (OPT_S | OPT_I))? 1 : 0;
    char *test = NULL;
    int err = 0;

    *fullname = '\0';

    if (fname == NULL || *fname == '\0') {
	return E_DATA;
    }

    strncat(fullname, fname, MAXLEN - 1);

    if (opt & OPT_W) {
	/* remote database: just use original name */
	return 0;
    }

    if (fullname[0] == '~' && fullname[1] == '/') {
	/* handle tilde == HOME */
	substitute_homedir(fullname);
    }

    if (g_path_is_absolute(fullname)) {
	goto test_open;
    }

    if (opt & OPT_I) {
	/* respect special "include" setting if present */
	char *ipath = getenv("GRETL_INCLUDE_PATH");

	if (ipath != NULL && *ipath != '\0') {
	    gretl_build_path(fullname, ipath, fname, NULL);
	    goto test_open;
	}
    }

    if (has_suffix(fullname, ".gfn") && get_gfn_special(fullname)) {
	/* checked for existence */
	return 0;
    }

    /* try a basic path search on this filename */
#if SEARCH_DEBUG
    fprintf(stderr, "get_full_filename: calling addpath on '%s'\n", fullname);
#endif
    test = gretl_addpath(fullname, script);
#if SEARCH_DEBUG
    fprintf(stderr, "get_full_filename: after: '%s'\n", fullname);
#endif

    /* If @test is non-NULL that means we actually found
       the file somewhere */

    if (test != NULL && (opt & OPT_S)) {
	/* set the 'current_dir' based on script name */
	int spos = gretl_slashpos(fullname);

	if (spos) {
	    *current_dir = '\0';
	    strncat(current_dir, fullname, spos + 1);
	} else {
	    current_dir[0] = '.';
	    current_dir[1] = SLASH;
	    current_dir[2] = '\0';
	}
    }

 test_open:

    if (!err && test == NULL) {
	err = gretl_test_fopen(fullname, "r");
	if (err) {
	    /* ensure we return a gretl error code */
	    err = E_FOPEN;
	}
    }

    return err;
}

int has_system_prefix (const char *fname, SearchType stype)
{
    const char *gretldir = gretl_home();
    int n = strlen(gretldir);
    int ret = 0;

    if (strlen(fname) < n) {
	return 0;
    }

    if (!strncmp(fname, gretldir, n)) {
	if (stype == DATA_SEARCH &&
	    !strncmp(fname + n, "data", 4)) {
	    ret = 1;
	} else if (stype == SCRIPT_SEARCH &&
		   !strncmp(fname + n, "scripts", 7)) {
	    ret = 1;
	}
    }

    return ret;
}

enum paths_status_flags {
    STRING_TABLE_WRITTEN = 1 << 0
};

static void set_gretl_plugpath (const char *path)
{
    *paths.plugpath = '\0';

#if defined(WIN32) && defined(PKGBUILD)
    strcpy(paths.plugpath, path);
    strcat(paths.plugpath, PLUGIN_SFX);
#else
# ifdef LIBDIR
    /* respect the libdir set at compile time, e.g. /usr/lib or
       /usr/lib64
    */
    gretl_build_path(paths.plugpath, LIBDIR, PLUGIN_SFX, NULL);
# else
#  ifdef WIN32
    char *p = strstr(path, "\\share");
#  else
    char *p = strstr(path, "/share");
#  endif

    if (p != NULL) {
	/* back up a level */
	size_t len = p - path;

	strncat(paths.plugpath, path, len);
    } else {
	strcpy(paths.plugpath, path);
    }
    slash_terminate(paths.plugpath);
    strcat(paths.plugpath, "lib");
    slash_terminate(paths.plugpath);
    strcat(paths.plugpath, PLUGIN_SFX);
# endif /* !LIBDIR */
#endif /* !WIN32 */

    slash_terminate(paths.plugpath);
}

static void set_gretl_binbase (const char *path)
{
    sprintf(paths.binbase, "%sdb", path);
}

/* This should be called after we're fairly confident that we
   have the dotdir setting right */

static int set_extra_dot_paths (void)
{
    char dirname[MAXLEN+128];
    size_t n;
    int err = 0;

    /* the personal function package directory */
    *dirname = '\0';
    gretl_build_path(dirname, paths.dotdir, "functions", NULL);
    gretl_mkdir(dirname);

    *paths.tramodir = '\0';
    *paths.x12adir = '\0';

#if !defined(HAVE_TRAMO) && !defined(HAVE_X12A)
    return 0;
#endif

    *dirname = '\0';
    strcpy(dirname, paths.dotdir);
    n = strlen(dirname);

    if (n > 0 && (dirname[n-1] == '\\' || dirname[n-1] == '/')) {
	dirname[n-1] = '\0';
    }

#ifdef HAVE_X12A
    gretl_build_path(paths.x12adir, paths.dotdir, "x12arima", NULL);
    err = gretl_mkdir(paths.x12adir);
    if (err) {
	*paths.x12adir = '\0';
    }
#endif

#ifdef HAVE_TRAMO
    gretl_build_path(paths.tramodir, paths.dotdir, "tramo", NULL);
    if (gretl_mkdir(paths.tramodir)) {
	*paths.tramodir = '\0';
	return 1;
    }

    sprintf(dirname, "%s%coutput", paths.tramodir, SLASH);
    gretl_mkdir(dirname);

    sprintf(dirname, "%s%cgraph", paths.tramodir, SLASH);
    if (gretl_mkdir(dirname)) {
	*paths.tramodir = '\0';
	return 1;
    }

    sprintf(dirname, "%s%cgraph%cacf", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cfilters", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cforecast", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cseries", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cspectra", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
#endif

    return err;
}

static void set_builtin_path_strings (int update)
{
    gretl_insert_builtin_string("gretldir", paths.gretldir);
    gretl_insert_builtin_string("gnuplot",  paths.gnuplot);
    gretl_insert_builtin_string("x12a",     paths.x12a);
    gretl_insert_builtin_string("tramo",    paths.tramo);

    if (!update) {
	/* these only have to be done once */
	gretl_insert_builtin_string("dotdir",   paths.dotdir);
	gretl_insert_builtin_string("workdir",  paths.workdir);
	gretl_insert_builtin_string("x12adir",  paths.x12adir);
	gretl_insert_builtin_string("tramodir", paths.tramodir);
    }

    if (*paths.tramo != '\0') {
	char s[MAXLEN];
	int n;

	*s = '\0';
	strncat(s, paths.tramo, MAXLEN - 1);
	n = strlen(s);
#ifdef WIN32
	if (n >= 9 && !strcmp(s + n - 9, "tramo.exe")) {
	    strcpy(s + n - 9, "seats.exe");
	    gretl_insert_builtin_string("seats", s);
	    return;
	}
#else
	if (n >= 5 && !strcmp(s + n - 5, "tramo")) {
	    strcpy(s + n - 5, "seats");
	    gretl_insert_builtin_string("seats", s);
	}
#endif
    }
}

const char *gretl_home (void)
{
    return paths.gretldir;
}

const char *gretl_bindir (void)
{
    static char bindir[MAXLEN];

    if (*bindir == '\0') {
	char *p;

	strcpy(bindir, paths.gretldir);
	p = strstr(bindir, "share/gretl");
	if (p != NULL) {
	    *p = '\0';
	    strcat(p, "bin/");
	}
#ifdef WIN32
	if (p == NULL) {
	    p = strstr(bindir, "share\\gretl");
	    if (p != NULL) {
		*p = '\0';
		strcat(p, "bin\\");
	    }
	}
#endif
    }

    return bindir;
}

const char *gretl_plugin_path (void)
{
    static int set;

    if (!set) {
	char *epath = getenv("GRETL_PLUGIN_PATH");

	if (epath != NULL) {
	    *paths.plugpath = '\0';
	    strncat(paths.plugpath, epath, MAXLEN - 2);
	    slash_terminate(paths.plugpath);
	}

#if defined(LIBDIR) || defined(GRETL_PREFIX)
	/* if blank, try drawing on compiled-in values */
	if (*paths.plugpath == '\0') {
# ifdef LIBDIR
	    strcat(paths.plugpath, LIBDIR);
# else
	    strcat(paths.plugpath, GRETL_PREFIX);
	    slash_terminate(paths.plugpath);
	    strcat(paths.plugpath, "lib");
# endif
	    slash_terminate(paths.plugpath);
	    strcat(paths.plugpath, PLUGIN_SFX);
	    slash_terminate(paths.plugpath);
	}
#endif /* LIBDIR or GRETL_PREFIX defined */
	set = 1;
    }

    return paths.plugpath;
}

const char *gretl_dotdir (void)
{
    return paths.dotdir;
}

char *gretl_make_dotpath (const char *basename)
{
    int n = strlen(paths.dotdir);
    int m = strlen(basename);
    char *path = NULL;

    if (paths.dotdir[n-1] != SLASH) {
	path = calloc(n + m + 2, 1);
	if (path != NULL) {
	    sprintf(path, "%s%c%s", paths.dotdir, SLASH, basename);
	}
    } else {
	path = calloc(n + m + 1, 1);
	if (path != NULL) {
	    sprintf(path, "%s%s", paths.dotdir, basename);
	}
    }

    return path;
}

const char *gretl_workdir (void)
{
    return paths.workdir;
}

#ifdef WIN32

static const char *win32_default_workdir (void)
{
    static char default_workdir[MAXLEN];
    char *base = mydocs_path();
    const char *retval = NULL;

    if (base != NULL) {
	sprintf(default_workdir, "%s\\gretl\\", base);
	if (strcmp(default_workdir, paths.workdir)) {
	    DIR *dir = gretl_opendir(default_workdir);

	    if (dir != NULL) {
		closedir(dir);
		retval = default_workdir;
	    }
	}
	free(base);
    }

    return retval;
}

#else /* !WIN32 */

static const char *regular_default_workdir (void)
{
    static char default_workdir[MAXLEN];
    char *home = getenv("HOME");
    const char *retval = NULL;

    if (home != NULL) {
	sprintf(default_workdir, "%s/gretl/", home);
	if (strcmp(default_workdir, paths.workdir)) {
	    DIR *dir = opendir(default_workdir);

	    if (dir != NULL) {
		closedir(dir);
		retval = default_workdir;
	    }
	}
    }

    return retval;
}

#endif /* WIN32 or not */

/**
 * maybe_get_default_workdir:
 *
 * Figures the full path to the default value of the
 * user's gretl working directory; call this "defdir".
 *
 * If this defdir turns out to be the same as the
 * current gretl working directory, as would be returned
 * by gretl_workdir(), this function returns NULL,
 * otherwise it returns the defdir value.
 *
 * Returns: a path, or NULL.
 */

const char *maybe_get_default_workdir (void)
{
#ifdef WIN32
    return win32_default_workdir();
#else
    return regular_default_workdir();
#endif
}

static int validate_writedir (const char *dirname)
{
    int err = 0;

    if (*dirname == '\0') {
	gretl_errmsg_set(_("User directory is not set"));
	return E_DATA;
    }

    err = gretl_mkdir(dirname);
    if (err) {
	gretl_errmsg_sprintf( _("Couldn't create directory '%s'"), dirname);
    }

    if (!err) {
	/* ensure the directory is writable */
	char testname[FILENAME_MAX];
	FILE *fp;

	gretl_build_path(testname, dirname, "write.chk", NULL);
	fp = gretl_fopen(testname, "w");
	if (fp == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't write to '%s': "
				   "gretl will not work properly!"),
				 dirname);
	    err = E_FOPEN;
	} else {
	    fclose(fp);
	    gretl_remove(testname);
	}
    }

    if (err) {
	set_gretl_alarm(1);
    }

    return err;
}

int set_gretl_work_dir (const char *path)
{
    DIR *test;
    int err = 0;

    errno = 0;

    test = gretl_opendir(path);

    if (test == NULL) {
	gretl_errmsg_set_from_errno(path, errno);
	fprintf(stderr, "set_gretl_work_dir: '%s': failed\n", path);
	err = E_FOPEN;
    } else {
	closedir(test);
	strcpy(paths.workdir, path);
	slash_terminate(paths.workdir);
	gretl_insert_builtin_string("workdir", paths.workdir);
    }

    return err;
}

const char *gretl_gnuplot_path (void)
{
    return paths.gnuplot;
}

const char *gretl_plotfile (void)
{
    return paths.plotfile;
}

char *set_gretl_plotfile (const char *fname)
{
    *paths.plotfile = 0;
    strncat(paths.plotfile, fname, MAXLEN - 1);

    return paths.plotfile;
}

void report_plot_written (PRN *prn)
{
    if (prn != NULL) {
	pprintf(prn, _("wrote %s\n"), paths.plotfile);
    }
}

const char *gretl_binbase (void)
{
    return paths.binbase;
}

const char *gretl_tramo (void)
{
    return paths.tramo;
}

const char *gretl_tramo_dir (void)
{
    return paths.tramodir;
}

const char *gretl_x12_arima (void)
{
    return paths.x12a;
}

const char *gretl_x12_arima_dir (void)
{
    return paths.x12adir;
}

int gretl_x12_is_x13 (void)
{
    return strstr(paths.x12a, "x13") != NULL;
}

#ifdef WIN32

/* try to avoid using a stale value saved to .gretl2rc */

static void R_path_try_registry (int which, char *targ)
{
    char tmp[MAX_PATH];
    int err;

    err = R_path_from_registry(tmp, which);

    if (!err) {
	*targ = '\0';
	strncat(targ, tmp, MAXLEN - 1);
    }
}

#endif

const char *gretl_rbin_path (void)
{
#ifdef WIN32
    static int checked;

    if (!checked) {
	R_path_try_registry(REXE, paths.rbinpath);
	checked = 1;
    }
#endif

#if 0
    fprintf(stderr, "gretl_rbin_path: '%s'\n", paths.rbinpath);
#endif

    return paths.rbinpath;
}

const char *gretl_rlib_path (void)
{
#ifdef WIN32
    static int checked;

    if (!checked) {
	R_path_try_registry(RLIB, paths.rlibpath);
	checked = 1;
    }
#endif

#if 0
    fprintf(stderr, "gretl_rlib_path: '%s'\n", paths.rlibpath);
#endif

    return paths.rlibpath;
}

const char *gretl_oxl_path (void)
{
    return paths.oxlpath;
}

const char *gretl_octave_path (void)
{
    return paths.octpath;
}

const char *gretl_stata_path (void)
{
    return paths.statapath;
}

const char *gretl_python_path (void)
{
    return paths.pypath;
}

const char *gretl_julia_path (void)
{
    return paths.jlpath;
}

const char *gretl_mpi_hosts (void)
{
    return paths.mpi_hosts;
}

const char *gretl_mpiexec (void)
{
    return paths.mpiexec;
}

const char *gretl_current_dir (void)
{
    return current_dir;
}

void gretl_set_current_dir (const char *s)
{
    int spos = gretl_slashpos(s);

    if (spos) {
	*current_dir = '\0';
	strncat(current_dir, s, spos + 1);
    }
}

const char *gretl_png_font (void)
{
    return paths.pngfont;
}

void set_gretl_png_font (const char *s)
{
    strcpy(paths.pngfont, s);
}

void set_string_table_written (void)
{
    paths.status |= STRING_TABLE_WRITTEN;
}

int gretl_string_table_written (void)
{
    int ret = 0;

    if (paths.status & STRING_TABLE_WRITTEN) ret = 1;

    paths.status &= ~STRING_TABLE_WRITTEN;

    return ret;
}

void show_paths (void)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", paths.gretldir);
    printf("workdir: %s\n", paths.workdir);
    printf("dotdir: %s\n", paths.dotdir);
    printf("gnuplot: %s\n", paths.gnuplot);
}

#ifdef WIN32

static char *rightmost (char *s1, char *s2)
{
    if (s1 == NULL) {
	return s2;
    } else if (s2 == NULL) {
	return s1;
    } else {
	return (s2 - s1 > 0)? s2 : s1;
    }
}

static char getsep (const char *s)
{
    int bak = 0, fwd = 0;

    while (*s) {
	if (*s == '\\') {
	    bak++;
	} else if (*s == '/') {
	    fwd++;
	}
	s++;
    }

    return fwd > bak ? '/' : '\\';
}

static char *rslashpos (const char *s)
{
    char *p1 = strrchr(s, '\\');
    char *p2 = strrchr(s, '/');

    return rightmost(p1, p2);
}

static int slash_terminated (const char *s, int n)
{
    return s[n-1] == '\\' || s[n-1] == '/';
}

/* note: @progname is argv[0] at startup; we use
   this apparatus only for MS Windows
*/

void win32_set_gretldir (const char *progname)
{
    int done = 0;

    *paths.gretldir = '\0';

    if (g_path_is_absolute(progname)) {
	strncat(paths.gretldir, progname, MAXLEN - 1);
    } else {
	char *test = getcwd(paths.gretldir, MAXLEN);

	if (test != NULL) {
	    int n = strlen(paths.gretldir);
	    int m = strlen(progname);
	    char sep = getsep(paths.gretldir);

	    if (n + m + 1 < MAXLEN) {
		slash_terminate(paths.gretldir);
		strncat(paths.gretldir, progname, m);
		if (!gretl_file_exists(paths.gretldir)) {
		    /* can't be right! */
		    *paths.gretldir = '\0';
		}
	    }
	}
    }

    if (*paths.gretldir != '\0') {
	char *s = rslashpos(paths.gretldir);

	if (s != NULL) {
	    /* chop off the program-name bit */
	    *(s+1) = '\0';
	    done = 1;
	}
    }

# ifdef PKGBUILD
    if (!done) {
	/* try the registry */
	char tmp[MAXLEN];
	int err;

	err = read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", tmp);
	if (!err) {
	    strcpy(paths.gretldir, tmp);
	    slash_terminate(paths.gretldir);
	}
    }
# else
    /* a non-pkgbuild Windows build */
    if (*paths.gretldir != '\0') {
	/* if we got "path-to\bin\" (with slashes of either sort)
	   we need to convert to "path-to\share\gretl\"
	*/
	int n = strlen(paths.gretldir);
	char *tail = paths.gretldir + n - 5;

	if (*tail == '\\' || *tail == '/' && !strncmp(tail+1, "bin", 3)) {
	    tail[1] = '\0';
	    strncat(tail + 1, "share", 5);
	    slash_terminate(paths.gretldir);
	    strncat(paths.gretldir, "gretl", 5);
	    slash_terminate(paths.gretldir);
	}
    } else {
	fprintf(stderr, "win32_set_gretldir: haven't got gretldir yet!\n");
    }
# endif
}

#else /* !WIN32 */

/* We have paths.gretldir in place: now test it by seeing if we can
   open the the GPL file "COPYING", which definitely should be in that
   directory.  If that doesn't work, try some remedial measures.
   Note, @config_path is the path garnered from the config file,
   which we may or may not have used to write paths.gretldir (and
   which may indeed be an empty string).
*/

static void check_gretldir (char *config_path)
{
    char testname[FILENAME_MAX];
    FILE *fp;
    int gotit = 0;

    sprintf(testname, "%sCOPYING", paths.gretldir);
    fp = gretl_fopen(testname, "r");

    if (fp != NULL) {
	/* should be fine as is */
	fclose(fp);
	gotit = 1;
    } else if (*config_path != '\0') {
	slash_terminate(config_path);
	if (strcmp(config_path, paths.gretldir)) {
	    /* we weren't using the config-file version: try it now */
	    sprintf(testname, "%sCOPYING", config_path);
	    fp = gretl_fopen(testname, "r");
	    if (fp != NULL) {
		strcpy(paths.gretldir, config_path);
		fclose(fp);
		gotit = 1;
	    }
	}
    }

    if (!gotit && !gretl_in_tool_mode()) {
	/* we're messed up; try to recover */
	pid_t pid = getpid();
	gchar *proc_exe;
	const char *s;
	ssize_t nr;

	proc_exe = g_strdup_printf("/proc/%d/exe", pid);
	nr = readlink(proc_exe, testname, FILENAME_MAX - 1);

	if (nr > 0) {
	    testname[nr] = '\0';
	    fprintf(stderr, "gretl is process %d, '%s'\n", (int) pid, testname);
	    /* should be something like /foo/bar/bin/gretl; we
	       want the /foo/bar bit to append to
	    */
	    s = strstr(testname, "bin/gretl");
	    if (s != NULL) {
		*paths.gretldir = '\0';
		strncat(paths.gretldir, testname, s - testname);
		strcat(paths.gretldir, "share/gretl/");
		fprintf(stderr, "gretldir is maybe '%s'?\n",
			paths.gretldir);
	    }
	}

	g_free(proc_exe);
    }
}

#endif

/* Called at start-up only: the @dirname argument is the value taken
   from the config file or registry.  In case we end up using a value
   other than the incoming one, sync back to @dirname.
*/

static void initialize_gretldir (char *dirname, gretlopt opt)
{
    char *ghome = getenv("GRETL_HOME");
    int err = 0;

    if (ghome != NULL) {
	/* environment setting, if any, takes precedence */
	strcpy(paths.gretldir, ghome);
	slash_terminate(paths.gretldir);
    } else if (dirname != NULL && *dirname != '\0' &&
	       *paths.gretldir == '\0') {
	/* use value from config/registry, unless we already got
	   a value somehow */
	strcpy(paths.gretldir, dirname);
	slash_terminate(paths.gretldir);
    }

    if (*paths.gretldir == '\0') {
#ifdef WIN32
	/* fall back on installation-time default */
	char *progfiles = program_files_path();

	sprintf(paths.gretldir, "%s\\gretl\\", progfiles);
	free(progfiles);
#else
	/* use the compile-time value */
	strcpy(paths.gretldir, GRETL_PREFIX);
	strcat(paths.gretldir, "/share/gretl/");
#endif
    }

#ifndef WIN32
    check_gretldir(dirname);
#endif

    if (!err) {
	set_helpfile_option(opt);
	set_gretl_plugpath(paths.gretldir);
	set_gretl_binbase(paths.gretldir);
    }

    strcpy(dirname, paths.gretldir);
}

/**
 * set_gretl_plugin_path:
 * @path: path to the gretl plugins.
 *
 * For use by third-party code: the purpose of this function
 * is to ensure that libgretl can find its plugins.
 *
 * @prefix, if given, should be the path under which the plugins
 * are installed. On a unix-type system this might be, for example,
 * /usr/local/lib/gretl-gtk2; on MS Windows it might be
 * c:\program files\gretl\plugins.
 **/

void set_gretl_plugin_path (const char *path)
{
    if (path != NULL) {
	*paths.plugpath = '\0';
	strncat(paths.plugpath, path, MAXLEN - 2);
	slash_terminate(paths.plugpath);
    }
}

/* Called at start-up only: set the "hidden" working dir,
   which is not user-configurable.
*/

static int initialize_dotdir (void)
{
    char *dirname;
    int err = 0;

    *paths.dotdir = '\0';

#ifdef WIN32
    dirname = appdata_path();
    if (dirname != NULL) {
	sprintf(paths.dotdir, "%s\\gretl\\", dirname);
	free(dirname);
    } else {
	sprintf(paths.dotdir, "%s\\user\\", paths.gretldir);
    }
#else
    dirname = getenv("HOME");
    if (dirname != NULL) {
	sprintf(paths.dotdir, "%s/.gretl/", dirname);
    }
#endif

    err = validate_writedir(paths.dotdir);

    if (err) {
	*paths.x12adir = '\0';
	*paths.tramodir = '\0';
    } else {
	/* these paths depend on dotdir */
	err = set_extra_dot_paths();
    }

    return err;
}

enum {
    PATH_NEEDS_SLASH = 1 << 0,
    PATH_BLANK_OK    = 1 << 1
};

/* Updating a gretl paths element: transcribe the new value unless it
   is unchanged; if it's a directory string that needs to be
   slash-terminated, check that; return 1 if any change was made to
   the internally recorded value, @targ, otherwise return 0.  Note
   that we ignore empty @src unless the PATH_BLANK_OK flag is given.
*/

static int maybe_transcribe_path (char *targ, char *src, int flags)
{
    int ret = 0;

    if (*src == '\0' && (flags & PATH_BLANK_OK)) {
	if (*targ != '\0') {
	    *targ = '\0';
	    ret = 1;
	}
    } else if (*src != '\0') {
	if (flags & PATH_NEEDS_SLASH) {
	    slash_terminate(src);
	}
	if (strcmp(src, targ)) {
	    strcpy(targ, src);
	    ret = 1;
	}
    } else {
	/* back-sync */
	strcpy(src, targ);
    }

    return ret;
}

#define CFG_DEBUG 0

/* gretl_update_paths is called from the GUI preferences dialog. The
   internal path elements that can be set in this way are:

   gretldir
   gnuplot (but not for MS Windows package)
   tramo, x12a, rbinpath, rlibpath, oxlpath, octpath, statapath,
     pypath, jlpath, dbhost

   * paths.workdir is updated via the separate working directory
     dialog

   * paths.pngfont is updated separately via the plot editing
     dialog

   The @opt argument can include OPT_N to force use of the English-
   language help file where this would not be the default.
*/

int gretl_update_paths (ConfigPaths *cpaths, gretlopt opt)
{
    int ndelta = 0;

    if (maybe_transcribe_path(paths.gretldir, cpaths->gretldir,
			      PATH_NEEDS_SLASH)) {
	set_helpfile_option(opt);
	set_gretl_plugpath(paths.gretldir);
	ndelta++;
    }

    /* native databases */
    maybe_transcribe_path(paths.dbhost, cpaths->dbhost,
			  PATH_BLANK_OK);

#if !defined(WIN32) || !defined(PKGBUILD)
    /* gnuplot path: this is set immutably at start-up in the
       gretl for Windows package */
    ndelta += maybe_transcribe_path(paths.gnuplot, cpaths->gnuplot, 0);
#endif

    /* other external programs */
    ndelta += maybe_transcribe_path(paths.x12a, cpaths->x12a, 0);
    ndelta += maybe_transcribe_path(paths.tramo, cpaths->tramo, 0);
    ndelta += maybe_transcribe_path(paths.rbinpath, cpaths->rbinpath, 0);
    ndelta += maybe_transcribe_path(paths.oxlpath, cpaths->oxlpath, 0);
    ndelta += maybe_transcribe_path(paths.octpath, cpaths->octpath, 0);
    ndelta += maybe_transcribe_path(paths.statapath, cpaths->statapath, 0);
    ndelta += maybe_transcribe_path(paths.pypath, cpaths->pypath, 0);
    ndelta += maybe_transcribe_path(paths.jlpath, cpaths->jlpath, 0);

#ifdef HAVE_MPI
    ndelta += maybe_transcribe_path(paths.mpiexec, cpaths->mpiexec, 0);
    ndelta += maybe_transcribe_path(paths.mpi_hosts, cpaths->mpi_hosts,
				    PATH_BLANK_OK);
#endif

#ifdef USE_RLIB
    if (maybe_transcribe_path(paths.rlibpath, cpaths->rlibpath, 0)) {
	gretl_R_reset_error();
	ndelta++;
    }
#endif

    if (ndelta > 0) {
	/* we changed at least one thing that should be
	   recorded in the builtin path strings */
	set_builtin_path_strings(1);
    }

#if CFG_DEBUG
    fprintf(stderr, "gretl_update_paths: ndelta = %d\n", ndelta);
#endif

    return 0;
}

#ifdef WIN32

/* MS Windows variants of defaults for any paths that
   we need that were not found in the Windows registry
   (or network config file).
*/

static void load_default_workdir (char *targ)
{
    char *home = mydocs_path();

    if (home != NULL) {
	sprintf(targ, "%s\\gretl\\", home);
	free(home);
    } else {
	sprintf(targ, "%suser\\", paths.gretldir);
    }
}

static void load_default_path (char *targ)
{
    char *progfiles = NULL;
    char *pfx86 = NULL;

#ifndef PKGBUILD
    if (targ == paths.gnuplot) {
	sprintf(targ, "%swgnuplot.exe", gretl_bindir());
	return;
    }
#endif

    progfiles = program_files_path();
    pfx86 = program_files_x86_path();

    if (targ == paths.workdir) {
	load_default_workdir(targ);
    } else if (targ == paths.dbhost) {
	strcpy(targ, "ricardo.ecn.wfu.edu");
    } else if (targ == paths.x12a) {
	sprintf(targ, "%s\\x13as\\x13as.exe", progfiles);
    } else if (targ == paths.tramo) {
	sprintf(targ, "%s\\tramo\\tramo.exe", pfx86);
    } else if (targ == paths.rbinpath) {
	R_path_from_registry(targ, REXE);
    } else if (targ == paths.rlibpath) {
	R_path_from_registry(targ, RLIB);
    } else if (targ == paths.oxlpath) {
	sprintf(targ, "%s\\OxMetrics6\\Ox\\bin\\oxl.exe", progfiles);
    } else if (targ == paths.octpath) {
	strcpy(targ, "C:\\Octave-3.6.4\\bin\\octave.exe");
    } else if (targ == paths.statapath) {
	sprintf(targ, "%s\\Stata\\stata.exe", progfiles);
    } else if (targ == paths.pypath) {
	strcpy(targ, "python.exe");
    } else if (targ == paths.jlpath) {
	strcpy(targ, "julia.exe");
    } else if (targ == paths.mpiexec) {
	strcpy(targ, "mpiexec.exe");
    } else if (targ == paths.mpi_hosts) {
	*targ = '\0';
    } else if (targ == paths.pngfont) {
	if (chinese_locale()) {
	    strcpy(targ, "SimSun 8");
	} else if (japanese_locale()) {
	    strcpy(targ, "Meiryo 8");
	} else {
	    strcpy(targ, "verdana 8");
	}
    }

    free(progfiles);
    free(pfx86);
}

# if CFG_DEBUG

static void show_paths_on_stderr (void)
{
    fprintf(stderr, "after gretl_set_paths:\n");
    fprintf(stderr, " gretldir = '%s'\n", paths.gretldir);
    fprintf(stderr, " workdir = '%s'\n", paths.workdir);
    fprintf(stderr, " dotdir = '%s'\n", paths.dotdir);
    fprintf(stderr, " gnuplot = '%s'\n", paths.gnuplot);
}

# endif

#else /* !WIN32 */

/* unix-type variants of defaults for any paths that we need
   that were not found in the gretl config file.
*/

static void load_default_workdir (char *targ)
{
    char *home = getenv("HOME");

    if (home != NULL) {
	sprintf(targ, "%s/gretl/", home);
    } else {
	home = getenv("GRETL_WORKDIR");
	if (home != NULL) {
	    strcpy(targ, home);
	} else {
	    gretl_path_compose(targ, MAXLEN, paths.gretldir, "user/");
	}
    }
}

static void load_default_path (char *targ)
{
#ifdef OS_OSX
    const char *app_paths[] = {
	"/Applications/OxMetrics6/ox/bin/oxl",
	"/Applications/Octave.app/Contents/Resources/bin/octave",
	"/Applications/Stata/Stata.app/Contents/MacOS/Stata"
    };
#else
    const char *app_paths[] = {
	"oxl",
	"octave",
	"stata"
    };
#endif

    if (targ == paths.workdir) {
	load_default_workdir(targ);
    } else if (targ == paths.dbhost) {
	strcpy(targ, "ricardo.ecn.wfu.edu");
    } else if (targ == paths.gnuplot) {
#if defined(OS_OSX) && defined(PKGBUILD)
	sprintf(targ, "%sgnuplot", gretl_bindir());
#else
	strcpy(targ, "gnuplot");
#endif
    } else if (targ == paths.x12a) {
#ifdef HAVE_X12A
	strcpy(targ, "x12a");
#else
	*targ = '\0';
#endif
    } else if (targ == paths.tramo) {
#ifdef HAVE_TRAMO
	strcpy(targ, "tramo");
#else
	*targ = '\0';
#endif
    } else if (targ == paths.rbinpath) {
	strcpy(paths.rbinpath, "R");
    } else if (targ == paths.rlibpath) {
#ifdef RLIBPATH
	strcpy(paths.rlibpath, RLIBPATH);
#else
	*paths.rlibpath = '\0';
#endif
    } else if (targ == paths.oxlpath) {
	strcpy(paths.oxlpath, app_paths[0]);
    } else if (targ == paths.octpath) {
	strcpy(paths.octpath, app_paths[1]);
    } else if (targ == paths.statapath) {
	strcpy(paths.statapath, app_paths[2]);
    } else if (targ == paths.pypath) {
	strcpy(paths.pypath, "python");
    } else if (targ == paths.jlpath) {
	strcpy(paths.jlpath, "julia");
    } else if (targ == paths.mpiexec) {
#if defined(MAC_NATIVE)
	strcpy(paths.mpiexec, "/opt/openmpi/bin/mpiexec");
#else
	strcpy(paths.mpiexec, "mpiexec");
#endif
    } else if (targ == paths.mpi_hosts) {
	*paths.mpi_hosts = '\0';
    } else if (targ == paths.pngfont) {
#if defined(MAC_NATIVE)
	strcpy(targ, "Sans 10"); /* was 13, why? */
#elif defined(OS_OSX)
	strcpy(targ, "Sans 9");
#else
	if (chinese_locale()) {
	    strcpy(targ, "NSimSun 10");
	} else {
	    strcpy(targ, "Vera 9");
	}
#endif
    }
}

#endif /* WIN32 or not */

int add_slash (char *s)
{
    if (s[strlen(s)-1] != SLASH) {
	strcat(s, SLASHSTR);
	return 1;
    }

    return 0;
}

static void path_init (char *targ, char *src, int needs_slash)
{
    if (*src) {
	strcpy(targ, src);
	if (needs_slash && slash_terminate(targ)) {
	    strcpy(src, targ);
	}
    } else {
	load_default_path(targ);
	strcpy(src, targ);
    }
}

/* Set paths, falling back to defaults if no value has been supplied.
   We do this only at startup.  If the path that we record differs
   from that given in @cpaths, sync the value back to @cpaths
   (via path_init, above).
*/

static void copy_paths_with_fallback (ConfigPaths *cpaths)
{
    /* working directory */
    path_init(paths.workdir, cpaths->workdir, 1);

    /* database server */
    path_init(paths.dbhost, cpaths->dbhost, 0);

    /* gnuplot */
#if defined(WIN32) && defined(PKGBUILD)
    /* "hard-wired" case for Windows package */
    sprintf(paths.gnuplot, "%swgnuplot.exe", paths.gretldir);
#else
    path_init(paths.gnuplot, cpaths->gnuplot, 0);
#endif

    /* other external programs */
    path_init(paths.x12a, cpaths->x12a, 0);
    path_init(paths.tramo, cpaths->tramo, 0);
    path_init(paths.rbinpath, cpaths->rbinpath, 0);
    path_init(paths.rlibpath, cpaths->rlibpath, 0);
    path_init(paths.oxlpath, cpaths->oxlpath, 0);
    path_init(paths.octpath, cpaths->octpath, 0);
    path_init(paths.statapath, cpaths->statapath, 0);
    path_init(paths.pypath, cpaths->pypath, 0);
    path_init(paths.jlpath, cpaths->jlpath, 0);
    path_init(paths.mpiexec, cpaths->mpiexec, 0);
    path_init(paths.mpi_hosts, cpaths->mpi_hosts, 0);

    /* graphing font */
    path_init(paths.pngfont, cpaths->pngfont, 0);
}

/* This is called after reading the gretl config file at startup
   (and only then).  Subsequent updates to paths via the GUI (if any)
   are handled by the function gretl_update_paths().

   The no_dotdir member of cpaths is used only when gretlcli is
   operating in "slave" mode (e.g. under a webserver). It forces gretl
   to use paths.workdir as the "dotdir" rather than using a directory
   under the executing user's HOME.  See
   http://gretl.sourceforge.net/slave/
*/

int gretl_set_paths (ConfigPaths *cpaths)
{
    int err0 = 0, err1 = 0;
    int retval = 0;

    *current_dir = '\0';
    *paths.workdir = '\0';
    *paths.plotfile = '\0';

    initialize_gretldir(cpaths->gretldir, OPT_NONE);

    if (!cpaths->no_dotdir) {
	err0 = initialize_dotdir();
    }

    copy_paths_with_fallback(cpaths);

    if (cpaths->no_dotdir) {
	strcpy(paths.dotdir, paths.workdir);
    }

    if (strcmp(paths.dotdir, paths.workdir)) {
	err1 = validate_writedir(paths.workdir);
	if (err1) {
	    /* try falling back on the default working dir */
	    const char *defpath = maybe_get_default_workdir();

	    if (defpath != NULL && *defpath != '\0' &&
		strcmp(defpath, paths.workdir)) {
		err1 = validate_writedir(defpath);
		if (err1 == 0) {
		    strcpy(paths.workdir, defpath);
		}
	    }
	}
    }

    set_builtin_path_strings(0);
    set_gretl_tex_preamble();

    retval = (err0)? err0 : err1;

#if CFG_DEBUG
    fprintf(stderr, "gretl_set_paths: returning %d\n", retval);
# ifdef WIN32
    show_paths_on_stderr();
# endif
#endif

    return retval;
}

/* For writing a file, name given by user: if the path is not
   absolute, switch to the gretl "workdir" unless @fname begins
   with '~' in which case we switch to the user's HOME.
*/

const char *gretl_maybe_switch_dir (const char *fname)
{
    if (fname[0] == '~' && fname[1] == '/') {
	char *home = getenv("HOME");

	if (home != NULL && chdir(home) == 0) {
	    fname += 2;
	}
    } else if (!g_path_is_absolute(fname)) {
	gretl_chdir(paths.workdir);
    }

    return fname;
}

/**
 * gretl_maybe_prepend_dir:
 * @fname: the original filename, which should be in a
 * location of length FILENAME_MAX.
 *
 * If @fname starts with the construction "~/" to indicate
 * the user's HOME, replace this with the full path to that
 * directory.  Otherwise, if @fname is not already an
 * absolute path, prepend the user's gretl working directory.
 * Otherwise do nothing.
 *
 * Returns: the possibly modified filename.
 */

char *gretl_maybe_prepend_dir (char *fname)
{
    char tmp[FILENAME_MAX];

    *tmp = '\0';

    if (fname[0] == '~' && fname[1] == '/') {
	char *home = getenv("HOME");

	if (home != NULL) {
	    gretl_build_path(tmp, home, fname + 2, NULL);
	}
    } else if (!g_path_is_absolute(fname)) {
	gretl_build_path(tmp, paths.workdir, fname, NULL);
    }

    if (*tmp != '\0') {
	strcpy(fname, tmp);
    }

    return fname;
}

/**
 * gretl_read_user_file:
 * @fname: name of file to open.
 *
 * Attempts to open @fname in read-only mode.  If the file
 * is not found when the name is used "as is", we use
 * gretl_maybe_prepend_dir() to prepend the user's gretl
 * working directory and try again.
 *
 * Returns: file pointer, or NULL on failure.
 */

FILE *gretl_read_user_file (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	char fullname[FILENAME_MAX];

	strcpy(fullname, fname);
	gretl_maybe_prepend_dir(fullname);
	if (*fullname != '\0') {
	    fp = gretl_fopen(fullname, "r");
	}
    }

    return fp;
}

/* remove '.' and '..' from @path */

int gretl_normalize_path (char *path)
{
    char tmp[FILENAME_MAX];
    char *pcpy, *pbit, *s = path;
    char **S, **P = NULL;
    int i, n;
    int err = 0;

    if (*path == '\0' || strstr(path, SLASHSTR) == NULL) {
	/* no-op */
	return 0;
    }

    if (*path == '.') {
	/* absolutize the path first, if necessary */
	char dirname[FILENAME_MAX];

	if (getcwd(dirname, FILENAME_MAX - 1) != NULL) {
	    char *tmp2 = gretl_strdup(path);

	    gretl_build_path(path, dirname, tmp2, NULL);
	    free(tmp2);
	}
    }

    pcpy = gretl_strdup(path);
    if (pcpy == NULL) {
	return E_ALLOC;
    }

    *tmp = '\0';
    s = pcpy;

#ifdef WIN32
    /* may be ok for a filename to start with a double backslash */
    if (!strncmp(path, "\\\\", 2)) {
	strcpy(tmp, SLASHSTR);
	s++;
    } else if (*path && path[1] == ':') {
	strncat(tmp, path, 2);
	s += 2;
    }
#endif

    /* split string s on the path separator and cumulate
       the pieces in array P, skipping any pieces which
       are just "." */

    n = 0;
    while ((pbit = strtok(s, SLASHSTR)) != NULL && !err) {
	if (strcmp(pbit, ".")) {
	    S = realloc(P, (n+1) * sizeof *P);
	    if (S == NULL) {
		err = E_ALLOC;
	    } else {
		P = S;
		P[n++] = pbit;
	    }
	}
	s = NULL; /* for next strtok call */
    }

    if (!err) {
	int j;

	/* let each ".." annihilate the preceding path chunk */

	for (i=n-1; i>0; i--) {
	    if (P[i] != NULL && !strcmp(P[i], "..")) {
		for (j=i-1; j>0; j--) {
		    if (P[j] != NULL && strcmp(P[j], "..")) {
			P[j] = NULL;
			break;
		    }
		}
	    }
	}

	/* re-assemble the path */

	for (i=0; i<n; i++) {
	    if (P[i] != NULL && strcmp(P[i], "..")) {
		strcat(tmp, SLASHSTR);
		strcat(tmp, P[i]);
	    }
	}

	strcpy(path, tmp);
    }

    free(P);
    free(pcpy);

    return err;
}

/**
 * slash_terminate:
 * @path: path string.
 *
 * Check whether @path is already slash-terminated, and if
 * not, append a #SLASH; @path should be a large enough
 * array to accept an extra byte.
 *
 * Returns: 1 if a slash was appended, otherwise 0.
 */

int slash_terminate (char *path)
{
#ifdef WIN32
    if (path != NULL && *path != '\0') {
	int n = strlen(path);

	if (path[n-1] != '\\' && path[n-1] != '/') {
	    char sep = getsep(path);

	    strcat(path, sep == '/' ? "/" : "\\");
	    return 1;
	}
    }
#else
    if (path != NULL && *path != '\0') {
	if (path[strlen(path) - 1] != '/') {
	    strcat(path, "/");
	    return 1;
	}
    }
#endif

    return 0;
}

static void rc_set_gp_colors (const char *gpcolors)
{
    const char *s = gpcolors;
    char cstr[N_GP_COLORS][8];
    int i, nc = 0;

    for (i=0; i<N_GP_COLORS; i++) {
	if (sscanf(s, "%7s", cstr[i]) == 1) {
	    nc++;
	    s += 7;
	    if (*s == ' ') {
		s++;
	    } else {
		break;
	    }
	} else {
	    *cstr[i] = '\0';
	    break;
	}
    }

    if (nc == 4) {
	/* old-style */
	for (i=0; i<3; i++) {
	    set_graph_palette_from_string(i, cstr[i]);
	}
	set_graph_palette_from_string(BOXCOLOR, cstr[3]);
    } else {
	for (i=0; i<nc; i++) {
	    set_graph_palette_from_string(i, cstr[i]);
	}
    }
}

static int rc_bool (const char *s)
{
    if (!strcmp(s, "true") || !strcmp(s, "1")) {
	return 1;
    } else {
	return 0;
    }
}

static void handle_use_cwd (int use_cwd, ConfigPaths *cpaths)
{
    libset_set_bool(USE_CWD, use_cwd);

    if (use_cwd) {
	char *s, cwd[MAXLEN];

	s = getcwd(cwd, MAXLEN);
	if (s != NULL) {
	    *cpaths->workdir = '\0';
	    strncat(cpaths->workdir, s, MAXLEN - 2);
	    slash_terminate(cpaths->workdir);
	}
    }
}

#define PROXLEN 64
#define GRETLCLI_USE_CWD 1

/* called only on behalf of gretlcli (for all platforms) */

void get_gretl_config_from_file (FILE *fp, ConfigPaths *cpaths,
				 char *dbproxy, int *use_proxy)
{
    char line[MAXLEN], key[32], val[MAXLEN];

    while (fgets(line, sizeof line, fp) != NULL) {
	if (*line == '#') {
	    continue;
	}
	if (!strncmp(line, "recent", 6)) {
	    /* reached the "recent files" section */
	    break;
	}
	if (sscanf(line, "%s", key) != 1) {
	    continue;
	}
	*val = '\0';
	/* get the string that follows " = " */
	strncat(val, line + strlen(key) + 3, MAXLEN - 1);
	gretl_strstrip(val);
	if (!strcmp(key, "gretldir")) {
	    strncat(cpaths->gretldir, val, MAXLEN - 1);
#ifndef WIN32
	} else if (!strcmp(key, "gnuplot")) {
	    strncat(cpaths->gnuplot, val, MAXLEN - 1);
#endif
	} else if (!strcmp(key, "workdir") || !strcmp(key, "userdir")) {
	    /* "userdir" is a legacy thing */
	    strncat(cpaths->workdir, val, MAXLEN - 1);
	} else if (!strcmp(key, "no_dotdir")) {
	    cpaths->no_dotdir = rc_bool(val);
	} else if (!strcmp(key, "shellok")) {
	    libset_set_bool(SHELL_OK, rc_bool(val));
	} else if (!strcmp(key, "usecwd")) {
#if GRETLCLI_USE_CWD
	    ; /* handled later */
#else
	    handle_use_cwd(rc_bool(val), cpaths);
#endif
	} else if (!strcmp(key, "lcnumeric")) {
	    libset_set_bool(FORCE_DECP, !rc_bool(val));
	} else if (!strcmp(key, "dbhost")) {
	    strncat(cpaths->dbhost, val, 32 - 1);
	} else if (!strcmp(key, "dbproxy")) {
	    strncat(dbproxy, val, PROXLEN - 1);
	} else if (!strcmp(key, "useproxy")) {
	    *use_proxy = rc_bool(val);
	} else if (!strcmp(key, "x12a")) {
	    strncat(cpaths->x12a, val, MAXLEN - 1);
	} else if (!strcmp(key, "tramo")) {
	    strncat(cpaths->tramo, val, MAXLEN - 1);
	} else if (!strcmp(key, "Rbin")) {
	    strncat(cpaths->rbinpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "Rlib")) {
	    strncat(cpaths->rlibpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "ox")) {
	    strncat(cpaths->oxlpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "octave")) {
	    strncat(cpaths->octpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "stata")) {
	    strncat(cpaths->statapath, val, MAXLEN - 1);
	} else if (!strcmp(key, "python")) {
	    strncat(cpaths->pypath, val, MAXLEN - 1);
	} else if (!strcmp(key, "julia")) {
	    strncat(cpaths->jlpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "mpiexec")) {
	    strncat(cpaths->mpiexec, val, MAXLEN - 1);
	} else if (!strcmp(key, "mpi_hosts")) {
	    strncat(cpaths->mpi_hosts, val, MAXLEN - 1);
	} else if (!strcmp(key, "mpi_pref")) {
#ifdef HAVE_MPI
	    set_mpi_variant(val);
#else
	    ;
#endif
	} else if (!strcmp(key, "Png_font")) {
	    strncat(cpaths->pngfont, val, 128 - 1);
	} else if (!strcmp(key, "Gp_colors")) {
	    rc_set_gp_colors(val);
	} else if (!strcmp(key, "HC_xsect")) {
	    set_xsect_hccme(val);
	} else if (!strcmp(key, "HC_tseri")) {
	    set_tseries_hccme(val);
	} else if (!strcmp(key, "HC_panel")) {
	    set_panel_hccme(val);
	} else if (!strcmp(key, "HC_garch")) {
	    set_garch_robust_vcv(val);
	}
    }

#if GRETLCLI_USE_CWD
    /* "workdir" is always the current directory */
    handle_use_cwd(1, cpaths);
#endif
}

#ifndef WIN32

void get_gretl_rc_path (char *rcfile)
{
    char *path = getenv("GRETL_CONFIG_FILE");

    if (path != NULL) {
	*rcfile = '\0';
	strncat(rcfile, path, FILENAME_MAX - 1);
#if 0
	fprintf(stderr, "rcfile from env: '%s'\n", rcfile);
#endif
    } else {
	path = getenv("HOME");
	if (path != NULL) {
	    sprintf(rcfile, "%s/.gretl2rc", path);
	} else {
	    strcpy(rcfile, ".gretl2rc");
	}
    }
}

/* non-Windows read of the gretl configuration file on behalf
   of the CLI program, gretlcli; the Windows variant of this,
   win32_cli_read_rc(), is in gretl_win32.c
*/

int cli_read_rc (void)
{
    ConfigPaths cpaths = {0};
    char rcfile[FILENAME_MAX];
    char dbproxy[PROXLEN] = {0};
    int use_proxy = 0;
    FILE *fp;
    int err = 0;

    get_gretl_rc_path(rcfile);
    fp = gretl_fopen(rcfile, "r");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	get_gretl_config_from_file(fp, &cpaths, dbproxy, &use_proxy);
	fclose(fp);
    }

    if (err) {
	gretl_set_paths(&cpaths);
    } else {
	err = gretl_set_paths(&cpaths);
    }

#ifdef USE_CURL
    gretl_www_init(cpaths.dbhost, dbproxy, use_proxy);
#endif

#if 0
    show_paths();
#endif

    return err;
}

#endif /* !WIN32 */

#ifdef OS_OSX

const char *gretl_app_support_dir (void)
{
    static char suppdir[FILENAME_MAX];

    if (*suppdir == '\0') {
	/* try to ensure that we have a per-user Application
	   Support dir, with appropriate subdirectories
	*/
	const char *home = getenv("HOME");

	if (home == NULL) {
	    fprintf(stderr, "problem: HOME is not defined\n");
	} else {
	    char *p;
	    int err;

	    sprintf(suppdir, "%s/Library/Application Support/gretl/functions",
		    home);
	    p = strrchr(suppdir, '/') + 1;
	    err = gretl_mkdir(suppdir);
	    if (!err) {
		strcpy(p, "data");
		err = gretl_mkdir(suppdir);
	    }
	    if (!err) {
		strcpy(p, "db");
		err = gretl_mkdir(suppdir);
	    }
	    if (err) {
		*suppdir = '\0';
	    } else {
		/* chop off subdir from name */
		*p = '\0';
	    }
	}
    }

    return suppdir;
}

#endif

static int dir_is_writable (const char *dirname)
{
    int ok = 0;

    if (gretl_mkdir(dirname) == 0) {
	gchar *test = g_strdup_printf("%s%c%s", dirname, SLASH, "wtest");

	if (test != NULL) {
	    ok = (gretl_test_fopen(test, "w") == 0);
	    g_free(test);
	}
    }

    return ok;
}

static int get_user_install_path (char *path, const char *subdir)
{
#ifdef OS_OSX
    const char *dirname = gretl_app_support_dir();
#else
    const char *dirname = gretl_dotdir();
#endif
    int err = 0;

    if (dirname == NULL || *dirname == '\0') {
	err = E_FOPEN;
    } else {
	sprintf(path, "%s%s", dirname, subdir);
	err = (dir_is_writable(path) == 0);
    }

    return err;
}

static int get_system_install_path (char *path, const char *subdir)
{
    sprintf(path, "%s%s", gretl_home(), subdir);

    if (dir_is_writable(path)) {
	return 0;
    } else {
	return E_FOPEN;
    }
}

/* get a path that's suitable for writing a function
   package on installation
*/

const char *gretl_function_package_path (void)
{
    static char path[FILENAME_MAX];

    if (*path == '\0') {
	int sys_first = 1;
	int err = 0;

#if defined(OS_OSX)
	/* we prefer writing to ~/Library/Application Support */
	sys_first = 0;
#elif defined(WIN32)
	if (win32_uses_virtual_store()) {
	    /* don't write to virtualized location */
	    sys_first = 0;
	}
#endif

	if (sys_first) {
	    err = get_system_install_path(path, "functions");
	    if (err) {
		err = get_user_install_path(path, "functions");
	    }
	} else {
	    err = get_user_install_path(path, "functions");
	}

	if (err) {
	    *path = '\0';
	} else {
	    slash_terminate(path);
	}
    }

    return path;
}

int gretl_path_compose (char *targ, int len,
			const char *s1,
			const char *s2)
{
    targ[0] = '\0';
    if (strlen(s1) + strlen(s2) >= len) {
	gretl_errmsg_set("filename is too long");
	return E_DATA;
    } else {
	strcpy(targ, s1);
	strcat(targ, s2);
	return 0;
    }
}

/* Code borrowed from Glib (gfileutils.c) and adapted to
   write to an input char * (@targ) instead of returning
   a newly allocated string. The code is also somewhat
   simplified by the assumption that if the platform is
   not MS Windows the directort separator is always '/':
   this is a safe assumptions for the platforms supported
   by gretl.
*/

#ifdef G_OS_WIN32

static void real_build_path_win32 (char *targ,
				   const gchar *first_element,
				   va_list *args)
{
    gboolean is_first = TRUE;
    gboolean have_leading = FALSE;
    const gchar *single_element = NULL;
    const gchar *next_element;
    const gchar *last_trailing = NULL;
    gchar current_separator = '\\';

    next_element = first_element;

    while (1) {
	const gchar *element;
	const gchar *start;
	const gchar *end;

	if (next_element) {
	    element = next_element;
	    next_element = va_arg(*args, gchar *);
	} else {
	    break;
	}

	/* ignore empty elements */
	if (*element == '\0') {
	    continue;
	}

	start = element;
	while (start && (*start == '\\' || *start == '/')) {
	    current_separator = *start;
	    start++;
	}

	end = start + strlen(start);
	while (end >= start + 1 && (end[-1] == '\\' || end[-1] == '/')) {
	    current_separator = end[-1];
	    end--;
	}

	last_trailing = end;
	while (last_trailing >= element + 1 &&
	       (last_trailing[-1] == '\\' || last_trailing[-1] == '/')) {
	    last_trailing--;
	}

	if (!have_leading) {
	    /* If the leading and trailing separator strings are in the
	       same element and overlap, the result is exactly that
	       element
	    */
	    if (last_trailing <= start) {
		single_element = element;
	    }
	    strncat(targ, element, start - element);
	    have_leading = TRUE;
	} else {
	    single_element = NULL;
	}

	if (end == start) {
	    continue;
	}

	if (!is_first) {
	    strncat(targ, &current_separator, 1);
	}
	strncat(targ, start, end - start);
	is_first = FALSE;
    }

    if (single_element) {
	*targ = '\0';
	strcat(targ, single_element);
    } else if (last_trailing) {
	strcat(targ, last_trailing);
    }
}

#else

static void real_build_path (char *targ,
			     const gchar *first_element,
			     va_list *args)
{
    gboolean is_first = TRUE;
    gboolean have_leading = FALSE;
    const gchar *single_element = NULL;
    const gchar *next_element;
    const gchar *last_trailing = NULL;

    next_element = first_element;

    while (1) {
	const gchar *element;
	const gchar *start;
	const gchar *end;

	if (next_element) {
	    element = next_element;
	    next_element = va_arg(*args, gchar *);
	} else {
	    break;
	}

	/* ignore empty elements */
	if (*element == '\0') {
	    continue;
	}

	start = element;
	while (*start == '/') {
	    start++;
	}

	end = start + strlen (start);
	while (end >= start + 1 && end[-1] == '/') {
	    end--;
	}

	last_trailing = end;
	while (last_trailing >= element + 1 && last_trailing[-1] == '/') {
	    last_trailing--;
	}

	if (!have_leading) {
	    /* If the leading and trailing separator strings are in the
	       same element and overlap, the result is exactly that
	       element
	    */
	    if (last_trailing <= start) {
		single_element = element;
	    }
	    strncat(targ, element, start - element);
	    have_leading = TRUE;
	} else {
	    single_element = NULL;
	}

	if (end == start) {
	    continue;
	}

	if (!is_first) {
	    strcat(targ, "/");
	}
	strncat(targ, start, end - start);
	is_first = FALSE;
    }

    if (single_element) {
	*targ = '\0';
	strcat(targ, single_element);
    } else if (last_trailing) {
	strcat(targ, last_trailing);
    }
}

#endif

/**
 * gretl_build_path:
 * @targ: target string to write to (must be pre-allocated).
 * @first_element: first component of path.
 *
 * Writes to @targ a path composed of @first_element
 * plus any additional string arguments supplied before
 * a terminating NULL. An appropriate separator is inserted
 * between the components of the path.
 *
 * Returns: the target string, @targ.
 */

char *gretl_build_path (char *targ, const gchar *first_element, ...)
{
    va_list args;

    *targ = '\0';

    va_start(args, first_element);
#ifdef G_OS_WIN32
    real_build_path_win32(targ, first_element, &args);
#else
    real_build_path(targ, first_element, &args);
#endif
    va_end(args);

    return targ;
}

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

/* gretl_win32.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "gretl_www.h"

#include <glib.h>

#include <windows.h>
#include <shlobj.h>
#include <aclapi.h>

static void win_print_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL); 

    if (buf != NULL) {
	fprintf(stderr, "Windows says: %s\n", (char *) buf);
    }

    LocalFree(buf);
}

/* returns 0 on success */

int read_reg_val (HKEY tree, const char *base, 
		  char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    char regpath[64];
    LONG ret;
    HKEY regkey;
    int err = 0;

    sprintf(regpath, "Software\\%s", base);

    ret = RegOpenKeyEx(tree,      /* handle to open key */
		       regpath,   /* subkey name */
		       0,         /* reserved */
		       KEY_READ,  /* access mask */
		       &regkey    /* key handle */
		       );

    if (ret != ERROR_SUCCESS) {
	fprintf(stderr, "Couldn't read registry path %s\n", regpath);
	win_print_last_error();
        return 1;
    }

    if (RegQueryValueEx(regkey,
			keyname,
			NULL,
			NULL,
			(LPBYTE) keyval,
			&datalen
			) != ERROR_SUCCESS) {
	*keyval = '\0';
	err = 1;
    }

    RegCloseKey(regkey);

    return err;
}

static char netfile[FILENAME_MAX];

const char *get_gretlnet_filename (void)
{
    return (*netfile != '\0')? netfile : NULL;
}

int set_gretlnet_filename (const char *prog)
{
    char *p;
    int i, n;

    strcpy(netfile, prog);
    n = strlen(netfile) - 1;
    p = netfile;

    for (i=n; i>0; i--) {
	if (p[i] == '\\' || p[i] == '/') {
	    strcpy(p + i,  "\\gretlnet.txt");
	    break;
	}
    }

    return 0;
}

static FILE *cli_gretlnet_open (const char *prog)
{
    FILE *fp = NULL;

    set_gretlnet_filename(prog);

    if (*netfile != '\0') {
	fp = gretl_fopen(netfile, "r");
    }

    return fp;
}

static FILE *cli_rcfile_open (void)
{
    char *appdata = appdata_path();
    FILE *fp = NULL;

    if (appdata != NULL) {
	char fname[FILENAME_MAX];

	sprintf(fname, "%s\\gretl\\.gretl2rc", appdata);
	free(appdata);
	fp = fopen(fname, "r");
    }

    return fp;
}

int read_rc_string (FILE *fp, const char *key, char *value)
{
    int ret = 0;

    if (fp != NULL) {
	char line[MAXLEN];
	char keystr[32];
	int n;

	rewind(fp);

	sprintf(keystr, "%s = ", key);
	n = strlen(keystr);

	while (fgets(line, sizeof line, fp) && !ret) {
	    gretl_strstrip(line);
	    if (!strncmp(line, keystr, n)) {
		strcpy(value, line + n);
		ret = 1;
	    }
	}
    }

    return ret;
}

static int cli_read_gretl_var (char *key, char *val,
			       FILE **fp, HKEY tree)
{
    int done = 0;

    *val = '\0';

    if (fp[0] != NULL) {
	done = read_rc_string(fp[0], key, val);
    }

    if (!done && fp[1] != NULL) {
	done = read_rc_string(fp[1], key, val);
    }

    if (!done) {
	read_reg_val(tree, "gretl", key, val);
    }

    return *val != '\0';
}

void cli_read_registry (char *callname)
{
    ConfigPaths cpaths = {0};
    char valstr[MAXLEN];
    char dbproxy[21];
    int done, use_proxy = 0;
    FILE *fp[2];

    fp[0] = cli_gretlnet_open(callname);
    fp[1] = cli_rcfile_open();

    /* gretl installation directory */
    cli_read_gretl_var("gretldir", cpaths.gretldir, fp, HKEY_LOCAL_MACHINE);

    /* user's working directory */
    cli_read_gretl_var("userdir", cpaths.workdir, fp, HKEY_CURRENT_USER);

    /* path to X-12-ARIMA */
    done = read_rc_string(fp[0], "x12a", cpaths.x12a);
    if (!done) {
	read_reg_val(HKEY_LOCAL_MACHINE, "x12arima", "x12a", cpaths.x12a);
    }

    /* path to tramo */
    done = read_rc_string(fp[0], "tramo", cpaths.tramo);
    if (!done) {
	read_reg_val(HKEY_LOCAL_MACHINE, "tramo", "tramo", cpaths.tramo);
    }

    /* path to R binary (non-interactive use) */
    read_rc_string(fp[0], "Rbin", cpaths.rbinpath);

    /* path to R shared library */
    read_rc_string(fp[0], "Rlib", cpaths.rlibpath);

    /* path to oxl */
    read_rc_string(fp[0], "ox", cpaths.oxlpath);

    /* remote database host */
    cli_read_gretl_var("dbhost", cpaths.dbhost, fp, HKEY_CURRENT_USER);

    /* www proxy for reading remote databases */
    cli_read_gretl_var("dbproxy", dbproxy, fp, HKEY_CURRENT_USER);

    /* should a proxy be used? */
    cli_read_gretl_var("useproxy", valstr, fp, HKEY_CURRENT_USER);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	use_proxy = 1;
    } 

    /* do we allow the shell command within gretl? */
    cli_read_gretl_var("shellok", valstr, fp, HKEY_CURRENT_USER);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	libset_set_bool(SHELL_OK, 1);
    } else {
	libset_set_bool(SHELL_OK, 0);
    }

    gretl_set_paths(&cpaths, OPT_NONE);
    gretl_www_init(cpaths.dbhost, dbproxy, use_proxy);

    if (fp[0] != NULL) {
	fclose(fp[0]);
    }

    if (fp[1] != NULL) {
	fclose(fp[1]);
    }
}

void win_show_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL); 

    MessageBox(NULL, (LPCTSTR) buf, "Error", MB_OK | MB_ICONERROR);
    LocalFree(buf);
}

void win_copy_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL);

    gretl_errmsg_set((const char *) buf);
    LocalFree(buf);
}

/* covers the cases of (a) execing a console application
   as "slave" (without opening a console window) and (b)
   execing a GUI app (in fact, just wgnuplot.exe) as slave
*/

static int real_win_run_sync (char *cmdline, const char *currdir,
			      int console_app) 
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 
    DWORD exitcode;
    DWORD flags;
    int ok, err = 0;

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);  

    si.cb = sizeof si;

    if (console_app) {
	flags = CREATE_NO_WINDOW | HIGH_PRIORITY_CLASS;
    } else {
	si.dwFlags = STARTF_USESHOWWINDOW;
	si.wShowWindow = SW_SHOWMINIMIZED;
	flags = HIGH_PRIORITY_CLASS;
    }

    /* zero return means failure */
    ok = CreateProcess(NULL, cmdline, 
		       NULL, NULL, FALSE,
		       flags,
		       NULL, currdir,
		       &si, &pi);

    if (!ok) {
	fprintf(stderr, "win_run_sync: failed command:\n%s\n", cmdline);
	win_copy_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pi.hProcess, INFINITE); 
	if (GetExitCodeProcess(pi.hProcess, &exitcode)) {
	    if (exitcode != 0) {
		gretl_errmsg_sprintf("%s: exit code %d\n", cmdline, 
				     exitcode);
		err = 1;
	    }
	} else {
	    fprintf(stderr, "win_run_sync: no exit code:\n%s\n", cmdline);
	    win_copy_last_error();
	    err = 1;
	}
    }
   
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return err;
}

/**
 * win_run_sync:
 * @cmdline: command line to execute.
 * @currdir: current directory for child process (or NULL to
 * inherit from parent)
 *
 * Run a command synchronously (i.e. block until it is
 * completed) under MS Windows. This is intended for use
 * with "slave" console applications such a latex, dvips,
 * tramo, x12a and so on.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int win_run_sync (char *cmdline, const char *currdir) 
{
    return real_win_run_sync(cmdline, currdir, 1);
}

int gretl_spawn (char *cmdline)
{
    return real_win_run_sync(cmdline, NULL, 0);
}

/* Retrieve various special paths from the bowels of MS
   Windows.  Note that these paths will be in the locale
   encoding, not UTF-8 
*/

static char *win_special_path (int folder)
{
    TCHAR dpath[MAX_PATH];
    LPITEMIDLIST id_list;
    DWORD result;
    LPMALLOC allocator;
    char *ret = NULL;

    if (SHGetSpecialFolderLocation(NULL, folder | CSIDL_FLAG_CREATE, 
				   &id_list) != S_OK) {
	return NULL;
    }

    result = SHGetPathFromIDList(id_list, dpath);

    if (result) {
	ret = gretl_strdup(dpath);
    }

    if (SHGetMalloc(&allocator) == S_OK) {
	allocator->lpVtbl->Free(allocator, id_list);
	allocator->lpVtbl->Release(allocator);
    }

    return ret;
}

char *desktop_path (void)
{
    return win_special_path(CSIDL_DESKTOPDIRECTORY);
}

char *appdata_path (void)
{
    return win_special_path(CSIDL_APPDATA);
}

#if 0 /* not yet? */

char *local_appdata_path (void)
{
    return win_special_path(CSIDL_LOCAL_APPDATA);
}

#endif

char *mydocs_path (void)
{
    return win_special_path(CSIDL_PERSONAL);
}

char *program_files_path (void)
{
    return win_special_path(CSIDL_PROGRAM_FILES);
}

static char *compose_command_line (const char *arg)
{
    CHAR cmddir[MAX_PATH];
    char *cmdline = NULL;
    
    GetSystemDirectory(cmddir, sizeof cmddir);

    if (getenv("SHELLDEBUG")) {
	cmdline = g_strdup_printf("%s\\cmd.exe /k %s", cmddir, arg);
    } else {
	cmdline = g_strdup_printf("%s\\cmd.exe /c %s", cmddir, arg);
    }

    return cmdline;
}

#define BUFSIZE 4096 
 
static int read_from_pipe (HANDLE hwrite, HANDLE hread, 
			   char **sout, PRN *inprn) 
{ 
    DWORD dwread;
    CHAR buf[BUFSIZE];
    PRN *prn;
    int ok;

    if (sout != NULL) {
	prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    } else {
	prn = inprn;
    }

    /* close the write end of the pipe */
    ok = CloseHandle(hwrite);
    
    if (!ok) {
	fputs("Closing handle failed\n", stderr); 
    } else {
	/* read output from the child process: note that the buffer
	   must be NUL-terminated for use with pputs() */
	while (1) { 
	    memset(buf, '\0', BUFSIZE);
	    ok = ReadFile(hread, buf, BUFSIZE-1, &dwread, NULL);
	    if (!ok || dwread == 0) {
		break;
	    }
	    pputs(prn, buf);
	} 
    }

    if (sout != NULL) {
	*sout = gretl_print_steal_buffer(prn);
	gretl_print_destroy(prn);
    }

    return ok;
} 

enum {
    SHELL_RUN,
    PROG_RUN
};

static int 
run_child_with_pipe (const char *arg, HANDLE hwrite, HANDLE hread,
		     int flag) 
{ 
    PROCESS_INFORMATION pinfo; 
    STARTUPINFO sinfo;
    char *cmdline = NULL;
    int ok;

    if (flag == SHELL_RUN) {
	cmdline = compose_command_line(arg);
    } else {
	cmdline = g_strdup(arg);
    }
 
    ZeroMemory(&pinfo, sizeof pinfo);
    ZeroMemory(&sinfo, sizeof sinfo);

    sinfo.cb = sizeof sinfo;
    sinfo.hStdError = hwrite;
    sinfo.hStdOutput = hwrite;
    sinfo.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    sinfo.dwFlags |= (STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW);
    sinfo.wShowWindow = SW_SHOWMINIMIZED;
 
    ok = CreateProcess(NULL, 
		       cmdline,
		       NULL,          /* process security attributes */
		       NULL,          /* primary thread security attributes */
		       TRUE,          /* handles are inherited */
		       CREATE_NO_WINDOW,
		       NULL,          /* use parent's environment */
		       get_shelldir(),          
		       &sinfo,
		       &pinfo);
   
    if (!ok) {
	win_show_last_error();
    } else {
	CloseHandle(pinfo.hProcess);
	CloseHandle(pinfo.hThread);
    }

    g_free(cmdline);

    return ok;
}

static int run_cmd_with_pipes (const char *arg, char **sout, PRN *prn,
			       int flag) 
{ 
    HANDLE hread, hwrite;
    SECURITY_ATTRIBUTES sa; 
    int ok; 
 
    /* set the bInheritHandle flag so pipe handles are inherited */
    sa.nLength = sizeof(SECURITY_ATTRIBUTES); 
    sa.bInheritHandle = TRUE; 
    sa.lpSecurityDescriptor = NULL; 

    /* create pipe for the child process's STDOUT */ 
    ok = CreatePipe(&hread, &hwrite, &sa, 0);

    if (!ok) {
	win_show_last_error();
    } else {
	/* ensure that the read handle to the child process's pipe for 
	   STDOUT is not inherited */
	SetHandleInformation(hread, HANDLE_FLAG_INHERIT, 0);
	ok = run_child_with_pipe(arg, hwrite, hread, flag);
	if (ok) {
	    /* read from child's output pipe */
	    read_from_pipe(hwrite, hread, sout, prn); 
	}
    }
 
    return 0; 
} 

static int run_cmd_wait (const char *cmd, PRN *prn)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    char *cmdline = NULL;
    int ok, err = 0;

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);

    si.cb = sizeof si;
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_SHOWMINIMIZED;

    cmdline = compose_command_line(cmd);

    ok = CreateProcess(NULL, cmdline, 
		       NULL, NULL, FALSE,
		       CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS,
		       NULL, get_shelldir(),
		       &si, &pi);

    if (!ok) {
	win_show_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pi.hProcess, INFINITE);
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
    }

    g_free(cmdline);

    return err;
}

int gretl_win32_grab_output (const char *cmdline, char **sout)
{
    return run_cmd_with_pipes(cmdline, sout, NULL, PROG_RUN);
}

int gretl_shell_grab (const char *arg, char **sout)
{
    return run_cmd_with_pipes(arg, sout, NULL, SHELL_RUN);
}

int gretl_shell (const char *arg, PRN *prn)
{
    UINT winret;
    int async = 0;
    int err = 0;

    if (arg == NULL || *arg == '\0') {
	return 0;
    }

    if (!libset_get_bool(SHELL_OK)) {
	gretl_errmsg_set(_("The shell command is not activated."));
	return 1;
    }

    if (!strncmp(arg, "launch ", 7)) {
	async = 1;
	arg += 7;
    } else if (*arg == '!') {
	arg++;
    }

    arg += strspn(arg, " \t");

    if (async) {
	winret = WinExec(arg, SW_SHOWNORMAL);
	if (winret <= 31) {
	    err = 1;
	}
    } else if (getenv("GRETL_SHELL_NEW")) {
	err = run_cmd_with_pipes(arg, NULL, prn, SHELL_RUN);
    } else {
	err = run_cmd_wait(arg, prn);
    } 

    return err;
}

/* unlike access(), returns 1 on success */

int win32_write_access (char *path)
{
    SID *sid = NULL;
    ACL *dacl = NULL;
    LPTSTR domain = NULL;
    SECURITY_DESCRIPTOR *sd = NULL;
    TRUSTEE t;
    DWORD sidsize = 0, dlen = 0;
    SID_NAME_USE stype;
    ACCESS_MASK amask;
    const char *username;
    int ret, ok = 0, err = 0;

    /* screen for the read-only attribute first */
    if (access(path, W_OK) != 0) {
	return 0;
    }

    username = g_get_user_name();

    /* get the size of the SID and domain */
    LookupAccountName(NULL, username, NULL, &sidsize, 
		      NULL, &dlen, &stype);

    sid = LocalAlloc(0, sidsize);
    domain = LocalAlloc(0, dlen * sizeof *domain);
    if (sid == NULL || domain == NULL) {
	err = 1;
    } 

    if (!err) {
	/* call the function for real */
	ret = LookupAccountName(NULL, username, sid, &sidsize, 
				domain, &dlen, &stype);
	err = (ret == 0);
    }

    if (!err) {
	/* build a trustee and get the file's DACL */
	BuildTrusteeWithSid(&t, sid);
	ret = GetNamedSecurityInfo(path, SE_FILE_OBJECT, 
				   DACL_SECURITY_INFORMATION, 
				   NULL, NULL, &dacl, NULL, &sd);
	err = (ret != ERROR_SUCCESS);
    }

    if (!err) {
	/* get the access mask for this trustee */
	ret = GetEffectiveRightsFromAcl(dacl, &t, &amask);
        if (ret != ERROR_SUCCESS) {
            fprintf(stderr, "GetEffectiveRights...: ret=%d\n", ret);   
            if (ret != RPC_S_SERVER_UNAVAILABLE && ret != ERROR_NO_SUCH_DOMAIN) {
                err = 1;
            }
        } else if (amask & STANDARD_RIGHTS_WRITE) {
	    ok = 1;
	}
    }

    if (dacl != NULL) {
	LocalFree(dacl);
    }
    if (sid != NULL) {
	LocalFree(sid);
    }    
    if (sd != NULL) {
	LocalFree(sd);
    }
    if (domain != NULL) {
	LocalFree(domain);
    }

    if (err) {
	win_show_last_error();
    }

    return ok;
}

int win32_delete_dir (const char *path)
{
    SHFILEOPSTRUCT op;
    char *from;
    int err = 0;

    from = calloc(strlen(path) + 2, 1);
    if (from == NULL) {
	return E_ALLOC;
    }

    strcpy(from, path);

    op.hwnd = NULL;
    op.wFunc = FO_DELETE;
    op.pFrom = from;
    op.pTo = NULL;
    op.fFlags = FOF_SILENT | FOF_NOCONFIRMATION | FOF_NOERRORUI;
    op.fAnyOperationsAborted = FALSE;
    op.hNameMappings = NULL;
    op.lpszProgressTitle = NULL;

    err = SHFileOperation(&op);

    free(from);

    return err;
}

char *slash_convert (char *str, int which)
{
    char *p;

    if (str == NULL) {
	return NULL;
    }

    p = str;
    while (*p) {
	if (which == FROM_BACKSLASH) {
	    if (*p == '\\') *p = '/';
	} else if (which == TO_BACKSLASH) {
	    if (*p == '/') *p = '\\';
	}
	p++;
    }

    return str;
}

static int try_for_R_path (HKEY tree, char *s)
{
    int err = 0;

    /* this used to work with R 2.9.1 */
    err = read_reg_val(tree, "R-core\\R", "InstallPath", s);

    if (err) {
	char version[8], path[32];

	/* new-style: path contains R version number */
	err = read_reg_val(tree, "R-core\\R", "Current Version", 
			   version);
	if (!err) {
	    sprintf(path, "R-core\\R\\%s", version);
	    err = read_reg_val(tree, path, "InstallPath", s);
	}
    }

    if (err) {
	/* did this variant work at one time? */
	err = read_reg_val(tree, "R", "InstallPath", s);
    }

    return err;
}

static void append_R_filename (char *s, int which)
{
    if (which == RGUI) {
	strcat(s, "Rgui.exe");
    } else if (which == RTERM) {
	strcat(s, "Rterm.exe");
    } else if (which == RLIB) {
	strcat(s, "R.dll");
    }
}

/* See if we can get the R installation path from the Windows
   registry. This is not a sure thing, since at least as of R
   2.11.1 recording the path in the registry on installation
   is an optional thing.

   To complicate matters, the path within the registry where
   we might find this information changed somewhere between
   R 2.9.1 and R 2.11.1.
*/

int R_path_from_registry (char *s, int which)
{
    int err;

    *s = '\0';

    err = try_for_R_path(HKEY_LOCAL_MACHINE, s);

    if (err) {
	/* maybe user is not an admin? */
	err = try_for_R_path(HKEY_CURRENT_USER, s);
    }

    if (!err && which != RBASE) {
	FILE *fp;

	strcat(s, "\\bin\\");
	append_R_filename(s, which);

	fp = fopen(s, "r");
	if (fp != NULL) {
	    fclose(fp);
	} else {
	    char *p = strrchr(s, 'R');

	    *p = '\0';
	    strcat(s, "i386\\");
	    append_R_filename(s, which);
	    fp = fopen(s, "r");
	    if (fp != NULL) {
		fclose(fp);
	    } else {
		err = E_FOPEN;
	    }
	}
    }

    return err;
}

/* for use in R, we need to form a version of the PATH with all
   backslashes doubled 
*/

static char *get_fixed_R_path (const char *path, const char *rpath)
{
    char *fixpath;
    int plen = (path != NULL)? strlen(path) : 0;
    int rlen = strlen(rpath);
    int i, ns = 0;

    for (i=0; i<plen; i++) {
	if (path[i] == '\\') ns++;
    }

    for (i=0; i<rlen; i++) {
	if (rpath[i] == '\\') ns++;
    }
 
    fixpath = malloc(plen + rlen + ns + 1);

    if (fixpath != NULL) {
	int j = 0;

	for (i=0; i<plen; i++) {
	    if (path[i] == '\\') {
		fixpath[j++] = '\\';
		fixpath[j++] = '\\';
	    } else {
		fixpath[j++] = path[i];
	    }
	}

	if (plen > 0) {
	    fixpath[j++] = ';';
	}

	for (i=0; i<rlen; i++) {
	    if (rpath[i] == '\\') {
		fixpath[j++] = '\\';
		fixpath[j++] = '\\';
	    } else {
		fixpath[j++] = rpath[i];
	    }
	}

	fixpath[j] = '\0';
    }

    return fixpath;
}

int maybe_print_R_path_addition (FILE *fp)
{
    static char *fixpath;
    static int ok;
    int err = 0;

    if (ok) {
	; /* no need to amend the path */
    } else if (fixpath != NULL) {
	/* revised path already built */
	fprintf(fp, "Sys.setenv(PATH=\"%s\")\n", fixpath);
    } else {
	char rpath[MAXLEN];

	strcpy(rpath, gretl_rlib_path());

	if (*rpath == '\0') {
	    err = 1;
	} else {
	    char *p = strrchr(rpath, '\\');
	    char *path = getenv("PATH");

	    if (p != NULL) {
		/* chop off "\R.dll" */
		*p = '\0';
	    }

	    if (path != NULL && strstr(path, rpath) != NULL) {
		ok = 1; /* nothing to be done */
	    } else {
		fixpath = get_fixed_R_path(path, rpath);
		if (fixpath == NULL) {
		    err = E_ALLOC;
		} else {
		    fprintf(fp, "Sys.setenv(PATH=\"%s\")\n", fixpath);
		}
	    }
	}
    }

    return err;
}

/*
  The original notice for the following functions
  (for which see http://ab-initio.mit.edu/~stevenj/align.c):

  _align_malloc and friends, implemented using Microsoft's public
  interfaces and with the help of the algorithm description provided
  by Wu Yongwei: http://sourceforge.net/mailarchive/message.php?msg_id=3847075

  I hereby place this implementation in the public domain.
               -- Steven G. Johnson (stevenj@alum.mit.edu)
*/

#include <errno.h>
#include <stddef.h> /* ptrdiff_t */

#ifdef HAVE_STDINT_H
# include <stdint.h> /* uintptr_t */
#else
# define uintptr_t size_t
#endif

#define NOT_POWER_OF_TWO(n) (((n) & ((n) - 1)))
#define UI(p) ((uintptr_t) (p))
#define CP(p) ((char *) p)

#define PTR_ALIGN(p0, alignment, offset)				\
            ((void *) (((UI(p0) + (alignment + sizeof(void*)) + offset)	\
			& (~UI(alignment - 1)))				\
		       - offset))

/* pointer must sometimes be aligned; assume sizeof(void*) is a power of two */
#define ORIG_PTR(p) (*(((void **) (UI(p) & (~UI(sizeof(void*) - 1)))) - 1))

static void *aligned_offset_malloc (size_t size, size_t alignment, 
				    size_t offset)
{
    void *p0, *p;

    if (NOT_POWER_OF_TWO(alignment)) {
	errno = EINVAL;
	return NULL;
    }

    if (size == 0)
	return NULL;

    if (alignment < sizeof(void *))
	alignment = sizeof(void *);

    /* including the extra sizeof(void*) is overkill on a 32-bit
       machine, since malloc is already 8-byte aligned, as long
       as we enforce alignment >= 8 ...but oh well */

    p0 = malloc(size + (alignment + sizeof(void *)));
    if (!p0) {
	return NULL;
    }

    p = PTR_ALIGN(p0, alignment, offset);
    ORIG_PTR(p) = p0;

    return p;
}

void *win32_memalign (size_t size, size_t alignment)
{
    return aligned_offset_malloc(size, alignment, 0);
}

void win32_aligned_free (void *mem)
{
    if (mem) {
	free(ORIG_PTR(mem));
    }
}

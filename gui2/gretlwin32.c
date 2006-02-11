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

#undef CHILD_DEBUG

#include "gretl.h"
#include "gretlwin32.h"
#include "textutil.h"
#include "guiprint.h"
#include "gpt_control.h"

#include <gdk/gdkwin32.h>
#include <dirent.h>
#include <mapi.h>
#include <shlobj.h>

#define HUSH_RUNTIME_WARNINGS

extern int wimp; /* settings.c */
extern int ws_startup (void);

int create_child_process (char *prog, char *env) 
{ 
    PROCESS_INFORMATION proc_info; 
    STARTUPINFO start_info; 
    int ret;

    ZeroMemory(&proc_info, sizeof proc_info);
    ZeroMemory(&start_info, sizeof start_info);
    start_info.cb = sizeof start_info; 

    ret = CreateProcess(NULL, 
			prog,          /* command line */
			NULL,          /* process security attributes  */
			NULL,          /* primary thread security attributes */ 
			FALSE,         /* handles are inherited?  */
			0,             /* creation flags  */
			(LPVOID) env,  /* NULL => use parent's environment  */
			NULL,          /* use parent's current directory  */
			&start_info,   /* receives STARTUPINFO */ 
			&proc_info);   /* receives PROCESS_INFORMATION  */

    if (ret == 0) {
	DWORD dw = GetLastError();
	win_show_error(dw);
    }

#ifdef CHILD_DEBUG
    if (1) {
	FILE *fp = fopen("gretlbug.txt", "w");

	if (fp != NULL) {
	    fprintf(fp, "gretl: create_child_process():\n"
		    " prog='%s'\n env='%s'\n", prog, env);
	    fprintf(fp, " return from CreateProcess() = %d\n", ret);
	    fclose(fp);
	}
    }	
#endif

    return ret;
}

void startR (const char *Rcommand)
{
    char Rprofile[MAXLEN], Rdata[MAXLEN], Rline[MAXLEN];
    const char *supp1 = "--no-init-file";
    const char *supp2 = "--no-restore-data";
    int *list;
    FILE *fp;
    int enverr;

    if (!data_status) {
	errbox(_("Please open a data file first"));
	return;
    }

    build_path(paths.userdir, "gretl.Rprofile", Rprofile, NULL);
    fp = gretl_fopen(Rprofile, "w");
    if (fp == NULL) {
	errbox(_("Couldn't write R startup file"));
	return;
    }

    enverr = ! SetEnvironmentVariable("R_PROFILE", Rprofile);
    if (enverr) {
	errbox(_("Couldn't set R_PROFILE environment variable"));
	fclose(fp);
	return;
    } 	

    build_path(paths.userdir, "Rdata.tmp", Rdata, NULL);

    sprintf(Rline, "store \"%s\" -r", Rdata);
    list = command_list_from_string(Rline);

    if (list == NULL ||
	write_data(Rdata, list, (const double **) Z, datainfo, 
		   OPT_R, NULL)) {
	errbox(_("Write of R data file failed"));
	fclose(fp);
	return; 
    }

    free(list);

    if (dataset_is_time_series(datainfo)) {
	fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fp);
	fputs("if (vnum > 1.89) library(stats) else library(ts)\n", fp);
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", 
		slash_convert(Rdata, FROM_BACKSLASH));
    } else {
	char Rtmp[MAXLEN];
	FILE *fq;

	build_path(paths.userdir, "Rtmp", Rtmp, NULL);
	fq = gretl_fopen(Rtmp, "w");
	fprintf(fq, "gretldata <- read.table(\"%s\")\n", 
		slash_convert(Rdata, FROM_BACKSLASH));
	fprintf(fq, "attach(gretldata)\n");
	fclose(fq);

	fprintf(fp, "source(\"%s\", echo=TRUE)\n", 
		slash_convert(Rtmp, FROM_BACKSLASH));
    }

    fclose(fp);

    sprintf(Rline, "\"%s\" %s %s", Rcommand, supp1, supp2);
    create_child_process(Rline, NULL);
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

static char *substr (char *targ, const char *src, const char *p)
{
    *targ = '\0';

    if (p == NULL) {
	strcpy(targ, src);
    } else {
	strncat(targ, src, p - src);
    }

    return targ;
}

int unmangle (const char *dosname, char *longname)
{
    HANDLE handle;
    WIN32_FIND_DATA fdata;
    char tmp[MAXLEN];
    const char *p;
    int err = 0, done = 0;
    char drive;

    *longname = '\0';
    
    if (sscanf(dosname, "%c:\\", &drive) != 1 ||
	dosname[strlen(dosname) - 1] == '\\') {
	strcpy(longname, dosname);
	return 0;
    }

    sprintf(longname, "%c:", tolower(drive));
    p = dosname + 2;

    while (!done) {
	p = strchr(p + 1, '\\');
	if (p == NULL) {
	    done = 1;
	} 
	substr(tmp, dosname, p);
	handle = FindFirstFile(tmp, &fdata);
	if (handle != INVALID_HANDLE_VALUE) {
	    strcat(longname, "\\");
	    strcat(longname, fdata.cFileName);
	    FindClose(handle);
	} else {
	    *longname = '\0';
	    err = 1;
	    break;
	}
    }

    return err;
}

void win_help (void)
{
    char hlpshow[MAXLEN];
    int found = 0;

    sprintf(hlpshow, "hh.exe \"%s\\%s\"", paths.gretldir, _("gretl.chm"));
    
    if (WinExec(hlpshow, SW_SHOWNORMAL) < 32) {
	if (strcmp("gretl.chm", _("gretl.chm"))) {
	    /* try falling back on untranslated helpfile */
	    sprintf(hlpshow, "hh.exe \"%s\\gretl.chm\"", paths.gretldir);
	    if (WinExec(hlpshow, SW_SHOWNORMAL) >= 32) found = 1;
	}
    } else {
	found = 1;
    }

    if (!found) {
	errbox(_("Couldn't access help file"));
    }
}

#ifdef HUSH_RUNTIME_WARNINGS

static void dummy_output_handler (const gchar *log_domain,
				  GLogLevelFlags log_level,
				  const gchar *message,
				  gpointer user_data)
{
    return;
}

static void hush_warnings (void)
{
    g_log_set_handler ("Gtk",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("Gdk",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("GLib",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("Pango",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
}

#endif /* HUSH_RUNTIME_WARNINGS */

char *default_windows_menu_fontspec (void)
{
    gchar *fontspec = NULL;
    NONCLIENTMETRICS ncm;

    memset(&ncm, 0, sizeof ncm);
    ncm.cbSize = sizeof ncm;

    if (SystemParametersInfo(SPI_GETNONCLIENTMETRICS, ncm.cbSize, &ncm, 0)) {
	HDC screen = GetDC(0);
	double y_scale = 72.0 / GetDeviceCaps(screen, LOGPIXELSY);
	int point_size = (int) (ncm.lfMenuFont.lfHeight * y_scale);

	if (point_size < 0) point_size = -point_size;
	fontspec = g_strdup_printf("%s %d", ncm.lfMenuFont.lfFaceName,
				   point_size);
	ReleaseDC(0, screen);
    }

    return fontspec;
}

static void try_to_get_windows_font (void)
{
    char regfont[MAXLEN];
    gchar *fontspec;

    /* don't override user's choice of font */
    if (read_reg_val(HKEY_CURRENT_USER, "gretl", "App_font", regfont) == 0) {
	if (*regfont != '\0') return;
    }

    fontspec = default_windows_menu_fontspec();

    if (fontspec != NULL) {
	int match = 0;
	PangoFontDescription *pfd;
	PangoFont *pfont;
	PangoContext *pc;
	GtkWidget *w;

	pfd = pango_font_description_from_string(fontspec);

	w = gtk_label_new(NULL);
	pc = gtk_widget_get_pango_context(w);
	pfont = pango_context_load_font(pc, pfd);
	match = (pfont != NULL);

	pango_font_description_free(pfd);
	g_object_unref(G_OBJECT(pc));
	gtk_widget_destroy(w);

	if (match) {
	    set_app_font(fontspec);
	}
	g_free(fontspec);
    }
}

void menu_font_option_off (void)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/File/Preferences/Menu font...", FALSE);
    }
}

void set_up_windows_look (void)
{
    if (wimp) { 
	/* "Windows Impersonator" wanted */
	size_t n = strlen(paths.gretldir);
	int needslash = (paths.gretldir[n-1] != SLASH);
	gchar *wimprc;

	wimprc = g_strdup_printf("%s%setc\\gtk-2.0\\gtkrc.wimp", 
				 paths.gretldir, (needslash)? "\\" : "");
	gtk_rc_parse(wimprc);
	g_free(wimprc);
    } else {
	try_to_get_windows_font();
    }
}

static char inifile[FILENAME_MAX];

const char *get_network_cfg_filename (void)
{
    return inifile;
}

static int set_network_cfg_filename (const char *prog)
{
    const char *p;
    int n;

    *inifile = '\0';
    
    n = strlen(prog);
    p = prog + n - 1;
    while (p - prog >= 0) {
	if (*p == '\\' || *p == '/') {
	    strncpy(inifile, prog, n - strlen(p));
	    strcat(inifile, "\\gretlnet.txt");
	    break;
	}
	p--;
    }

    return 0;
}

static int set_gd_fontpath (void)
{
    char fpath[MAX_PATH];
    char *envstr;
    UINT usz;

    usz = GetWindowsDirectory(fpath, MAX_PATH);
    if (usz == 0 || usz > MAX_PATH - 6) {
        return 1;
    }

    fpath[2] = '/';
    strcat(fpath, "/fonts");

    SetEnvironmentVariable("GDFONTPATH", fpath);

    envstr = malloc(strlen(fpath) + 12);
    if (envstr != NULL) {
	sprintf(envstr, "GDFONTPATH=%s", fpath);
	putenv(envstr);
	free(envstr);
    }

    return 0;        
}

void gretl_win32_init (const char *progname)
{
    set_network_cfg_filename(progname);

    read_rc(); /* get config info from registry */
    set_gd_fontpath();

# ifdef HUSH_RUNTIME_WARNINGS
    hush_warnings();
# endif 

    ws_startup(); 
    atexit(write_rc);
}

static int win_mkdir (const char *path)
{
    DIR *test;
    int done;

    test = opendir(path);
    if (test != NULL) {
	closedir(test);
	return 0;
    }

    done = CreateDirectory(path, NULL);
    
    return !done;
}

void win32_make_user_dirs (void)
{
    char dirname[MAXLEN];
    extern char *tramodir;
    size_t n;

    strcpy(dirname, paths.userdir);
    n = strlen(dirname);

    if (n > 0 && (dirname[n-1] == '\\' || dirname[n-1] == '/')) {
	dirname[n-1] = '\0';
    }

    if (win_mkdir(dirname)) {
	gchar *msg;

	msg = g_strdup_printf("Couldn't open or create gretl "
			      "user directory\n%s", paths.userdir);
	errbox(msg);
	g_free(msg);
	return;
    }

    build_path(paths.userdir, "x12arima", paths.x12adir, NULL);
    win_mkdir(paths.x12adir);

    build_path(paths.userdir, "tramo", tramodir, NULL);
    if (win_mkdir(tramodir)) return;

    sprintf(dirname, "%s\\output", tramodir);
    win_mkdir(dirname);

    sprintf(dirname, "%s\\graph", tramodir);
    if (win_mkdir(dirname)) return;

    sprintf(dirname, "%s\\graph\\acf", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\filters", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\forecast", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\series", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\spectra", tramodir);
    win_mkdir(dirname);
}

static int win_copy_buf (const char *buf, int fmt, size_t buflen)
{
    HGLOBAL winclip;
    LPTSTR ptr;
    char *winbuf;
    unsigned clip_format;
    size_t len;
    gchar *tr = NULL;

    if (buf == NULL) {
	errbox(_("Copy buffer was empty"));
	return 0;
    }

    if (!OpenClipboard(NULL)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    if (doing_nls() && (fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT)) { 
	tr = my_locale_from_utf8(buf);
	if (tr == NULL) {
	    CloseClipboard();
	    return 1;
	}
	winbuf = dosify_buffer(tr, fmt);
    } else {
	winbuf = dosify_buffer(buf, fmt);
    }

    if (winbuf == NULL) {
	CloseClipboard();
	return 1;
    }

    if (buflen == 0) {
	len = strlen(winbuf) + 1; 
    } else {
	len = buflen;
    }
        
    winclip = GlobalAlloc(GMEM_MOVEABLE, len * sizeof(TCHAR));        

    ptr = GlobalLock(winclip);
    memcpy(ptr, winbuf, len);
    GlobalUnlock(winclip); 

    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) { 
	clip_format = RegisterClipboardFormat("Rich Text Format");
    } else if (fmt == GRETL_FORMAT_CSV) {
	clip_format = RegisterClipboardFormat("CSV");
    } else {
	clip_format = CF_TEXT;
    }

    SetClipboardData(clip_format, winclip);

    CloseClipboard();

    if (tr != NULL) {
	free(tr);
    }

    free(winbuf);

    return 0;
}

int prn_to_clipboard (PRN *prn, int fmt)
{
    const char *buf = gretl_print_get_buffer(prn);

    return win_copy_buf(buf, fmt, 0);
}

int win_buf_to_clipboard (const char *buf)
{
    HGLOBAL winclip;
    LPTSTR ptr;
    size_t len;

    if (!OpenClipboard(NULL)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    len = strlen(buf);
    winclip = GlobalAlloc(GMEM_MOVEABLE, (len + 1) * sizeof(TCHAR));        

    ptr = GlobalLock(winclip);
    memcpy(ptr, buf, len + 1);
    GlobalUnlock(winclip); 

    SetClipboardData(CF_TEXT, winclip);

    CloseClipboard();

    return 0;
}

int fnamecmp_win32 (const char *f1, const char *f2)
{
    gchar *c1, *c2;
    int ret;

    c1 = g_utf8_casefold(f1, -1);
    c2 = g_utf8_casefold(f2, -1);

    ret = strcmp(c1, c2);

    g_free(c1);
    g_free(c2);

    return ret;
}

static char *fname_from_fullname (char *fullname)
{
    char *fname = fullname;
    char *p;

    p = strrchr(fullname, '\\');
    if (p == NULL) {
	p = strrchr(fullname, '/');
    }

    if (p != NULL) {
	fname = p + 1;
    }

    return fname;
}

int send_file (char *fullname)
{
    HINSTANCE mapilib = NULL;
    LPMAPISENDMAIL send_mail = NULL;
    int err = 0;

    mapilib = LoadLibrary("MAPI32.DLL");
    if (mapilib == NULL) {
	err = 1;
    } else {
	send_mail = (LPMAPISENDMAIL) GetProcAddress(mapilib, "MAPISendMail");
	if (send_mail == NULL) {
	    err = 1;
	} 
    }

    if (err) {
	errbox("Couldn't access Windows MAPI system");
    } else {
	char *fname = fname_from_fullname(fullname);
	gchar *note = NULL;
	MapiFileDesc mfd;
	MapiMessage msg;
	ULONG sd;

	memset(&mfd, 0, sizeof mfd);
	memset(&msg, 0, sizeof msg);

	mfd.lpszPathName = fullname;
	mfd.lpszFileName = fname;
	mfd.nPosition = 1; /* ? */

	if (strstr(fname, ".gdt") != NULL) {
	    note = g_strdup_printf("Please find the gretl data file %s attached.\n", fname);
	    msg.lpszSubject  = "dataset";
	} else {
	    note = g_strdup_printf("Please find the gretl script %s attached.\n", fname);
	    msg.lpszSubject  = "script";
	}

	msg.lpszNoteText = note;
	msg.nFileCount = 1;
	msg.lpFiles = &mfd;

	sd = send_mail(0L, 0, &msg, MAPI_DIALOG, 0L);

	if (sd != SUCCESS_SUCCESS && sd != MAPI_E_USER_ABORT) {
	    errbox("MAPI error sending message");
	}

	g_free(note);
    }

    if (mapilib != NULL) {
	FreeLibrary(mapilib);
    }

    return err;
}

/* copy plot to clipboard by generating an EMF file (enhanced
   metafile), reading it into a buffer, and putting it on the
   clipboard.

   Weirdness: when an emf is put on the clipboard as below, Word 2000
   behaves thus: a straight "Paste" puts in a version of the graph
   with squashed up numbers on the axes and no legend text; but a
   "Paste special" (where one accepts the default of pasting it as an
   enhanced metafile) puts in an accurate version with correct text.
   Go figure.  (This is on win98)
*/

static int emf_to_clip (char *emfname)
{
    HWND mainw;
    HENHMETAFILE hemf, hemfclip;
    HANDLE htest;

    mainw = GDK_WINDOW_HWND(mdata->w->window);
    if (mainw == NULL) {
	errbox("Got NULL HWND");
	return 1;
    }	

    if (!OpenClipboard(mainw)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    hemf = GetEnhMetaFile(emfname);
    if (hemf == NULL) {
	errbox("Couldn't get handle to graphic metafile");
	return 1;
    }

    hemfclip = CopyEnhMetaFile(hemf, NULL);
    if (hemfclip == NULL) {
	errbox("Couldn't copy graphic metafile");
	return 1;
    }    

    htest = SetClipboardData(CF_ENHMETAFILE, hemfclip);
    if (htest == NULL) {
	errbox("Failed to put data on clipboard");
	return 1;
    }  	

    CloseClipboard();

    DeleteEnhMetaFile(hemf);

    return 0;
}

void win32_process_graph (GPT_SPEC *spec, int color, int dest)
{
    FILE *fq;
    PRN *prn;
    char plottmp[MAXLEN], plotline[MAXLEN];
    gchar *plotcmd = NULL;
    gchar *emfname = NULL;
    int err, done_pt2 = 0;

    /* create temporary file to hold the special gnuplot commands */
    if (user_fopen("gptout.tmp", plottmp, &prn)) {
	return;
    }

    /* open the gnuplot source file for the graph */
    fq = gretl_fopen(spec->fname, "r");
    if (fq == NULL) {
	errbox(_("Couldn't access graph info"));
	gretl_print_destroy(prn);
	return;
    }

    /* generate gnuplot source file to make emf */
    pprintf(prn, "%s\n", get_gretl_emf_term_line(spec->code, color));
    emfname = g_strdup_printf("%sgpttmp.emf", paths.userdir);
    pprintf(prn, "set output '%s'\n", emfname);
    while (fgets(plotline, MAXLEN-1, fq)) {
	if (!done_pt2 && strstr(plotline, "using 1:2")) {
	    done_pt2 = maybe_switch_emf_point_style(plotline, prn);
	} else if (strncmp(plotline, "set term", 8) && 
	    strncmp(plotline, "set output", 10)) {
	    pputs(prn, plotline);
	}
    }

    gretl_print_destroy(prn);
    fclose(fq);

    /* get gnuplot to create the emf file */
    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
			      plottmp);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
    g_free(plotcmd);
    remove(plottmp);
    
    if (err) {
        errbox(_("Gnuplot error creating graph"));
    } else if (dest == WIN32_TO_CLIPBOARD) {
	err = emf_to_clip(emfname);
	if (!err) {
	    infobox(_("To paste, use Edit/Paste special.../Enhanced metafile"));
	}
    } else if (dest == WIN32_TO_PRINTER) {
	err = winprint_graph(emfname);
    }

    remove(emfname);
    g_free(emfname);
}

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

#include "gretl.h"
#include "gpt_control.h"
#include "gpt_dialog.h"
#include "session.h"
#include "textbuf.h"
#include "textutil.h"
#include "filelists.h"
#include "fileselect.h"
#include "fnsave.h"
#include "database.h"
#include "datafiles.h"
#include "graphics.h"
#include "bootstrap.h"

#include <sys/stat.h>
#include <unistd.h>

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#define IS_DAT_ACTION(i) (i == SAVE_DATA || \
                          i == SAVE_DATA_AS || \
                          i == OPEN_DATA)

#define OPEN_DATA_ACTION(i)  (i == OPEN_DATA || \
                              i == OPEN_CSV || \
                              i == OPEN_OCTAVE || \
                              i == OPEN_GNUMERIC || \
	                      i == OPEN_XLS || \
                              i == OPEN_WF1 || \
                              i == OPEN_DTA || \
	                      i == OPEN_SAV || \
			      i == OPEN_SAS || \
                              i == OPEN_JMULTI || \
                              i == OPEN_ODS)

#define APPEND_DATA_ACTION(i) (i == APPEND_DATA || \
                               i == APPEND_CSV || \
                               i == APPEND_OCTAVE || \
                               i == APPEND_GNUMERIC || \
                               i == APPEND_XLS || \
                               i == APPEND_WF1 || \
                               i == APPEND_DTA || \
                               i == APPEND_SAV || \
			       i == APPEND_SAS || \
                               i == APPEND_JMULTI || \
                               i == APPEND_ODS)

#define EXPORT_ACTION(a,s) ((a == EXPORT_OCTAVE || \
                             a == EXPORT_R || \
                             a == EXPORT_CSV || \
                             a == EXPORT_DAT || \
                             a == EXPORT_JM) && s != FSEL_DATA_PRN)

#define GDT_ACTION(i) (i == SAVE_DATA || \
                       i == SAVE_DATA_AS || \
                       i == SAVE_BOOT_DATA || \
                       i == OPEN_DATA || \
                       i == APPEND_DATA)

#define SET_DIR_ACTION(i) (i == SET_DIR || i == SET_WDIR || \
                           i == SET_FDIR || i == SET_DBDIR)

struct extmap {
    int action;
    char *ext;
};

static struct extmap action_map[] = {
    { SAVE_DBDATA,       ".bin" },
    { SAVE_SCRIPT,       ".inp" },
    { SAVE_FUNCTIONS_AS, ".inp" },
    { SAVE_CMD_LOG,      ".inp" },
    { SAVE_CONSOLE,      ".txt" },
    { SAVE_SESSION,      ".gretl" },
    { SAVE_GP_CMDS,      ".plt" },
    { SAVE_R_CMDS,       ".R" },
    { SAVE_OX_CMDS,      ".ox" },
    { SAVE_OCTAVE_CMDS,  ".m" },    
    { SAVE_FUNCTIONS,    ".gfn" },
    { SAVE_MARKERS,      ".txt" },
    { SAVE_LABELS,       ".txt" },
    { SAVE_GFN_SPEC,     ".spec" },
    { EXPORT_CSV,        ".csv" },
    { EXPORT_R,          ".R" },
    { OPEN_OCTAVE,       ".m" },
    { APPEND_OCTAVE,     ".m" },
    { EXPORT_OCTAVE,     ".m" },
    { EXPORT_DAT,        ".dat" },
    { SAVE_OUTPUT,       ".txt" },
    { SAVE_TEX,          ".tex" },
    { SAVE_RTF,          ".rtf" },
    { SAVE_TEXT,         ".txt" },
    { OPEN_SCRIPT,       ".inp" },
    { OPEN_SESSION,      ".gretl" },
    { OPEN_CSV,          ".csv" },
    { APPEND_CSV,        ".csv" },
    { OPEN_GNUMERIC,     ".gnumeric" },
    { APPEND_GNUMERIC,   ".gnumeric" },
    { OPEN_XLS,          ".xls" },
    { APPEND_XLS,        ".xls" },
    { OPEN_WF1,          ".wf1" },
    { APPEND_WF1,        ".wf1" },
    { OPEN_DTA,          ".dta" },
    { APPEND_DTA,        ".dta" },
    { OPEN_SAV,          ".sav" },
    { APPEND_SAV,        ".sav" },
    { OPEN_SAS,          ".xpt" },
    { APPEND_SAS,        ".xpt" },
    { OPEN_JMULTI,       ".dat" },
    { APPEND_JMULTI,     ".dat" },
    { OPEN_ODS,          ".ods" },
    { APPEND_ODS,        ".ods" },
    { OPEN_RATS_DB,      ".rat" },
    { OPEN_PCGIVE_DB,    ".bn7" },
    { OPEN_GFN,          ".gfn" },
    { OPEN_BARS,         ".txt" },
    { FILE_OP_MAX,       NULL }
};

static gretlopt save_action_to_opt (int action, gpointer p)
{
    gretlopt opt = OPT_NONE;

    switch (action) {
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_BOOT_DATA: opt = OPT_Z; break;
    case SAVE_DBDATA:    opt = OPT_D; break;
    case EXPORT_OCTAVE:  opt = OPT_M; break;
    case EXPORT_R:       opt = OPT_R; break;
    case EXPORT_CSV:     opt = OPT_C; break;
    case EXPORT_DAT:     opt = OPT_G; break; /* PcGive */
    case EXPORT_JM:      opt = OPT_J; break; /* JMulti */
    default: break;
    }

    if (action == SAVE_DATA_AS && session_file_is_open()) {
	opt |= OPT_X; /* "exporting" to GDT */
    }

    if (p != NULL) {
	opt |= GPOINTER_TO_INT(p);
    }

    return opt;
}

static const char *get_extension_for_action (int action, gpointer data)
{
    const char *s = NULL;

    if (GDT_ACTION(action)) {
	return ".gdt";
    } else if (action == SAVE_GNUPLOT || action == SAVE_GRAPHIC) {
	int ttype = gp_term_code(data, action);

	if (ttype == GP_TERM_EPS)      s = ".eps";
	else if (ttype == GP_TERM_PDF) s = ".pdf";
	else if (ttype == GP_TERM_FIG) s = ".fig";
	else if (ttype == GP_TERM_TEX) s = ".tex";
	else if (ttype == GP_TERM_PNG) s = ".png";
	else if (ttype == GP_TERM_EMF) s = ".emf";
	else if (ttype == GP_TERM_SVG) s = ".svg";
	else if (ttype == GP_TERM_PLT) s = ".plt";
    } else {
	size_t i;

	for (i=0; i < sizeof action_map / sizeof *action_map; i++) {
	    if (action == action_map[i].action) {
		s = action_map[i].ext;
		break;
	    }
	}
    }

    return s;
}

static int fname_has_extension (const char *fname)
{
    const char *p;
    int ret = 0;

#ifdef G_OS_WIN32
    p = strrchr(fname, '\\');
    if (p == NULL) {
	p = strrchr(fname, '/');
    }
#else
    p = strrchr(fname, '/');
#endif

    if (p == NULL) {
	p = fname;
    }

    p = strrchr(p, '.');

    if (p != NULL) {
	/* got a '.' in the basename: we'll count fname
	   as having an extension unless the trailing
	   portion is an integer
	*/
	if (!integer_string(p+1)) {
	    ret = 1;
	}
    }

    return ret;
}

static int post_process_savename (char *fname, int action, gpointer data)
{
    int err = 0;

    if (fname == NULL) {
	err = 1;
    } else if (gretl_isdir(fname)) {
	/* fname is actually a directory */
	err = !SET_DIR_ACTION(action);
    } else {
	/* get the recommended extension for this action */
	const char *ext = get_extension_for_action(action, data);

	if (ext != NULL && !has_suffix(fname, ext)) {
	    /* We found a recommended extension, and it's not
	       already present on fname; so we should probably
	       append the extension -- unless perhaps we're
	       looking at an action where there isn't exactly
	       a canonical extension and fname already has some
	       sort of extension in place.
	    */
	    if (strcmp(ext, ".txt") && strcmp(ext, ".plt")) {
		/* append canonical extension */
		strcat(fname, ext);
	    } else if (!fname_has_extension(fname)) {
		strcat(fname, ext);
	    }
	}
    }    

    return err;
}

static void script_window_update (windata_t *vwin, const char *fname)
{
    gchar *trfname, *title;
    const char *p = strrchr(fname, SLASH);

    /* ensure UTF-8 filename for display */
    if (p != NULL) {
	trfname = my_filename_to_utf8(p + 1);
    } else {
	trfname = my_filename_to_utf8(fname);
    }

    /* update the window title */
    title = g_strdup_printf("gretl: %s", trfname);
    gtk_window_set_title(GTK_WINDOW(vwin->main), title);
    g_free(trfname);
    g_free(title);

    /* and update internal filename record */
    strcpy(vwin->fname, fname);

    if (vwin->role == VIEW_LOG || vwin->role == VIEW_SCRIPT ||
	vwin->role == VIEW_PKG_CODE) {
	/* change role of window for editing */
	vwin->role = EDIT_SCRIPT;
    } else if (vwin->role == EDIT_GP) {
	/* plot file no longer under session control */
	vwin->flags &= ~VWIN_SESSION_GRAPH;
    }

    mark_vwin_content_saved(vwin);

    /* make the window editable */
    if (!gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->text))) {
	view_window_set_editable(vwin);
    }
}

static void 
save_editable_content (int action, const char *fname, windata_t *vwin)
{
    FILE *fp;
    gchar *buf;
    gchar *trbuf;

    buf = textview_get_text(vwin->text);
    if (buf == NULL) {
	errbox("Couldn't retrieve buffer");
	return;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	file_write_errbox(fname);
	g_free(buf);
	return;
    }

    trbuf = my_locale_from_utf8(buf);
    if (trbuf != NULL) {
	/* UTF minuses? */
	system_print_buf(trbuf, fp);
	g_free(trbuf);
    }

    g_free(buf);
    fclose(fp);
    
    if (action == SAVE_SCRIPT) {
	strcpy(scriptfile, fname);
	mkfilelist(FILE_LIST_SCRIPT, scriptfile);
	script_window_update(vwin, fname);
    } else if (action == SAVE_GP_CMDS || 
	       action == SAVE_R_CMDS ||
	       action == SAVE_OX_CMDS ||
	       action == SAVE_OCTAVE_CMDS) {
	script_window_update(vwin, fname);
    } 
}

static void filesel_save_prn_buffer (PRN *prn, const char *fname)
{
    FILE *fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	file_write_errbox(fname);
    } else {
	const char *buf = gretl_print_get_buffer(prn);

	fputs(buf, fp);
	fclose(fp);
    }
}

static void filesel_open_script (const char *fname, windata_t *vwin)
{
    if (vwin != NULL) {
	strcpy(tryfile, fname);
	sourceview_insert_file(vwin, fname);
    } else if (has_suffix(fname, ".R")) {
	view_file(fname, 1, 0, 78, 370, EDIT_R);
    } else if (has_suffix(fname, ".plt")) {
	view_file(fname, 1, 0, 78, 370, EDIT_GP);
    } else if (ox_support && has_suffix(fname, ".ox")) {
	view_file(fname, 1, 0, 78, 370, EDIT_OX);
    } else if (has_suffix(fname, ".m")) {
	view_file(fname, 1, 0, 78, 370, EDIT_OCTAVE);
    } else {
	strcpy(tryfile, fname);
	if (view_file(tryfile, 1, 0, 78, 370, EDIT_SCRIPT) != NULL) {
	    strcpy(scriptfile, tryfile);
	    mkfilelist(FILE_LIST_SCRIPT, scriptfile);
	    gretl_set_current_dir(scriptfile);
	}
    }
}

static void filesel_open_session (const char *fname)
{
    strcpy(tryfile, fname);

    if (gretl_is_pkzip_file(fname)) {
	verify_open_session();
    } else {
	/* old script-style session file? */
	windata_t *vwin;

	if (has_system_prefix(tryfile, SCRIPT_SEARCH)) {
	    vwin = view_file(tryfile, 0, 0, 78, 370, VIEW_SCRIPT);
	} else {
	    vwin = view_file(tryfile, 1, 0, 78, 370, EDIT_SCRIPT);
	}

	if (vwin != NULL) {
	    strcpy(scriptfile, tryfile);
	}
    }
}

static char *suggested_savename (const char *fname)
{
    const char *ss = strrchr(fname, SLASH);
    char *s, *sfx;

    if (ss == NULL) {
	s = g_strdup(fname);
    } else {
	s = g_strdup(ss + 1);
    }

    sfx = strrchr(s, '.');

    if (sfx != NULL && (strlen(sfx) == 4 || !strcmp(sfx, ".gnumeric"))) {
	if (strcmp(sfx, ".gdt")) {
	    strcpy(sfx, ".gdt");
	}
    }

    return s;
}

static char *suggested_exportname (const char *fname, int action)
{
    const char *ss = strrchr(fname, SLASH);
    char *s, *sfx;

    if (ss == NULL) {
	s = g_strdup(fname);
    } else {
	s = g_strdup(ss + 1);
    }

    sfx = strrchr(s, '.');

    if (sfx != NULL && (strlen(sfx) == 4 || !strcmp(sfx, ".gnumeric"))) {
	const char *test;

	switch(action) {
	case EXPORT_OCTAVE:
	    test = ".m";
	    break;
	case EXPORT_R:
	    test = ".R";
	    break;
	case EXPORT_CSV:
	    test = ".csv";
	    break;
	case EXPORT_DAT:
	case EXPORT_JM:
	    test = ".dat";
	    break;
	default:
	    test = NULL;
	    break;
	}

	if (test != NULL && strcmp(test, sfx)) {
	    strcpy(sfx, test);
	}
    }

    return s;
}

static void bootstrap_save_callback (const char *fname)
{
    int err = bootstrap_save_data(fname);

    if (err) {
	gui_errmsg(err);
    } 
}

static void maybe_set_fsel_status (int action, FselDataSrc src, 
				   gpointer p, int val)
{
    if (src == FSEL_DATA_STATUS && p != NULL) {
	* (int *) p = val;
    } else if (action == OPEN_BARS && val) {
	/* error or cancel */
	set_plotbars_filename(NULL, p);
    }
}

static int overwrite_stop (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");

    if (fp != NULL) {
	fclose(fp);
	if (yes_no_dialog(NULL, 
			  _("There is already a file of this name.\n"
			    "OK to overwrite it?"), 
			  0) == GRETL_NO) {
	    return 1;
	}
    }

    return 0;
}

static void
file_selector_process_result (const char *in_fname, int action, FselDataSrc src,
			      gpointer data)
{
    char fname[FILENAME_MAX];
    int err = 0;

    *fname = 0;
    strncat(fname, in_fname, FILENAME_MAX - 1);

    if (action < END_OPEN) {
	/* opening a file: check that the file is accessible */
	FILE *fp = gretl_fopen(fname, "r");

	if (fp == NULL) {
	    file_read_errbox(fname);
	    maybe_set_fsel_status(action, src, data, E_FOPEN);
	    return;
	} else {
	    fclose(fp);
	}
    } 

    if (OPEN_DATA_ACTION(action)) {
	strcpy(tryfile, fname);
	verify_open_data(NULL, action);
    } else if (APPEND_DATA_ACTION(action)) {
	strcpy(tryfile, fname);
	do_open_data(NULL, action);
    } else if (action == OPEN_SCRIPT) {
	if (src == FSEL_DATA_VWIN) {
	    filesel_open_script(fname, data);
	} else {
	    filesel_open_script(fname, NULL);
	}
    } else if (action == OPEN_SESSION) {
	filesel_open_session(fname);
    } else if (action == OPEN_MARKERS) {
	do_add_markers(fname);
    } else if (action == OPEN_LABELS) {
	do_add_labels(fname);
    } else if (action == OPEN_BARS) {
	set_plotbars_filename(fname, data);
    } else if (action == OPEN_GFN) {
	edit_function_package(fname);
    } else if (action == OPEN_RATS_DB) {
	open_rats_window(fname);
    } else if (action == OPEN_PCGIVE_DB) {
	open_bn7_window(fname);
    }

    if (action < END_OPEN) {
	return;
    }

    /* now for the save/export options */

    if (post_process_savename(fname, action, data)) {
	return;
    }

    if (EXPORT_ACTION(action, src) || 
	(action > END_SAVE_DATA && action < END_SAVE_OTHER)) {
	/* saving CSV, graphs etc.; check overwrite */
	if (overwrite_stop(fname)) {
	    return;
	}
    }

    if (action == SAVE_TEX) {
	if (src == FSEL_DATA_PRN) {
	    save_latex((PRN *) data, fname);
	} else {
	    save_latex(NULL, fname);
	}
    } else if (src == FSEL_DATA_PRN) {
	filesel_save_prn_buffer((PRN *) data, fname);
    } else if (SAVE_DATA_ACTION(action)) {
	err = do_store(fname, save_action_to_opt(action, data));
    } else if (action == SAVE_GNUPLOT) {
	save_graph_to_file(data, fname);
    } else if (action == SAVE_GRAPHIC) {
	save_graphic_to_file(data, fname);
    } else if (action == SAVE_SESSION) {
	err = save_session(fname);
    } else if (action == SAVE_CMD_LOG) {
	err = save_session_commands(fname);
    } else if (action == SAVE_FUNCTIONS) {
	err = save_function_package(fname, data);
    } else if (action == SAVE_FUNCTIONS_AS) {
	err = save_function_package_as_script(fname, data);
    } else if (action == SAVE_GFN_SPEC) {
	err = save_function_package_spec(fname, data);
    } else if (action == SAVE_BOOT_DATA) {
	bootstrap_save_callback(fname);
    } else if (action == SAVE_MARKERS) {
	err = do_save_markers(fname);
    } else if (action == SAVE_LABELS) {
	err = do_save_labels(fname);
    } else if (action == SET_PROG || action == SET_DIR) {
	set_path_callback(data, fname);
    } else if (action == SET_WDIR) {
	set_working_dir_callback(data, fname);
    } else if (action == SET_FDIR) {
	set_funcs_dir_callback(data, fname);
    } else if (action == SET_DBDIR) {
	set_db_dir_callback(data, fname);
    } else {
	windata_t *vwin = (windata_t *) data;

	save_editable_content(action, fname, vwin);
    }

    maybe_set_fsel_status(action, src, data, err);
}

static char *get_filter_suffix (int action, gpointer data, char *suffix)
{
    
    const char *ext = get_extension_for_action(action, data);

    if (ext == NULL) { 
	strcpy(suffix, "*");
    } else {
	sprintf(suffix, "*%s", ext);
    }

    return suffix;
}

static GtkFileFilter *get_file_filter (int action, gpointer data)
{
    GtkFileFilter *filter;
    char suffix[16];
    
    filter = gtk_file_filter_new();
    get_filter_suffix(action, data, suffix);
    gtk_file_filter_add_pattern(filter, suffix);

    return filter;
}

#ifdef G_OS_WIN32

static gchar *inplace_windows_filename (gchar *s)
{
    char tmp[MAXLEN];
    int err;

    err = filename_to_win32(tmp, s);

    if (!err) {
	g_free(s);
	s = g_strdup(tmp);
    }

    return s;
}

static char *win32_correct_path (char *s)
{
    if (!g_utf8_validate(s, -1, NULL)) {
	gchar *tmp = my_locale_to_utf8(s);

	if (tmp != NULL) {
	    strcpy(s, tmp);
	    g_free(tmp);
	}
    }

    return s;
}

#endif

static void set_default_progs_path (GtkFileChooser *fsel)
{
#ifdef G_OS_WIN32
    char *progs = program_files_path();

    if (progs != NULL) {
	gchar *path = my_filename_to_utf8(progs);

	gtk_file_chooser_set_current_folder(fsel, path);
	g_free(path);
	free(progs);
    }
#else
    gtk_file_chooser_set_current_folder(fsel, "/usr/bin");
#endif
}

static void filesel_maybe_set_current_name (GtkFileChooser *filesel,
					    int action, FselDataSrc src,
					    gpointer data)
{
    gchar *currname;

    if (action == SAVE_DATA && *datafile != '\0') {
	currname = suggested_savename(datafile);
	gtk_file_chooser_set_current_name(filesel, currname);
	g_free(currname);
    } else if (EXPORT_ACTION(action, src) && *datafile != '\0') {
	currname = suggested_exportname(datafile, action);
	gtk_file_chooser_set_current_name(filesel, currname);
	g_free(currname);
    } else if (action == SAVE_CMD_LOG) {
	gtk_file_chooser_set_current_name(filesel, "session.inp");
    } else if (action == SET_PROG || action == SET_DIR) {
	currname = (gchar *) data;
	if (currname != NULL && g_path_is_absolute(currname)) {
	    gtk_file_chooser_set_filename(filesel, currname);
	} else if (action == SET_PROG) {
	    set_default_progs_path(filesel);
	}	    
    } else if (action == SAVE_FUNCTIONS || 
	       action == SAVE_FUNCTIONS_AS ||
	       action == SAVE_GFN_SPEC) {
	char fname[MAXLEN];

	*fname = '\0';
	get_default_package_name(fname, data, action);
	gtk_file_chooser_set_current_name(filesel, fname);
    } else if (action == SAVE_MARKERS) {
	gtk_file_chooser_set_current_name(filesel, "markers.txt");
    } else if (action == SAVE_LABELS) {
	gtk_file_chooser_set_current_name(filesel, "labels.txt");
    }
	
}

static void filesel_add_filter (GtkWidget *filesel,
				const char *desc, 
				const char *pat)
{
    GtkFileFilter *filt = gtk_file_filter_new();

    gtk_file_filter_set_name(filt, desc);
    gtk_file_filter_add_pattern(filt, pat);
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(filesel), filt);
}

static void filesel_set_filters (GtkWidget *filesel, int action,
				 FselDataSrc src, gpointer data)
{
    if (action == OPEN_CSV || action == APPEND_CSV) {
	filesel_add_filter(filesel, _("CSV files (*.csv)"), "*.csv");
	filesel_add_filter(filesel, _("ASCII files (*.txt)"), "*.txt");
	filesel_add_filter(filesel, _("all files (*.*)"), "*");
    } else if (action == OPEN_SCRIPT) {
	filesel_add_filter(filesel, _("gretl script files (*.inp)"), "*.inp");
	if (src != FSEL_DATA_FNPKG) {
	    filesel_add_filter(filesel, _("GNU R files (*.R)"), "*.R");
	    filesel_add_filter(filesel, _("gnuplot files (*.plt)"), "*.plt");
	    filesel_add_filter(filesel, _("GNU Octave files (*.m)"), "*.m");
	    if (ox_support) {
		filesel_add_filter(filesel, _("Ox files (*.ox)"), "*.ox");
	    }
	}
    } else if (action == OPEN_LABELS) {
	filesel_add_filter(filesel, _("ASCII files (*.txt)"), "*.txt");
	filesel_add_filter(filesel, _("all files (*.*)"), "*");
    } else {
	GtkFileFilter *filter = get_file_filter(action, data);

	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), filter);
    }
}

static void remember_folder (GtkFileChooser *chooser, char *savedir)
{
    gchar *folder = gtk_file_chooser_get_current_folder(chooser);

    strcpy(savedir, folder);
    g_free(folder);
}

static void gtk_file_selector (int action, FselDataSrc src, 
			       gpointer data, GtkWidget *parent) 
{
    static char savedir[MAXLEN];
    char startdir[MAXLEN];
    GtkWidget *filesel;
    GtkFileChooserAction fa;
    const gchar *okstr;
    int remember = get_keep_folder();
    gint response;

    if (SET_DIR_ACTION(action)) {
	fa = GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER;
	okstr = GTK_STOCK_OK;
	remember = 0;
    } else if (action == SET_PROG) {
	fa = GTK_FILE_CHOOSER_ACTION_OPEN;
	okstr = GTK_STOCK_OK;
	remember = 0;
    } else if (action < END_OPEN) {
	fa = GTK_FILE_CHOOSER_ACTION_OPEN;
	okstr = GTK_STOCK_OPEN;
    } else {
	fa = GTK_FILE_CHOOSER_ACTION_SAVE;
	okstr = GTK_STOCK_SAVE;
	if (action == SAVE_FUNCTIONS || 
	    action == SAVE_DATA_PKG ||
	    action == SAVE_GFN_SPEC) {
	    remember = 0;
	}
    }

    if (remember && *savedir != '\0') {
	strcpy(startdir, savedir);
    } else {
	get_default_dir(startdir, action);
    }

#ifdef G_OS_WIN32
    win32_correct_path(startdir);
#endif

    if (parent == NULL) {
	parent = mdata->main;
    }

    filesel = gtk_file_chooser_dialog_new(NULL, GTK_WINDOW(parent), fa,
					  GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
					  okstr, GTK_RESPONSE_ACCEPT,
					  NULL);

    gtk_dialog_set_default_response(GTK_DIALOG(filesel), GTK_RESPONSE_ACCEPT);

    if (SET_DIR_ACTION(action)) {
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(filesel), startdir);
    } else {
	filesel_set_filters(filesel, action, src, data);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filesel), startdir);
	filesel_maybe_set_current_name(GTK_FILE_CHOOSER(filesel), action,
				       src, data);
    }

    response = gtk_dialog_run(GTK_DIALOG(filesel));

    if (response == GTK_RESPONSE_ACCEPT) {
	gchar *fname;
	
	fname = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(filesel));
	if (remember) {
	    remember_folder(GTK_FILE_CHOOSER(filesel), savedir);
	}
	gtk_widget_destroy(filesel);
	filesel = NULL;

	if (fname != NULL) {
#ifdef G_OS_WIN32
	    fname = inplace_windows_filename(fname);
#endif
	    file_selector_process_result(fname, action, src, data);
	    g_free(fname);
	}
    } else {
	maybe_set_fsel_status(action, src, data, GRETL_CANCEL);
    }

    if (filesel != NULL) {
	gtk_widget_destroy(filesel);
    }
}

void file_selector (int action, FselDataSrc src, gpointer data)
{
    GtkWidget *w = NULL;

    if (src == FSEL_DATA_VWIN) {
	windata_t *vwin = (windata_t *) data;

	w = vwin->main;
    }

    gtk_file_selector(action, src, data, w);
}

void file_selector_with_parent (int action, FselDataSrc src, 
				gpointer data, GtkWidget *w)
{
    gtk_file_selector(action, src, data, w);
}

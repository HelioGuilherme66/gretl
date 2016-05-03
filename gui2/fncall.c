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
#include "version.h"
#include "dlgutils.h"
#include "selector.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "uservar.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "gretl_xml.h"
#include "gretl_typemap.h"
#include "gretl_zip.h"
#include "database.h"
#include "guiprint.h"
#include "ssheet.h"
#include "datafiles.h"
#include "toolbar.h"
#include "obsbutton.h"
#include "cmdstack.h"
#include "winstack.h"
#include "treeutils.h"
#include "gfn_arglists.h"
#include "fncall.h"

#include <errno.h>

#define FCDEBUG 0

enum {
    SHOW_GUI_MAIN = 1 << 0,
    MODEL_CALL    = 1 << 1
};

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;      /* main dialog */
    GtkWidget *top_hbox; /* upper hbox in dialog */
    windata_t *vwin;     /* gretl caller window */
    GList *vsels;        /* series argument selectors */
    GList *lsels;        /* list argument selectors */
    GList *msels;        /* matrix arg selectors */
    GList *ssels;        /* scalar arg selectors */
    GList *bsels;        /* bundle arg selectors */
    GList *asels;        /* array arg selectors */
    fnpkg *pkg;          /* the active function package */
    gchar *pkgname;      /* and its name */
    gchar *pkgver;       /* plus its version */
    int *publist;        /* list of public interfaces */
    int iface;           /* selected interface */
    int flags;           /* misc. info on package */
    const ufunc *func;   /* the function we're calling */
    FuncDataReq dreq;    /* the function's data requirement */
    int minver;          /* minimum gretl version for pkg */
    int n_params;        /* its number of parameters */
    char rettype;        /* its return type */
    gchar **args;        /* its arguments */
    gchar *ret;          /* return assignment name */
    gchar *label;        /* the function's label */
};

#define scalar_arg(t) (t == GRETL_TYPE_DOUBLE || t == GRETL_TYPE_SCALAR_REF)
#define series_arg(t) (t == GRETL_TYPE_SERIES || t == GRETL_TYPE_SERIES_REF)
#define matrix_arg(t) (t == GRETL_TYPE_MATRIX || t == GRETL_TYPE_MATRIX_REF)
#define bundle_arg(t) (t == GRETL_TYPE_BUNDLE || t == GRETL_TYPE_BUNDLE_REF)
#define array_arg(t)  (t == GRETL_TYPE_ARRAY  || t == GRETL_TYPE_ARRAY_REF)

#define AUTOLIST "LTmp___"
#define SELNAME "selected series"

static GtkWidget *open_fncall_dlg;
static gboolean close_on_OK = TRUE;

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo);
static void maybe_record_include (const char *pkgname, int model_id);

static gchar **glib_str_array_new (int n)
{
    gchar **S = g_malloc0(n * sizeof *S);

    return S;
}

static void glib_str_array_free (gchar **S, int n)
{
    if (S != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    g_free(S[i]);
	}
	g_free(S);
    }
}

static int caller_is_model_window (windata_t *vwin)
{
    if (vwin != NULL && 
	(vwin->role == VIEW_MODEL ||
	 vwin->role == VAR ||
	 vwin->role == VECM ||
	 vwin->role == SYSTEM) &&
	vwin->data != NULL) {
	return 1;
    }

    return 0;
}

static call_info *cinfo_new (fnpkg *pkg, windata_t *vwin)
{
    call_info *cinfo = mymalloc(sizeof *cinfo);

    if (cinfo == NULL) {
	return NULL;
    }

    cinfo->pkg = pkg;
    cinfo->pkgname = NULL;
    cinfo->pkgver = NULL;

    cinfo->vwin = vwin;
    cinfo->dlg = NULL;
    cinfo->top_hbox = NULL;

    cinfo->publist = NULL;
    cinfo->iface = -1;
    cinfo->flags = 0;

    if (caller_is_model_window(vwin)) {
	cinfo->flags |= MODEL_CALL;
    }

    cinfo->vsels = NULL;
    cinfo->lsels = NULL;
    cinfo->msels = NULL;
    cinfo->bsels = NULL;
    cinfo->asels = NULL;
    cinfo->ssels = NULL;

    cinfo->func = NULL;
    cinfo->n_params = 0;

    cinfo->rettype = GRETL_TYPE_NONE;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    cinfo->dreq = 0;
    cinfo->label = NULL;

    return cinfo;
}

static int cinfo_args_init (call_info *cinfo)
{
    int err = 0;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    if (cinfo->n_params > 0) {
	cinfo->args = glib_str_array_new(cinfo->n_params);
	if (cinfo->args == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static void cinfo_free (call_info *cinfo)
{
    if (cinfo->n_params > 0) {
	glib_str_array_free(cinfo->args, cinfo->n_params);
    }
    if (cinfo->ret != NULL) {
	g_free(cinfo->ret);
    }
    if (cinfo->vsels != NULL) {
	g_list_free(cinfo->vsels);
    }
    if (cinfo->lsels != NULL) {
	g_list_free(cinfo->lsels);
    }
    if (cinfo->msels != NULL) {
	g_list_free(cinfo->msels);
    }
    if (cinfo->bsels != NULL) {
	g_list_free(cinfo->bsels);
    }
    if (cinfo->asels != NULL) {
	g_list_free(cinfo->asels);
    }
    if (cinfo->ssels != NULL) {
	g_list_free(cinfo->ssels);
    }

    g_free(cinfo->pkgname);
    g_free(cinfo->pkgver);

    g_free(cinfo->label);
    free(cinfo->publist);
    free(cinfo);
}

static int check_args (call_info *cinfo)
{
    int i;

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    if (cinfo->args[i] == NULL) {
		if (fn_param_optional(cinfo->func, i)) {
		    cinfo->args[i] = g_strdup("null");
		} else {
		    errbox_printf(_("Argument %d (%s) is missing"), i + 1,
				  fn_param_name(cinfo->func, i));
		    return 1;
		}
	    }
	}
    }

    return 0;
}

static void fncall_dialog_destruction (GtkWidget *w, call_info *cinfo)
{
    cinfo_free(cinfo);
    open_fncall_dlg = NULL;
}

static void fncall_close (GtkWidget *w, call_info *cinfo)
{
    gtk_widget_destroy(cinfo->dlg);
}

static GtkWidget *label_hbox (call_info *cinfo, GtkWidget *w)
{
    GtkWidget *hbox, *lbl;
    gchar *buf = NULL;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 5);

    if (cinfo->label != NULL) {
	buf = g_markup_printf_escaped("<span weight=\"bold\">%s</span>",
				      _(cinfo->label));
    } else {
	const char *funcname;
	
	funcname = user_function_name_by_index(cinfo->iface);
	buf = g_markup_printf_escaped("<span weight=\"bold\">%s</span>",
				      funcname);
    }
    
    lbl = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(lbl), buf);
    g_free(buf);

    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_widget_show(lbl);

    return hbox;
}

static gboolean update_double_arg (GtkWidget *w, call_info *cinfo)
{
    double val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%g", val);

    return FALSE;
}

static gboolean update_int_arg (GtkWidget *w, call_info *cinfo)
{
    int val = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);

    return FALSE;
}

static gboolean update_bool_arg (GtkWidget *w, call_info *cinfo)
{
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	cinfo->args[i] = g_strdup("1");
    } else {
	cinfo->args[i] = g_strdup("0");
    }

    return FALSE;
}

static gchar *combo_box_get_trimmed_text (GtkComboBox *combo)
{
    gchar *s = combo_box_get_active_text(combo);
    gchar *ret = NULL;

    if (s != NULL && *s != '\0') {
	while (isspace(*s)) s++;
	if (*s != '\0') {
	    int i, len = strlen(s);

	    for (i=len-1; i>0; i--) {
		if (!isspace(s[i])) break;
		len--;
	    }

	    if (len > 0) {
		ret = g_strndup(s, len);
	    }
	}
    }

    g_free(s);

    return ret;
}

static gboolean update_arg (GtkComboBox *combo, 
			    call_info *cinfo)
{
    int i = widget_get_int(combo, "argnum");
    char *s;

    g_free(cinfo->args[i]);
    s = cinfo->args[i] = combo_box_get_trimmed_text(combo);

    if (s != NULL && fn_param_type(cinfo->func, i) == GRETL_TYPE_DOUBLE) {
	if (isdigit(*s) || *s == '-' || *s == '+' || *s == ',') {
	    gretl_charsub(s, ',', '.');
	}
    }

    return FALSE;
}

static gboolean update_return (GtkComboBox *combo, 
			       call_info *cinfo)
{
    g_free(cinfo->ret);
    cinfo->ret = combo_box_get_trimmed_text(combo);

    return FALSE;
}

/* simple heuristic for whether or not a series probably 
   represents a stochastic variable
*/

static int probably_stochastic (int v)
{
    int ret = 1;

    if (sample_size(dataset) >= 3) {
	/* rule our vars that seem to be integer-valued with
	   a constant increment */
	int t = dataset->t1;
	double d1 = dataset->Z[v][t+1] - dataset->Z[v][t];
	double d2 = dataset->Z[v][t+2] - dataset->Z[v][t+1];

	if (d1 == floor(d1) && d2 == d1) {
	    ret = 0;
	}
    }

    return ret;
}

static GList *add_names_for_type (GList *list, GretlType type)
{
    GList *tlist = user_var_names_for_type(type);
    GList *tail = tlist;

    while (tail != NULL) {
	list = g_list_append(list, tail->data);
	tail = tail->next;
    }

    if (type == GRETL_TYPE_LIST && mdata_selection_count() > 1) {
	list = g_list_append(list, SELNAME);
    }

    g_list_free(tlist);

    return list;
}

static GList *add_series_names (GList *list)
{
    int i;

    for (i=1; i<dataset->v; i++) {
	if (!series_is_hidden(dataset, i)) {
	    list = g_list_append(list, (gpointer) dataset->varname[i]);
	} 
    }

    list = g_list_append(list, (gpointer) dataset->varname[0]);

    return list;
}

static GList *get_selection_list (int type)
{
    GList *list = NULL;

    if (series_arg(type)) {
	list = add_series_names(list);
    } else if (scalar_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_DOUBLE);
    } else if (type == GRETL_TYPE_LIST) {
	list = add_names_for_type(list, GRETL_TYPE_LIST);
    } else if (matrix_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_MATRIX);
    } else if (bundle_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_BUNDLE);
    } else if (array_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_ARRAY);
    }

    return list;
}

static windata_t *make_help_viewer (const char *fnname,
				    const char *pdfname,
				    PRN *prn)
{
    windata_t *vwin;
    gchar *title;
    char *buf;

    if (pdfname != NULL) {
	/* append a link to the PDF file */
	gchar *localpdf = g_strdup(pdfname);
	gchar *p = strrchr(localpdf, '.');

	*p = '\0';
	strncat(p, ".pdf", 4);
	pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n", localpdf);
	g_free(localpdf);	
    }

    buf = gretl_print_steal_buffer(prn);
    title = g_strdup_printf(_("help on %s"), fnname);
    vwin = view_formatted_text_buffer(title, buf, 76, 350, VIEW_PKG_INFO);
    g_free(title);
    free(buf);

    return vwin;
}

static void fncall_help (GtkWidget *w, call_info *cinfo)
{
    char *pdfname = NULL;
    int show_ghlp = 0;
    int have_pdf = 0;

    if ((cinfo->flags & SHOW_GUI_MAIN) &&
	function_package_has_gui_help(cinfo->pkg)) {
	show_ghlp = 1;
    }

    have_pdf = function_package_has_PDF_doc(cinfo->pkg, &pdfname);

    if (have_pdf && !show_ghlp) {
	/* simple: just show PDF doc */
	FILE *fp = gretl_fopen(pdfname, "r");

	if (fp != NULL) {
	    fclose(fp);
	    gretl_show_pdf(pdfname, NULL);
	} else {
	    gui_errmsg(E_FOPEN);
	}
    } else {
	/* show help text, either "plain" or GUI */
	const char *fnname;
	PRN *prn = NULL;
	gretlopt opt = OPT_M;
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	if (show_ghlp) {
	    opt |= OPT_G;
	}

	fnname = user_function_name_by_index(cinfo->iface);
	err = user_function_help(fnname, opt, prn);

	if (err) {
	    gretl_print_destroy(prn);
	    errbox("Couldn't find any help");
	} else {
	    make_help_viewer(fnname, pdfname, prn);
	    gretl_print_destroy(prn);
	}
    }

    free(pdfname);
}

static int combo_list_index (const gchar *s, GList *list)
{
    GList *mylist = list;
    int i;

    for (i=0; mylist != NULL; i++) {
	if (!strcmp(s, (gchar *) mylist->data)) {
	    return i;
	}
	mylist = mylist->next;
    }
    
    return -1;
}

/* Update the combo argument selector(s) for series, matrices
   lists or scalars after defining a new variable of one of 
   these types.
*/

static void update_combo_selectors (call_info *cinfo, 
				    GtkWidget *refsel,
				    int ptype)
{
    GList *sellist, *newlist;
    int llen;

    /* get the list of relevant selectors and the
       list of relevant variables to put into the
       selectors
    */

    if (ptype == GRETL_TYPE_MATRIX) {
	sellist = g_list_first(cinfo->msels);
    } else if (ptype == GRETL_TYPE_BUNDLE) {
	sellist = g_list_first(cinfo->bsels);
    } else if (ptype == GRETL_TYPE_ARRAY) {
	sellist = g_list_first(cinfo->asels);
    } else if (ptype == GRETL_TYPE_LIST) {
	sellist = g_list_first(cinfo->lsels);
    } else if (ptype == GRETL_TYPE_DOUBLE) {
	sellist = g_list_first(cinfo->ssels);
    } else {
	sellist = g_list_first(cinfo->vsels);
    }

    newlist = get_selection_list(ptype);
    llen = g_list_length(newlist);

    while (sellist != NULL) {
	/* iterate over the affected selectors */
	GtkComboBox *sel = GTK_COMBO_BOX(sellist->data);
	int target = GTK_WIDGET(sel) == refsel;
	int null_OK, selpos;
	gchar *saved = NULL;

	/* target == 1 means that we're looking at the
	   selector whose button was clicked to add a
	   variable: for this selector the newly added
	   variable should be marked as selected;
	   otherwise we modify the list of choices but
	   preserve the previous selection.
	*/

	if (!target) {
	    /* make a record of the old selected item */
	    saved = combo_box_get_active_text(sel);
	} 

	depopulate_combo_box(sel);
	set_combo_box_strings_from_list(GTK_WIDGET(sel), newlist);
	null_OK = widget_get_int(sel, "null_OK");
	if (null_OK) {
	    combo_box_append_text(sel, "null");
	}

	if (target) {
	    /* select the newly added var, which will be at the 
	       end, or thereabouts */
	    selpos = llen - 1;
	    if (series_arg(ptype)) {
		selpos--; /* the const is always in last place */
	    } 
	    gtk_combo_box_set_active(sel, selpos);
	} else if (saved != NULL) {
	    /* reinstate the previous selection */
	    selpos = combo_list_index(saved, newlist);
	    if (selpos < 0) {
		if (*saved == '\0') {
		    combo_box_prepend_text(sel, "");
		    selpos = 0;
		} else if (!strcmp(saved, "null")) {
		    selpos = llen;
		}
	    }
	    gtk_combo_box_set_active(sel, selpos);
	    g_free(saved);
	} else {
	    /* reinstate empty selection */
	    gtk_combo_box_set_active(sel, -1);
	} 

	sellist = sellist->next;
    }

    g_list_free(newlist);
}

static int do_make_list (selector *sr)
{
    const char *buf = selector_list(sr);
    const char *lname = selector_entry_text(sr);
    gpointer data = selector_get_data(sr);
    call_info *cinfo = NULL;
    GtkWidget *aux = NULL;
    const char *msg = NULL;
    PRN *prn = NULL;
    int *list = NULL;
    int empty = 0;
    int nl, err = 0;

    if (lname == NULL || *lname == '\0') {
	errbox(_("No name was given for the list"));
	return 1;
    }

    if (data != NULL) {
	/* called from elsewhere in fncall.c */
	GtkWidget *entry = GTK_WIDGET(data);

	cinfo = g_object_get_data(G_OBJECT(entry), "cinfo");
	aux = gtk_widget_get_parent(entry);
    }

    /* record initial status */
    nl = n_user_lists();

    if (buf == NULL || *buf == '\0') {
	int resp;

	resp = yes_no_dialog("gretl", _("Really create an empty list?"),
			     selector_get_window(sr));
	if (resp == GRETL_YES) {
	    list = gretl_null_list();
	    if (list == NULL) {
		err = E_ALLOC;
	    } else {
		empty = 1;
	    }
	} else {
	    /* canceled */
	    return 0;
	}
    } else {
	list = command_list_from_string(buf, &err);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    if (cinfo != NULL) {
	/* don't bother with "Added list..." message */
	err = remember_list(list, lname, NULL);
	if (err) {
	    gui_errmsg(err);
	}
    } else {
	if (bufopen(&prn)) {
	    free(list);
	    return 1;
	}
	err = remember_list(list, lname, prn);
	msg = gretl_print_get_buffer(prn);
	if (err) {
	    errbox(msg);
	}
    }

    if (!err) {
	lib_command_sprintf("list %s =%s", lname, empty ? " null" : buf);
	record_command_verbatim();
	gtk_widget_hide(selector_get_window(sr));
	if (cinfo != NULL) {
	    if (n_user_lists() > nl) {
		update_combo_selectors(cinfo, aux, GRETL_TYPE_LIST);
	    } 
	} else {
	    infobox(msg);
	}
    }

    free(list);

    if (prn != NULL) {
	gretl_print_destroy(prn);
    }

    return err;
} 

static void launch_list_maker (GtkWidget *button, GtkWidget *entry)
{
    call_info *cinfo;

    cinfo = g_object_get_data(G_OBJECT(button), "cinfo");
    g_object_set_data(G_OBJECT(cinfo->dlg), "button", button);
    simple_selection_with_data(DEFINE_LIST, _("Define list"), 
			       do_make_list, cinfo->dlg, 
			       entry);
}

void gui_define_list (void)
{
    simple_selection_with_data(DEFINE_LIST, _("Define list"), 
			       do_make_list, NULL, NULL);
}

static void launch_matrix_maker (GtkWidget *button, call_info *cinfo)
{
    int n = n_user_matrices();

    if (!strcmp(cinfo->pkgname, "SVAR")) {
	widget_set_int(cinfo->dlg, "matrix-no-series", 1);
    }

    g_object_set_data(G_OBJECT(cinfo->dlg), "button", button);
    fncall_add_matrix(cinfo->dlg);

    if (n_user_matrices() > n) {
	GtkWidget *sel = g_object_get_data(G_OBJECT(button), "combo");

	update_combo_selectors(cinfo, sel, GRETL_TYPE_MATRIX);
    }

    gtk_window_present(GTK_WINDOW(cinfo->dlg));
}

/* callback after invoking "genr" via the "+" button
   beside a combo argument selector */

void fncall_register_genr (int addv, gpointer p)
{
    GtkWidget *combo = p;
    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(combo));
    call_info *cinfo = g_object_get_data(G_OBJECT(entry), "cinfo");
    int ptype = widget_get_int(combo, "ptype");

    if (addv > 0) {
	update_combo_selectors(cinfo, combo, ptype);
    }

    gtk_window_present(GTK_WINDOW(cinfo->dlg));
}

static void launch_series_maker (GtkWidget *button, call_info *cinfo)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(button), "combo");

    edit_dialog(GENR, _("add series"), 
		_("Enter name=formula for new series"), NULL,
		do_fncall_genr, combo, 
		VARCLICK_INSERT_NAME, cinfo->dlg);  
}

static void launch_scalar_maker (GtkWidget *button, call_info *cinfo)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(button), "combo");

    edit_dialog(GENR, _("add scalar"), 
		_("Enter name=formula for new scalar"), NULL,
		do_fncall_genr, combo, 
		VARCLICK_INSERT_NAME, cinfo->dlg);  
}

static GtkWidget *bool_arg_selector (call_info *cinfo, int i,
				     const char *prior_val)
{
    GtkWidget *button;
    int active;

    if (prior_val != NULL) {
	active = *prior_val == '1';
    } else {
	double deflt = fn_param_default(cinfo->func, i);
	
	active = !na(deflt) && deflt != 0.0;
    }

    button = gtk_check_button_new();
    widget_set_int(button, "argnum", i);
    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(update_bool_arg), cinfo);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), active);
    cinfo->args[i] = g_strdup((active)? "1" : "0");

    return button;
}

static void update_xlist_arg (GtkComboBox *combo,
			      call_info *cinfo)
{
    MODEL *pmod = cinfo->vwin->data;
    int i = widget_get_int(combo, "argnum");
    int k = gtk_combo_box_get_active(combo);

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", k + 1 + pmod->ifc);
}

static GtkWidget *xlist_int_selector (call_info *cinfo, int i)
{
    MODEL *pmod;
    int *xlist;
    GtkWidget *combo;

    if (cinfo->vwin == NULL || cinfo->vwin->data == NULL) {
	return NULL;
    }

    pmod = cinfo->vwin->data;

    combo = gtk_combo_box_text_new();
    widget_set_int(combo, "argnum", i);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_xlist_arg), cinfo);

    xlist = gretl_model_get_x_list(pmod);

    if (xlist != NULL) {
	const char *s;
	int i, vi;

	for (i=1; i<=xlist[0]; i++) {
	    vi = xlist[i];
	    if (vi > 0) {
		s = dataset->varname[xlist[i]];
		combo_box_append_text(combo, s);
	    }
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	free(xlist);
    } 

    return combo;
}

static void update_enum_arg (GtkComboBox *combo, call_info *cinfo)
{
    int val = gtk_combo_box_get_active(combo);
    int i = widget_get_int(combo, "argnum");
    
    val += widget_get_int(combo, "minv");
    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);
}

static GtkWidget *enum_arg_selector (call_info *cinfo, int i,
				     const char **S, int nvals,
				     int minv, int initv)
{
    GtkWidget *combo;
    int j;

    combo = gtk_combo_box_text_new();
    widget_set_int(combo, "argnum", i);
    widget_set_int(combo, "minv", minv);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_enum_arg), cinfo);
    for (j=0; j<nvals; j++) {
	combo_box_append_text(combo, (const char *) S[j]);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), initv - minv);    

    return combo;
}

static GtkWidget *spin_arg_selector (call_info *cinfo, int i,
				     int minv, int maxv, int initv, 
				     GretlType type)
{
    GtkAdjustment *adj;
    GtkWidget *spin;

    adj = (GtkAdjustment *) gtk_adjustment_new(initv, minv, maxv, 
					       1, 1, 0);
    if (type == GRETL_TYPE_OBS) {
	spin = obs_button_new(adj, dataset, 0);
    } else {
	spin = gtk_spin_button_new(adj, 1, 0);
    }
    widget_set_int(spin, "argnum", i);
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed", 
		     G_CALLBACK(update_int_arg), cinfo);

    cinfo->args[i] = g_strdup_printf("%d", initv);

    return spin;
}

static GtkWidget *int_arg_selector (call_info *cinfo,
				    int i, GretlType type,
				    const char *prior_val)
{
    double dminv = fn_param_minval(cinfo->func, i);
    double dmaxv = fn_param_maxval(cinfo->func, i);
    double deflt = fn_param_default(cinfo->func, i);
    int minv, maxv, initv = 0;

    if (type == GRETL_TYPE_OBS) {
	/* the incoming vals will be 1-based */
	minv = (na(dminv) || dminv < 1)? 0 : (int) dminv - 1;
	maxv = (na(dmaxv) || dmaxv > dataset->n)? 
	    (dataset->n - 1) : (int) dmaxv - 1;
    } else {
	minv = (na(dminv))? INT_MIN : (int) dminv;
	maxv = (na(dmaxv))? INT_MAX : (int) dmaxv;
    }

    if (prior_val != NULL) {
	initv = atoi(prior_val);
    } else if (!na(deflt)) {
	initv = (int) deflt;
    } else if (!na(dminv)) {
	initv = (int) dminv;
    } 

    if (type == GRETL_TYPE_INT && !na(dminv) && !na(dmaxv)) {
	const char **S;
	int nvals;

	S = fn_param_value_labels(cinfo->func, i, &nvals);
	if (S != NULL) {
	    return enum_arg_selector(cinfo, i, S, nvals, minv, initv);
	}
    }

    return spin_arg_selector(cinfo, i, minv, maxv, initv, type);
}

static GtkWidget *double_arg_selector (call_info *cinfo, int i,
				       const char *prior_val)
{
    double minv  = fn_param_minval(cinfo->func, i);
    double maxv  = fn_param_maxval(cinfo->func, i);
    double deflt = fn_param_default(cinfo->func, i);
    double step  = fn_param_step(cinfo->func, i);
    GtkAdjustment *adj;
    GtkWidget *spin;
    gchar *p, *tmp;
    int ndec = 0;

    tmp = g_strdup_printf("%g", maxv - step);
    p = strchr(tmp, '.');
    if (p == NULL) {
	p = strchr(tmp, ',');
    }
    if (p != NULL) {
	ndec = strlen(p + 1);
    }
    g_free(tmp);

    if (prior_val != NULL) {
	/* locale? */
	deflt = atof(prior_val);
    }

    if (deflt > maxv) {
	/* note that default may be NADBL */
	deflt = minv;
    }

    adj = (GtkAdjustment *) gtk_adjustment_new(deflt, minv, maxv, 
					       step, step, 0);
    spin = gtk_spin_button_new(adj, 1, ndec);
    widget_set_int(spin, "argnum", i);
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed", 
		     G_CALLBACK(update_double_arg), cinfo);
    cinfo->args[i] = g_strdup_printf("%g", deflt);

    return spin;
}

/* see if the variable named @name of type @ptype
   has already been set as the default argument in
   an "upstream" combo argument selector
*/

static int already_set_as_default (call_info *cinfo,
				   const char *name,
				   int ptype)
{
    GList *slist;
    int ret = 0;

    if (series_arg(ptype)) {
	slist = g_list_first(cinfo->vsels);
    } else if (matrix_arg(ptype)) {
	slist = g_list_first(cinfo->msels);
    } else if (bundle_arg(ptype)) {
	slist = g_list_first(cinfo->bsels);
    } else if (array_arg(ptype)) {
	slist = g_list_first(cinfo->asels);
    } else if (scalar_arg(ptype)) {
	slist = g_list_first(cinfo->ssels);
    } else if (ptype == GRETL_TYPE_LIST) {
	slist = g_list_first(cinfo->lsels);
    } else {
	return 0;
    }

    while (slist != NULL && !ret) {
	GtkComboBox *sel = GTK_COMBO_BOX(slist->data);
	gchar *s = combo_box_get_active_text(sel);

	if (!strcmp(s, name)) {
	    ret = 1;
	} else {
	    slist = g_list_next(slist);
	}
	g_free(s);
    }

    return ret;
}

static int has_single_arg_of_type (call_info *cinfo,
				   GretlType type)
{
    int i, n = 0;

    for (i=0; i<cinfo->n_params; i++) {
	if (fn_param_type(cinfo->func, i) == type) {
	    n++;
	}
    }

    return n == 1;
}

/* Try to be somewhat clever in selecting the default values to show
   in function-argument drop-down "combo" selectors.

   Heuristics: (a) when a series is wanted, it's more likely to be a
   stochastic series rather than (e.g.) a time trend or panel group
   variable, so we try to avoid the latter as defaults; and (b) it's
   unlikely that the user wants to select the same named variable in
   more than one argument slot, so we try to avoid setting duplicate
   defaults.

   Special case: the function has exactly one series argument, and
   a single series is selected in the main gretl window: in that
   case we pre-select that series.
*/

static void arg_combo_set_default (call_info *cinfo,
				   GtkComboBox *combo, 
				   GList *list,
				   int ptype)
{
    GList *mylist = g_list_first(list);
    const char *targname = NULL;
    int i, v, k = 0;

    if (ptype == GRETL_TYPE_SERIES) {
	if (has_single_arg_of_type(cinfo, ptype)) {
	    v = mdata_active_var();
	    if (v > 0 && probably_stochastic(v)) {
		targname = dataset->varname[v];
	    }
	}
    } else if (ptype == GRETL_TYPE_LIST) {
	if (has_single_arg_of_type(cinfo, ptype) &&
	    mdata_selection_count() > 1) {
	    targname = SELNAME;
	    mylist = g_list_prepend(mylist, SELNAME);
	}
    }    

    for (i=0; mylist != NULL; i++) {
	gchar *name = mylist->data;
	int ok = 0;

	if (targname != NULL) {
	    ok = strcmp(name, targname) == 0;
	} else if (series_arg(ptype)) {
	    v = current_series_index(dataset, name);
	    if (v > 0 && probably_stochastic(v)) {
		ok = !already_set_as_default(cinfo, name, ptype);
	    }
	} else {
	    ok = !already_set_as_default(cinfo, name, ptype);
	}

	if (ok) {
	    k = i;
	    break;
	} else {
	    mylist = g_list_next(mylist);
	}
    }
 
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), k);
}

/* create an argument selector widget in the form of a
   GtkComboBox, with an entry field plus a drop-down
   list (which may initally be empty)
*/

static GtkWidget *combo_arg_selector (call_info *cinfo, int ptype, int i,
				      const char *prior_val)
{
    GList *list = NULL;
    GtkWidget *combo;
    GtkWidget *entry;
    int k = 0, null_OK = 0;

    combo = combo_box_text_new_with_entry();
    entry = gtk_bin_get_child(GTK_BIN(combo));
    g_object_set_data(G_OBJECT(entry), "cinfo", cinfo);
    widget_set_int(combo, "argnum", i);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_arg), cinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    if (fn_param_optional(cinfo->func, i)) {
	null_OK = 1;
	widget_set_int(combo, "null_OK", 1);
    }

    list = get_selection_list(ptype);
    if (list != NULL) {
	set_combo_box_strings_from_list(combo, list);
	arg_combo_set_default(cinfo, GTK_COMBO_BOX(combo), 
			      list, ptype);
	k = g_list_length(list);
	g_list_free(list);
    } 

    if (null_OK) {
	combo_box_append_text(combo, "null");
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), k);	
    }

    if (prior_val != NULL) {
	gtk_entry_set_text(GTK_ENTRY(entry), prior_val);
    } else if (ptype == GRETL_TYPE_INT) {
	double x = fn_param_default(cinfo->func, i);

	if (!na(x)) {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	    g_free(tmp);
	} 
    } else if (ptype == GRETL_TYPE_DOUBLE &&
	       fn_param_has_default(cinfo->func, i)) {
	double x = fn_param_default(cinfo->func, i);

	if (na(x)) {
	    gtk_entry_set_text(GTK_ENTRY(entry), "NA");
	} else {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	    g_free(tmp);
	}
    }

    return combo;
}

static void add_table_header (GtkWidget *tbl, gchar *txt,
			      int cols, int r0, int ypad)
{
    GtkWidget *label = gtk_label_new(txt);
    GtkWidget *align = gtk_alignment_new(0.0, 0.5, 0.0, 0.0);

    gtk_container_add(GTK_CONTAINER(align), label);
    gtk_table_attach(GTK_TABLE(tbl), align, 0, cols, r0, r0 + 1,
		     GTK_FILL, GTK_FILL, 5, ypad);
}

static void add_table_cell (GtkWidget *tbl, GtkWidget *w,
			    int c0, int c1, int r0)
{
    gtk_table_attach(GTK_TABLE(tbl), w, c0, c1, r0, r0 + 1,
		     GTK_FILL, GTK_FILL, 5, 3);
}

static GtkWidget *add_object_button (int ptype, GtkWidget *combo,
				     const char *parname)
{
    GtkWidget *img = gtk_image_new_from_stock(GTK_STOCK_ADD, 
					      GTK_ICON_SIZE_MENU);
    GtkWidget *button = gtk_button_new();

    gtk_container_add(GTK_CONTAINER(button), img);
    g_object_set_data(G_OBJECT(button), "combo", combo);
    if (parname != NULL) {
	/* FIXME is the cast OK here? */
	g_object_set_data(G_OBJECT(button), "parname", (char *) parname);
    }

    if (series_arg(ptype)) {
	gretl_tooltips_add(button, _("New variable"));
    } else if (matrix_arg(ptype)) {
	gretl_tooltips_add(button, _("Define matrix"));
    } else if (ptype == GRETL_TYPE_LIST) {
	gretl_tooltips_add(button, _("Define list"));
    } 

    return button;
}

static int spinnable_scalar_arg (call_info *cinfo, int i)
{
    const ufunc *func = cinfo->func;
    double mi = fn_param_minval(func, i);
    double ma = fn_param_maxval(func, i);
    double s = fn_param_step(func, i);

    return !na(mi) && !na(ma) && !na(s);
}

static void set_close_on_OK (GtkWidget *b, gpointer p)
{
    close_on_OK = button_is_active(b);
}

static int cinfo_show_return (call_info *c)
{
    if (c->rettype == GRETL_TYPE_NONE ||
	c->rettype == GRETL_TYPE_VOID) {
	return 0;
    } else if (c->rettype == GRETL_TYPE_BUNDLE &&
	       (c->flags & SHOW_GUI_MAIN)) {
	return 0;
    } else {
	return 1;
    }
}

static gchar *cinfo_pkg_title (call_info *cinfo)
{
    return g_strdup_printf("gretl: %s %s", cinfo->pkgname,
			   cinfo->pkgver);
}

static void function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *sel, *tbl = NULL;
    GtkWidget *vbox, *hbox, *bbox;
    arglist *alist = NULL;
    gchar *txt;
    int trows = 0, tcols = 0;
    int show_ret;
    int i, row;
    int err;

    if (open_fncall_dlg != NULL) {
	gtk_window_present(GTK_WINDOW(open_fncall_dlg));
	return;
    }

    err = cinfo_args_init(cinfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    cinfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    txt = cinfo_pkg_title(cinfo);
    gtk_window_set_title(GTK_WINDOW(cinfo->dlg), txt);
    g_free(txt);
    gretl_emulated_dialog_add_structure(cinfo->dlg, &vbox, &bbox);
    open_fncall_dlg = cinfo->dlg;
    g_signal_connect(G_OBJECT(cinfo->dlg), "destroy",
		     G_CALLBACK(fncall_dialog_destruction), cinfo);

    /* above table: label or name of function being called */
    cinfo->top_hbox = hbox = label_hbox(cinfo, vbox);

    show_ret = cinfo_show_return(cinfo);

    if (cinfo->n_params > 0) {
	tcols = 3; /* label, selector, add-button */
	trows = cinfo->n_params + 1;
	if (show_ret) { 
	    trows += 4;
	}
	alist = arglist_lookup(cinfo->pkgname);
    } else if (show_ret) {
	tcols = 2;
	trows = 3;
    }

    if (trows > 0 && tcols > 0) {
	tbl = gtk_table_new(trows, tcols, FALSE);
    }

    row = 0; /* initialize writing row */

    for (i=0; i<cinfo->n_params; i++) {
	const char *desc = fn_param_descrip(cinfo->func, i);
	const char *parname = fn_param_name(cinfo->func, i);
	const char *prior_val = NULL;
	int ptype = fn_param_type(cinfo->func, i);
	int spinnable = 0;
	gchar *argtxt;

	if (i == 0) {
	    add_table_header(tbl, _("Select arguments:"), tcols, row, 5);
	}

	if (alist != NULL) {
	    prior_val = arglist_lookup_val(alist, i);
	}

	row++;

	if (ptype == GRETL_TYPE_DOUBLE) {
	    spinnable = spinnable_scalar_arg(cinfo, i);
	}

	/* label for name (and maybe type) of argument, using
	   descriptive string if available */
	if (ptype == GRETL_TYPE_INT ||
	    ptype == GRETL_TYPE_BOOL ||
	    ptype == GRETL_TYPE_OBS ||
	    spinnable) {
	    argtxt = g_strdup_printf("%s",
				     (desc != NULL)? _(desc) :
				     parname);
	} else {
	    const char *astr = gretl_type_get_name(ptype);

	    if (desc != NULL && strstr(desc, astr)) {
		argtxt = g_strdup_printf("%s", _(desc));
	    } else {
		argtxt = g_strdup_printf("%s (%s)",
					 (desc != NULL)? _(desc) :
					 parname, astr);
	    }
	}

	label = gtk_label_new(argtxt);
	g_free(argtxt);
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, row);

	/* make the appropriate type of selector widget */

	if (fn_param_uses_xlist(cinfo->func, i)) {
	    sel = xlist_int_selector(cinfo, i);
	} else if (ptype == GRETL_TYPE_BOOL) {
	    sel = bool_arg_selector(cinfo, i, prior_val);
	} else if (ptype == GRETL_TYPE_INT ||
		   ptype == GRETL_TYPE_OBS) {
	    sel = int_arg_selector(cinfo, i, ptype, prior_val);
	} else if (spinnable) {
	    sel = double_arg_selector(cinfo, i, prior_val);
	} else {
	    sel = combo_arg_selector(cinfo, ptype, i, prior_val);
	}

	add_table_cell(tbl, sel, 1, 2, row);

	/* hook up signals and "+" add buttons for the
	   selectors for most types of arguments (though
	   not for bool and spinner-type args)
	*/

	if (series_arg(ptype)) {
	    cinfo->vsels = g_list_append(cinfo->vsels, sel);
	    widget_set_int(sel, "ptype", GRETL_TYPE_SERIES);
	    button = add_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_series_maker), 
			     cinfo);
	} else if (scalar_arg(ptype) && !spinnable) {
	    cinfo->ssels = g_list_append(cinfo->ssels, sel);
	    widget_set_int(sel, "ptype", GRETL_TYPE_DOUBLE);
	    button = add_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_scalar_maker), 
			     cinfo);
	} else if (matrix_arg(ptype)) {
	    cinfo->msels = g_list_append(cinfo->msels, sel);
	    button = add_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_matrix_maker), 
			     cinfo);
	} else if (bundle_arg(ptype)) {
	    cinfo->bsels = g_list_append(cinfo->bsels, sel);
	} else if (array_arg(ptype)) {
	    cinfo->asels = g_list_append(cinfo->asels, sel);
	} else if (ptype == GRETL_TYPE_LIST) {
	    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(sel));
	    
	    cinfo->lsels = g_list_append(cinfo->lsels, sel);
	    button = add_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    widget_set_int(entry, "argnum", i);
	    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_list_maker),
			     entry);
	} 
    }
	
    if (show_ret) {
	/* selector/entry for return value */
	GtkWidget *child;
	GList *list = NULL;

	if (cinfo->n_params > 0) {
	    /* separator row */
	    add_table_header(tbl, "", tcols, ++row, 0);
	}	    

	add_table_header(tbl, _("Assign return value (optional):"), 
			 tcols, ++row, 5);

	label = gtk_label_new(_("selection (or new variable)"));
	add_table_cell(tbl, label, 1, 2, ++row);

	label = gtk_label_new(gretl_type_get_name(cinfo->rettype));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, ++row);

	sel = combo_box_text_new_with_entry();
	g_signal_connect(G_OBJECT(sel), "changed",
			 G_CALLBACK(update_return), cinfo);
	list = get_selection_list(cinfo->rettype);
	if (list != NULL) {
	    set_combo_box_strings_from_list(sel, list);
	    g_list_free(list);
	}

	/* prepend blank option and select it */
	combo_box_prepend_text(sel, "");
	gtk_combo_box_set_active(GTK_COMBO_BOX(sel), 0);
	child = gtk_bin_get_child(GTK_BIN(sel));
	gtk_entry_set_activates_default(GTK_ENTRY(child), TRUE);
	add_table_cell(tbl, sel, 1, 2, row); /* same row as above */
    }

    if (tbl != NULL) {
	/* the table is complete: pack it now */
	gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    }

    /* option button */
    hbox = gtk_hbox_new(FALSE, 5);
    button = gtk_check_button_new_with_label(_("close this dialog on \"OK\""));
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(set_close_on_OK), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 close_on_OK);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* Close button */
    button = close_button(bbox);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(fncall_close), cinfo);

    /* "OK" button */
    button = ok_button(bbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(fncall_exec_callback), cinfo);
    gtk_widget_grab_default(button);

    /* Help button */
    button = context_help_button(bbox, -1);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(fncall_help), cinfo);

    if (cinfo->vwin != NULL) {
	gtk_window_set_transient_for(GTK_WINDOW(cinfo->dlg), 
				     GTK_WINDOW(cinfo->vwin->main));
	gtk_window_set_destroy_with_parent(GTK_WINDOW(cinfo->dlg), TRUE);
    }

    gtk_widget_show_all(cinfo->dlg);
}

/* called when defining a matrix for use as an argument:
   @dlg will be the function call dialog
*/

void get_fncall_param_info (GtkWidget *dlg, int *series_ok,
			    char **pname)
{
    if (dlg == NULL) {
	return;
    }
    
    if (series_ok != NULL) {
	*series_ok = !widget_get_int(dlg, "matrix-no-series");
    }

    if (pname != NULL) {
	GtkWidget *button =
	    g_object_get_data(G_OBJECT(dlg), "button");

	if (button != NULL) {
	    char *name = g_object_get_data(G_OBJECT(button),
					   "parname");

	    if (name != NULL &&
		current_series_index(dataset, name) < 0 &&
		gretl_get_object_by_name(name) == NULL) {
		/* OK, not the name of current object, so
		   offer it as default */
		*pname = g_strdup(name);
	    }
	}
    }	
}

static int function_data_check (call_info *cinfo)
{
    int i, err = 0;

    if (cinfo->dreq != FN_NODATA_OK) {
	if (dataset == NULL || dataset->v == 0) {
	    warnbox(_("Please open a data file first"));
	    return 1;
	}
    }

    for (i=0; i<cinfo->n_params; i++) {
	int type = fn_param_type(cinfo->func, i);

	if (type == GRETL_TYPE_SERIES || type == GRETL_TYPE_LIST ||
	    type == GRETL_TYPE_SERIES_REF) {
	    if (dataset == NULL || dataset->v == 0) {
		warnbox(_("Please open a data file first"));
		err = 1;
		break;
	    }
	}
    }

    return err;
}

/* detect the case where we need a "pointer" variable but
   have been given a scalar or matrix constant 
*/

static int should_addressify_var (call_info *cinfo, int i)
{
    char *numchars = "0123456789+-.,";
    int t = fn_param_type(cinfo->func, i);
    gchar *s = cinfo->args[i];

    return (t == GRETL_TYPE_SCALAR_REF && strchr(numchars, *s)) ||
	(t == GRETL_TYPE_MATRIX_REF && *s == '{');
}

static int maybe_add_amp (call_info *cinfo, int i, PRN *prn, int *add)
{
    int t = fn_param_type(cinfo->func, i);
    gchar *name = cinfo->args[i];
    int err = 0;

    *add = 0;

    if (!gretl_ref_type(t)) {
	return 0;
    }

    if (*name == '&' || !strcmp(name, "null")) {
	return 0;
    }

    /* Handle cases where an "indirect return" variable 
       does not yet exist: we need to declare it.
    */

    if (t == GRETL_TYPE_MATRIX_REF) {
	if (get_matrix_by_name(name) == NULL) {
	    gretl_matrix *m = gretl_null_matrix_new();

	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = user_var_add_or_replace(name,
					      GRETL_TYPE_MATRIX,
					      m);
	    }
	    if (!err) {
		pprintf(prn, "? matrix %s\n", name);
	    }
	}
    } else if (t == GRETL_TYPE_SERIES_REF) {
	if (current_series_index(dataset, name) < 0) {
	    err = generate(name, dataset, GRETL_TYPE_SERIES,
			   OPT_Q, NULL);
	    if (!err) {
		pprintf(prn, "? series %s\n", name);
	    }
	}
    }	

    if (!err) {
	*add = 1;
    }

    return err;
}

static int needs_quoting (call_info *cinfo, int i)
{
    int t = fn_param_type(cinfo->func, i);
    gchar *s = cinfo->args[i];

    return (t == GRETL_TYPE_STRING && 
	    strcmp(s, "null") && 
	    get_string_by_name(s) == NULL &&
	    *s != '"');
}

static int pre_process_args (call_info *cinfo, int *autolist,
			     PRN *prn)
{
    char auxline[MAXLINE];
    char auxname[VNAMELEN+2];
    int i, add = 0, err = 0;

    for (i=0; i<cinfo->n_params && !err; i++) {
	if (should_addressify_var(cinfo, i)) {
	    sprintf(auxname, "FNARG%d", i + 1);
	    sprintf(auxline, "%s=%s", auxname, cinfo->args[i]);
	    err = generate(auxline, dataset, GRETL_TYPE_ANY,
			   OPT_NONE, NULL);
	    if (!err) {
		g_free(cinfo->args[i]);
		cinfo->args[i] = g_strdup(auxname);
		pprintf(prn, "? %s\n", auxline);
	    } 
	}
	
	err = maybe_add_amp(cinfo, i, prn, &add);
	
	if (add) {
	    strcpy(auxname, "&");
	    strncat(auxname, cinfo->args[i], VNAMELEN);
	    g_free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup(auxname);
	} else if (needs_quoting(cinfo, i)) {
	    sprintf(auxname, "\"%s\"", cinfo->args[i]);
	    g_free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup(auxname);
	}
	
	if (fn_param_type(cinfo->func, i) == GRETL_TYPE_OBS) {
	    /* convert integer value from 0- to 1-based */
	    int val = atoi(cinfo->args[i]) + 1;

	    g_free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup_printf("%d", val);
	} else if (fn_param_type(cinfo->func, i) == GRETL_TYPE_LIST) {
	    /* do we have an automatic list arg? */
	    if (!strcmp(cinfo->args[i], SELNAME)) {
		user_var_add(AUTOLIST, GRETL_TYPE_LIST,
			     main_window_selection_as_list());
		g_free(cinfo->args[i]);
		cinfo->args[i] = g_strdup(AUTOLIST);
		*autolist = i;
	    }
	}	
    }

    return err;
}

static void set_genr_model_from_vwin (windata_t *vwin)
{
    GretlObjType type = GRETL_OBJ_EQN;

    if (vwin->role == VAR || vwin->role == VECM) {
	type = GRETL_OBJ_VAR;
    } else if (vwin->role == SYSTEM) {
	type = GRETL_OBJ_SYS;
    }

    set_genr_model(vwin->data, type);
}

/* Compose the command line that calls a packaged function with the
   appropriate arguments, if any, and possible assignment of the
   return value, if any.
*/

static void compose_fncall_line (char *line, 
				 call_info *cinfo,
				 const char *funname,
				 char **tmpname,
				 int *grab_bundle)
{
    arglist *alist = arglist_lookup(cinfo->pkgname);

    if (alist == NULL) {
	alist = arglist_new(cinfo->pkgname, cinfo->n_params);
    }
    
    *line = '\0';
    
    if (cinfo->ret != NULL) {
	strcat(line, cinfo->ret);
	strcat(line, " = ");
    } else if (cinfo->rettype == GRETL_TYPE_BUNDLE) {
	/* the function offers a bundle return but this has not been
	   assigned by the user; make a special arrangement to grab
	   the bundle for GUI purposes
	*/
	*tmpname = temp_name_for_bundle();
	strcat(line, *tmpname);
	strcat(line, " = ");
	*grab_bundle = 1;
    }	

    strcat(line, funname);
    strcat(line, "(");

    if (cinfo->args != NULL) {
	int i;

	for (i=0; i<cinfo->n_params; i++) {
	    strcat(line, cinfo->args[i]);
	    if (alist != NULL) {
		arglist_record_arg(alist, i, cinfo->args[i]);
	    }
	    if (i < cinfo->n_params - 1) {
		strcat(line, ", ");
	    }
	}
    }

    strcat(line, ")");
}

static int real_GUI_function_call (call_info *cinfo, PRN *prn)
{
    windata_t *outwin = NULL;
    ExecState state;
    char fnline[MAXLINE];
    char *tmpname = NULL;
    const char *funname;
    const char *title;
    gretl_bundle *bundle = NULL;
    int grab_bundle = 0;
    int show = 1;
    int orig_v = dataset->v;
    int err = 0;

    funname = user_function_name_by_index(cinfo->iface);
    title = cinfo->label != NULL ? cinfo->label : funname;

    compose_fncall_line(fnline, cinfo, funname,
			&tmpname, &grab_bundle);

    /* FIXME: the following conditionality may be wrong? */

    if (!grab_bundle && strncmp(funname, "GUI", 3)) {
	pprintf(prn, "? %s\n", fnline);
    }

#if FCDEBUG
    fprintf(stderr, "fnline: %s\n", fnline);
#endif

    /* note: gretl_exec_state_init zeros the first byte of its
       'line' member
    */

    gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			  model, prn);
    state.line = fnline;

    if (cinfo->flags & MODEL_CALL) {
	set_genr_model_from_vwin(cinfo->vwin);
    }

    show = !user_func_is_noprint(cinfo->func);

#if FCDEBUG
    fprintf(stderr, "show = %d, grab_bundle = %d\n", show, grab_bundle);
#endif    

    if (close_on_OK && cinfo->dlg != NULL) {
	gtk_widget_hide(cinfo->dlg);
    }

    if (show) {
	/* allow the "flush" mechanism to operate */
	err = exec_line_with_output_handler(&state, dataset,
					    title, &outwin);
    } else {
	/* execute "invisibly" */
	err = gui_exec_line(&state, dataset, NULL);
    }

    if (!err && strstr(fnline, AUTOLIST) == NULL) {
	int ID = 0;

	if (cinfo->flags & MODEL_CALL) {
	    ID = get_genr_model_ID();
	}

	maybe_record_include(cinfo->pkgname, ID);
	
	if (ID > 0) {
	    lib_command_sprintf("# %s", fnline);
	    record_model_command_verbatim(ID);
	} else {
	    lib_command_strcpy(fnline);
	    record_command_verbatim();
	}
    }

    if (cinfo->flags & MODEL_CALL) {
	unset_genr_model();
    }

    if (!err && cinfo->rettype == GRETL_TYPE_BUNDLE) {
	if (grab_bundle) {
	    bundle = get_bundle_by_name(tmpname);
	    if (bundle != NULL && !gretl_bundle_has_content(bundle)) {
		/* we got a useless empty bundle */
		gretl_bundle_pull_from_stack(tmpname, &err);
		gretl_bundle_destroy(bundle);
		bundle = NULL;
	    }
	} else if (cinfo->ret != NULL) {
	    bundle = get_bundle_by_name(cinfo->ret);
	}
    }

    if (!err && bundle != NULL && !show) {
	gretl_print_reset_buffer(prn);
	if (try_exec_bundle_print_function(bundle, prn)) {
	    /* flag the fact that we do have something to show */
	    show = 1;
	}
    }

    if (grab_bundle && bundle != NULL) {
	/* If the user specified an assignment of a returned bundle,
	   we should leave it in the user_vars stack; but if we added
	   the assignment automatically, we should pull it out of the
	   stack, leaving the saving (or not) up to the user.
	*/
	gretl_bundle_pull_from_stack(tmpname, &err);
    }

    if (!err && !show) {
	gretl_print_destroy(prn);
    } else if (outwin == NULL) {
	/* output window not already in place */
	view_buffer(prn, 80, 400, title, 
		    (bundle == NULL)? PRINT : VIEW_BUNDLE,
		    bundle);
    } else {
	/* an output window has already been opened
	   via "flush" */
	if (bundle != NULL) {
	    finalize_script_output_window(VIEW_BUNDLE, bundle);
	} else {
	    finalize_script_output_window(0, NULL); 
	}
    }

    free(tmpname);

    if (err && !show) {
	gui_errmsg(err);
    } 

    if (dataset->v != orig_v) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    return err;
}

/* For interface selection via the GUI, when no gui-main
   is set: we'll suppose that if there's an interface with
   the same name as the package itself, it should probably
   be the first option on the list
*/

static void maybe_reshuffle_iface_order (call_info *cinfo)
{
    int *list = cinfo->publist;
    const char *s;
    int i, pref = 1;

    for (i=1; i<=list[0]; i++) {
	s = user_function_name_by_index(list[i]);
	if (!strcmp(s, cinfo->pkgname)) {
	    pref = i;
	    break;
	}
    }

    if (pref > 1) {
	int tmp = list[1];

	list[1] = list[pref];
	list[pref] = tmp;
    }
}

/* In case a function package offers more than one public
   interface, give the user a selector: for four or fewer
   options we use radio buttons, otherwise we use a pull-down
   list.
*/

static void pkg_select_interface (call_info *cinfo, int npub)
{
    const char *funname;
    char **opts = NULL;
    GList *ilist = NULL;
    int radios = (npub < 5);
    int i, nopts = 0;
    int err = 0;

    maybe_reshuffle_iface_order(cinfo);

    for (i=1; i<=npub && !err; i++) {
	funname = user_function_name_by_index(cinfo->publist[i]);
	if (funname == NULL) {
	    err = E_DATA;
	} else if (radios) {
	    err = strings_array_add(&opts, &nopts, funname);
	} else {
	    ilist = g_list_append(ilist, (gpointer) funname);
	}
    }

    if (err) {
	cinfo->iface = -1;
	gui_errmsg(err);
    } else {
	GtkWidget *parent = vwin_toplevel(cinfo->vwin);
	int resp;

	if (radios) {
	    gchar *title = g_strdup_printf("gretl: %s\n", cinfo->pkgname);

	    resp = radio_dialog(title, _("Select function"), 
				(const char **) opts, 
				nopts, 0, 0, parent);
	    if (resp >= 0) {
		cinfo->iface = cinfo->publist[resp+1];
	    } else {
		cinfo->iface = -1;
	    }
	    strings_array_free(opts, nopts);
	    g_free(title);
	} else {
	    resp = combo_selector_dialog(ilist, _("Select function"), 
					 0, parent);
	    if (resp >= 0) {
		cinfo->iface = cinfo->publist[resp+1];
	    } else {
		cinfo->iface = -1;
	    }	
	    g_list_free(ilist);
	}
    }
}

/* Callback from "OK" button in function call GUI: if there's a
   problem with the argument selection just return so the dialog stays
   in place and the user can correct matters; otherwise really execute
   the function.
*/

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo)
{
    if (check_args(cinfo)) {
	return;
    } else {
	PRN *prn = NULL;
	int autopos = -1;
	int err;

	err = bufopen(&prn);

	if (!err && cinfo->args != NULL) {
	    err = pre_process_args(cinfo, &autopos, prn);
	    if (err) {
		gui_errmsg(err);
	    }
	}

	if (!err) {
	    err = real_GUI_function_call(cinfo, prn);
	} else {
	    gretl_print_destroy(prn);
	}

	if (autopos >= 0) {
	    user_var_delete_by_name(AUTOLIST, NULL);
	    g_free(cinfo->args[autopos]);
	    cinfo->args[autopos] = g_strdup(SELNAME);
	}

	if (cinfo->dlg != NULL && close_on_OK) {
	    gtk_widget_destroy(cinfo->dlg);
	} else if (cinfo->dlg == NULL) {
	    cinfo_free(cinfo);
	}
    }
}

/* Here we're testing whether @pkg has a "gui-main" function
   that should be displayed as the default interface.
*/

static void maybe_set_gui_interface (call_info *cinfo,
				     int from_browser)
{
    int gmid = -1, fid = -1;

    function_package_get_properties(cinfo->pkg, 
				    "gui-main-id", &gmid,
				    NULL);

    if (gmid >= 0 && !from_browser) {
	fid = gmid;
    } else if (cinfo->publist[0] == 1) {
	/* single suitable interface: implicit gui-main */
	fid = cinfo->publist[1];
    } else if (gmid >= 0) {
	/* called from browser: check for masking */
	const ufunc *u = get_user_function_by_index(gmid);

	if (!user_func_is_menu_only(u)) {
	    fid = gmid;
	}
    }

    if (fid >= 0) {
	/* we found a usable gui-main */
	gchar *name = NULL, *label = NULL;
	
	cinfo->iface = fid;
	cinfo->flags |= SHOW_GUI_MAIN;
	function_package_get_properties(cinfo->pkg,
					"name", &name,
					"label", &label,
					NULL);
	if (label != NULL) {
	    cinfo->label = label;
	    g_free(name);
	} else {
	    cinfo->label = name;
	}
    }
}

static int need_model_check (call_info *cinfo)
{
    int i, err = 0;

    for (i=0; i<cinfo->n_params; i++) {
	if (fn_param_uses_xlist(cinfo->func, i)) {
	    if (cinfo->vwin == NULL || cinfo->vwin->role != VIEW_MODEL) {
		err = E_DATA;
		errbox(_("This function needs a model in place"));
		break;
	    }
	}
    }

    return err;
}

static call_info *start_cinfo_for_package (const char *pkgname,
					   const char *fname,
					   windata_t *vwin,
					   int *err)
{
    call_info *cinfo;
    fnpkg *pkg;

    pkg = get_function_package_by_name(pkgname);

    if (pkg == NULL) {
	/* not already loaded */
	pkg = get_function_package_by_filename(fname, err);
	if (*err) {
	    gui_errmsg(*err);
	    return NULL;
	}
    }

    cinfo = cinfo_new(pkg, vwin);
    if (cinfo == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* get the interface list and other basic info for package */

    *err = function_package_get_properties(pkg,
					   "name", &cinfo->pkgname,
					   "version", &cinfo->pkgver,
					   "gui-publist", &cinfo->publist,
					   "data-requirement", &cinfo->dreq,
					   "min-version", &cinfo->minver,
					   NULL);

    if (*err) {
	gui_errmsg(*err);
    } else if (cinfo->publist == NULL) {
	/* no available interfaces */
	errbox(_("Function package is broken"));
	*err = E_DATA;
    }

    if (*err) {
	cinfo_free(cinfo);
	cinfo = NULL;
    }

    return cinfo;
}

/* Call to execute a function from the package pre-attached to
   @cinfo. We may or may not end up offering a list of
   interfaces. This is called both from menu items (see below in this
   file) and (indirectly) from the "browser" window that lists
   installed function packages -- see open_function_package() below.
*/

static int call_function_package (call_info *cinfo, windata_t *vwin,
				  int from_browser)
{
    int err = 0;

    if (!from_browser) {
	/* Do we have suitable data in place? (This is already
	   checked if @from_browser is non-zero).
	*/	
	err = check_function_needs(dataset, cinfo->dreq, cinfo->minver);
	if (err) {
	    gui_errmsg(err);
	}
    }

    if (!err) {
	maybe_set_gui_interface(cinfo, from_browser);
    }

    if (!err && cinfo->iface < 0) {
	pkg_select_interface(cinfo, cinfo->publist[0]);
	if (cinfo->iface < 0) {
	    /* failed, or cancelled */
	    cinfo_free(cinfo);
	    return 0; /* note: handled */
	}
    }

    if (!err) {
	cinfo->func = get_user_function_by_index(cinfo->iface);
	if (cinfo->func == NULL) {
	    fprintf(stderr, "get_user_function_by_index: failed\n");
	    errbox(_("Couldn't get function package information"));
	}
    }	
    
    if (!err) {
	cinfo->n_params = fn_n_params(cinfo->func);
	err = function_data_check(cinfo);
    }

    if (!err) {
	cinfo->rettype = user_func_get_return_type(cinfo->func);
	if (err) {
	    fprintf(stderr, "user_func_get_return_type: failed\n");
	    errbox(_("Couldn't get function package information"));
	}
    }

    if (!err) {
	/* Should this check come earlier? */
	err = need_model_check(cinfo);
    }

    if (!err) {
	if (fn_n_params(cinfo->func) == 0) {
	    /* no arguments to be gathered */
	    fncall_exec_callback(NULL, cinfo);
	} else {
	    /* put up a dialog to collect arguments */
	    function_call_dialog(cinfo);
	}
    } else {
	cinfo_free(cinfo);
    }

    return err;
}

/* Called from the function-package browser: unless the
   package can't be loaded we should return 0 to signal
   that loading happened.
*/

int open_function_package (const char *pkgname,
			   const char *fname,
			   windata_t *vwin)
{
    call_info *cinfo;
    int can_call = 1;
    int free_cinfo = 1;
    int err = 0;

    /* note: this ensures the package gets loaded */
    cinfo = start_cinfo_for_package(pkgname, fname, vwin, &err);

    if (err) {
	return err;
    }

    /* do we have suitable data in place? */
    err = check_function_needs(dataset, cinfo->dreq, cinfo->minver);
    
    if (err == E_DATA) {
	/* we might still run the sample script */
	can_call = 0;
	err = 0;
	gretl_error_clear();
    } else if (err) {
	/* fatal error */
	gui_errmsg(err);
	return err;
    }

    if (can_call) {
	/* actually call the package: preserve @cinfo! */
	free_cinfo = 0;
	call_function_package(cinfo, vwin, 1);
    } else {
	/* notify and give choice of running sample */
	const char *ts_msg = N_("This package needs time series data.");
	const char *qm_msg = N_("This package needs quarterly or monthly data.");
	const char *pn_msg = N_("This package needs panel data.");
	const char *ds_msg = N_("This package needs a dataset in place.");
	const char *query = N_("Would you like to open its sample script?");
	gchar *msg, *title = cinfo_pkg_title(cinfo);
	const char *req;
	int resp;

	if (cinfo->dreq == FN_NEEDS_TS) {
	    req = ts_msg;
	} else if (cinfo->dreq == FN_NEEDS_QM) {
	    req = qm_msg;
	} else if (cinfo->dreq == FN_NEEDS_PANEL) {
	    req = pn_msg;
	} else {
	    req = ds_msg;
	}

	msg = g_strdup_printf("%s\n%s", _(req), _(query));
	resp = yes_no_dialog(title, msg, vwin_toplevel(vwin));
	if (resp == GRETL_YES) {
	    display_function_package_data(cinfo->pkgname, fname,
					  VIEW_PKG_SAMPLE);
	}
	g_free(msg);
    }

    if (free_cinfo) {
	cinfo_free(cinfo);
    }

    return 0;
}

void function_call_cleanup (void)
{
    if (open_fncall_dlg != NULL) {
	gtk_widget_destroy(open_fncall_dlg);
    }

    arglist_cleanup();
}

/* Execute the plotting function made available by the function
   package that produced bundle @b, possibly inflected by an 
   integer option -- if an option is present it's packed into 
   @aname, after a colon. 
*/

int exec_bundle_plot_function (gretl_bundle *b, const char *aname)
{
    ufunc *func;
    char funname[32];
    int iopt = -1;
    int err = 0;

    if (aname != NULL) {
	if (strchr(aname, ':') != NULL) {
	    /* extract option */
	    sscanf(aname, "%31[^:]:%d", funname, &iopt);
	} else {
	    /* name but no option present */
	    strcpy(funname, aname);
	}
    } else {
	gchar *pf = get_bundle_plot_function(b);

	if (pf == NULL) {
	    return E_DATA;
	} else {
	    strcpy(funname, pf);
	    g_free(pf);
	}
    }

    func = get_user_function_by_name(funname);

    if (func == NULL) {
	err = E_DATA;
    } else {
	const char *bname = user_var_get_name_by_data(b);
#if 1
	fprintf(stderr, "bundle plot: using bundle %p (%s)\n",
		(void *) b, bname);
#endif
	err = push_function_arg(func, bname, GRETL_TYPE_BUNDLE_REF, b);

	if (!err && iopt >= 0) {
	    /* add the option flag, if any, to args */
	    double minv = fn_param_minval(func, 1);

	    if (!na(minv)) {
		iopt += (int) minv;
		err = push_function_arg(func, NULL, GRETL_TYPE_INT, &iopt);
	    }
	}
    }

    if (!err) {
	/* Note that the function may need a non-NULL prn for
	   use with printing redirection (outfile).
	*/
	PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, &err);

	err = gretl_function_exec(func, GRETL_TYPE_NONE,
				  dataset, NULL, NULL, prn);
	gretl_print_destroy(prn);
    }

    if (err) {
	function_clear_args(func);
	gui_errmsg(err);
    } 

    return err;
}

/* See if a bundle has the name of a "creator" function package
   recorded on it. If so, see whether that package is already loaded,
   or can be loaded.  And if that works, see if the package has a
   default function for @task (e.g. BUNDLE_PRINT).
*/

static gchar *get_bundle_special_function (gretl_bundle *b,
					   const char *task)
{
    const char *pkgname = gretl_bundle_get_creator(b);
    gchar *ret = NULL;

    if (pkgname != NULL && *pkgname != '\0') {
	fnpkg *pkg = get_function_package_by_name(pkgname);

	if (pkg == NULL) {
	    char *fname = 
		gretl_function_package_get_path(pkgname, PKG_ALL);
	    int err = 0;

	    if (fname != NULL) {
		pkg = get_function_package_by_filename(fname, &err);
		free(fname);
	    }
	}

	if (pkg != NULL) {
	    function_package_get_properties(pkg, task, &ret, NULL);
	}
    }

    return ret;
}

gchar *get_bundle_plot_function (gretl_bundle *b)
{
    return get_bundle_special_function(b, BUNDLE_PLOT);
}

/* See if we can find a "native" printing function for a
   gretl bundle. If we can find this, try executing it.

   Notice that this function returns 1 on success, 0 on
   failure.
*/

int try_exec_bundle_print_function (gretl_bundle *b, PRN *prn)
{
    gchar *funname;
    int ret = 0;

    funname = get_bundle_special_function(b, BUNDLE_PRINT);

    if (funname != NULL) {
	const char *name = user_var_get_name_by_data(b);
	char fnline[MAXLINE];
	ExecState state;
	int err;

	sprintf(fnline, "%s(&%s)", funname, name);
	g_free(funname);
	gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			      NULL, prn);
	state.line = fnline;
	err = gui_exec_line(&state, dataset, NULL);

	if (err) {
	    gui_errmsg(err);
	} else {
	    ret = 1;
	}
    }

    return ret;
}

/* get a listing of available "official" addons along with (a) the
   name of the versioned subdirectory containing the most recent
   usable version (given the gretl version) and (b) the date of that
   package version, taken from its spec file. E.g.

   gig 1.9.5 2011-04-22
   ivpanel 1.9.4 2011-02-10
*/

int query_addons (void)
{
    gchar *query;
    char *buf = NULL;
    int err = 0;

    query = g_strdup_printf("/addons-data/pkginfo.php?gretl_version=%s",
			    GRETL_VERSION);
    err = query_sourceforge(query, &buf);
    g_free(query);

    if (!err && buf == NULL) {
	/* shouldn't happen */
	err = E_DATA;
    }

    if (!err) {
	infobox(buf);
    }

    free(buf);

    return err;
}

#define PKG_DEBUG 0

static int query_addons_dir (const char *pkgname, char *sfdir)
{
    gchar *query;
    char *buf = NULL;
    int err = 0;

    *sfdir = '\0';

    query = g_strdup_printf("/addons-data/pkgdir.php?gretl_version=%s"
			    "&pkg=%s", GRETL_VERSION, pkgname);
    err = query_sourceforge(query, &buf);
    g_free(query);

    if (!err && buf == NULL) {
	/* shouldn't happen */
	err = E_DATA;
    }

    if (!err && buf != NULL && strstr(buf, "<head>")) {
	/* got some sort of garbage */
	err = E_DATA;
    }
    
    if (!err) {
	char *p = strchr(buf, ':');

	if (p == NULL || 
	    sscanf(p + 2, "%15s", sfdir) != 1 ||
	    strcmp(sfdir, "none") == 0) {
	    err = E_DATA;
	}
    }

    free(buf);

    if (err) {
	errbox_printf("Couldn't find %s for gretl %s",
		      pkgname, GRETL_VERSION);
    } 

    return err;
}

int download_addon (const char *pkgname, char **local_path)
{
    const char *SF = "http://downloads.sourceforge.net/"
	"project/gretl/addons";
    char pkgdir[16];
    int err;

    err = query_addons_dir(pkgname, pkgdir);

    if (!err) {
	const char *path = gretl_function_package_path();
	gchar *uri = g_strdup_printf("%s/%s/%s.zip", SF, pkgdir, pkgname);
	gchar *fullname = g_strdup_printf("%s%s.zip", path, pkgname);

#if 1
	fprintf(stderr, "uri   = '%s'\n", uri);
	fprintf(stderr, "fname = '%s'\n", fullname);
#endif

	err = retrieve_public_file(uri, fullname);
	fprintf(stderr, "retrieve_public_file: err = %d\n", err);
	if (!err) {
	    err = gretl_unzip_into(fullname, path);
	    fprintf(stderr, "gretl_unzip_into: err = %d\n", err);
	    gretl_remove(fullname);
	}
	if (err) {
	    gui_errmsg(err);
	} else if (local_path != NULL) {
	    /* return local path to gfn file */
	    *local_path =
		gretl_strdup_printf("%s%s%s%s.gfn", path,
				    pkgname, SLASHSTR,
				    pkgname);
	    fprintf(stderr, "local_path: '%s'\n", *local_path);
	}
	g_free(uri);
	g_free(fullname);
    }

    return err;
}

#if 0 

/* for a package that's not a (zipfile) addon
   work in progress, 2015-07-21
*/

static int download_package (const char *pkgname, char **local_path)
{
    int err;

    err = script_install_function_package(pkgname, OPT_NONE,
					  NULL, NULL,
					  local_path);
    fprintf(stderr, "d/l package: err=%d, local_path '%s'\n",
	    err, *local_path);
    return err;
}

#endif

/* information about a function package that offers a
   menu attachment point */

typedef enum {
    GPI_MODELWIN = 1 << 0,
    GPI_SUBDIR   = 1 << 1,
    GPI_SYSFILE  = 1 << 2,
    GPI_INCLUDED = 1 << 3
} GpiFlags;

enum {
    SYS_PACKAGES,
    USER_PACKAGES
};

#define gpi_included(g) (g->flags & GPI_INCLUDED)
#define gpi_modelwin(g) (g->flags & GPI_MODELWIN)
#define gpi_sysfile(g) (g->flags & GPI_SYSFILE)
#define gpi_subdir(g) (g->flags & GPI_SUBDIR)
#define gpi_ptype(g) ((g->flags & GPI_SUBDIR)? PKG_SUBDIR : PKG_TOPLEV)

struct gui_package_info_ {
    char *pkgname;  /* @name element from packages.xml */
    char *label;    /* @label element from packages.xml */
    char *menupath; /* @path element from packages.xml */
    char *filepath; /* actual filesystem path */
    GpiFlags flags; /* state flags */
    guint merge_id; /* created at runtime when added to GUI */
    GtkActionGroup *ag; /* run-time UI thing */
};

typedef struct gui_package_info_ gui_package_info;

static void add_package_to_menu (gui_package_info *gpi,
				 windata_t *vwin);

static gui_package_info *gpkgs;
static int n_gpkgs;
static int gpkgs_changed;

/* Return the total number of slots for function packages
   "registered" for use via menus. Note that this number
   may include some slots that are actually empty, if the
   user has removed a dynamic menu item during the
   current gretl session.
*/

int n_registered_packages (void)
{
    return n_gpkgs;
}

int n_user_handled_packages (void)
{
    int i, n = 0;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkgs[i].pkgname == NULL) {
	    /* a vacant slot */
	    continue;
	} else if (gpkgs[i].flags & GPI_SYSFILE) {
	    /* a file under control of the "system" 
	       packages.xml file */
	    continue;
	} else {
	    n++;
	}
    }    

    return n;
}

/* On adding a new entry to the gui package info
   array, zero it appropriately */

static void gpi_entry_init (gui_package_info *gpi)
{
    gpi->pkgname = NULL;
    gpi->label = NULL;
    gpi->menupath = NULL;
    gpi->filepath = NULL;
    gpi->flags = 0;
    gpi->merge_id = 0;
    gpi->ag = NULL;
}

static gchar *packages_xml_path (int which)
{
    if (which == SYS_PACKAGES) {
	return g_strdup_printf("%sfunctions%cpackages.xml", 
			       gretl_home(), SLASH);
    } else {
	return g_strdup_printf("%sfunctions%cpackages.xml", 
			       gretl_dotdir(), SLASH);
    }
}

static void write_packages_xml (void)
{
    char *fname = packages_xml_path(USER_PACKAGES);
    int i, n_write = n_user_handled_packages();
	
    if (n_write == 0) {
	gretl_remove(fname);
    } else {
	FILE *fp = gretl_fopen(fname, "w");

	if (fp == NULL) {
	    fprintf(stderr, "Couldn't write to %s\n", fname);
	    return;
	}

	fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", fp);
	fputs("<gretl-package-info>\n", fp);

	for (i=0; i<n_gpkgs; i++) {
	    if (gpkgs[i].pkgname != NULL &&
		!(gpkgs[i].flags & GPI_SYSFILE)) {
		fprintf(fp, "<package name=\"%s\"", gpkgs[i].pkgname);
		fprintf(fp, " label=\"%s\"", gpkgs[i].label);
		if (gpkgs[i].flags & GPI_MODELWIN) {
		    fputs(" model-window=\"true\"", fp);
		}
		fprintf(fp, " path=\"%s\"", gpkgs[i].menupath);
		if (gpkgs[i].flags & GPI_SUBDIR) {
		    fputs("/>\n", fp);
		} else {
		    fputs(" toplev=\"true\"/>\n", fp);
		}
	    }
	}

	fputs("</gretl-package-info>\n", fp);
	fclose(fp);
    }

    g_free(fname);
}

/* At exit, clean up the entire array of gui package
   info, after first saving the info to file if anything
   relevant has changed.
*/

void destroy_gui_package_info (void)
{
    if (gpkgs_changed) {
	/* save info to packages.xml first */
	write_packages_xml();
    }
    
    if (gpkgs != NULL && n_gpkgs > 0) {
	int i;

	for (i=0; i<n_gpkgs; i++) {
	    free(gpkgs[i].pkgname);
	    free(gpkgs[i].label);
	    free(gpkgs[i].menupath);
	    free(gpkgs[i].filepath);
	}
	free(gpkgs);
    }

    gpkgs = NULL;
    n_gpkgs = 0;
    gpkgs_changed = 0;
}

static int gpkg_name_match (char *s1, const char *s2)
{
    return s1 != NULL && strcmp(s1, s2) == 0;
}

static gui_package_info *get_gpi_entry (const gchar *pkgname)
{
    int i;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkg_name_match(gpkgs[i].pkgname, pkgname)) {
	    return &gpkgs[i];
	}
    }

    return NULL;
}

void get_registered_pkg_info (int i, char **name, char **path,
			      char **label, int *modelwin)
{
    if (i < 0 || i >= n_gpkgs || 
	gpkgs[i].pkgname == NULL ||
	(gpkgs[i].flags & GPI_SYSFILE)) {
	*name = *path = *label = NULL;
    } else {
	*name = gpkgs[i].pkgname;
	*path = gpkgs[i].menupath;
	*label = gpkgs[i].label;
	*modelwin = (gpkgs[i].flags & GPI_MODELWIN)? 1 : 0;
    }
}

static void maybe_record_include (const char *pkgname,
				  int model_id)
{
    gui_package_info *gpi;
    int i;

    for (i=0; i<n_gpkgs; i++) {
	gpi = &gpkgs[i];
	if (gpkg_name_match(gpi->pkgname, pkgname)) {
	    if (!gpi_included(gpi)) {
		lib_command_sprintf("include %s.gfn", pkgname);
		if (model_id > 0) {
		    record_model_command_verbatim(model_id);
		} else {
		    record_command_verbatim();
		}
		gpi->flags |= GPI_INCLUDED;
	    }
	    break;
	}
    }
}

int package_is_available_for_menu (const gchar *pkgname,
				   const char *fname)
{
    int present = 0;
    int i, ret = 0;

    for (i=0; i<n_gpkgs && !present; i++) {
	if (gpkg_name_match(gpkgs[i].pkgname, pkgname)) {
	    present = 1;
	}
    }

    if (!present) {
	/* not already present in menus: can it be added? */
	ret = package_has_menu_attachment(fname, NULL, NULL, NULL);
    }

    return ret;
}

static int official_addon (const char *pkgname)
{
    /* All of the following are available in zip format
       from sourceforge */
    if (!strcmp(pkgname, "SVAR") ||
	!strcmp(pkgname, "gig") ||
	!strcmp(pkgname, "HIP") ||
	!strcmp(pkgname, "ivpanel")) {
	return 1;
    } else {
	return 0;
    }
}

/* Callback for a menu item representing a function package whose
   name is attached to @action. We first see if we can find the full
   path to the corresponding gfn file; if so we initiate a GUI call to
   the package.
*/

static void gfn_menu_callback (GtkAction *action, windata_t *vwin)
{
    const gchar *pkgname = gtk_action_get_name(action);
    gui_package_info *gpi;

    gpi = get_gpi_entry(pkgname);
    if (gpi == NULL) {
	/* "can't happen" */
	return;
    }

    if (gpi->filepath == NULL) {
	gpi->filepath =
	    gretl_function_package_get_path(pkgname, gpi_ptype(gpi));
    }

    if (gpi->filepath == NULL && official_addon(pkgname)) {
	gchar *msg = g_strdup_printf(_("The %s package was not found, or is not "
				       "up to date.\nWould you like to try "
				       "downloading it now?"), pkgname);
	int resp = yes_no_dialog(NULL, msg, vwin_toplevel(vwin));

	g_free(msg);
	if (resp == GRETL_YES) {
	    int err;
	    
	    if (official_addon(pkgname)) {
		/* FIXME maybe record download in command log? */
		err = download_addon(pkgname, &gpi->filepath);
	    } else {
		/* ?? err = download_package(pkgname, &gpi->filepath); */
		err = 1;
		errbox_printf("Sorry, could not find %s", pkgname);
	    }
	    if (err) {
		return;
	    }
	} else {
	    /* canceled, effectively */
	    return;
	}
    }

    if (gpi->filepath != NULL) {
	call_info *cinfo;
	int err = 0;

	cinfo = start_cinfo_for_package(pkgname, gpi->filepath, vwin, &err);
	if (cinfo != NULL) {
	    call_function_package(cinfo, vwin, 0);
	}
    } else {
	errbox_printf("Sorry, could not find %s", pkgname);
	/* ? */
	gui_function_pkg_unregister(pkgname);
    }
}

static int package_is_unseen (const char *name, int n)
{
    int i;

    for (i=0; i<n; i++) {
	if (gpkg_name_match(gpkgs[i].pkgname, name)) {
	    return 0;
	}
    }

    return 1;
}

static void gpi_set_flags (gui_package_info *gpi,
			   int which,
			   int subdir,
			   int modelwin)
{
    gpi->flags = 0;
    
    if (which == SYS_PACKAGES) {
	gpi->flags |= GPI_SYSFILE;
    }
    if (subdir) {
	gpi->flags |= GPI_SUBDIR;
    }
    if (modelwin) {
	gpi->flags |= GPI_MODELWIN;
    }
}

static int read_packages_file (const char *fname, int *pn, int which)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    int err, n = *pn;

    err = gretl_xml_open_doc_root(fname, "gretl-package-info", &doc, &cur);

    if (err) {
	return (which == SYS_PACKAGES)? err : 0;
    } 

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "package")) {
	    int mw = gretl_xml_get_prop_as_bool(cur, "model-window");
	    int top = gretl_xml_get_prop_as_bool(cur, "toplev");
	    xmlChar *name, *desc, *path;
	    int freeit = 1;

	    name = xmlGetProp(cur, (XUC) "name");
	    desc = xmlGetProp(cur, (XUC) "label");
	    path = xmlGetProp(cur, (XUC) "path");

	    if (name == NULL || desc == NULL || path == NULL) {
		err = E_DATA;
	    } else if (package_is_unseen((const char *) name, n)) {
		gpkgs = myrealloc(gpkgs, (n+1) * sizeof *gpkgs);
		if (gpkgs == NULL) {
		    err = E_ALLOC;
		} else {
		    freeit = 0;
		    gpkgs[n].pkgname = (char *) name;
		    gpkgs[n].label = (char *) desc;
		    gpkgs[n].menupath = (char *) path;
		    gpkgs[n].filepath = NULL;
		    gpi_set_flags(&gpkgs[n], which, !top, mw);
		    gpkgs[n].merge_id = 0;
		    gpkgs[n].ag = NULL;
		    n++;
		} 
	    }
	    if (freeit) {
		free(name);
		free(desc);
		free(path);
	    }		
	}
	if (!err) {
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }    

    *pn = n;

    return err;
}

static void destroy_gpi_ui (gui_package_info *gpi)
{
    fprintf(stderr, "removing UI for %s\n", gpi->pkgname);
    gtk_ui_manager_remove_ui(mdata->ui, gpi->merge_id);
    gtk_ui_manager_remove_action_group(mdata->ui, gpi->ag);
    g_object_unref(gpi->ag);
    
    gpi->merge_id = 0;
    gpi->ag = NULL;
}

/* On the call to remove a dynamic menu item, tear down
   its UI and blank out @gpi so that it may be safely
   reused by another package.
*/

static void clear_gpi_entry (gui_package_info *gpi)
{
    if (gpi->merge_id > 0) {
	destroy_gpi_ui(gpi);
    }    
    
    free(gpi->pkgname);
    free(gpi->label);
    free(gpi->menupath);
    free(gpi->filepath);

    gpi->pkgname = NULL;
    gpi->label = NULL;
    gpi->menupath = NULL;
    gpi->filepath = NULL;

    gpi->flags = 0;

    gpkgs_changed = 1;
}

static int fill_gpi_entry (gui_package_info *gpi,
			   const char *pkgname,
			   const char *fname,
			   const char *label,
			   const char *relpath,
			   int modelwin,
			   int uses_subdir,
			   int replace)
{
    int err = 0;

    if (replace) {
	free(gpi->pkgname);
	free(gpi->label);
	free(gpi->menupath);
	free(gpi->filepath);
    }
	
    gpi->pkgname = gretl_strdup(pkgname);
    gpi->label = gretl_strdup(label);
    gpi->menupath = malloc(9 + strlen(relpath));
    if (gpi->menupath != NULL) {
	sprintf(gpi->menupath, "/menubar%s", relpath);
    }
    gpi->filepath = gretl_strdup(fname);
    gpi_set_flags(gpi, USER_PACKAGES, uses_subdir, modelwin);

    if (gpi->pkgname == NULL || gpi->label == NULL ||
	gpi->menupath == NULL || gpi->filepath == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static gui_package_info *get_new_package_info (int *err)
{
    int i, pos = -1;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkgs[i].pkgname == NULL) {
	    /* found an unused slot */
	    pos = i;
	    break;
	}
    }

    if (pos < 0) {
	/* we need to extend the array */
	gui_package_info *gpkgs_new;
	int n = n_gpkgs + 1;
	
	gpkgs_new = realloc(gpkgs, n * sizeof *gpkgs);
	if (gpkgs_new == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	} else {
	    gpkgs = gpkgs_new;
	    n_gpkgs = n;
	    pos = n - 1;
	    gpi_entry_init(&gpkgs[pos]);
	}
    }

    return &gpkgs[pos];
}

static int menu_paths_differ (const char *rp, const char *mp)
{
    if ((rp == NULL && mp != NULL) ||
	(rp != NULL && mp == NULL)) {
	return 1;
    } else if (rp == NULL && mp == NULL) {
	return 0;
    } else {
	if (!strncmp(mp, "/menubar", 8)) {
	    mp += 8;
	}

	return strcmp(rp, mp);
    }
}

static int gpi_strings_differ (const char *s1, const char *s2)
{
    if ((s1 == NULL && s2 != NULL) ||
	(s1 != NULL && s2 == NULL)) {
	return 1;
    } else if (s1 == NULL && s2 == NULL) {
	return 0;
    } else {
	return strcmp(s1, s2);
    }
}

/* This is called when creating a new entry in the
   gui package registry, and also when updating
   the entry for a previously registered package.
*/

static int update_gui_package_info (const char *pkgname,
				    const char *fname,
				    const char *label,
				    const char *relpath,
				    int modelwin,
				    int uses_subdir)
{
    gui_package_info *gpi;
    int menu_update = 0;
    int update = 0;
    int replace = 1;
    int err = 0;

    gpi = get_gpi_entry(pkgname);

    if (gpi != NULL) {
	/* found a pre-existing entry for @pkgname:
	   let's see if there are really any changes
	*/
	if (gpi_strings_differ(label, gpi->label)) {
	    menu_update = update = 1;
	} else if (menu_paths_differ(relpath, gpi->menupath)) {
	    menu_update = update = 1;
	} else if (gpi_strings_differ(fname, gpi->filepath)) {
	    update = 1;
	} else if (modelwin && !gpi_modelwin(gpi)) {
	    menu_update = update = 1;
	} else if (uses_subdir && !gpi_subdir(gpi)) {
	    update = 1;
	}
    } else {
	gpi = get_new_package_info(&err);
	menu_update = update = 1;
	replace = 0;
    }

    fprintf(stderr, "update_gui_package_info: update=%d, menu_update=%d, "
	    "err=%d\n", update, menu_update, err);

    if (!err && menu_update && gpi->merge_id > 0) {
	/* trash stale UI info */
	destroy_gpi_ui(gpi);
    }	

    if (!err && update) {
	err = fill_gpi_entry(gpi,
			     pkgname,
			     fname,
			     label,
			     relpath,
			     modelwin,
			     uses_subdir,
			     replace);
    }

    if (update && !err) {
	gpkgs_changed = 1;
	if (menu_update) {
	    if (!modelwin) {
		/* add to main-window menu */
		add_package_to_menu(gpi, mdata);
	    }
	}
	/* sync GUI */
	maybe_update_pkg_registry_window(pkgname, MENU_ADD_FN_PKG);
    }
    
    return err;
}

/* read "packages.xml" to find out what's what among packages
   that offer to place themselves in the gretl menu system
*/

static int gui_package_info_init (void)
{
    gchar *fname;
    int err, n = 0;

    /* start with the "system" packages.xml, which should
       always be present */
    fname = packages_xml_path(SYS_PACKAGES);
    err = read_packages_file(fname, &n, SYS_PACKAGES);
    g_free(fname);

    if (!err) {
	/* then read the per-user packages.xml, if present */
	fname = packages_xml_path(USER_PACKAGES);
	if (gretl_file_exists(fname)) {
	    err = read_packages_file(fname, &n, USER_PACKAGES);
	}
	g_free(fname);
    }

    if (err) {
	destroy_gui_package_info();
    } else {
	n_gpkgs = n;
    }

    return err;
}

/* For function packages offering a menu attachment point: given the
   internal package name (e.g. "gig") and a menu path where we'd like
   it to appear (e.g. "/menubar/Model/TSModels" -- see the ui
   definition file gui2/gretlmain.xml), construct the appropriate menu
   item and connect it to gfn_menu_callback(), for which see above.

   This is called for both main-window and model-window menu
   attachments.
*/

static void add_package_to_menu (gui_package_info *gpi,
				 windata_t *vwin)
{
    static GtkActionEntry item = {
	NULL, NULL, NULL, NULL, NULL, G_CALLBACK(gfn_menu_callback)
    };
    char *fixed_label = NULL;
    guint merge_id;

#if PKG_DEBUG
    fprintf(stderr, "add_package_to_menu:\n pkgname='%s', menupath='%s', label='%s'\n",
	    gpi->pkgname, gpi->menupath, gpi->label);
#endif

    item.name = gpi->pkgname;
    item.label = gpi->label != NULL ? gpi->label : gpi->pkgname;
    
    if (strchr(item.label, '_')) {
	const char *s = item.label;
	int n = 0;

	while (*s && (s = strchr(s, '_')) != NULL) {
	    n++;
	    s++;
	}
	fixed_label = malloc(strlen(item.label) + n + 1);
	double_underscores(fixed_label, item.label);
	item.label = fixed_label;
    }

    merge_id = gtk_ui_manager_new_merge_id(vwin->ui);

    gtk_ui_manager_add_ui(vwin->ui, merge_id, gpi->menupath, 
			  _(item.label), item.name,
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);

    if (vwin == mdata) {
	gpi->merge_id = merge_id;
    }

    gpi->ag = gtk_action_group_new(item.name);
    gtk_action_group_set_translation_domain(gpi->ag, "gretl");
    gtk_action_group_add_actions(gpi->ag, &item, 1, vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, gpi->ag, 0);
    // g_object_unref(gpi->ag);

    free(fixed_label);
}

/* run a package's gui-precheck function to determine if
   it's OK to add its GUI interface to a gretl
   model-window menu
*/

static int precheck_error (ufunc *func, windata_t *vwin)
{
    PRN *prn;
    double check_err = 0;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
    set_genr_model_from_vwin(vwin);
    err = gretl_function_exec(func, GRETL_TYPE_DOUBLE,
			      dataset, &check_err, NULL, prn);
    unset_genr_model();
    gretl_print_destroy(prn);

    if (err == 0 && check_err != 0) {
	err = 1;
    }
    
    return err;
}

static int maybe_add_model_pkg (gui_package_info *gpi,
				windata_t *vwin)
{
    int dreq, mreq, minver = 0;
    gchar *precheck = NULL;
    fnpkg *pkg;
    int err = 0;

    if (gpi->filepath == NULL) {
	gpi->filepath =
	    gretl_function_package_get_path(gpi->pkgname, gpi_ptype(gpi));
    }
	
    if (gpi->filepath == NULL) {
	fprintf(stderr, "%s: couldn't find it\n", gpi->pkgname);
	return E_FOPEN;
    }

    pkg = get_function_package_by_filename(gpi->filepath, &err);

    if (!err) {
	err = function_package_get_properties(pkg,
					      "data-requirement", &dreq,
					      "model-requirement", &mreq,
					      "min-version", &minver,
					      "gui-precheck", &precheck,
					      NULL);
    }

    if (!err) {
	/* "skip" = skip this package 'cos it won't work
	   with the current model */
	int skip = 0;
	
	if (mreq > 0) {
	    MODEL *pmod = vwin->data;

	    skip = pmod->ci != mreq;
	}
	if (!skip) {
	    skip = check_function_needs(dataset, dreq, minver);
	}
	if (!skip && precheck != NULL) {
	    ufunc *func = get_function_from_package(precheck, pkg);
	    
	    if (func == NULL || precheck_error(func, vwin)) {
		skip = 1;
	    }
	} 
	if (!skip) {
	    add_package_to_menu(gpi, vwin);
	}
    }

    g_free(precheck);

    return err;
}

/* Called from gretl.c on initializing the GUI: put suitable
   function packages (other than those that are designed to
   appear in model-window menus) into the appropriate main-
   window menus.
*/

void maybe_add_packages_to_menus (windata_t *vwin)
{
    gui_package_info *gpi;
    int i;

    if (gpkgs == NULL) {
	gui_package_info_init();
    }

    for (i=0; i<n_gpkgs; i++) {
	gpi = &gpkgs[i];
	if (gpi->pkgname != NULL && !gpi_modelwin(gpi)) {
	    add_package_to_menu(gpi, vwin);
	}
    }
}

/* Called from gui_utils.c when setting up the UI for a
   model window.
*/

void maybe_add_packages_to_model_menus (windata_t *vwin)
{
    gui_package_info *gpi;
    int i, err;

    if (gpkgs == NULL) {
	gui_package_info_init();
    }

    for (i=0; i<n_gpkgs; i++) {
	gpi = &gpkgs[i];
	if (gpi->pkgname != NULL && gpi_modelwin(gpi)) {
	    err = maybe_add_model_pkg(gpi, vwin);
	    if (err) {
		if (err == E_FOPEN && gpi_sysfile(gpi)) {
		    ;
		} else {
		    /* delete entry from registry */
		    gui_function_pkg_unregister(gpi->pkgname);
		}
	    }
	}
    }
}

/* Below: apparatus for activation when the user installs a function
   package (other than an official addon) from the gretl server.

   We check to see if (a) the package offers a menu attachment, and if
   so (b) that the package is not already "registered" in the user's
   packages.xml file. If both of these conditions are met we put up a
   dialog asking if the user wants to add the package to the menu
   system. If the answer is Yes we write an appropriate entry into
   packages.xml, or write this file from scratch if it doesn't yet
   exist.
*/

/* Find out where the package is supposed to attach: the
   @mpath string (which gets into the gfn file from its
   associated spec file or the GUI) should look something like

   MODELWIN/Analysis or
   MAINWIN/Model

   The first portion just tells us in which window it should
   appear.
*/

static const gchar *pkg_get_attachment (const gchar *mpath,
					int *modelwin)
{
    const gchar *relpath = mpath;

#if PKG_DEBUG
    fprintf(stderr, "pkg_get_attachment: mpath = '%s'\n", mpath);
#endif    

    if (!strncmp(mpath, "MAINWIN/", 8)) {
	relpath = mpath + 7;
    } else if (!strncmp(mpath, "MODELWIN/", 9)) {
	relpath = mpath + 8;
	*modelwin = 1;
    }

    return relpath;
}

static int pkg_attach_query (const gchar *name,
			     const gchar *label,
			     const gchar *relpath,
			     int modelwin,
			     GtkWidget *parent)
{
    int resp = -1;

    if (relpath != NULL && *relpath != '\0') {
	const gchar *window_names[] = {
	    N_("main window"),
	    N_("model window")
	};
	gchar *msg, *ustr = NULL;

	ustr = user_friendly_menu_path(relpath, modelwin);
	msg = g_strdup_printf(_("The package %s can be attached to the "
				"gretl menus\n"
				"as \"%s/%s\" in the %s.\n"
				"Do you want to do this?"),
			      name, ustr ? ustr : relpath, _(label),
			      modelwin ? _(window_names[1]) :
			      _(window_names[0]));
	resp = yes_no_dialog(NULL, msg, parent);
	g_free(msg);
	g_free(ustr);
    }

    return resp;
}

/* Called from fnsave.c, when edits to a function package
   are being saved.
*/

int gui_function_pkg_revise_status (const gchar *pkgname,
				    const gchar *fname,
				    const gchar *label,
				    const gchar *mpath,
				    gboolean uses_subdir)
{
    gui_package_info *gpi;
    int has_attachment = 0;
    int do_update = 0;
    int err = 0;

    /* In this context we may, in principle, be adding,
       removing, or revising the gui-menu status of
       a package.
    */

    if (label != NULL && mpath != NULL) {
	has_attachment = 1;
    }

    gpi = get_gpi_entry(pkgname);

    if (gpi == NULL) {
	/* not in registry yet */
	if (has_attachment) {
	    do_update = 1;
	}
    } else {
	/* already in registry */
	if (has_attachment) {
	    do_update = 1;
	} else {
	    gui_function_pkg_unregister(pkgname);
	}
    }

    if (do_update) {
	const char *relpath;
	int modelwin = 0;
	    
	relpath = pkg_get_attachment(mpath, &modelwin);
	err = update_gui_package_info(pkgname,
				      fname,
				      label,
				      relpath,
				      modelwin,
				      uses_subdir);
	if (err) {
	    gui_errmsg(err);
	}
    }

    return err;
}

/* Remove a package from the in-memory representation of 
   menu-attached packages.
*/

void gui_function_pkg_unregister (const gchar *pkgname)
{
    gui_package_info *gpi;

    gpi = get_gpi_entry(pkgname);
    if (gpi != NULL) {
	clear_gpi_entry(gpi);
    }     

    /* sync the package registry window, if it happens to
       be open currently */
    maybe_update_pkg_registry_window(pkgname, DELETE_FN_PKG);
}

static int in_own_subdir (const char *pkgname, const char *path)
{
    gchar *test;
    int ret = 0;

    /* e.g. "mypkg/mypkg" */
    test = g_strdup_printf("%s%c%s", pkgname, SLASH, pkgname);
    
    if (strstr(path, test) != NULL) {
	ret = 1;
    } else {
	gretl_errmsg_sprintf("The function file %s.gfn is not installed correctly:\n"
			     "it should be in a subdirectory named '%s'.",
			     pkgname, pkgname);
    }

    g_free(test);

    return ret;
}

/* Actually do the business of registering a function package
   to appear in a menu */

static int gui_function_pkg_register (const char *fname,
				      const char *pkgname,
				      const char *label,
				      const char *relpath,
				      int modelwin)
{
    fnpkg *pkg;
    int err = 0;

    pkg = get_function_package_by_filename(fname, &err);

#if PKG_DEBUG
    fprintf(stderr, "add_gfn_to_registry: %s: err = %d\n", fname, err);
#endif       

    if (!err) {
	int uses_subdir = 0;

	err = function_package_get_properties(pkg, "lives-in-subdir",
					      &uses_subdir, NULL);

	if (!err && !uses_subdir) {
	    /* fallback detection: packages that have PDF doc
	       must be in their own subdir */
	    uses_subdir = function_package_has_PDF_doc(pkg, NULL);
	}

	if (!err && uses_subdir && !in_own_subdir(pkgname, fname)) {
	    /* detect mis-installed package: should be in
	       own subdirectory but is not */
	    err = E_DATA;
	}

	if (!err) {
	    err = update_gui_package_info(pkgname,
					  fname,
					  label,
					  relpath,
					  modelwin,
					  uses_subdir);
	}
    }

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

/* The following is called in two contexts:

   (1) From the handler for installing a function package from the
   gretl server.

   (2) From the popup menu-item or button "Add to menu" in the window
   displaying installed function packages.

   We return non-zero if we show a dialog here: that's for the
   benefit of the installation handler, to tell it not to put
   up a second, redundant confirmation dialog.
*/

int gui_function_pkg_query_register (const char *fname,
				     GtkWidget *parent)
{
    char *pkgname = NULL;
    char *menupath = NULL;
    char *label = NULL;
    int notified = 0;

    if (package_has_menu_attachment(fname, &pkgname, &menupath,
				    &label)) {
	if (menupath == NULL || label == NULL) {
	    msgbox(_("This package lacks a label or menu-path"),
		   GTK_MESSAGE_ERROR, parent);
	} else {
	    const gchar *relpath;
	    int resp, modelwin = 0;
	    
	    relpath = pkg_get_attachment(menupath, &modelwin);
	    resp = pkg_attach_query(pkgname, label, relpath,
				    modelwin, parent);
	    if (resp == GRETL_YES) {
		gui_function_pkg_register(fname, pkgname,
					  label, relpath,
					  modelwin);
	    }
	    notified = 1;
	}
	
	free(pkgname);
	free(menupath);
	free(label);
    }

    return notified;
}

char *installed_addon_status_string (const char *path,
				     const char *svstr)
{
    fnpkg *pkg;
    int err = 0;
    char *ret = NULL;

    pkg = get_function_package_by_filename(path, &err);

    if (pkg != NULL) {
	int current = 0;
	int update_ok = 0;
	char reqstr[8] = {0};
	gchar *ivstr = NULL;
	int minver = 0;

	/* @ivstr = installed package version string
	   @svstr = package version string from server
	*/
	err = function_package_get_properties(pkg, 
					      "version", &ivstr, 
					      "min-version", &minver,
					      NULL);
	if (!err) {
	    double svnum = dot_atof(svstr);
	    double ivnum = dot_atof(ivstr);

	    current = ivnum >= svnum;
	    g_free(ivstr);

	    if (!current) {
		/* Not current, but can the addon be updated?  It may
		   be that the running instance of gretl is too old.
		*/
		update_ok = package_version_ok(minver, reqstr);
	    }

	    if (current) {
		ret = gretl_strdup(_("Up to date"));
	    } else if (update_ok) {
		ret = gretl_strdup(_("Not up to date"));
	    } else if (*reqstr != '\0') {
		ret = gretl_strdup_printf(_("Requires gretl %s"), reqstr);
	    }
	}
    }

    if (ret == NULL) {
	ret = gretl_strdup(_("Error reading package"));
    }

    return ret;
}

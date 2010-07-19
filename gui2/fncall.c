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
#include "dlgutils.h"
#include "selector.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "usermat.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "gretl_string_table.h"
#include "gretl_scalar.h"
#include "gretl_bundle.h"
#include "database.h"
#include "guiprint.h"
#include "ssheet.h"
#include "datafiles.h"

#define FCDEBUG 0

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;      /* main dialog */
    GtkWidget *top_hbox; /* upper hbox in dialog */
    GList *lsels;        /* list argument selectors */
    GList *msels;        /* matrix arg selectors */
    GList *ssels;        /* string arg selectors */
    int *publist;        /* list of public interfaces */
    int iface;           /* selected interface */
    int extracol;
    const ufunc *func;
    FuncDataReq dreq;
    int n_params;
    char rettype;
    gchar **args;
    gchar *ret;
};

#define scalar_arg(t) (t == GRETL_TYPE_DOUBLE || t == GRETL_TYPE_SCALAR_REF)
#define series_arg(t) (t == GRETL_TYPE_SERIES || t == GRETL_TYPE_SERIES_REF)
#define matrix_arg(t) (t == GRETL_TYPE_MATRIX || t == GRETL_TYPE_MATRIX_REF)

static GtkWidget *open_fncall_dlg;
static gboolean close_on_OK = TRUE;

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo);

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

static call_info *cinfo_new (void)
{
    call_info *cinfo = mymalloc(sizeof *cinfo);

    if (cinfo == NULL) {
	return NULL;
    }

    cinfo->publist = NULL;
    cinfo->iface = -1;

    cinfo->lsels = NULL;
    cinfo->msels = NULL;

    cinfo->func = NULL;
    cinfo->n_params = 0;

    cinfo->rettype = GRETL_TYPE_NONE;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    cinfo->extracol = 0;
    cinfo->dreq = 0;

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
    if (cinfo->lsels != NULL) {
	g_list_free(cinfo->lsels);
    }
    if (cinfo->msels != NULL) {
	g_list_free(cinfo->msels);
    }
    
    free(cinfo->publist);
    free(cinfo);
}

static const char *arg_type_string (int t)
{
    if (t == GRETL_TYPE_BOOL)   return "boolean";
    if (t == GRETL_TYPE_INT)    return "int";
    if (t == GRETL_TYPE_LIST)   return "list";
    if (t == GRETL_TYPE_DOUBLE) return "scalar";
    if (t == GRETL_TYPE_SERIES) return "series";
    if (t == GRETL_TYPE_MATRIX) return "matrix";
    if (t == GRETL_TYPE_STRING) return "string";
    if (t == GRETL_TYPE_BUNDLE) return "bundle";
    
    if (t == GRETL_TYPE_SCALAR_REF) return "scalar *";
    if (t == GRETL_TYPE_SERIES_REF) return "series *";
    if (t == GRETL_TYPE_MATRIX_REF) return "matrix *";
    if (t == GRETL_TYPE_BUNDLE_REF) return "bundle *";

    return "";
}

static int check_args (call_info *cinfo)
{
    int i;

    /* FIXME optional args? */

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    if (cinfo->args[i] == NULL) {
		errbox(_("Argument %d (%s) is missing"), i + 1,
		       fn_param_name(cinfo->func, i));
		return 1;
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

static GtkWidget *label_hbox (GtkWidget *w, const char *txt, 
			      int vspace, int center)
{
    GtkWidget *hbox, *label;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, vspace);

    label = gtk_label_new(txt);

    if (center) {
	gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    } else {
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    }

    gtk_widget_show(label);

    return hbox;
}

static gboolean update_int_arg (GtkWidget *w, call_info *cinfo)
{
    int val = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "argnum"));

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);

    return FALSE;
}

static gboolean update_bool_arg (GtkWidget *w, call_info *cinfo)
{
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "argnum"));

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
    gchar *s = gtk_combo_box_get_active_text(combo);
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
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(combo), "argnum"));
    char *s;

    g_free(cinfo->args[i]);
    s = cinfo->args[i] = combo_box_get_trimmed_text(combo);

    if (s != NULL && fn_param_type(cinfo->func, i) == GRETL_TYPE_DOUBLE) {
	if (isdigit(*s) || *s == '-' || *s == '+' || *s == ',') {
	    charsub(s, ',', '.');
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

static GList *get_selection_list (call_info *cinfo, int i, int type,
				  int set_default)
{
    GList *list = NULL;
    const char *name;
    int n, optional = 0;

    if (i >= 0) {
	optional = fn_param_optional(cinfo->func, i);
    }

    if (!set_default) {
	list = g_list_append(list, "");
    }

    if (scalar_arg(type)) {
	n = n_saved_scalars();
	for (i=0; i<n; i++) {
	    name = gretl_scalar_get_name(i);
	    list = g_list_append(list, (gpointer) name);
	}	
    } else if (series_arg(type)) {
	for (i=1; i<datainfo->v; i++) {
	    if (!var_is_hidden(datainfo, i)) {
		list = g_list_append(list, (gpointer) datainfo->varname[i]);
	    } 
	}
	list = g_list_append(list, (gpointer) datainfo->varname[0]);
    } else if (type == GRETL_TYPE_LIST) {
	n = n_saved_lists();

	if (optional) {
	    list = g_list_append(list, "null");
	}

	for (i=0; i<n; i++) {
	    name = get_list_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}
    } else if (matrix_arg(type)) {
	n = n_user_matrices();
	for (i=0; i<n; i++) {
	    name = get_matrix_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}	
    } 

    if (optional && type != GRETL_TYPE_LIST) {
	list = g_list_append(list, "null");
    }

    return list;
}

static windata_t *make_help_viewer (const char *fnname, char *buf)
{
    windata_t *vwin;
    gchar *title;

    title = g_strdup_printf(_("help on %s"), fnname);
    vwin = view_formatted_text_buffer(title, buf, 70, 350);
    g_free(title);

    return vwin;
}

static void fncall_help (GtkWidget *w, call_info *cinfo)
{
    const char *fnname = user_function_name_by_index(cinfo->iface);
    char *pdfname = NULL;
    PRN *prn;
    int err;

    /* FIXME incorpoate this into markup */
    if (user_function_has_PDF_doc(fnname, &pdfname)) {
	err = display_gfn_help(pdfname);
	free(pdfname);
	if (!err) {
	    return;
	}
    }

    if (bufopen(&prn)) {
	return;
    }
    
    err = user_function_help(fnname, OPT_M, prn);

    if (err) {
	gretl_print_destroy(prn);
	errbox("Couldn't find any help");
    } else {
	char *buf = gretl_print_steal_buffer(prn);

	make_help_viewer(fnname, buf);
	free(buf);
	gretl_print_destroy(prn);
    }
}

static int combo_list_index (const gchar *s, GList *list)
{
    GList *mylist = list;
    int i = 0;

    while (mylist != NULL) {
	if (!strcmp(s, (gchar *) mylist->data)) {
	    return i;
	}
	mylist = mylist->next;
	i++;
    }
    
    return -1;
}

static void update_matrix_selectors (call_info *cinfo)
{
    GList *slist = cinfo->msels;
    GList *mlist = NULL;
    GtkComboBox *sel;
    const char *mname;
    gchar *saved;
    int nm = n_user_matrices();
    int i, old;

    for (i=0; i<nm; i++) {
	mname = get_matrix_name_by_index(i);
	mlist = g_list_append(mlist, (gpointer) mname);
    }

    while (slist != NULL) {
	sel = GTK_COMBO_BOX(slist->data);
	saved = gtk_combo_box_get_active_text(sel);
	depopulate_combo_box(sel);
	set_combo_box_strings_from_list(sel, mlist);
	if (saved != NULL) {
	    old = combo_list_index(saved, mlist);
	    gtk_combo_box_set_active(sel, (old >= 0)? old : 0);
	    g_free(saved);
	}
	slist = slist->next;
    }

    g_list_free(mlist);
}

static int combo_accepts_null (GtkComboBox *c)
{
    gpointer p = g_object_get_data(G_OBJECT(c), "null_OK");

    return (p != NULL);
}

static void update_list_selectors (call_info *cinfo)
{
    GList *slist = cinfo->lsels;
    GList *llist1 = NULL;
    GList *llist2 = NULL;
    GtkComboBox *sel;
    const char *lname;
    gchar *saved;
    int nl = n_saved_lists();
    int null_OK;
    int i, old;

    llist2 = g_list_append(llist2, "null");

    for (i=0; i<nl; i++) {
	lname = get_list_name_by_index(i);
	llist1 = g_list_append(llist1, (gpointer) lname);
	llist2 = g_list_append(llist2, (gpointer) lname);
    }

    while (slist != NULL) {
	GList *llist;

	sel = GTK_COMBO_BOX(slist->data);
	null_OK = combo_accepts_null(sel);
	llist = (null_OK)? llist2 : llist1; 
	saved = gtk_combo_box_get_active_text(sel);
	depopulate_combo_box(sel);
	set_combo_box_strings_from_list(sel, llist);
	if (saved != NULL) {
	    old = combo_list_index(saved, llist);
	    gtk_combo_box_set_active(sel, (old >= 0)? old : 0);
	    g_free(saved);
	} else if (null_OK) {
	    gtk_combo_box_set_active(sel, 0);
	}
	slist = slist->next;
    }

    g_list_free(llist1);
    g_list_free(llist2);
}

int do_make_list (selector *sr)
{
    const char *buf = selector_list(sr);
    const char *lname = selector_entry_text(sr);
    gpointer data = selector_get_data(sr);
    call_info *cinfo = NULL;
    const char *msg = NULL;
    PRN *prn = NULL;
    int *list = NULL;
    int err = 0;

    if (data != NULL) {
	GtkWidget *w = GTK_WIDGET(data);
	
	cinfo = g_object_get_data(G_OBJECT(w), "cinfo");
    }

    if (lname == NULL || *lname == 0) {
	errbox(_("No name was given for the list"));
	return 1;
    }   

    if (buf == NULL || *buf == '\0') {
	int resp;

	resp = yes_no_dialog("gretl", _("Really create an empty list?"), 0);
	if (resp == GRETL_YES) {
	    list = gretl_null_list();
	    if (list == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    return 0;
	}
    } else {
	list = gretl_list_from_string(buf, &err);
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
	gtk_widget_hide(selector_get_window(sr));
	if (cinfo != NULL) {
	    update_list_selectors(cinfo);
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

static void launch_list_maker (GtkWidget *w, GtkWidget *entry)
{
    selector *sr;
    GtkWidget *dlg;

    sr = simple_selection(_("Define list"), do_make_list, DEFINE_LIST, 
			  entry);
    dlg = selector_get_window(sr);
    gtk_window_set_keep_above(GTK_WINDOW(dlg), TRUE);
}

void gui_define_list (void)
{
    launch_list_maker(NULL, NULL);
}

static void launch_matrix_maker (GtkWidget *w, call_info *cinfo)
{
    int n = n_user_matrices();

    gui_new_matrix();

    if (n_user_matrices() > n) {
	update_matrix_selectors(cinfo);
    }

    gtk_window_present(GTK_WINDOW(cinfo->dlg));
}

static int spinner_arg (call_info *cinfo, int i)
{
    double x = fn_param_minval(cinfo->func, i);
    double y = fn_param_maxval(cinfo->func, i);

    return !na(x) && !na(y);
}

static GtkWidget *bool_arg_selector (call_info *cinfo, int i)
{
    double deflt = fn_param_default(cinfo->func, i);
    int active = !na(deflt) && deflt != 0.0;
    GtkWidget *button;

    button = gtk_check_button_new();
    g_object_set_data(G_OBJECT(button), "argnum", GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(update_bool_arg), cinfo);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), active);
    cinfo->args[i] = g_strdup((active)? "1" : "0");

    return button;
}

static GtkWidget *spin_arg_selector (call_info *cinfo, int i)
{
    int minv = (int) fn_param_minval(cinfo->func, i);
    int maxv = (int) fn_param_maxval(cinfo->func, i);
    double deflt = fn_param_default(cinfo->func, i);
    int initv = (na(deflt))? minv : (int) deflt;
    GtkObject *adj;
    GtkWidget *spin;

    adj = gtk_adjustment_new(initv, minv, maxv, 1, 1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    g_object_set_data(G_OBJECT(spin), "argnum", GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed", 
		     G_CALLBACK(update_int_arg), cinfo);

    cinfo->args[i] = g_strdup_printf("%d", (na(deflt))? minv : 
				     (int) deflt);

    return spin;
}

static GtkWidget *combo_arg_selector (call_info *cinfo, int ptype, int i)
{
    GList *list = NULL;
    GtkWidget *combo;
    GtkWidget *entry;

    combo = gtk_combo_box_entry_new_text();
    entry = gtk_bin_get_child(GTK_BIN(combo));
    g_object_set_data(G_OBJECT(entry), "cinfo", cinfo);
    g_object_set_data(G_OBJECT(combo), "argnum", GINT_TO_POINTER(i));
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_arg), cinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    if (i >= 0 && fn_param_optional(cinfo->func, i)) {
	g_object_set_data(G_OBJECT(combo), "null_OK", GINT_TO_POINTER(1));
    }

    list = get_selection_list(cinfo, i, ptype, 1);
    if (list != NULL) {
	set_combo_box_strings_from_list(GTK_COMBO_BOX(combo), list);
	g_list_free(list);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    } 

    /* FIXME bool etc */

    if (ptype == GRETL_TYPE_INT || ptype == GRETL_TYPE_DOUBLE) {
	double x = fn_param_default(cinfo->func, i);

	if (!na(x)) {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	    g_free(tmp);
	}
    }

    return combo;
}

static void add_table_header (GtkWidget *tbl, gchar *txt,
			      int cols, int r0)
{
    GtkWidget *label = gtk_label_new(txt);
    GtkWidget *align = gtk_alignment_new(0.0, 0.5, 0.0, 0.0);

    gtk_container_add(GTK_CONTAINER(align), label);
    gtk_table_attach(GTK_TABLE(tbl), align, 0, cols, r0, r0 + 1,
		     GTK_FILL, GTK_FILL, 5, 5);
}

static void add_table_cell (GtkWidget *tbl, GtkWidget *w,
			    int c0, int c1, int r0)
{
    gtk_table_attach(GTK_TABLE(tbl), w, c0, c1, r0, r0 + 1,
		     GTK_FILL, GTK_FILL, 5, 3);
}

static void set_close_on_OK (GtkWidget *b, gpointer p)
{
    close_on_OK = button_is_active(b);
}

#define cinfo_has_return(c) (c->rettype != GRETL_TYPE_NONE && \
			     c->rettype != GRETL_TYPE_VOID)

#define cinfo_offer_return(c) (c->rettype != GRETL_TYPE_NONE && \
			       c->rettype != GRETL_TYPE_VOID && \
			       c->rettype != GRETL_TYPE_BUNDLE)

static void function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *sel, *tbl = NULL;
    GtkWidget *vbox, *hbox, *bbox;
    gchar *txt;
    const char *fnname;
    int trows = 0, tcols = 0;
    int i, row = 0;
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

    fnname = user_function_name_by_index(cinfo->iface);

    cinfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    txt = g_strdup_printf("gretl: %s", fnname);
    gtk_window_set_title(GTK_WINDOW(cinfo->dlg), txt);
    g_free(txt);
    gretl_emulated_dialog_add_structure(cinfo->dlg, &vbox, &bbox);
    open_fncall_dlg = cinfo->dlg;
    g_signal_connect(G_OBJECT(cinfo->dlg), "destroy",
		     G_CALLBACK(fncall_dialog_destruction), cinfo);

    cinfo->top_hbox = hbox = label_hbox(vbox, fnname, 5, 0);

    if (cinfo->n_params > 0) {
	tcols = (cinfo->extracol)? 3 : 2;
	trows = cinfo->n_params + 1;
	if (cinfo_offer_return(cinfo)) { 
	    trows += 4;
	}
    } else if (cinfo_offer_return(cinfo)) {
	tcols = 2;
	trows = 3;
    }

    if (trows > 0 && tcols > 0) {
	tbl = gtk_table_new(trows, tcols, FALSE);
    }

    if (cinfo->n_params > 0) {
	add_table_header(tbl, _("Select arguments:"), tcols, row);
    }

    for (i=0; i<cinfo->n_params; i++) {
	const char *desc = fn_param_descrip(cinfo->func, i);
	int ptype = fn_param_type(cinfo->func, i);

	if (desc != NULL) {
	    label = gtk_label_new(desc);
	} else {
	    txt = g_strdup_printf("%s (%s)",
				  fn_param_name(cinfo->func, i), 
				  arg_type_string(ptype));
	    label = gtk_label_new(txt);
	    g_free(txt);			     
	}

	row++;
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, row);

	if (ptype == GRETL_TYPE_BOOL) {
	    sel = bool_arg_selector(cinfo, i);
	} else if (ptype == GRETL_TYPE_INT && spinner_arg(cinfo, i)) {
	    sel = spin_arg_selector(cinfo, i);
	} else {
	    sel = combo_arg_selector(cinfo, ptype, i);
	}

	add_table_cell(tbl, sel, 1, 2, row);

	if (ptype == GRETL_TYPE_LIST) {
	    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(sel));
	    
	    cinfo->lsels = g_list_append(cinfo->lsels, sel);
	    button = gtk_button_new_with_label(_("More..."));
	    add_table_cell(tbl, button, 2, 3, row);
	    g_object_set_data(G_OBJECT(entry), "argnum", 
			      GINT_TO_POINTER(i+1));
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_list_maker),
			     entry);
	} else if (ptype == GRETL_TYPE_MATRIX) {
	    cinfo->msels = g_list_append(cinfo->msels, sel);
	    button = gtk_button_new_with_label(_("New..."));
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_matrix_maker), 
			     cinfo);
	} 
    }
	
    if (cinfo_offer_return(cinfo)) {
	GtkWidget *child;
	GList *list = NULL;

	if (cinfo->n_params > 0) {
	    row++;
	}	    

	add_table_header(tbl, _("Assign return value (optional):"), tcols, row);
	row++;

	label = gtk_label_new(_("selection (or new variable)"));
	add_table_cell(tbl, label, 1, 2, row);
	row++;

	label = gtk_label_new(arg_type_string(cinfo->rettype));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, row);

	sel = gtk_combo_box_entry_new_text();
	g_signal_connect(G_OBJECT(sel), "changed",
			 G_CALLBACK(update_return), cinfo);
	list = get_selection_list(cinfo, -1, cinfo->rettype, 0);
	if (list != NULL) {
	    set_combo_box_strings_from_list(GTK_COMBO_BOX(sel), list);
	    g_list_free(list);
	}
	child = gtk_bin_get_child(GTK_BIN(sel));
	gtk_entry_set_activates_default(GTK_ENTRY(child), TRUE);
	add_table_cell(tbl, sel, 1, 2, row);
    }

    if (tbl != NULL) {
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

    gtk_widget_show_all(cinfo->dlg);
}

static int function_data_check (call_info *cinfo)
{
    int i, err = 0;

    if (cinfo->dreq != FN_NODATA_OK) {
	if (datainfo == NULL || datainfo->v == 0) {
	    warnbox(_("Please open a data file first"));
	    return 1;
	}
    }

    for (i=0; i<cinfo->n_params; i++) {
	int type = fn_param_type(cinfo->func, i);

	if (type == GRETL_TYPE_SERIES || type == GRETL_TYPE_LIST ||
	    type == GRETL_TYPE_SERIES_REF) {
	    if (datainfo == NULL || datainfo->v == 0) {
		warnbox(_("Please open a data file first"));
		err = 1;
		break;
	    }
	}
	if (type == GRETL_TYPE_LIST) {
	    cinfo->extracol = 1;
	} else if (type == GRETL_TYPE_MATRIX || type == GRETL_TYPE_MATRIX_REF) {
	    cinfo->extracol = 1;
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
    gchar *s = cinfo->args[i];
    int err = 0;

    *add = 0;

    if (!gretl_ref_type(t)) {
	return 0;
    }

    if (*s == '&' || !strcmp(s, "null")) {
	return 0;
    }

    if (t == GRETL_TYPE_MATRIX_REF) {
	/* handle case where indirect return matrix does not yet exist */
	if (get_matrix_by_name(s) == NULL) {
	    gretl_matrix *m = gretl_null_matrix_new();

	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = add_or_replace_user_matrix(m, s);
	    }
	    if (!err) {
		pprintf(prn, "? matrix %s\n", s);
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

static int pre_process_args (call_info *cinfo, PRN *prn)
{
    char auxline[MAXLINE];
    char auxname[VNAMELEN+2];
    int i, add = 0, err = 0;

    for (i=0; i<cinfo->n_params && !err; i++) {
	if (should_addressify_var(cinfo, i)) {
	    sprintf(auxname, "FNARG%d", i + 1);
	    sprintf(auxline, "genr %s=%s", auxname, cinfo->args[i]);
	    err = generate(auxline, &Z, datainfo, OPT_NONE, NULL);
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
    }

    return err;
}

#define AUTO_BUNDLE "BUNDLE_RET__"

static int real_GUI_function_call (call_info *cinfo, PRN *prn)
{
    ExecState state;
    char fnline[MAXLINE];
    const char *funname;
    gretl_bundle *bundle = NULL;
    int attach_bundle = 0;
    int orig_v = datainfo->v;
    int i, err = 0;

    funname = user_function_name_by_index(cinfo->iface);
    *fnline = 0;

    /* compose the function command-line */

    if (cinfo->ret != NULL) {
	strcat(fnline, cinfo->ret);
	strcat(fnline, " = ");
    } else if (cinfo->rettype == GRETL_TYPE_BUNDLE) {
	strcat(fnline, AUTO_BUNDLE);
	strcat(fnline, " = ");
	attach_bundle = 1;
    }	

    strcat(fnline, funname);
    strcat(fnline, "(");

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    strcat(fnline, cinfo->args[i]);
	    if (i < cinfo->n_params - 1) {
		strcat(fnline, ", ");
	    }
	}
    }

    strcat(fnline, ")");

    if (!attach_bundle) {
	pprintf(prn, "? %s\n", fnline);
    }

    /* note: gretl_exec_state_init zeros the first byte of the
       supplied 'line' */

    gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			  models, prn);
    state.line = fnline;

#if USE_GTK_SPINNER
    start_wait_for_output(cinfo->top_hbox, 0); 
#endif
    err = gui_exec_line(&state, &Z, datainfo);
#if USE_GTK_SPINNER
    stop_wait_for_output(cinfo->top_hbox);
#endif

    /* destroy any "ARG" vars or matrices that were created? */

    if (!err && attach_bundle) {
	bundle = gretl_bundle_pull_from_stack(AUTO_BUNDLE, &err);
    }

    view_buffer(prn, 80, 400, funname, 
		(bundle == NULL)? PRINT : VIEW_BUNDLE,
		bundle);

    if (err) {
	gui_errmsg(err);
    } 

    if (datainfo->v > orig_v) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    return err;
}

/* In case a function package offers more than one public
   interface, give the user a selector: for four or fewer
   options we use radio buttons, otherwise we use a pull-down
   list.
*/

static void select_interface (call_info *cinfo)
{
    const char *funname;
    char **opts = NULL;
    GList *ilist = NULL;
    int npub = cinfo->publist[0];
    int radios = (npub < 5);
    int i, nopts = 0;
    int resp, err = 0;

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
    } else if (radios) {
	resp = radio_dialog("gretl", "select function", 
			    (const char **) opts, 
			    nopts, 0, 0);
	if (resp >= 0) {
	    cinfo->iface = cinfo->publist[resp+1];
	} else {
	    cinfo->iface = -1;
	}
	free_strings_array(opts, nopts);
    } else {
	resp = combo_selector_dialog(ilist, "select function", 0);
	if (resp >= 0) {
	    cinfo->iface = cinfo->publist[resp+1];
	} else {
	    cinfo->iface = -1;
	}	
	g_list_free(ilist);
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
	int err;

	err = bufopen(&prn);

	if (!err && cinfo->args != NULL) {
	    err = pre_process_args(cinfo, prn);
	    if (err) {
		gui_errmsg(err);
	    }
	}

	if (!err) {
	    err = real_GUI_function_call(cinfo, prn);
	} else {
	    gretl_print_destroy(prn);
	}

	if (close_on_OK) {
	    gtk_widget_destroy(cinfo->dlg);
	}
    }
}

/* call to execute a function from the specified package: we do
   this only for locally installed packages */

void call_function_package (const char *fname, GtkWidget *w,
			    int *loaderr)
{
    int minver = 0;
    call_info *cinfo;
    fnpkg *pkg;
    int err = 0;

    pkg = get_function_package_by_filename(fname, &err);

    if (err) {
	gui_errmsg(err); /* FIXME error not very informative? */
	if (loaderr != NULL) {
	    *loaderr = 1;
	}
	return;
    }

    cinfo = cinfo_new();
    if (cinfo == NULL) {
	return;
    }

    /* get the interface list and other info for package */

    err = function_package_get_properties(pkg,
					  "publist", &cinfo->publist,
					  "data-requirement", &cinfo->dreq,
					  "min-version", &minver,
					  NULL);

    if (err) {
	gui_errmsg(err);
    } else if (cinfo->publist == NULL) {
	/* no available interfaces */
	err = E_DATA;
	errbox(_("Function package is broken"));
    }

    if (!err) {
	/* do we have suitable data in place? */
	err = check_function_needs(datainfo, cinfo->dreq, minver);
	if (err) {
	    gui_errmsg(err);
	}
    }

    if (!err) {
	if (cinfo->publist[0] > 1) {
	    select_interface(cinfo);
	    if (cinfo->iface < 0) {
		/* failed, or cancelled */
		cinfo_free(cinfo);
		return; /* note: handled */
	    }
	} else {
	    /* only one interface available */
	    cinfo->iface = cinfo->publist[1];
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
	function_call_dialog(cinfo);
    } else {
	cinfo_free(cinfo);
    }
}

void function_call_cleanup (void)
{
    if (open_fncall_dlg != NULL) {
	gtk_widget_destroy(open_fncall_dlg);
    }
}

#if 1

/* This is quite complicated, because we're bypassing the "standard"
   mechanism for a function call, and instead using output from the
   gretl model selection dialog and constructing a "direct" call to
   the function, as is done inside geneval.c. There's probably a
   simpler way of doing this.
*/

static int ols_bundle_callback (selector *sr)
{
    const char *funname = "ols_pack";
    const char *pkgname = "olsbundle";
    gchar *path = NULL;
    ufunc *uf = NULL;
    fnargs *args = NULL;
    int yno, err = 0;

    if (selector_error(sr)) {
	return 1;
    }

    /* get path to package; load package if not already loaded;
       get specific function from package */

    path = gretl_function_package_get_path(pkgname);
    if (path == NULL) {
	gretl_errmsg_sprintf("Couldn't find package '%s'", pkgname);
	err = E_DATA;
    } else {
	fnpkg *pkg = get_function_package_by_filename(path, &err);

	if (!err) {
	    uf = get_function_from_package(funname, pkg);
	}
	if (uf == NULL) {
	    gretl_errmsg_sprintf("Couldn't find function '%s'", funname);
	    err = E_DATA;
	}
	g_free(path);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    } else {
	/* args: get ID of dependent var plus X list */
	const char *buf = selector_list(sr);

	if (buf == NULL) {
	    err = E_DATA;
	} else {
	    int *list = gretl_list_from_string(buf, &err);

	    if (!err) {
		yno = list[1];
		gretl_list_delete_at_pos(list, 1);
		remember_list(list, "arg1temp", NULL);
		free(list);
	    }
	}
    }

    if (!err) {
	/* push args to function */
	args = fn_args_new();
	if (args == NULL) {
	    err = E_ALLOC;
	} else {
	    err = push_fn_arg(args, GRETL_TYPE_USERIES, &yno);
	    if (!err) {
		err = push_fn_arg(args, GRETL_TYPE_LIST, "arg1temp");
	    }
	}
    }

    if (!err) {
	/* execute function and retrieve bundle */
	gretl_bundle *bundle = NULL;
	PRN *prn = NULL;

	err = bufopen(&prn);

	if (!err) {
	    err = gretl_function_exec(uf, args, GRETL_TYPE_BUNDLE, &Z, datainfo, 
				      &bundle, NULL, prn);
	}
	if (!err) {
	    view_buffer(prn, 80, 400, funname, VIEW_BUNDLE, bundle);
	}
    }

    if (err) {
	gui_errmsg(err);
    }

    fn_args_free(args);
    delete_list_by_name("arg1temp");

    return err;
}

#endif

/* Callback for a menu item representing a function package, whose
   name (e.g. "gig") is attached to @action. We first see if we can
   find the full path to the corresponding gfn file; if so we
   initiate a GUI call to the package.
*/

static void gfn_menu_callback (GtkAction *action, windata_t *vwin)
{
    const gchar *name = gtk_action_get_name(action);
    gchar *path;

    path = gretl_function_package_get_path(name);

    if (path == NULL) {
	errbox("Couldn't find package %s\n", name);
    } else if (!strcmp(name, "olsbundle")) {
	/* see if we can get cosy with olsbundle */
	selection_dialog(_("gretl: specify model"), 
			 ols_bundle_callback, OLS);
    } else {
	call_function_package(path, vwin->main, NULL);
	g_free(path);
    }
}

/* For "privileged" function packages: given the internal
   package name (e.g. "gig") and a menu path where we'd like
   it to appear (e.g. "/menubar/Model/TSModels" -- see the
   ui definition file gui2/gretlmain.xml), construct the
   appropriate menu item and connect it to gfn_menu_callback(), 
   which see above.
*/

static void maybe_add_package_to_menu (const char *pkgname, 
				       const char *menupath,
				       windata_t *vwin)
{
    static GtkActionEntry pkg_item = {
	NULL, NULL, NULL, NULL, NULL, G_CALLBACK(gfn_menu_callback)
    };
    GtkActionGroup *actions;
    gchar *label;

    pkg_item.name = pkgname;
    label = g_strdup_printf("_%s...", pkgname);
    pkg_item.label = label;

    gtk_ui_manager_add_ui(vwin->ui, gtk_ui_manager_new_merge_id(vwin->ui),
			  menupath, label, pkgname,
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);

    actions = gtk_action_group_new(pkgname);
    gtk_action_group_add_actions(actions, &pkg_item, 1, vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);
    g_free(label);
}

void maybe_add_packages_to_menus (windata_t *vwin)
{
    maybe_add_package_to_menu("gig", "/menubar/Model/TSModels", vwin);
#if 0
    maybe_add_package_to_menu("olsbundle", "/menubar/Model", vwin);
#endif
}


/*
 *  Copyright (c) 2003 by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"

#include <gtk/gtk.h>

#undef PCA_DEBUG

static double *standardize (const double *x, int n)
{
    double *sx;
    double xbar, sd;
    int i, err;

    err = moments(0, n-1, x, &xbar, &sd, NULL, NULL, 1);
    if (err) return NULL;

    sx = malloc(n * sizeof *sx);
    if (sx == NULL) return NULL;

    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    sx[i] = NADBL;
	} else {
	    sx[i] = (x[i] - xbar) / sd;
	}
    }

    return sx;
}

struct flag_info {
    GtkWidget *dialog;
    gint *flag;
};

enum pca_flags {
    PCA_SAVE_NONE,
    PCA_SAVE_MAIN,
    PCA_SAVE_ALL
};

static gboolean destroy_pca_dialog (GtkWidget *w, struct flag_info *finfo)
{
    free(finfo);
    gtk_main_quit();
    return FALSE;
}

static gboolean set_pca_flag (GtkWidget *w, struct flag_info *finfo)
{
#if GTK_MAJOR_VERSION >= 2
    gint opt = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
#else
    gint opt = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "opt"));
#endif

    *(finfo->flag) = opt;
    return FALSE;
}

static gboolean cancel_set_flag (GtkWidget *w, struct flag_info *finfo)
{
    *(finfo->flag) = PCA_SAVE_NONE;
    gtk_widget_destroy(finfo->dialog);
    return FALSE;
}

static gboolean pca_dialog_finalize (GtkWidget *w, struct flag_info *finfo)
{
    gtk_widget_destroy(finfo->dialog);
    return FALSE;
}

static unsigned long pca_flag_dialog (void)
{
    struct flag_info *finfo;
    GtkWidget *dialog, *tempwid, *button, *hbox;
    GtkWidget *internal_vbox;
    GSList *group;
    gint flag = PCA_SAVE_MAIN;

    finfo = malloc(sizeof *finfo);
    if (finfo == NULL) return 0;

    dialog = gtk_dialog_new();

    finfo->dialog = dialog;
    finfo->flag = &flag;
    
    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: save data")); 
#if GTK_MAJOR_VERSION >= 2
    gtk_window_set_resizable (GTK_WINDOW (dialog), FALSE);
#endif
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);

    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

#if GTK_MAJOR_VERSION >= 2
    g_signal_connect (G_OBJECT(dialog), "destroy", 
		      G_CALLBACK(destroy_pca_dialog), finfo);
#else
    gtk_signal_connect (GTK_OBJECT(dialog), "destroy", 
			GTK_SIGNAL_FUNC(destroy_pca_dialog), finfo);
#endif

    internal_vbox = gtk_vbox_new (FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("Variables to save:"));
    gtk_box_pack_start (GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);
    gtk_box_pack_start (GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    /* Only those with eigenvalues > 1.0 */
    button = gtk_radio_button_new_with_label(NULL, 
					     _("Components with eigenvalues > 1.0"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#if GTK_MAJOR_VERSION >= 2
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_pca_flag), finfo);
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_MAIN)); 
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_pca_flag), finfo);
    gtk_object_set_data(GTK_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_MAIN)); 
#endif   
    gtk_widget_show (button);   

    /* All components */
#if GTK_MAJOR_VERSION >= 2
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
#else
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
#endif
    button = gtk_radio_button_new_with_label(group, _("All components"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
#if GTK_MAJOR_VERSION >= 2
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_pca_flag), finfo);
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_ALL)); 
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_pca_flag), finfo);
    gtk_object_set_data(GTK_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_ALL)); 
#endif
   
    gtk_widget_show (button);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    gtk_widget_show (internal_vbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    /* Create the "OK" button */
#if GTK_MAJOR_VERSION >= 2
    tempwid = gtk_button_new_from_stock (GTK_STOCK_OK);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(pca_dialog_finalize), finfo);
#else
    tempwid = gtk_button_new_with_label(_("OK"));
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(pca_dialog_finalize), finfo);
#endif
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* "Cancel" button */
#if GTK_MAJOR_VERSION >= 2
    tempwid = gtk_button_new_from_stock (GTK_STOCK_CANCEL);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(cancel_set_flag), finfo);
#else
    tempwid = gtk_button_new_with_label(_("Cancel"));
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(cancel_set_flag), finfo);
#endif    
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_widget_show (tempwid);

    gtk_widget_show(dialog);

    gtk_main();

    if (flag == PCA_SAVE_MAIN) return OPT_O;
    if (flag == PCA_SAVE_ALL) return OPT_A;

    return 0L;
}

static void pca_print (CORRMAT *corrmat, gretl_matrix *m,
		       double *evals, DATAINFO *pdinfo, 
		       PRN *prn)
{
    double x, y;
    int n = corrmat->list[0];
    int i, j, cols;

    pprintf(prn, "%s\n\n", _("Principal Components Analysis"));
    pprintf(prn, "%s\n\n", _("Eigenanalysis of the Correlation Matrix"));

    pputs(prn, _("Component  Eigenvalue  Proportion   Cumulative\n"));

    x = 0.0;
    y = 0.0;
    for (i=n-1; i>=0; i--) {
	y += evals[i] / n;
	pprintf(prn, "%5d%13.4f%13.4f%13.4f\n", n - i,
		evals[i], evals[i] / n, y);
	x += evals[i];
    }
    pputc(prn, '\n');

#ifdef PCA_DEBUG
    fprintf(stderr, "check: sum of evals = %g\n", x);
#endif

    pprintf(prn, "%s\n\n", _("Eigenvectors (component loadings)"));

    cols = n;
    while (cols > 0) {
	int colsdone = 0;

	pputs(prn, "Variable  ");
	for (i=n-cols; i<n-cols+7 && i<n; i++) {
	    char pcname[8];

	    sprintf(pcname, "PC%d", i + 1);
	    pprintf(prn, "%9s", pcname);
	    colsdone++;
	}
	pputc(prn, '\n');
	for (i=0; i<n; i++) {
	    pprintf(prn, "%-10s", pdinfo->varname[corrmat->list[i+1]]);
	    for (j=cols-1; j>cols-8 && j>=0; j--) {
		pprintf(prn, "%9.3f", gretl_matrix_get(m, i, j));
	    }
	    pputc(prn, '\n');
	}
	cols -= colsdone;
	pputc(prn, '\n');
    }
}

int pca_from_corrmat (CORRMAT *corrmat, double ***pZ,
		      DATAINFO *pdinfo, unsigned long *pflag,
		      PRN *prn)
{
    gretl_matrix *m;
    double x;
    int i, j, idx, n = corrmat->list[0];
    double *evals;
    unsigned long oflag = 0L;

    if (pflag != NULL) oflag = *pflag;

    if (oflag & OPT_D) { 
	oflag = pca_flag_dialog();
	if (!oflag) {
	    /* canceled */
	    *pflag = 0L;
	    return 0; 
	}
    }    

    m = gretl_matrix_alloc(n, n);
    if (m == NULL) return E_ALLOC;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    idx = ijton(i, j, n);
	    x = corrmat->xpx[idx];
	    gretl_matrix_set(m, i, j, x);
	}
    }

    evals = gretl_symmetric_matrix_eigenvals(m, 1);
    if (evals == NULL) {
	gretl_matrix_free(m);
	return 1;
    }

    if (prn != NULL) {
	pca_print(corrmat, m, evals, pdinfo, prn);
    }

    if (oflag) {
	/* add components with eigenvalues > 1 to the dataset */
	int v = pdinfo->v;
	int nc = 0, err = 0;
	double **sZ = NULL;
	int add_all = (oflag == OPT_A);
	int *list;

	if (add_all) {
	    nc = n;
	} else {
	    for (i=0; i<n; i++) {
		if (evals[i] > 1.0) nc++;
	    }
	}

	list = malloc((nc + 1) * sizeof *list);
	if (list == NULL) err = E_ALLOC;

	if (!err) {
	    /* build list of PCs (with eigenvals > 1?) */
	    list[0] = nc;
	    j = 1;
	    for (i=n-1; i>=0; i--) {
		if (add_all || evals[i] > 1.0) list[j++] = i;
	    }
#ifdef PCA_DEBUG
	    printlist(list, "pclist");
#endif
	    err = dataset_add_vars(nc, pZ, pdinfo);
	}

	if (!err) {
	    /* construct standardized versions of variables */
	    sZ = malloc(n * sizeof *sZ);
	    if (sZ == NULL) err = E_ALLOC;
	    else {
		for (i=0; i<n; i++) sZ[i] = NULL;
		for (i=0; i<n; i++) {
		    int oldv = corrmat->list[i+1];

#ifdef PCA_DEBUG
		    fprintf(stderr, "Getting standardized version of "
			    "var %d\n", oldv);
#endif
		    sZ[i] = standardize((const double *) (*pZ)[oldv], 
					pdinfo->n);
		    if (sZ[i] == NULL) {
			err = E_ALLOC;
			break;
		    }
		}
		if (err) {
		    for (i=0; i<n; i++) free(sZ[i]);
		    free(sZ);
		    sZ = NULL;
		}
	    }
	}

	if (!err) {
	    for (i=1; i<=list[0]; i++) {
		int newv = v + i - 1;
		int pcnum = list[i];
		int t;

		sprintf(pdinfo->varname[newv], "PC%d", i);
		make_varname_unique(pdinfo->varname[newv], newv, pdinfo);
		sprintf(VARLABEL(pdinfo, newv), "Component with "
			"eigenvalue = %.4f", evals[pcnum]);
		for (t=0; t<pdinfo->n; t++) {
#ifdef PCA_DEBUG
		    fprintf(stderr, "Obs %d\n", t);
#endif
		    (*pZ)[newv][t] = 0.0;
		    for (j=0; j<n; j++) {
			double load = gretl_matrix_get(m, j, pcnum);
			double val = sZ[j][t];

#ifdef PCA_DEBUG
			fprintf(stderr, "j=%d,pcnum=%d,load=%g,val=%g\n",
				j,pcnum,load,val);
#endif

			if (na(val)) {
			    (*pZ)[newv][t] = NADBL;
			    break;
			} else {
			    (*pZ)[newv][t] += load * val;
			}
		    }
		} /* end loop over observations */
	    } /* end loop over components */
	} /* end !err conditional */

	free(list);
	if (sZ != NULL) {
	    for (i=0; i<n; i++) free(sZ[i]);
	    free(sZ);
	}

    } /* end oflag conditional */

    free(evals);
    gretl_matrix_free(m);

    if (pflag != NULL) *pflag = oflag;

    return 0;
}

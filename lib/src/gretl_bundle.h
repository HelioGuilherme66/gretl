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

#ifndef GRETL_BUNDLE_H_
#define GRETL_BUNDLE_H_

#define BUNDLE_PRINT "bundle-print"
#define BUNDLE_PLOT  "bundle-plot"
#define BUNDLE_TEST  "bundle-test"
#define BUNDLE_FCAST "bundle-fcast"
#define BUNDLE_EXTRA "bundle-extra"
#define GUI_MAIN     "gui-main"
#define GUI_PRECHECK "gui-precheck"

typedef enum {
    BUNDLE_PLAIN,
    BUNDLE_KALMAN
} BundleType;

typedef struct gretl_bundle_ gretl_bundle;
typedef struct bundled_item_ bundled_item;

gretl_bundle *gretl_bundle_new (void);

int gretl_is_bundle (const char *name);

int type_can_be_bundled (GretlType type);

gretl_bundle *get_bundle_by_name (const char *name);

BundleType gretl_bundle_get_type (gretl_bundle *bundle);

void *gretl_bundle_get_content (gretl_bundle *bundle);

void *gretl_bundle_get_data (gretl_bundle *bundle, const char *key,
			     GretlType *type, int *size, int *err);

void *gretl_bundle_steal_data (gretl_bundle *bundle, const char *key,
			       GretlType *type, int *size, int *err);

void *gretl_bundle_get_private_data (gretl_bundle *bundle);

GretlType gretl_bundle_get_member_type (gretl_bundle *bundle,
					const char *key,
					int *err);

gretl_matrix *gretl_bundle_get_matrix (gretl_bundle *bundle,
				       const char *key,
				       int *err);

double *gretl_bundle_get_series (gretl_bundle *bundle,
				 const char *key,
				 int *n, int *err);

double gretl_bundle_get_scalar (gretl_bundle *bundle,
				const char *key,
				int *err);

int gretl_bundle_get_int (gretl_bundle *bundle,
			  const char *key,
			  int *err);

const char *gretl_bundle_get_string (gretl_bundle *bundle,
				     const char *key,
				     int *err);

void *gretl_bundle_get_array (gretl_bundle *bundle,
			      const char *key,
			      int *err);

const char *gretl_bundle_get_note (gretl_bundle *bundle, const char *key);

const char *gretl_bundle_get_creator (gretl_bundle *bundle);

void *bundled_item_get_data (bundled_item *item, GretlType *type,
			     int *size);

const char *bundled_item_get_note (bundled_item *item);

int gretl_bundle_donate_data (gretl_bundle *bundle, const char *key,
			      void *ptr, GretlType type, int size);

int gretl_bundle_set_data (gretl_bundle *bundle, const char *key,
			   void *ptr, GretlType type, int size);

int gretl_bundle_set_string (gretl_bundle *bundle, const char *key,
			     const char *str);

int gretl_bundle_set_scalar (gretl_bundle *bundle, const char *key,
			     double val);

int gretl_bundle_set_int (gretl_bundle *bundle, const char *key,
			  int val);

int gretl_bundle_set_series (gretl_bundle *bundle, const char *key,
			     const double *x, int n);

int gretl_bundle_set_matrix (gretl_bundle *bundle, const char *key,
			     const gretl_matrix *m);

int gretl_bundle_set_note (gretl_bundle *bundle, const char *key,
			   const char *note);

int gretl_bundle_delete_data (gretl_bundle *bundle, const char *key);

int gretl_bundle_add_or_replace (gretl_bundle *bundle, const char *name);

int gretl_bundle_copy_as (const char *name, const char *copyname);

gretl_bundle *gretl_bundle_copy (const gretl_bundle *bundle, int *err);

int gretl_bundle_set_name (gretl_bundle *b, const char *name);

int gretl_bundle_set_creator (gretl_bundle *b, const char *name);

int gretl_bundle_print (gretl_bundle *bundle, PRN *prn);

int gretl_bundle_is_stacked (gretl_bundle *b);

int gretl_bundle_get_n_keys (gretl_bundle *b);

int gretl_bundle_get_n_members (gretl_bundle *b);

int gretl_bundle_has_content (gretl_bundle *b);

gretl_bundle *gretl_bundle_pull_from_stack (const char *name,
					    int *err);

gretl_bundle *gretl_bundle_union (const gretl_bundle *bundle1,
				  const gretl_bundle *bundle2,
				  int *err);

void gretl_bundle_destroy (gretl_bundle *bundle);

void gretl_bundle_void_content (gretl_bundle *bundle);

void gretl_bundle_serialize (gretl_bundle *b, const char *name, 
			     FILE *fp);

gretl_bundle *gretl_bundle_deserialize (void *p1, void *p2, 
					int *err);

int gretl_bundle_write_to_file (gretl_bundle *b, 
				const char *fname,
				int to_dotdir);

gretl_bundle *gretl_bundle_read_from_file (const char *fname, 
					   int from_dotdir,
					   int *err);

gretl_bundle *gretl_bundle_read_from_buffer (const char *buf, 
					     int len,
					     int *err);

int bundle_contains_data (gretl_bundle *b, void *data);

gretl_bundle *get_sysinfo_bundle (int *err);

void *sysinfo_bundle_get_data (const char *key,
			       GretlType *type, 
			       int *err);

gretl_bundle *bundle_from_model (MODEL *pmod,
				 DATASET *dset,
				 int *err);

gretl_bundle *kalman_bundle_new (gretl_matrix *M[],
				 int copy[], int nmat,
				 int *err);

void gretl_bundle_cleanup (void);

#endif /* GRETL_BUNDLE_H_ */




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

typedef struct gretl_bundle_ gretl_bundle;
typedef struct bundled_item_ bundled_item;

int gretl_is_bundle (const char *name);

gretl_bundle *get_gretl_bundle_by_name (const char *name);

void *gretl_bundle_get_content (gretl_bundle *bundle);

void *gretl_bundle_get_data (gretl_bundle *bundle, const char *key,
			     GretlType *type, int *size);

void *bundled_item_get_data (bundled_item *item, GretlType *type,
			     int *size);

int gretl_bundle_set_data (gretl_bundle *bundle, const char *key,
			   void *ptr, GretlType type, int size);

int gretl_bundle_add_or_replace (gretl_bundle *bundle, const char *name);

int save_named_bundle (const char *name);

int gretl_bundle_copy_as (const char *name, const char *copyname);

gretl_bundle *gretl_bundle_copy (gretl_bundle *bundle, int *err);

int gretl_bundle_delete_by_name (const char *name, PRN *prn);

gretl_bundle *gretl_bundle_pull_from_stack (const char *name,
					    int *err);

int gretl_bundle_localize (const char *origname,
			   const char *localname);

int gretl_bundle_unlocalize (const char *localname,
			     const char *origname);

gretl_bundle *gretl_bundle_union (const gretl_bundle *bundle1,
				  const gretl_bundle *bundle2,
				  int *err);

void gretl_bundle_destroy (gretl_bundle *bundle);

int destroy_saved_bundles_at_level (int level);

void destroy_user_bundles (void);




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

#define FULL_XML_HEADERS 1

#include "libgretl.h"
#include "gretl_func.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_array.h"
#include "gretl_xml.h"
#include "gretl_foreign.h"
#include "gretl_typemap.h"
#include "genparse.h"
#include "kalman.h"
#include "gretl_bundle.h"

#define BDEBUG 0

/**
 * gretl_bundle:
 *
 * An opaque type; use the relevant accessor functions.
 */

struct gretl_bundle_ {
    BundleType type; /* see enum above */
    GHashTable *ht;  /* holds key/value pairs */
    char *creator;   /* name of function that built the bundle */
    void *data;      /* holds pointer to struct for some uses */
};

/**
 * bundled_item:
 *
 * An item of data within a gretl_bundle. This is an
 * opaque type; use the relevant accessor functions.
 */

struct bundled_item_ {
    GretlType type;
    int size;
    gpointer data;
    char *note;
};

static gretl_bundle *sysinfo_bundle;

static int real_bundle_set_data (gretl_bundle *b, const char *key,
				 void *ptr, GretlType type,
				 int size, int copy,
				 const char *note);

int gretl_bundle_set_name (gretl_bundle *b, const char *name)
{
    user_var *u = get_user_var_by_data(b);
    
    if (u == NULL) {
	return E_UNKVAR;
    } else {
	user_var_set_name(u, name);
	return 0;
    }
}

int gretl_bundle_is_stacked (gretl_bundle *b)
{
    user_var *u = get_user_var_by_data(b);

    return u != NULL;
}

/* gets the number of keys in the bundle's hash table */

int gretl_bundle_get_n_keys (gretl_bundle *b)
{
    int n_items = 0;

    if (b != NULL && b->ht != NULL) {
	n_items = g_hash_table_size(b->ht);
    }

    return n_items;
}

/* gets total number of members including any "special"
   contents outside of the hash table */

int gretl_bundle_get_n_members (gretl_bundle *b)
{
    int nmemb = 0;

    if (b != NULL) {
	if (b->type == BUNDLE_KALMAN) {
	    nmemb += kalman_bundle_n_members(b);
	}
	if (b->ht != NULL) {
	    nmemb += g_hash_table_size(b->ht);
	}
    }

    return nmemb;
}

int gretl_bundle_has_content (gretl_bundle *b)
{
    int ret = 0;
    
    if (b != NULL && b->ht != NULL &&
	(b->type == BUNDLE_KALMAN ||
	 g_hash_table_size(b->ht) > 0)) {
	ret = 1;
    }    

    return ret;
}

int type_can_be_bundled (GretlType type)
{
    if (type == GRETL_TYPE_INT ||
	type == GRETL_TYPE_BOOL) {
	type = GRETL_TYPE_DOUBLE;
    }

    return (type == GRETL_TYPE_DOUBLE ||
	    type == GRETL_TYPE_STRING ||
	    type == GRETL_TYPE_MATRIX ||
	    type == GRETL_TYPE_MATRIX_REF ||
	    type == GRETL_TYPE_SERIES ||
	    type == GRETL_TYPE_BUNDLE ||
	    type == GRETL_TYPE_ARRAY);
}

/* allocate and fill out a 'value' (type plus data pointer) that will
   be inserted into a bundle's hash table */

static bundled_item *bundled_item_new (GretlType type, void *ptr, 
				       int size, int copy,
				       const char *note,
				       int *err)
{
    bundled_item *item = malloc(sizeof *item);

    if (item == NULL) {
	*err = E_ALLOC;
    } else {
	item->type = type;
	item->size = 0;
	item->note = NULL;

	switch (item->type) {
	case GRETL_TYPE_DOUBLE:
	    item->data = malloc(sizeof(double));
	    if (item->data != NULL) {
		double *dp = item->data;

		*dp = *(double *) ptr;
	    }
	    break;
	case GRETL_TYPE_STRING:	
	    if (copy) {
		item->data = gretl_strdup((char *) ptr);
	    } else {
		item->data = ptr;
	    }
	    break;
	case GRETL_TYPE_MATRIX:
	    if (copy) {
		item->data = gretl_matrix_copy((gretl_matrix *) ptr);
	    } else {
		item->data = ptr;
	    }
	    break;
	case GRETL_TYPE_MATRIX_REF:
	    item->data = ptr;
	    break;
	case GRETL_TYPE_SERIES:
	    if (copy) {
		item->data = copyvec((const double *) ptr, size);
	    } else {
		item->data = ptr;
	    }
	    item->size = size;
	    break;
	case GRETL_TYPE_BUNDLE:
	    if (copy) {
		item->data = gretl_bundle_copy((gretl_bundle *) ptr, err);
	    } else {
		item->data = ptr;
	    }
	    break;
	case GRETL_TYPE_ARRAY:
	    if (copy) {
		item->data = gretl_array_copy((gretl_array *) ptr, err);
	    } else {
		item->data = ptr;
	    }
	    break;
	default:
	    *err = E_TYPES;
	    break;
	}

	if (!*err && item->data == NULL) {
	    free(item);
	    item = NULL;
	    *err = E_ALLOC;
	}

	if (item != NULL && note != NULL) {
	    item->note = gretl_strdup(note);
	}	
    }

    return item;
}

static void release_matrix_pointer (gretl_matrix **pm)
{
    const void *data = *pm;

    if (get_user_var_by_data(data) == NULL) {
	/* the bundle now has the only pointer to
	   this matrix */
	gretl_matrix_free(*pm);
    }

    *pm = NULL;
}

static int bundled_item_replace_data (bundled_item *item,
				      GretlType type, void *ptr, 
				      int size, int copy)
{
    int err = 0;

    if (ptr == item->data) {
	return 0;
    }

    if (item->type == GRETL_TYPE_DOUBLE) {
	double *dp = item->data;

	*dp = *(double *) ptr;
    } else if (item->type == GRETL_TYPE_STRING) {
	free(item->data);
	if (copy) {
	    item->data = gretl_strdup((char *) ptr);
	} else {
	    item->data = ptr;
	}
    } else if (item->type == GRETL_TYPE_MATRIX) {
	gretl_matrix_free(item->data);
	if (copy) {
	    item->data = gretl_matrix_copy((gretl_matrix *) ptr);
	} else {
	    item->data = ptr;
	}
    } else if (item->type == GRETL_TYPE_MATRIX_REF) {
	release_matrix_pointer((gretl_matrix **) &item->data);
	item->data = ptr;
    } else if (item->type == GRETL_TYPE_SERIES) {
	free(item->data);
	if (copy) {
	    item->data = copyvec((const double *) ptr, size);
	} else {
	    item->data = ptr;
	}
	item->size = size;
    } else if (item->type == GRETL_TYPE_BUNDLE) {
	gretl_bundle_destroy((gretl_bundle *) item->data);
	if (copy) {
	    item->data = gretl_bundle_copy((gretl_bundle *) ptr, &err);
	} else {
	    item->data = ptr;
	}
    } else if (item->type == GRETL_TYPE_ARRAY) {
	gretl_array_destroy((gretl_array*) item->data);
	if (copy) {
	    item->data = gretl_array_copy((gretl_array *) ptr, &err);
	} else {
	    item->data = ptr;
	}	
    } else {
	return E_DATA;
    }

    if (!err && item->data == NULL) {
	err = E_ALLOC;
    }

    if (item->note != NULL) {
	free(item->note);
	item->note = NULL;
    }

    return err;
}

/* callback invoked when a bundle's hash table is destroyed */

static void bundled_item_destroy (gpointer data)
{
    bundled_item *item = data;

#if BDEBUG
    fprintf(stderr, "bundled_item_destroy: type %d\n", item->type);
    if (item->type == GRETL_TYPE_STRING) {
	fprintf(stderr, " string: '%s'\n", (char *) item->data);
    }
#endif

    switch (item->type) {
    case GRETL_TYPE_DOUBLE:
    case GRETL_TYPE_STRING:
    case GRETL_TYPE_SERIES:
	free(item->data);
	break;
    case GRETL_TYPE_MATRIX:
	gretl_matrix_free((gretl_matrix *) item->data);
	break;
    case GRETL_TYPE_MATRIX_REF:
	release_matrix_pointer((gretl_matrix **) &item->data);
	break;
    case GRETL_TYPE_BUNDLE:
	gretl_bundle_destroy((gretl_bundle *) item->data);
	break;
    case GRETL_TYPE_ARRAY:
	gretl_array_destroy((gretl_array*) item->data);
	break;
    default:
	break;
    }

    free(item->note);
    free(item);
}

static void bundle_key_destroy (gpointer data)
{
#if BDEBUG
    fprintf(stderr, "freeing key '%s'\n", (char *) data);
#endif
    free(data);
}

/**
 * gretl_bundle_destroy:
 * @bundle: bundle to destroy.
 *
 * Frees all contents of @bundle as well as the pointer itself.
 */

void gretl_bundle_destroy (gretl_bundle *bundle)
{
    if (bundle != NULL) {
	if (bundle->ht != NULL) {
	    g_hash_table_destroy(bundle->ht);
	}
	free(bundle->creator);
	if (bundle->type == BUNDLE_KALMAN) {
	    kalman_free(bundle->data);
	}
	free(bundle);
    }
}

/**
 * gretl_bundle_void_content:
 * @bundle: target bundle.
 *
 * Frees all contents of @bundle.
 */

void gretl_bundle_void_content (gretl_bundle *bundle)
{
    if (bundle == NULL) {
	return;
    }

    if (bundle->creator != NULL) {
	free(bundle->creator);
	bundle->creator = NULL;
    }
    
    if (bundle->ht != NULL && g_hash_table_size(bundle->ht) > 0) {
	g_hash_table_destroy(bundle->ht);
	bundle->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
					   bundle_key_destroy, 
					   bundled_item_destroy);
    }

    if (bundle->type == BUNDLE_KALMAN) {
	kalman_free(bundle->data);
	bundle->data = NULL;
	bundle->type = BUNDLE_PLAIN;
    }
}

/**
 * gretl_bundle_new:
 *
 * Returns: a newly allocated, empty gretl bundle.
 */

gretl_bundle *gretl_bundle_new (void)
{
    gretl_bundle *b = malloc(sizeof *b);

    if (b != NULL) {
	b->type = BUNDLE_PLAIN;
	b->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
				      bundle_key_destroy, 
				      bundled_item_destroy);
	b->creator = NULL;
	b->data = NULL;
    }

    return b;
}

/* Determine whether @name is the name of a saved bundle. */

int gretl_is_bundle (const char *name)
{
    return get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE) != NULL;
}

/**
 * get_bundle_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved bundle, if found, else NULL.
 */

gretl_bundle *get_bundle_by_name (const char *name)
{
    gretl_bundle *b = NULL;

    if (name != NULL && *name != '\0') {
	user_var *u = get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE);

	if (u != NULL) {
	    b = user_var_get_value(u);
	}
    }

    return b;
}

static int gretl_bundle_has_data (gretl_bundle *b, const char *key)
{
    gpointer p = g_hash_table_lookup(b->ht, key);

    return (p != NULL);
}

/**
 * gretl_bundle_get_data:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0), or NULL.
 * @err: location to receive error code.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 * If @err is non-NULL, its content is set to a non-zero
 * value if @bundle contains no item with key @key. If
 * the intent is simply to determine whether @bundle contains 
 * an item under the specified @key, @err should generally be 
 * left NULL.
 *
 * Note that the value returned is the actual data pointer from
 * within the bundle, not a copy of the data; so the pointer
 * must not be freed, and in general its content should not
 * be modified.
 */

void *gretl_bundle_get_data (gretl_bundle *bundle, const char *key,
			     GretlType *type, int *size, int *err)
{
    void *ret = NULL;
    int reserved = 0;
    int myerr = 0;

    if (bundle == NULL) {
	myerr = E_DATA;
	goto finish;
    }

    if (bundle->type == BUNDLE_KALMAN) {
	ret = maybe_retrieve_kalman_element(bundle->data, key,
					    type, &reserved,
					    &myerr);
    }

    if (!myerr && ret == NULL && !reserved) {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;

	    ret = item->data;
	    if (type != NULL) {
		*type = item->type;
	    }
	    if (size != NULL) {
		*size = item->size;
	    }
	} else {
	    if (err != NULL) {
		gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	    }
	    myerr = E_DATA;
	}
    }

 finish:

    if (err != NULL) {
	*err = myerr;
    }

    return ret;
}

/**
 * gretl_bundle_steal_data:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0), or NULL.
 * @err: location to receive error code, or NULL.
 *
 * Works like gretl_bundle_get_data() except that the data
 * pointer in question is removed from @bundle before it is
 * returned to the caller; so in effect the caller assumes
 * ownership of the item.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 */

void *gretl_bundle_steal_data (gretl_bundle *bundle, const char *key,
			       GretlType *type, int *size, int *err)
{
    void *ret = NULL;

    if (bundle == NULL) {
	if (err != NULL) {
	    *err = E_DATA;
	}
    } else {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    GList *keys = g_hash_table_get_keys(bundle->ht);
	    bundled_item *item = p;
	    gchar *keycpy = NULL;

	    ret = item->data;
	    if (type != NULL) {
		*type = item->type;
	    }
	    if (size != NULL) {
		*size = item->size;
	    }
	    while (keys) {
		if (!strcmp(keys->data, key)) {
		    keycpy = keys->data;
		    break;
		} else if (keys->next) {
		    keys = keys->next;
		} else {
		    break;
		}
	    }
	    g_hash_table_steal(bundle->ht, key);
	    g_free(keycpy);
	    g_list_free(keys);
	    free(item);
	} else if (err != NULL) {
	    gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	    *err = E_DATA;
	}	    
    }

    return ret;
}

void *gretl_bundle_get_private_data (gretl_bundle *bundle)
{
    return bundle->data;
}

BundleType gretl_bundle_get_type (gretl_bundle *bundle)
{
    return bundle->type;
}

/**
 * gretl_bundle_get_member_type:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err:location to receive error code.
 *
 * Returns: the data type associated with @key in the
 * specified @bundle, or 0 on failure.
 */

GretlType gretl_bundle_get_member_type (gretl_bundle *bundle,
					const char *key,
					int *err)
{
    GretlType ret = GRETL_TYPE_NONE;
    int reserved = 0;

    if (bundle == NULL) {
	*err = E_DATA;
	return GRETL_TYPE_NONE;
    }

    if (bundle->type == BUNDLE_KALMAN) {
	maybe_retrieve_kalman_element(bundle->data, key,
				      &ret, &reserved,
				      err);
    }
    
    if (!*err && ret == GRETL_TYPE_NONE && !reserved) {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;
	    
	    ret = item->type;
	} else {
	    gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	    *err = E_DATA;
	}
				 
    }

    return ret;
}

/**
 * gretl_bundle_get_matrix:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the matrix associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

gretl_matrix *gretl_bundle_get_matrix (gretl_bundle *bundle,
				       const char *key,
				       int *err)
{
    gretl_matrix *m = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_MATRIX && 
	type != GRETL_TYPE_MATRIX_REF) {
	myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
	m = (gretl_matrix *) ptr;
    }

    if (err != NULL) {
	*err = myerr;
    }

    return m;
}

/**
 * gretl_bundle_get_array:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the array associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

void *gretl_bundle_get_array (gretl_bundle *bundle,
			      const char *key,
			      int *err)
{
    gretl_array *a = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_ARRAY) { 
	myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
	a = (gretl_array *) ptr;
    }

    if (err != NULL) {
	*err = myerr;
    }

    return a;
}

/**
 * gretl_bundle_get_series:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @n: location to receive length of series.
 * @err: location to receive error code.
 *
 * Returns: the series associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

double *gretl_bundle_get_series (gretl_bundle *bundle,
				 const char *key,
				 int *n, int *err)
{
    double *x = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, n, err);
    if (ptr != NULL && type != GRETL_TYPE_SERIES) {
	myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
	x = (double *) ptr;
    }

    if (err != NULL) {
	*err = myerr;
    }    

    return x;
}

/**
 * gretl_bundle_get_scalar:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the scalar value associated with @key in the
 * specified @bundle, if any; otherwise #NADBL.
 */

double gretl_bundle_get_scalar (gretl_bundle *bundle,
				const char *key,
				int *err)
{
    double x = NADBL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_DOUBLE) {
	myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
	double *px = (double *) ptr;

	x = *px;
    }

    if (err != NULL) {
	*err = myerr;
    }

    return x;
}

/**
 * gretl_bundle_get_string:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the string value associated with @key in the
 * specified @bundle, if any; otherwise #NADBL.
 */

const char *gretl_bundle_get_string (gretl_bundle *bundle,
				     const char *key,
				     int *err)
{
    const char *ret = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_STRING) {
	myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
	ret = (const char *) ptr;
    }

    if (err != NULL) {
	*err = myerr;
    }    

    return ret;
}

/**
 * gretl_bundle_get_note:
 * @bundle: bundle to access.
 * @key: name of key to access.
 *
 * Returns: the note associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const char *gretl_bundle_get_note (gretl_bundle *bundle, 
				   const char *key)
{
    const char *ret = NULL;

    if (bundle != NULL) {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;

	    ret = item->note;
	}
    }

    return ret;
}

/**
 * gretl_bundle_get_creator:
 * @bundle: bundle to access.
 *
 * Returns: the name of the package that created @bundle, if any, 
 * otherwise NULL.
 */

const char *gretl_bundle_get_creator (gretl_bundle *bundle)
{
    return (bundle != NULL)? bundle->creator : NULL;
}

/**
 * bundled_item_get_data:
 * @item: bundled item to access.
 * @type: location to receive data type.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0).
 *
 * Returns: the data pointer associated with @item, or
 * NULL on failure.
 */

void *bundled_item_get_data (bundled_item *item, GretlType *type,
			     int *size)
{
    *type = item->type;

    if (size != NULL) {
	*size = item->size;
    }

    return item->data;
}

/**
 * bundled_item_get_note:
 * @item: bundled item.
 *
 * Returns: the note associated with @item (may be NULL).
 */

const char *bundled_item_get_note (bundled_item *item)
{
    return item->note;
}

/**
 * gretl_bundle_get_content:
 * @bundle: bundle to access.
 *
 * Returns: the content of @bundle, which is in fact
 * a GHashTable object.
 */

void *gretl_bundle_get_content (gretl_bundle *bundle)
{
    return (bundle == NULL)? NULL : (void *) bundle->ht;
}

static int real_bundle_set_data (gretl_bundle *b, const char *key,
				 void *ptr, GretlType type,
				 int size, int copy, 
				 const char *note)
{
    int err, done = 0;

    err = strlen(key) >= VNAMELEN ? E_DATA : 0;

    if (err) {
	gretl_errmsg_sprintf("'%s': invalid key string", key);
	return err;
    }

    if (b->type == BUNDLE_KALMAN) {
	done = maybe_set_kalman_element(b->data, key,
					ptr, type, copy,
					&err);
    }

    if (!done && !err) {
	bundled_item *item = g_hash_table_lookup(b->ht, key);
	int replace = 0;
	
	if (item != NULL) {
	    replace = 1;
	    if (item->type == type) {
		return bundled_item_replace_data(item, type, ptr, 
						 size, copy);
	    }
	}
	
	item = bundled_item_new(type, ptr, size, copy, note, &err);

	if (!err) {
	    gchar *k = g_strdup(key);

	    if (replace) {
		g_hash_table_replace(b->ht, k, item);
	    } else {
		g_hash_table_insert(b->ht, k, item);
	    }
	}
    }

    return err;
}

/**
 * gretl_bundle_donate_data:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: type of data.
 * @size: if @type == GRETL_TYPE_SERIES, the length of
 * the series, otherwise 0.
 * 
 * Sets the data type and pointer to be associated with @key in 
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed. The value of @ptr is transcribed into the
 * bundle, which therefore "takes ownership" of the data;
 * compare gretl_bundle_set_data().
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_donate_data (gretl_bundle *bundle, const char *key,
			      void *ptr, GretlType type, int size)
{
    return real_bundle_set_data(bundle, key, ptr, type, size, 0, NULL);
}

/**
 * gretl_bundle_set_data:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: type of data.
 * @size: if @type == GRETL_TYPE_SERIES, the length of
 * the series, otherwise 0.
 * 
 * Sets the data type and pointer to be associated with @key in 
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed. The content of @ptr is copied into the
 * bundle; compare gretl_bundle_donate_data().
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_data (gretl_bundle *bundle, const char *key,
			   void *ptr, GretlType type, int size)
{
    int err;

    if (bundle == NULL) {
	err = E_UNKVAR;
    } else {
	err = real_bundle_set_data(bundle, key, ptr, type, 
				   size, 1, NULL);
    }

    return err;
}

/**
 * gretl_bundle_set_string:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @str: the string to set.
 * 
 * Sets @str as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_string (gretl_bundle *bundle, const char *key,
			     const char *str)
{
    return gretl_bundle_set_data(bundle, key, (void *) str, 
				 GRETL_TYPE_STRING, 0);
}

/**
 * gretl_bundle_set_scalar:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the value to set.
 * 
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_scalar (gretl_bundle *bundle, const char *key,
			     double val)
{
    return gretl_bundle_set_data(bundle, key, &val, 
				 GRETL_TYPE_DOUBLE, 0);
}

/**
 * gretl_bundle_set_series:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @x: array of doubles.
 * @n: the length of @x.
 * 
 * Sets @x as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_series (gretl_bundle *bundle, const char *key,
			     const double *x, int n)
{
    return gretl_bundle_set_data(bundle, key, (void *) x, 
				 GRETL_TYPE_SERIES, n);
}

/**
 * gretl_bundle_set_matrix:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @m: gretl matrix.
 * 
 * Sets @m as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_matrix (gretl_bundle *bundle, const char *key,
			     const gretl_matrix *m)
{
    return gretl_bundle_set_data(bundle, key, (void *) m, 
				 GRETL_TYPE_MATRIX, 0);
}

/**
 * gretl_bundle_delete_data:
 * @bundle: target bundle.
 * @key: name of key to delete.
 * 
 * Deletes the data item under @key from @bundle, if
 * such an item is present.
 */

int gretl_bundle_delete_data (gretl_bundle *bundle, const char *key)
{
    int done = 0;
    int err = 0;

    if (bundle == NULL) {
	return E_DATA;
    }

    if (bundle->type == BUNDLE_KALMAN) {
	done = maybe_delete_kalman_element(bundle->data, key, &err);
    }

    if (!done && !err) {
	done = g_hash_table_remove(bundle->ht, key);
	if (!done) {
	    err = E_DATA;
	}
    }
    
    return err;
}

/**
 * gretl_bundle_set_note:
 * @bundle: target bundle.
 * @key: name of key to access.
 * @note: note to add.
 * 
 * Adds a descriptive note to the item under @key in @bundle.
 * If a note is already present it is replaced by the new one.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_note (gretl_bundle *bundle, const char *key,
			   const char *note)
{
    int err = 0;

    if (bundle == NULL) {
	err = E_UNKVAR;
    } else {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p == NULL) {
	    err = E_DATA;
	} else {
	    bundled_item *item = p;

	    free(item->note);
	    item->note = gretl_strdup(note);
	}
    }

    return err;
}

/**
 * gretl_bundle_set_creator:
 * @bundle: target bundle.
 * @name: name of function package that built the bundle.
 * 
 * Sets the "creator" attribute of @bundle. This is called 
 * automatically when a bundle is returned to top-level 
 * userspace by a packaged function.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_creator (gretl_bundle *bundle, const char *name)
{
    int err = 0;

    if (bundle == NULL) {
	err = E_DATA;
    } else {
	free(bundle->creator);
	if (name == NULL) {
	    bundle->creator = NULL;
	} else {
	    bundle->creator = gretl_strdup(name);
	    if (bundle->creator == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

/* replicate on a target bundle a bundled_item from some other
   other bundle, provided the target bundle does not already
   have a bundled_item under the same key
*/

static void copy_new_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = (bundled_item *) value;
    gretl_bundle *targ = (gretl_bundle *) p;

    if (!gretl_bundle_has_data(targ, (const char *) key)) {
	real_bundle_set_data(targ, (const char *) key,
			     item->data, item->type,
			     item->size, 1, item->note);
    }
}

/* replicate on a target bundle a bundled_item from some other
   bundle */

static void copy_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = (bundled_item *) value;
    gretl_bundle *targ = (gretl_bundle *) p;

    real_bundle_set_data(targ, (const char *) key,
			 item->data, item->type,
			 item->size, 1, item->note);
}

/* Create a new bundle as the union of two existing bundles:
   we first copy bundle1 in its entirety, then append any elements
   of bundle2 whose keys that are not already present in the
   copy-target. In case bundle1 and bundle2 share any keys, the
   value copied to the target is therefore that from bundle1.
*/

gretl_bundle *gretl_bundle_union (const gretl_bundle *bundle1,
				  const gretl_bundle *bundle2,
				  int *err)
{
    gretl_bundle *b = NULL;

    if (bundle2->type == BUNDLE_KALMAN) {
	gretl_errmsg_set("bundle union: the right-hand operand cannot "
			 "be a kalman bundle");
	*err = E_DATA;
    } else {
	b = gretl_bundle_copy(bundle1, err);
    }

    if (!*err) {
	g_hash_table_foreach(bundle2->ht, copy_new_bundled_item, b);
    }
    
    return b;
}

/**
 * gretl_bundle_copy:
 * @bundle: gretl bundle to be copied.
 * @err: location to receive error code.
 *
 * Returns: a "deep copy" of @bundle (all the items in @bundle
 * are themselves copied), or NULL on failure.
 */

gretl_bundle *gretl_bundle_copy (const gretl_bundle *bundle, int *err)
{
    gretl_bundle *bcpy = NULL;

    if (bundle == NULL) {
	*err = E_DATA;
    } else {
	if (bundle->type == BUNDLE_KALMAN) {
	    bcpy = kalman_bundle_copy(bundle, err);
	} else {
	    bcpy = gretl_bundle_new();
	    if (bcpy == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (!*err) {
	    g_hash_table_foreach(bundle->ht, copy_bundled_item, bcpy);
	}
    }

    return bcpy;
}

/**
 * gretl_bundle_copy_as:
 * @name: name of source bundle.
 * @copyname: name for copy.
 *
 * Look for a saved bundle specified by @name, and if found,
 * make a full copy and save it under @copyname. This is
 * called from geneval.c on completion of assignment to a
 * bundle named @copyname, where the returned value on the
 * right-hand side is a pre-existing saved bundle.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_bundle_copy_as (const char *name, const char *copyname)
{
    gretl_bundle *b0, *b1 = NULL;
    user_var *u;
    int prev = 0;
    int err = 0;

    if (!strcmp(name, "$sysinfo")) {
	b0 = sysinfo_bundle;
    } else {
	u = get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE);
	if (u == NULL) {
	    return E_DATA;
	} else {
	    b0 = user_var_get_value(u);
	}
    }

    u = get_user_var_of_type_by_name(copyname, GRETL_TYPE_BUNDLE);
    
    if (u != NULL) {
	b1 = user_var_steal_value(u);
	if (b1 != NULL) {
	    gretl_bundle_destroy(b1);
	}
	prev = 1;
    }

    b1 = gretl_bundle_copy(b0, &err);

    if (!err) {
	if (prev) {
	    err = user_var_replace_value(u, b1);
	} else {
	    err = user_var_add(copyname, GRETL_TYPE_BUNDLE, b1);
	}
    }

    return err;
}

static void print_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = value;
    const gchar *kstr = key;
    gretl_array *a;
    gretl_matrix *m;
    double x;
    char *s;
    PRN *prn = p;

    switch (item->type) {
    case GRETL_TYPE_DOUBLE:
	x = *(double *) item->data;
	if (na(x)) {
	    pprintf(prn, " %s = NA", kstr);
	} else {
	    pprintf(prn, " %s = %g", kstr, x);
	}	
	break;
    case GRETL_TYPE_STRING:
	s = (char *) item->data;
	if (strlen(s) < 64) {
	    pprintf(prn, " %s = %s", kstr, s);
	} else {
	    pprintf(prn, " %s (%s)", kstr, 
		    gretl_type_get_name(item->type));
	}	
	break;
    case GRETL_TYPE_BUNDLE:
	pprintf(prn, " %s (%s)", kstr, 
		gretl_type_get_name(item->type));
	break;
    case GRETL_TYPE_MATRIX:
    case GRETL_TYPE_MATRIX_REF:
	m = item->data;
	if (m->rows == 1 && m->cols == 1) {
	    pprintf(prn, " %s = %g", kstr, m->val[0]);
	} else {
	    pprintf(prn, " %s (%s: %d x %d)", kstr, 
		    gretl_type_get_name(item->type), 
		    m->rows, m->cols);
	}
	break;
    case GRETL_TYPE_SERIES:
	pprintf(prn, " %s (%s: length %d)", kstr, 
		gretl_type_get_name(item->type), item->size);
	break;
    case GRETL_TYPE_ARRAY:
	a = item->data;
	{
	    GretlType t = gretl_array_get_type(a);
	    int n = gretl_array_get_length(a);

	    pprintf(prn, " %s = array of %s, length %d", kstr, 
		    gretl_type_get_name(t), n);
	}
	break;
    default:
	break;
    }

    if (item->note != NULL) {
	pprintf(prn, " %s\n", item->note);
    } else {
	pputc(prn, '\n');
    }
}

/**
 * gretl_bundle_print:
 * @bundle: gretl bundle.
 * @prn: gretl printer.
 *
 * Prints to @prn a list of the keys defined in @bundle, along
 * with descriptive notes, if any.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int gretl_bundle_print (gretl_bundle *bundle, PRN *prn)
{
    if (bundle == NULL) {
	return E_DATA;
    } else {
	int n_items = g_hash_table_size(bundle->ht);
	user_var *u = get_user_var_by_data(bundle);
	const char *name = NULL;

	if (u != NULL) {
	    name = user_var_get_name(u);
	} else {
	    name = "anonymous";
	}

	if (bundle->type == BUNDLE_PLAIN && n_items == 0) {
	    pprintf(prn, "bundle %s: empty\n", name);
	} else {
	    if (bundle->creator != NULL) {
		pprintf(prn, "bundle %s, created by %s:\n", 
			name, bundle->creator);
	    } else {
		pprintf(prn, "bundle %s:\n", name);
	    }
	    if (bundle->type == BUNDLE_KALMAN) {
		print_kalman_bundle_info(bundle->data, prn);
		if (n_items > 0) {
		    pputs(prn, "\nOther content\n");
		    g_hash_table_foreach(bundle->ht, print_bundled_item, prn);
		}
	    } else if (n_items > 0) {
		g_hash_table_foreach(bundle->ht, print_bundled_item, prn);
	    }
	    pputc(prn, '\n');
	}
	
	return 0;
    }
}

static gboolean match_by_data (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = value;

    return item->data == p;
}

int bundle_contains_data (gretl_bundle *b, void *data)
{
    return g_hash_table_find(b->ht, match_by_data, data) != NULL;
}

/* Called from gretl_func.c on return, to remove
   a given bundle from the stack of named bundles in
   preparation for handing it over to the caller,
   who will take ownership of it.
*/

gretl_bundle *gretl_bundle_pull_from_stack (const char *name,
					    int *err)
{
    gretl_bundle *b = NULL;
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE);

    if (u != NULL) {
	b = user_var_unstack_value(u);
    }

    if (b == NULL) {
	*err = E_DATA;
    } 

    return b;
}

/* serialize a gretl bundled item as XML */

static void xml_put_bundled_item (gpointer keyp, gpointer value, gpointer p)
{
    const char *key = keyp;
    bundled_item *item = value;
    FILE *fp = p;
    int j;

    if (item->type == GRETL_TYPE_STRING) {
	char *s = item->data;

	if (s == NULL || *s == '\0') {
	    fprintf(stderr, "bundle -> XML: skipping empty string %s\n", key);
	    return;
	} 
    }

    fprintf(fp, "<bundled-item key=\"%s\" type=\"%s\"", key,
	    gretl_type_get_name(item->type));

    if (item->note != NULL) {
	fprintf(fp, " note=\"%s\"", item->note);
    } 

    if (item->size > 0) {
	fprintf(fp, " size=\"%d\">\n", item->size);
    } else if (item->type == GRETL_TYPE_STRING) {
	fputc('>', fp);
    } else {
	fputs(">\n", fp);
    }

    if (item->type == GRETL_TYPE_DOUBLE) {
	double x = *(double *) item->data;
	
	if (na(x)) {
	    fputs("NA", fp);
	} else {
	    fprintf(fp, "%.16g", x);
	}
    } else if (item->type == GRETL_TYPE_SERIES) {
	double *vals = (double *) item->data;

	for (j=0; j<item->size; j++) {
	    if (na(vals[j])) {
		fputs("NA ", fp);
	    } else {
		fprintf(fp, "%.16g ", vals[j]);
	    }
	}	    
    } else if (item->type == GRETL_TYPE_STRING) {
	gretl_xml_put_string((char *) item->data, fp);
    } else if (item->type == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = (gretl_matrix *) item->data;

	gretl_matrix_serialize(m, NULL, fp);
    } else if (item->type == GRETL_TYPE_BUNDLE) {
	gretl_bundle *b = (gretl_bundle *) item->data;

	gretl_bundle_serialize(b, NULL, fp);
    } else if (item->type == GRETL_TYPE_ARRAY) {
	gretl_array *a = (gretl_array *) item->data;

	gretl_array_serialize(a, fp);
    } else {
	fprintf(stderr, "bundle -> XML: skipping %s\n", key);
    }

    fputs("</bundled-item>\n", fp);    
}

void gretl_bundle_serialize (gretl_bundle *b, const char *name, 
			     FILE *fp)
{
    fputs("<gretl-bundle", fp);
    if (name != NULL) {
	fprintf(fp, " name=\"%s\"", name);
    }
    if (b->creator != NULL && *b->creator != '\0') {
	fprintf(fp, " creator=\"%s\"", b->creator);
    }
    if (b->type == BUNDLE_KALMAN) {
	fputs(" type=\"kalman\"", fp);
    }
    fputs(">\n", fp);

    if (b->type == BUNDLE_KALMAN) {
	kalman_serialize(b->data, fp);
    }

    if (b->ht != NULL) {
	g_hash_table_foreach(b->ht, xml_put_bundled_item, fp);
    }

    fputs("</gretl-bundle>\n", fp); 
}

static int load_bundled_items (gretl_bundle *b, xmlNodePtr cur, xmlDocPtr doc)
{
    GretlType type;
    char *key;
    int err = 0;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "bundled-item")) {
	    key = (char *) xmlGetProp(cur, (XUC) "key");
	    type = gretl_xml_get_type_property(cur);
	    if (key == NULL || type == 0) {
		err = E_DATA;
	    } else {
		int size = 0;

		if (type == GRETL_TYPE_DOUBLE) {
		    double x;

		    if (!gretl_xml_node_get_double(cur, doc, &x)) {
			err = E_DATA;
		    } else {
			err = gretl_bundle_set_data(b, key, &x, type, size);
		    }
		} else if (type == GRETL_TYPE_STRING) {
		    char *s;

		    if (!gretl_xml_node_get_trimmed_string(cur, doc, &s)) {
			err = E_DATA;
		    } else {
			err = gretl_bundle_donate_data(b, key, s, type, size);
		    }
		} else if (type == GRETL_TYPE_SERIES) {
		    double *xvec = gretl_xml_get_double_array(cur, doc, &size, &err);

		    if (!err) {
			err = gretl_bundle_donate_data(b, key, xvec, type, size);
		    }
		} else if (type == GRETL_TYPE_MATRIX) {
		    xmlNodePtr child = cur->xmlChildrenNode;
		    gretl_matrix *m;

		    if (child == NULL) {
			err = E_DATA;
		    } else {
			m = gretl_xml_get_matrix(child, doc, &err);
			if (!err) {
			    err = gretl_bundle_donate_data(b, key, m, type, size);
			}
		    }
		} else if (type == GRETL_TYPE_BUNDLE) {
		    xmlNodePtr child = cur->xmlChildrenNode;
		    gretl_bundle *baby;

		    if (child == NULL) {
			err = E_DATA;
		    } else {
			baby = gretl_bundle_deserialize(child, doc, &err);
			if (!err) {
			    err = gretl_bundle_donate_data(b, key, baby, type, size);
			}
		    }
		} else if (type == GRETL_TYPE_ARRAY) {
		    xmlNodePtr child = cur->xmlChildrenNode;
		    gretl_array *a;

		    if (child == NULL) {
			err = E_DATA;
		    } else {
			a = gretl_array_deserialize(child, doc, &err);
			if (!err) {
			    err = gretl_bundle_donate_data(b, key, a, type, size);
			}
		    }
		} else {
		    fprintf(stderr, "bundle: ignoring unhandled type %d\n", type);
		}

		if (!err) {
		    char *note = (char *) xmlGetProp(cur, (XUC) "note");

		    if (note != NULL) {
			gretl_bundle_set_note(b, key, note);
			xmlFree(note);
		    }
		}

		xmlFree(key);
	    }
	}
	cur = cur->next;
    }

    return err;
}

/* For internal use only: @p1 should be of type xmlNodePtr and @p2
   should be an xmlDocPtr. We suppress the actual pointer types in the
   prototype so that it's possible for a module to include
   gretl_bundle.h without including the full libxml headers.
*/

gretl_bundle *gretl_bundle_deserialize (void *p1, void *p2, 
					int *err)
{
    xmlNodePtr cur, node = p1;
    xmlDocPtr doc = p2;
    gretl_bundle *b = NULL;
    char *btype = NULL;

    btype = (char *) xmlGetProp(node, (XUC) "type");

    cur = node->xmlChildrenNode;

    if (btype != NULL && !strcmp(btype, "kalman")) {
	b = kalman_deserialize(cur, doc, err);
    } else {
	b = gretl_bundle_new();
	if (b == NULL) {
	    *err = E_ALLOC;
	}
    }

    free(btype);

    if (b != NULL) {
	*err = load_bundled_items(b, cur, doc);
	if (*err) {
	    fprintf(stderr, "deserialize bundle: "
		    "bundle is broken (err = %d)\n", *err);
	    gretl_bundle_destroy(b);
	    b = NULL;
	}
    }

    return b;
}

int gretl_bundle_write_to_file (gretl_bundle *b, 
				const char *fname,
				int to_dotdir)
{
    char fullname[FILENAME_MAX];
    FILE *fp;
    int err = 0;

    if (to_dotdir) {
	build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
	strcpy(fullname, fname);
	gretl_maybe_switch_dir(fullname);
    }

    fp = gretl_fopen(fullname, "wb");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	gretl_push_c_numeric_locale();
	gretl_xml_header(fp);
	gretl_bundle_serialize(b, NULL, fp);
	fclose(fp);
	gretl_pop_c_numeric_locale();
    }

    return err;
}

gretl_bundle *gretl_bundle_read_from_file (const char *fname, 
					   int from_dotdir,
					   int *err)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    char fullname[FILENAME_MAX];
    gretl_bundle *b = NULL;

    if (from_dotdir) {
	build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
	strcpy(fullname, fname);
	gretl_maybe_prepend_dir(fullname);
    }

    *err = gretl_xml_open_doc_root(fullname, "gretl-bundle", &doc, &cur);

    if (!*err) {
	gretl_push_c_numeric_locale();
	b = gretl_bundle_deserialize(cur, doc, err);
	gretl_pop_c_numeric_locale();
	xmlFreeDoc(doc);
    }

    if (*err && b != NULL) {
	gretl_bundle_destroy(b);
	b = NULL;
    }

    return b;
}

gretl_bundle *gretl_bundle_read_from_buffer (const char *buf, int len,
					     int *err)
{
    xmlDocPtr doc = NULL;
    gretl_bundle *b = NULL;

    xmlKeepBlanksDefault(0);
    doc = xmlParseMemory(buf, len);

    if (doc == NULL) {
	gretl_errmsg_set(_("xmlParseMemory failed"));
	*err = 1;
    } else {
	xmlNodePtr cur = xmlDocGetRootElement(doc);

	if (cur == NULL) {
	    gretl_errmsg_set(_("xmlDocGetRootElement failed"));
	    *err = 1;
	} else {
	    gretl_push_c_numeric_locale();
	    b = gretl_bundle_deserialize(cur, doc, err);
	    gretl_pop_c_numeric_locale();
	}
	xmlFreeDoc(doc);
    }

    if (*err && b != NULL) {
	gretl_bundle_destroy(b);
	b = NULL;
    }

    return b;
}

gretl_bundle *get_sysinfo_bundle (int *err)
{
    if (sysinfo_bundle == NULL) {
	gretl_bundle *b = gretl_bundle_new();

	if (b == NULL) {
	    *err = E_ALLOC;
	} else {
	    char *s1, *s2;
	    int ival = 0;

#if HAVE_MPI
	    ival = check_for_mpiexec();
#endif	    
	    gretl_bundle_set_scalar(b, "mpi", (double) ival);
	    ival = gretl_max_mpi_processes();
	    gretl_bundle_set_scalar(b, "mpimax", (double) ival);
	    ival = gretl_n_processors();
	    gretl_bundle_set_scalar(b, "nproc", (double) ival);
	    ival = 0;
#ifdef _OPENMP
	    ival = 1;
#endif
	    gretl_bundle_set_scalar(b, "omp", (double) ival);
	    ival = sizeof(void*) == 8 ? 64 : 32;
	    gretl_bundle_set_scalar(b, "wordlen", (double) ival);
	    gretl_bundle_set_scalar(b, "omp_num_threads", get_omp_n_threads());
#if defined(G_OS_WIN32)
	    gretl_bundle_set_string(b, "os", "windows");
#elif defined(OS_OSX) 
	    gretl_bundle_set_string(b, "os", "osx");
#elif defined(linux)
	    gretl_bundle_set_string(b, "os", "linux");
#else
	    gretl_bundle_set_string(b, "os", "other");
#endif
	    gretl_bundle_set_string(b, "hostname", g_get_host_name());
	    gretl_bundle_set_string(b, "blas", blas_variant_string());
	    if (get_openblas_details(&s1, &s2)) {
		gretl_bundle_set_string(b, "blascore", s1);
		gretl_bundle_set_string(b, "blas_parallel", s2);
	    }
	}
	sysinfo_bundle = b;
    }

    return sysinfo_bundle;
}

void *sysinfo_bundle_get_data (const char *key, GretlType *type, 
			       int *err)
{
    gretl_bundle *b = get_sysinfo_bundle(err);
    void *ret = NULL;

    if (b != NULL) {
	ret = gretl_bundle_get_data(b, key, type, NULL, err);
    }

    return ret;
}

static gretl_matrix *matrix_from_list (const int *list)
{
    gretl_matrix *m = NULL;
    int i;
    
    if (list != NULL && list[0] > 0) {
	m = gretl_matrix_alloc(1, list[0]);
	if (m != NULL) {
	    for (i=0; i<list[0]; i++) {
		m->val[i] = list[i+1];
	    }
	}
    }

    return m;
}

/* For a single-equation model, create a bundle containing
   all the data available via $-accessors.
*/

gretl_bundle *bundle_from_model (MODEL *pmod,
				 DATASET *dset,
				 int *err)
{
    gretl_bundle *b = NULL;
    gretl_matrix *m;
    double *x;
    double val;
    int *list;
    const char *s;
    const char *key;
    int i, t, berr;

    if (pmod == NULL) {
	/* get the "last model" */
	GretlObjType type = 0;
	void *p = get_last_model(&type);

	if (p == NULL || type != GRETL_OBJ_EQN) {
	    *err = E_DATA;
	    return NULL;
	} else {
	    pmod = p;
	}
    }

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    b = gretl_bundle_new();
    if (b == NULL) {
	free(x);
	*err = E_ALLOC;
	return NULL;
    }

    for (i=M_ESS; i<M_SCALAR_MAX && !*err; i++) {
	berr = 0;
	val = gretl_model_get_scalar(pmod, i, dset, &berr);
	if (!berr) {
	    key = mvarname(i) + 1;
	    *err = gretl_bundle_set_scalar(b, key, val);	    
	}
    }

    for (i=M_SCALAR_MAX+1; i<M_SERIES_MAX && !*err; i++) {
	for (t=0; t<dset->n; t++) {
	    x[t] = NADBL;
	}
	berr = gretl_model_get_series(x, pmod, dset, i);
	if (!berr) {
	    key = mvarname(i) + 1;
	    *err = gretl_bundle_set_series(b, key, x, dset->n);
	}	
    }    

    for (i=M_SERIES_MAX+1; i<M_MATRIX_MAX && !*err; i++) {
	berr = 0;
	m = gretl_model_get_matrix(pmod, i, &berr);
	if (!berr) {
	    key = mvarname(i) + 1;
	    *err = gretl_bundle_donate_data(b, key, m,
					    GRETL_TYPE_MATRIX,
					    0);
	}
    }

    for (i=M_MBUILD_MAX+1; i<M_LIST_MAX && !*err; i++) {
	list = NULL;
	if (i == M_XLIST) {
	    list = gretl_model_get_x_list(pmod);
	} else if (i == M_YLIST) {
	    list = gretl_model_get_y_list(pmod);
	}
	if (list != NULL) {
	    /* convert list to matrix for bundling */
	    m = matrix_from_list(list);
	    if (m != NULL) {
		key = mvarname(i) + 1;
		*err = gretl_bundle_donate_data(b, key, m,
						GRETL_TYPE_MATRIX,
						0);
	    }
	    free(list);
	}
    }

    for (i=M_LIST_MAX+1; i<M_MAX && !*err; i++) {
	s = NULL;
	if (i == M_DEPVAR) {
	    s = gretl_model_get_depvar_name(pmod, dset);
	} else if (i == M_COMMAND) {
	    s = gretl_command_word(pmod->ci);
	}
	if (s != NULL && *s != '\0') {
	    key = mvarname(i) + 1;
	    *err = gretl_bundle_set_string(b, key, s);	    
	}	    
    }

    if (!*err) {
	*err = bundlize_model_data_scalars(pmod, b);
    }

    free(x);

    /* don't return a broken bundle */
    if (*err && b != NULL) {
	gretl_bundle_destroy(b);
	b = NULL;
    }

    return b;
}

gretl_bundle *kalman_bundle_new (gretl_matrix *M[],
				 int copy[], int nmat,
				 int *err)
{
    gretl_bundle *b = gretl_bundle_new();

    if (b == NULL) {
	*err = E_ALLOC;
    } else {
	b->type = BUNDLE_KALMAN;
	b->data = kalman_new_minimal(M, copy, nmat, err);
    }

    /* don't return a broken bundle */
    if (*err && b != NULL) {
	gretl_bundle_destroy(b);
	b = NULL;
    }

    return b;
}

void gretl_bundle_cleanup (void)
{
    if (sysinfo_bundle != NULL) {
	gretl_bundle_destroy(sysinfo_bundle);
	sysinfo_bundle = NULL;
    }
}


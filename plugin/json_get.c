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

/* parsing of JSON buffer using the json-glib library */

#include "libgretl.h"
#include "version.h"

#include <glib-object.h>
#include <json-glib/json-glib.h>

#define handled_type(t) (t == G_TYPE_STRING || \
			 t == G_TYPE_DOUBLE || \
			 t == G_TYPE_INT64)

static int real_json_get (JsonParser *parser, const char *pathstr,
			  int *n_objects, PRN *prn)
{
    GError *gerr = NULL;
    JsonNode *match, *node;
    JsonPath *path;
    GType ntype;
    double x;
    int err = 0;

    *n_objects = 0;

    node = json_parser_get_root(parser);
    path = json_path_new();

    if (!json_path_compile(path, pathstr, &gerr)) {
	if (gerr != NULL) {
	    gretl_errmsg_sprintf("Failed to compile JsonPath: %s",
				 gerr->message);
	    g_error_free(gerr);
	} else {
	    gretl_errmsg_set("Failed to compile JsonPath");
	}	    
	g_object_unref(path);
	return E_DATA;
    }

    match = json_path_match(path, node);
    if (match == NULL) {
	/* FIXME : maybe return empty string? */
	g_object_unref(path);
	return E_DATA;
    }

    /* in case we get floating-point output */
    gretl_push_c_numeric_locale();

    if (JSON_NODE_HOLDS_ARRAY(match)) {
	JsonArray *array;

	array = json_node_get_array(match);
	node = json_array_get_element(array, 0);

    repeat:

	if (node == NULL) {
	    gretl_errmsg_set("Failed to match JsonPath");
	    ntype = 0;
	} else {
	    ntype = json_node_get_value_type(node);
	}

	if (!handled_type(ntype)) {
	    if (JSON_NODE_HOLDS_ARRAY(node)) {
		array = json_node_get_array(node);
		node = json_array_get_element(array, 0);
		goto repeat;
	    } else {
		gretl_errmsg_sprintf("Unhandled array type '%s'", 
				     g_type_name(ntype));
		err = E_DATA;
	    }
	} else {
	    int i, n = json_array_get_length(array);

	    for (i=0; i<n; i++) {
		node = json_array_get_element(array, i);
		if (ntype == G_TYPE_STRING) {
		    pputs(prn, json_node_get_string(node));
		} else {
		    x = json_node_get_double(node);
		    pprintf(prn, "%.15g", x);
		}
		if (n > 1) {
		    pputc(prn, '\n');
		}
	    }
	    *n_objects = n;
	}
    } else {
	ntype = json_node_get_value_type(match);
	if (!handled_type(ntype)) {
	    gretl_errmsg_sprintf("Unhandled object type '%s'", 
				 g_type_name(ntype));
	    err = E_DATA;
	} else {
	    if (ntype == G_TYPE_STRING) {
		pputs(prn, json_node_get_string(match));
	    } else {
		x = json_node_get_double(match);
		pprintf(prn, "%.15g", x);
	    }
	    *n_objects = 1;
	}
    }

    gretl_pop_c_numeric_locale();

    json_node_free(match);
    g_object_unref(path);

    return err;
}

/*
  @data: JSON buffer.
  @path: the JsonPath to the target info.
  @n_objects: location to receive the number of pieces
  of information retrieved, or NULL.
  @err: location to receive error code.

  On success, returns an allocated string. If the "target"
  is an array, the members are printed one per line. This
  function handles target types of double, int or string;
  in the case of doubles or ints, their string representation
  is returned (using the C locale for doubles).
*/

char *json_get (const char *data, const char *path, int *n_objects,
		int *err)
{
    GError *gerr = NULL;
    JsonParser *parser;
    char *ret = NULL;
    int n = 0;

    parser = json_parser_new();
    if (parser == NULL) {
	gretl_errmsg_set("json_parser_new returned NULL!\n");
	*err = 1;
	return NULL;
    }

    json_parser_load_from_data(parser, data, -1, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_sprintf("Couldn't parse JSON input: %s",
			     gerr->message);
	g_error_free(gerr);
	*err = E_DATA;
    } else {
	PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER, err);

	if (!*err) {
	    *err = real_json_get(parser, path, &n, prn);
	    if (!*err) {
		ret = gretl_print_steal_buffer(prn);
	    }
	    gretl_print_destroy(prn);
	}
    }

    if (*err) {
	fprintf(stderr, "json_get: err = %d\n", *err);
    }

    if (n_objects != NULL) {
	*n_objects = n;
    }    

    g_object_unref(parser);

    return ret;
}

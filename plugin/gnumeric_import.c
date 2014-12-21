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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "version.h"
#include "gretl_xml.h"
#include "importer.h"
#include "csvdata.h"

#include <gtk/gtk.h>

#define IDEBUG 0

/* from gnumeric's value.h */
typedef enum {
    VALUE_EMPTY     = 10,
    VALUE_BOOLEAN   = 20, 
    VALUE_INTEGER   = 30,
    VALUE_FLOAT     = 40,
    VALUE_ERROR     = 50,
    VALUE_STRING    = 60,
    VALUE_CELLRANGE = 70,
    VALUE_ARRAY     = 80
} ValueType;

#include "import_common.c"

static void wsheet_init (wsheet *sheet)
{
    sheet->col_offset = sheet->row_offset = 0;
    sheet->maxcol = sheet->maxrow = 0;
    sheet->text_cols = sheet->text_rows = 0;
    sheet->colheads = 0;
    sheet->ID = 0;
    sheet->flags = 0;
    sheet->name = NULL;
    sheet->Z = NULL;
    sheet->varname = NULL;
    sheet->label = NULL;
}

static void wsheet_free (wsheet *sheet)
{
    int rows = sheet->maxrow + 1 - sheet->row_offset;
    int cols = sheet->maxcol + 1 - sheet->col_offset;
    int i;

    for (i=0; i<cols; i++) {
	if (sheet->varname != NULL) {
	    free(sheet->varname[i]);
	}
	if (sheet->Z != NULL) {
	    free(sheet->Z[i]);
	}
    }

    free(sheet->varname);
    free(sheet->Z);

    if (sheet->label != NULL) { 
	for (i=0; i<rows; i++) {
	    free(sheet->label[i]);
	}
	free(sheet->label);
    }

    free(sheet->name);

    wsheet_init(sheet);
}

static void wsheet_print_info (wsheet *sheet)
{
    int cols = sheet->maxcol + 1 - sheet->col_offset;
    int i;
#if IDEBUG
    int rows = sheet->maxrow + 1 - sheet->row_offset;
    int t;
#endif

    fputs("*** wsheet info after reading cells ***\n", stderr);
    fprintf(stderr, " maxcol = %d\n", sheet->maxcol);
    fprintf(stderr, " maxrow = %d\n", sheet->maxrow);
    fprintf(stderr, " text_cols = %d\n", sheet->text_cols);
    fprintf(stderr, " text rows = %d\n", sheet->text_rows);
    fprintf(stderr, " col_offset = %d\n", sheet->col_offset);
    fprintf(stderr, " row_offset = %d\n", sheet->row_offset);
    fputs(" varnames?\n", stderr);

    for (i=0; i<cols; i++) {
	fprintf(stderr, "  %d: '%s'\n", i, sheet->varname[i]);
    }	

#if IDEBUG
    fputs(" observations?\n", stderr);
    for (t=0; t<rows; t++) {
	fprintf(stderr, "  %d: ", t);
	if (sheet->text_cols) {
	    fprintf(stderr, "label='%s', ", sheet->label[t]);
	}
	for (i=0; i<cols; i++) {
	    if (na(sheet->Z[i][t])) {
		fprintf(stderr, "NA%c", (i == cols - 1)? '\n' : ' ');
	    } else {
		fprintf(stderr, "%g%c", sheet->Z[i][t],
			(i == cols - 1)? '\n' : ' ');
	    }
	}
    }
#endif

    fputs("*** end wsheet info ***\n", stderr);
}

#define VTYPE_IS_NUMERIC(v) ((v) == VALUE_BOOLEAN || \
                             (v) == VALUE_INTEGER || \
                             (v) == VALUE_FLOAT)

static int wsheet_allocate (wsheet *sheet, int cols, int rows)
{
    int i, t;

#if 1
    fprintf(stderr, "wsheet_allocate: allocating %d variables, each %d obs\n",
	    cols, rows);
#endif

    sheet->Z = doubles_array_new(cols, rows);
    if (sheet->Z == NULL) {
	return 1;
    }

    for (i=0; i<cols; i++) {
	for (t=0; t<rows; t++) {
	    sheet->Z[i][t] = NADBL;
	}
    }

    sheet->varname = strings_array_new_with_length(cols, VNAMELEN);
    if (sheet->varname == NULL) {
	return 1;
    }

    sheet->label = strings_array_new_with_length(rows, VNAMELEN);
    if (sheet->label == NULL) {
	return 1;
    }

    return 0;
}

static void check_for_date_format (wsheet *sheet, const char *fmt)
{
#if 0
    fprintf(stderr, "check_for_date_format: fmt = '%s'\n", fmt);
#endif

    if (strchr(fmt, '/') || 
	(strstr(fmt, "mm") && !(strchr(fmt, ':'))) || 
	strstr(fmt, "yy")) {
	book_set_numeric_dates(sheet);
    }
}

static int stray_numeric (int vtype, char *tmp, double *x)
{
    if (vtype == VALUE_STRING) {
	if (string_is_blank(tmp)) {
	    *x = NADBL;
	    return 1;
	} else if (import_na_string(tmp)) {
	    *x = NADBL;
	    return 1;
	} else if (numeric_string(tmp)) {
	    *x = atof(tmp);
	    return 1;
	}
    }

    return 0;
}

static int node_get_vtype_and_content (xmlNodePtr p, int *vtype,
				       char **content)
{
    char *tmp;
    int err = 0;

    tmp = (char *) xmlGetProp(p, (XUC) "ValueType");

    if (tmp != NULL) {
	*vtype = atoi(tmp);
	free(tmp);
	*content = (char *) xmlNodeGetContent(p);
    } else { 
	err = E_DATA;
    }

    return err;
}

static int inspect_top_left (xmlNodePtr p, int *obscol)
{
    char *content = NULL;
    int err, vtype = 0;

    err = node_get_vtype_and_content(p, &vtype, &content);

    if (!err) {
	if (vtype == VALUE_EMPTY) {
	    *obscol = 1;
	} else if (vtype == VALUE_STRING) {
	    if (import_obs_label(content)) {
		*obscol = 1;
	    }
	}
    }
	
    free(content);

    return err;
}

/* Crawl over all the cells and determine the maximum row and column
   indices. While we're at it, inspect the top left cell.
*/

static int wsheet_get_real_size_etc (xmlNodePtr node, wsheet *sheet,
				     int *obscol)
{
    xmlNodePtr p = node->xmlChildrenNode;
    char *tmp;
    int err = 0;

    sheet->maxrow = 0;
    sheet->maxcol = 0;

    while (p != NULL && !err) {
	if (!xmlStrcmp(p->name, (XUC) "Cell")) {
	    int i = -1, j = -1;

	    tmp = (char *) xmlGetProp(p, (XUC) "Row");
	    if (tmp) {
		i = atoi(tmp);
		free(tmp);
		if (i > sheet->maxrow) {
		    sheet->maxrow = i;
		}
	    }
	    tmp = (char *) xmlGetProp(p, (XUC) "Col");
	    if (tmp) {
		j = atoi(tmp);
		free(tmp);
		if (j > sheet->maxcol) {
		    sheet->maxcol = j;
		}
	    }
	    if (i == sheet->row_offset && j == sheet->col_offset) {
		err = inspect_top_left(p, obscol);
	    }
	}
	p = p->next;
    }

    if (!err) {
	fprintf(stderr, "wsheet_get_real_size: maxrow=%d, maxcol=%d\n",
		sheet->maxrow, sheet->maxcol);
    }

    return err;
}

/* Note that below we're being agnostic regarding the presence/absence
   of observation labels in the first column. We're writing first
   column values into sheet->labels if they're of string type and
   entering them into sheet->Z if they're numeric. Once we're finished
   we can decide what to do with the first column and the labels.
*/

static int wsheet_parse_cells (xmlNodePtr node, wsheet *sheet, 
			       int obscol, PRN *prn)
{
    xmlNodePtr p = node->xmlChildrenNode;
    char *tmp;
    double x;
    int vtype = 0;
    int gotlabels = 0;
    int cols, rows;
    int i, t, r, c;
    int err = 0;

    cols = sheet->maxcol + 1 - sheet->col_offset;
    rows = sheet->maxrow + 1 - sheet->row_offset;

    if (rows < 1) {
	pputs(prn, _("Starting row is out of bounds.\n"));
	return 1;
    }
    
    if (cols < 1) {
	pputs(prn, _("Starting column is out of bounds.\n"));
	return 1;
    }	

    if (wsheet_allocate(sheet, cols, rows)) {
	return 1;
    }

    sheet->colheads = 0;

    while (p != NULL && !err) {
	if (!xmlStrcmp(p->name, (XUC) "Cell")) {
	    x = NADBL;
	    c = r = 0;
	    i = t = -1;

	    /* what column are we in? */
	    tmp = (char *) xmlGetProp(p, (XUC) "Col");
	    if (tmp) {
		c = atoi(tmp);
		i = c - sheet->col_offset;
		free(tmp);
	    }

	    /* what row are we on? */
	    tmp = (char *) xmlGetProp(p, (XUC) "Row");
	    if (tmp) {
		r = atoi(tmp);
		t = r - sheet->row_offset;
		free(tmp);
	    }

	    if (i < 0 || t < 0) {
		/* we're not in the user-specified reading area */
		p = p->next;
		continue;
	    }

	    /* get cell type and content */
	    err = node_get_vtype_and_content(p, &vtype, &tmp);
	    if (err) {
		/* a formula perhaps? */
		pprintf(prn, _("Couldn't get value for col %d, row %d.\n"
			       "Maybe there's a formula in the sheet?"),
			c+1, r+1);
		break;
	    }

	    if (tmp != NULL) {
		if (VTYPE_IS_NUMERIC(vtype) || vtype == VALUE_STRING) {
		    if (i == 0) {
			/* first column: write content to labels */
			gretl_utf8_strncat_trim(sheet->label[t], tmp, OBSLEN - 1);
		    }
		}

		if (i == 0 && t == 1 && VTYPE_IS_NUMERIC(vtype)) {
		    char *fmt = (char *) xmlGetProp(p, (XUC) "ValueFormat");

		    if (fmt) {
			check_for_date_format(sheet, fmt);
			free(fmt);
		    }
		}

		if (VTYPE_IS_NUMERIC(vtype)) {
		    x = atof(tmp);
		    sheet->Z[i][t] = x;
		} else if (i > 0 && stray_numeric(vtype, tmp, &x)) {
		    sheet->Z[i][t] = x;
		} else if (vtype == VALUE_STRING) {
		    if (t == 0) {
			/* first row: look for varnames */
			strncat(sheet->varname[i], tmp, VNAMELEN - 1);
			sheet->colheads += 1;
			if (i == 0 && obscol) {
			    ; /* keep going */
			} else {
			    err = check_imported_varname(sheet->varname[i],
							 i, r, c, prn);
			}
		    } else if (i == 0 && obscol) {
			/* first column, not first row */
			if (!gotlabels) {
			    gotlabels = 1;
			}
			sheet->text_cols = 1;
		    } else {
			pprintf(prn, _("Expected numeric data, found string:\n"
				       "'%s' at row %d, column %d\n"), 
				tmp, r+1, c+1);
			err = 1;
		    }
		}
		free(tmp);
	    }
	}
	p = p->next;
    }

    if (gotlabels && sheet->colheads == 1) {
	/* rough notion here: if there's only one heading, it's
	   probably not really a variable name, but rather
	   a first observation label 
	*/
	sheet->colheads = 0;
    }

    return err;
}

static int wsheet_get_data (const char *fname, wsheet *sheet, 
			    int *obscol, PRN *prn) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int got_sheet = 0;
    int err;

    err = gretl_xml_open_doc_root(fname, "Workbook", &doc, &cur);
    if (err) {
	return err;
    }

    cur = cur->xmlChildrenNode;

    /* Now walk the tree */

    while (!err && cur != NULL && !got_sheet) {
	if (!xmlStrcmp(cur->name, (XUC) "Sheets")) {
	    int sheetcount = 0;

	    sub = cur->xmlChildrenNode;

	    while (sub != NULL && !got_sheet && !err) {
		if (!xmlStrcmp(sub->name, (XUC) "Sheet")) {
		    xmlNodePtr snode = sub->xmlChildrenNode;

		    while (snode != NULL && !err) {
			if (!xmlStrcmp(snode->name, (XUC) "Name")) {
			    sheetcount++;
			    tmp = (char *) xmlNodeGetContent(snode);
			    if (tmp) {
				tailstrip(tmp);
				if (!strcmp(tmp, sheet->name) &&
				    sheetcount == sheet->ID + 1) {
				    got_sheet = 1;
				}
				free(tmp);
			    }
			} else if (got_sheet && !xmlStrcmp(snode->name, (XUC) "Cells")) {
			    err = wsheet_get_real_size_etc(snode, sheet, obscol);
			    if (!err) {
				err = wsheet_parse_cells(snode, sheet, *obscol, prn);
			    }
			}
			snode = snode->next;
		    }
		}
		sub = sub->next;
	    }
	}
	cur = cur->next;
    }

    xmlFreeDoc(doc);

    if (!got_sheet) {
	fprintf(stderr, "'%s': couldn't find the requested sheet\n", sheet->name);
	err = 1;
    }

    return err;
}

static int wbook_record_name (char *name, wbook *book)
{
    char **sheetnames;
    int ns = book->nsheets + 1;

    sheetnames = realloc(book->sheetnames, ns * sizeof *sheetnames);
    if (sheetnames == NULL) {
	return 1;
    }

    book->sheetnames = sheetnames;
    book->nsheets = ns;
    book->sheetnames[ns - 1] = name;
    tailstrip(name);

    return 0;
}

static int wbook_get_info (const char *fname, const int *list,
			   char *sheetname, wbook *book, 
			   PRN *prn) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int got_index = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "Workbook", 
				  &doc, &cur);
    if (err) {
	return err;
    }

    wbook_init(book, list, sheetname);

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !got_index && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "SheetNameIndex")) {
	    got_index = 1;
	    sub = cur->xmlChildrenNode;
	    while (sub != NULL && !err) {
		if (!xmlStrcmp(sub->name, (XUC) "SheetName")) {
		    tmp = (char *) xmlNodeGetContent(sub);
		    if (tmp != NULL) {
			if (wbook_record_name(tmp, book)) {
			    err = 1;
			    free(tmp);
			}
		    }
		}
		sub = sub->next;
	    }
        }
	cur = cur->next;
    }

    xmlFreeDoc(doc);

    return err;
}

static int wsheet_setup (wsheet *sheet, wbook *book, int n)
{
    int err = 0;

    sheet->name = gretl_strdup(book->sheetnames[n]);

    if (sheet->name == NULL) {
	err = E_ALLOC;
    } else {
	sheet->ID = n;
	sheet->col_offset = book->col_offset;
	sheet->row_offset = book->row_offset;
    }   
    
    return err;
}

static int wsheet_labels_complete (wsheet *sheet)
{
    int rmin = (sheet->colheads)? 1 : 0;
    int rmax = sheet->maxrow + 1 - sheet->row_offset;
    int i, ret = 1;
    
    for (i=rmin; i<rmax; i++) {
	if (sheet->label[i][0] == '\0') {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

static void 
sheet_time_series_setup (wsheet *sheet, wbook *book, DATASET *newinfo, int pd)
{
    newinfo->pd = pd;
    newinfo->structure = TIME_SERIES;

    fprintf(stderr, "stobs='%s'\n", newinfo->stobs);
    newinfo->sd0 = get_date_x(newinfo->pd, newinfo->stobs);
    fprintf(stderr, "sd0=%g\n", newinfo->sd0);

    sheet->text_cols = 1;
    book_set_time_series(book);
    book_unset_obs_labels(book);
}

/* check that a given column doesn't contain all NAs */

static int column_is_blank (wsheet *sheet, int k, int n)
{
    int t, s = (sheet->colheads)? 1 : 0;

    for (t=0; t<n; t++) {
	if (!na(sheet->Z[k][s++])) {
	    return 0;
	}
    }

    return 1;
}

static int labels_are_index (char **S, int n)
{
    int t;

    /* note: the first label will be blank or "obs" or similar */

    for (t=1; t<=n; t++) {
	if (!integer_string(S[t])) {
	    return 0;
	} else if (atoi(S[t]) != t) {
	    return 0;
	}
    }  

    return 1;
}

int gnumeric_get_data (const char *fname, int *list, char *sheetname,
		       DATASET *dset, gretlopt opt, PRN *prn)
{
    int gui = (opt & OPT_G);
    wbook gbook;
    wbook *book = &gbook;
    wsheet gsheet;
    wsheet *sheet = &gsheet;
    int sheetnum = -1;
    int obscol = 0;
    DATASET *newset;
    int err = 0;

    newset = datainfo_new();
    if (newset == NULL) {
	pputs(prn, _("Out of memory\n"));
	return 1;
    }

    wsheet_init(sheet);

    gretl_push_c_numeric_locale();

    if (wbook_get_info(fname, list, sheetname, book, prn)) {
	pputs(prn, _("Failed to get workbook info"));
	err = 1;
	goto getout;
    } 

    wbook_print_info(book);

    if (book->nsheets == 0) {
	pputs(prn, _("No worksheets found"));
	err = 1;
	goto getout;
    }

    if (gui) {
	if (book->nsheets > 1) {
	    wsheet_menu(book, 1);
	    sheetnum = book->selected;
	} else {
	    wsheet_menu(book, 0);
	    sheetnum = 0;
	}
    } else {
	err = wbook_check_params(book);
	if (err) {
	    gretl_errmsg_set(_("Invalid argument for worksheet import"));
	} else if (book->selected >= 0) {
	    sheetnum = book->selected;
	} else {
	    sheetnum = 0;
	}
    }

    if (book->selected == -1) {
	/* canceled */
	err = -1;
    }

    if (!err && sheetnum >= 0) {
	fprintf(stderr, "Getting data...\n");
	err = wsheet_setup(sheet, book, sheetnum);
	if (!err) {
	    err = wsheet_get_data(fname, sheet, &obscol, prn);
	    if (err) {
		fprintf(stderr, "wsheet_get_data returned %d\n", err);
	    } else {
		wsheet_print_info(sheet);
		book->flags |= sheet->flags;
	    } 
	}
    } 

    if (err) {
	goto getout;
    } else {
	int r0 = 1; /* the first data row */
	int i, j, t;
	int ts_markers = 0;
	int merge = (dset->Z != NULL);
	char **ts_S = NULL;
	int blank_cols = 0;
	int missvals = 0;
	int pd = 0;

	if (obscol) {
	    book_set_obs_labels(book);
	    if (sheet->text_cols == 0) {
		sheet->text_cols = 1;
	    }
	} else if (sheet->text_cols > 0) {
	    /* string-valued variable? */
	    fprintf(stderr, "Problem: sheet->text_cols = %d\n", sheet->text_cols);
	}

	if (sheet->colheads == 0) {
	    book_set_auto_varnames(book);
	    r0 = 0;
	}

	if (book_numeric_dates(book)) {
	    fputs("found calendar dates in first imported column\n", stderr);
	} else if (obscol) {
	    fprintf(stderr, "found label strings in first imported column (text_cols = %d)\n",
		    sheet->text_cols);
	} else if (sheet->text_cols > 0) {
	    fputs("found string-valued variable in first imported column?\n", stderr);
	} else {
	    fputs("check for label strings in first imported column: not found\n", stderr);
	}

	newset->n = sheet->maxrow - sheet->row_offset;

	if (!sheet->colheads) {
	    pputs(prn, _("it seems there are no variable names\n"));
	    newset->n += 1;
	}

	if (book_numeric_dates(book) || obscol) {
	    pd = importer_dates_check(sheet->label + r0, &book->flags,
				      newset, prn, &err);
	    if (pd > 0) {
		/* got time-series info from dates/labels */
		sheet_time_series_setup(sheet, book, newset, pd);
		ts_markers = newset->markers;
		ts_S = newset->S;
	    } else if (!book_numeric_dates(book)) {
		if (labels_are_index(sheet->label, newset->n)) {
		    /* trash the labels */
		    book_unset_obs_labels(book);
		}
	    }
	}

	newset->v = sheet->maxcol + 2 - sheet->col_offset - sheet->text_cols;
	fprintf(stderr, "newset->v = %d, newset->n = %d\n",
		newset->v, newset->n);

	/* create import dataset */
	err = worksheet_start_dataset(newset);
	if (err) {
	    goto getout;
	}

	if (book_time_series(book)) {
	    newset->markers = ts_markers;
	    newset->S = ts_S;
	} else {
	    dataset_obs_info_default(newset);
	} 

	j = 1;
	for (i=1; i<newset->v; i++) {
	    int s = (sheet->colheads)? 1 : 0;
	    int k = i - 1 + sheet->text_cols;
	    double zkt;

	    if (column_is_blank(sheet, k, newset->n)) {
		blank_cols++;
		continue;
	    } 

	    if (sheet->colheads && *sheet->varname[k] != '\0') {
		strcpy(newset->varname[j], sheet->varname[k]);
	    } else {
		sprintf(newset->varname[j], "v%d", j);
	    }
	    for (t=0; t<newset->n; t++) {
		zkt = sheet->Z[k][s++];
		if (zkt == -999 || zkt == -9999) {
		    newset->Z[j][t] = NADBL;
		} else {
		    newset->Z[j][t] = zkt;
		}
		if (na(newset->Z[j][t])) {
		    missvals = 1;
		}
	    }
	    j++;
	}

	if (blank_cols > 0) {
	    fprintf(stderr, "Dropping %d apparently blank column(s)\n", 
		    blank_cols);
	    dataset_drop_last_variables(newset, blank_cols);
	}

	if (missvals) {
	    pputs(prn, _("Warning: there were missing values\n"));
	}

	if (fix_varname_duplicates(newset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (book_obs_labels(book) && wsheet_labels_complete(sheet)) {
	    int offset = (sheet->colheads)? 1 : 0;

	    dataset_allocate_obs_markers(newset);
	    if (newset->S != NULL) {
		for (t=0; t<newset->n; t++) {
		    strcpy(newset->S[t], sheet->label[t+offset]);
		}
	    }
	}

	if (book->flags & BOOK_DATA_REVERSED) {
	    reverse_data(newset, prn);
	}

	if (!err && !dataset_is_time_series(newset) && newset->S != NULL) {
	    /* we didn't time series info above, but it's possible
	       the observation strings carry such info
	    */
	    import_ts_check(newset);
	}

	err = merge_or_replace_data(dset, &newset, opt, prn);

	if (!err && !merge) {
	    dataset_add_import_info(dset, fname, GRETL_GNUMERIC);
	}

	if (!err && gui) {
	    wbook_record_params(book, list);
	}
    } 

 getout:

    wbook_free(book);
    wsheet_free(sheet);

    gretl_pop_c_numeric_locale();

    if (err && newset != NULL) {
	destroy_dataset(newset);
    }

    return err;
}


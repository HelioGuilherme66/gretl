/*
 *  Copyright (c) by Allin Cottrell 2002-2004
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

/*
  Based on the Gnumeric excel plugin by Michael Meeks.
*/

#include <gtk/gtk.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "libgretl.h"
#include "version.h"

#include "importer.h"
#include "biff.h"

static void free_sheet (void);
static int allocate_row_col (int row, int col, wbook *book);
static char *copy_unicode_string (unsigned char *src, int remlen, 
				  int *skip, int *slop);
static char *convert8to7 (const unsigned char *s, int count);
static char *convert16to7 (const unsigned char *s, int count);
static char *make_string (char *str);

int debug_print;

#ifdef WIN32
# include "build.h"
char debug_fname[FILENAME_MAX];
FILE *fdb;
static void make_debug_fname (void);
#endif

#define cell_record(r) (r == BIFF_LABEL || \
                        r == BIFF_STRING || \
                        r == BIFF_NUMBER || \
                        r == BIFF_RK || \
                        r == BIFF_MULRK || \
                        r == BIFF_FORMULA || \
                        r == BIFF_LABELSST)

struct sheetrow {
    int last, end;
    char **cells;
};	

enum {
    VARNAMES_OK = 0,
    VARNAMES_NULL,
    VARNAMES_NOTSTR,
    VARNAMES_INVALID
} varname_errors;

char **sst = NULL;
int sstsize = 0, sstnext = 0;
struct sheetrow *rows = NULL;
int nrows = 0;

#define EXCEL_IMPORTER
#include "import_common.c"

const char *adjust_rc = N_("Perhaps you need to adjust the "
			   "starting column or row?");

#ifdef WIN32
static void make_debug_fname (void)
{
    read_reg_val(HKEY_CURRENT_USER, "gretl", "userdir", debug_fname);
    strcat(debug_fname, "xls.log");
}

static void open_debug_stream (void)
{
    if (debug_fname[0] != '\0') {
	fdb = fopen(debug_fname, "w");
    }
}
#endif

static int dprintf (const char *format, ...)
{
    va_list args;
    int len = 0;

#ifdef WIN32
    if (debug_print) {
	if (fdb == NULL) {
	    open_debug_stream();
	}
	if (fdb != NULL) {
	    va_start(args, format);
	    len = vfprintf(fdb, format, args);
	    va_end(args);
	    fflush(fdb);
	}
    }
#else
    if (debug_print) {
	va_start(args, format);
	len = vfprintf(stderr, format, args);
	va_end(args);
	fflush(stderr);
    }
#endif

    return len;
}

static void print_version (void)
{
#ifdef WIN32
    dprintf("gretl, version %s, %s\n", GRETL_VERSION, BUILD_DATE);
#else
    dprintf("gretl, version %s\n", GRETL_VERSION); 
#endif
}

static double get_le_double (const unsigned char *rec) 
{
    union { 
        unsigned char cc[8];
        double d;
    } dconv;

    unsigned char *d;
    const unsigned char *s;
    int i;

    if (sizeof(double) != 8) {
	fputs("Size of double != 8; this won't work!\n", stderr);
	return NADBL;
    }

#if G_BYTE_ORDER == G_BIG_ENDIAN
    for (s=rec+8, d=dconv.cc, i=0; i<8; i++) *(d++) = *(--s);
#else       
    for (s=rec, d=dconv.cc, i=0; i<8; i++) *(d++) = *(s++);
#endif     

    return dconv.d;
}

static double biff_get_rk (const unsigned char *ptr)
{
    gint32 number;
    enum eType {
	eIEEE = 0, 
	eIEEEx100, 
	eInt, 
	eIntx100
    } type;

    number = MS_OLE_GET_GUINT32(ptr);
    type = (number & 0x3);

    switch (type) {
    case eIEEE:
    case eIEEEx100:
        {
	    guint8 tmp[8];
	    double answer;
	    int lp;

	    for (lp = 0; lp < 4; lp++) {
		tmp[lp + 4] = (lp > 0) ? ptr[lp]: (ptr[lp] & 0xfc);
		tmp[lp] = 0;
	    }
	    answer = get_le_double(tmp);
	    return (type == eIEEEx100)? answer / 100 : answer;
        }
    case eInt:
	return (double) (number >> 2);
    case eIntx100:
	number >>= 2;
	if ((number % 100) == 0) {
	    return (double) (number/100);
	} else {
	    return (double) (number/100.0);
	}
    }

    return NADBL;
}

static int row_col_err (int row, int col, PRN *prn) 
{
    static int prevrow = -1, prevcol = -1;
    int err = 0;
    
    if (row < 0 || col < 0) {
	fprintf(stderr, "Error: got row=%d, col=%d\n", row, col);
	err = 1;
    } else if (row == prevrow && col == prevcol) {
	pprintf(prn, "Error: found a second cell entry for cell (%d, %d)\n",
		prevrow, prevcol);
	err = 1;
    }
	
    prevrow = row;
    prevcol = col;

    return err;
}

static int check_copy_string (struct sheetrow *prow, int row, int col, 
			      int idx, const char *s)
{
    if (row > 0 && col > 0) {
	int i, len = strlen(s);
	int commas = 0;
	static int warned = 0;

	for (i=0; i<len; i++) {
	    if (!isdigit(s[i]) && s[i] != ' ' && s[i] != '-' &&
		s[i] != ',' && s[i] != '.') {
		len = 0;
		break;
	    }
	    if (s[i] == ',') commas++;
	}

	if (len > 0) {
	    char *p, *q = malloc(len + 1);

	    if (q == NULL) return 1;

	    p = q;
	    for (i=0; i<len; i++) {
		if (s[i] != ' ' && s[i] != ',') {
		    *p++ = s[i];
		} 
		if (commas == 1 && s[i] == ',') {
		    /* single comma could be for 000s, or decimal */
		    if (!warned) {
			fprintf(stderr, "Warning: found ambiguous comma in '%s'\n", s);
			warned = 1;
		    }
		    if (len - i != 4) {
			/* probably decimal? */
			*p++ = '.';
		    }
		}
	    }
	    *p = '\0';
	    dprintf("converting sst[%d] '%s' to numeric as %s\n", idx, s, q);
	    prow->cells[col] = g_strdup_printf("%s", q);
	    free(q);
	    return 0;
	}
    }

    dprintf("copying sst[%d] '%s' into place\n", idx, s);
    prow->cells[col] = g_strdup_printf("\"%s", s);

    return 0;
}

static int wbook_find_format (wbook *book, int xfref)
{
    int fmt = -1;

    if (book->xf_list != NULL && xfref < book->xf_list[0]) {
	fmt = book->xf_list[xfref + 1];
    }

    return fmt;
}

#undef FORMAT_INFO

static int process_item (BiffQuery *q, wbook *book, PRN *prn) 
{
    struct sheetrow *prow = NULL;
    static char **string_targ;
    static int slop; /* SST overslop */
    unsigned char *ptr = NULL;
    int i = 0, j = 0;
    double val;

    if (cell_record(q->ls_op)) {
	i = EX_GETROW(q);
	j = EX_GETCOL(q);
	if (row_col_err(i, j, prn)) {
	    return 1;
	}
	if (q->ls_op == BIFF_NUMBER || q->ls_op == BIFF_RK) {
	    guint16 xfref = EX_GETXF(q);
	    int fmt = wbook_find_format(book, xfref);

	    fprintf(stderr, "Numeric cell (%d, %d), XF index = %d, fmt = %d\n", 
		    i, j, (int) xfref, fmt);
	}
    }

    switch (q->ls_op) {

    case BIFF_SST: {
	int k, skip, remlen, oldsz = sstsize;
	guint16 sz;

	if (sst != NULL) {
	    fprintf(stderr, "Got a second string table: this is nonsense\n");
	    return 1;
	}

	sz = MS_OLE_GET_GUINT16(q->data + 4);
	sstsize += sz;
	sst = realloc(sst, sstsize * sizeof *sst);
	if (sst == NULL) {
	    return 1;
	}

	dprintf("Got SST: allocated for %d strings (%d bytes), %p\n", 
		sstsize, sstsize * sizeof *sst, (void *) sst);

	for (k=oldsz; k<sstsize; k++) {
	    /* careful: initialize all pointers to NULL */
	    sst[k] = NULL;
	}

	ptr = q->data + 8;

	for (k=oldsz; k<sstsize; k++) {
	    remlen = q->length - (ptr - q->data);
	    dprintf("Working on sst[%d], data offset=%d, remlen=%d\n", 
		    k, (int) (ptr - q->data), remlen);
	    if (remlen <= 0) {
		break;
	    }
	    sst[k] = copy_unicode_string(ptr, remlen, &skip, &slop);
	    ptr += skip;
	}

	if (k < sstsize) {
	    sstnext = k;
	}

	break;
    }	

    case BIFF_CONTINUE: 
	dprintf("Got CONTINUE, sstnext = %d, len = %d\n", 
		sstnext, (int) q->length);
	if (sstnext > 0) {
	    int k, skip, remlen;

	    ptr = q->data;
	    if (slop > 0) {
		unsigned char flags = *ptr;
		int csize = (flags & 0x01)? 2 : 1;

		ptr += 1 + csize * slop;
	    }
	    for (k=sstnext; k<sstsize; k++) {
		remlen = q->length - (ptr - q->data);
		if (remlen <= 0) {
		    break;
		}
		dprintf("Working on sst[%d]\n", k);
		sst[k] = copy_unicode_string(ptr, remlen, &skip, &slop);
		ptr += skip;
	    }
	    if (k < sstsize) {
		sstnext = k;
	    }
	}
	break;
			   
    case BIFF_LABEL: 
	dprintf("Got LABEL, row=%d, col=%d\n", i, j);
	if (allocate_row_col(i, j, book)) {
	    return 1;
	} else {
	    unsigned int len = MS_OLE_GET_GUINT16(q->data + 6);
	
	    prow = rows + i;
	    ptr = q->data + 8;
	    prow->cells[j] = make_string(convert8to7(ptr, len));
	}
	break;
  
    case BIFF_LABELSST:
	dprintf("Got LABELSST, row=%d, col=%d\n", i, j);
	if (allocate_row_col(i, j, book)) {
	    return 1;
	} else {
	    unsigned int sidx = MS_OLE_GET_GUINT16(q->data + 6);

	    prow = rows + i;
	    if (sidx >= sstsize) {
		pprintf(prn, _("String index too large"));
		pputc(prn, '\n');
	    } else if (sst[sidx] != NULL) {
		check_copy_string(prow, i, j, sidx, sst[sidx]);
	    } else {
		prow->cells[j] = malloc(2);
		if (prow->cells[j] != NULL) {
		    prow->cells[j][0] = '\0';
		}
	    }
	}	
	break;

    case BIFF_NUMBER: 
	if (allocate_row_col(i, j, book)) {
	    return 1;
	} else {
	    val = get_le_double(q->data + 6);
	    prow = rows + i;
	    prow->cells[j] = g_strdup_printf("%.10g", val);
	    dprintf("Got NUMBER (%g), row=%d, col=%d\n", val, i, j);
	}
	break;

    case BIFF_RK: 
	if (allocate_row_col(i, j, book)) {
	    return 1;
	} else {
	    val = biff_get_rk(q->data + 6);
	    prow = rows + i;
	    prow->cells[j] = g_strdup_printf("%.10g", val);
	    dprintf("Got RK (%g), row=%d, col=%d\n", val, i, j);
	}
	break;

    case BIFF_MULRK: {
	int k, ncols = (q->length - 6) / 6;

	dprintf("Got MULRK, row=%d, first_col=%d, ncols=%d\n", i, j, ncols);
	for (k=0; k<ncols; k++) {
	    if (allocate_row_col(i, j, book)) {
		return 1;
	    }
	    val = biff_get_rk(q->data + 6 + 6 * k);
	    prow = rows + i; /* might have moved */
	    prow->cells[j] = g_strdup_printf("%.10g", val);
	    dprintf(" MULRK[col=%d] = %g\n", j, val);
	    j++;
	}
	break;
    }

    case BIFF_FORMULA:  
	dprintf("Got FORMULA, row=%d, col=%d\n", i, j);
	if (allocate_row_col(i, j, book)) {
	    return 1;
	} else {
	    ptr = q->data + 6;
	    prow = rows + i;
	    if (ptr[6] == 0xff && ptr[7] == 0xff) {
		unsigned char fcode = ptr[0];

		if (fcode == 0x0) {
		    /* string formula: record target for following 
		       STRING record */
		    string_targ = prow->cells + j;
		} else if (fcode == 0x1) {
		    /* boolean value */
		    prow->cells[j] = g_strdup((ptr[2])? "1" : "0");
		} else if (fcode == 0x2 || fcode == 0x3) {
		    /* error code or empty */
		    prow->cells[j] = g_strdup("-999.0");
		} else {
		    fprintf(stderr, "Bad formula code 0x%u\n", 
			    (unsigned) fcode);
		    prow->cells[j] = g_strdup("-999.0");
		}
	    } else {
		/* floating-point */
		val = get_le_double(ptr);
		dprintf(" floating-point value = %g\n", val);
		if (isnan(val)) {
		    fprintf(stderr, "Got a NaN\n");
		    prow->cells[j] = g_strdup("-999.0");
		} else {
		    prow->cells[j] = g_strdup_printf("%.10g", val);
		}
	    }
	}
	break;

    case BIFF_STRING: 
	if (string_targ == NULL) {
	    dprintf("String record without preceding string formula\n");
	} else {
	    char *tmp = copy_unicode_string(q->data, 0, NULL, NULL);

	    *string_targ = make_string(tmp);
	    dprintf("Filled out string formula with '%s'\n", *string_targ);	
	    string_targ = NULL;
	}
	break;

    case BIFF_BOF: 
	if (rows != NULL) {
	    fprintf(stderr, "BOF when current sheet is not flushed\n");
	    return 1;
	}
	if (1) {
	    unsigned version, boftype;

	    version = MS_OLE_GET_GUINT16(q->data + 0);
	    boftype = MS_OLE_GET_GUINT16(q->data + 2);
	    dprintf("Got BOF: version=%x, type=%x\n", version, boftype);
	}
	break;

    case BIFF_FORMAT: {
	int idx = MS_OLE_GET_GUINT16(q->data + 0);

	if (idx >= 14 && idx <= 17) {
	    fprintf(stderr, "Got date format: index %d\n", idx);
	}
	break;
    }

#ifdef FORMAT_INFO
    case BIFF_COLINFO:
	fprintf(stderr, "Got BIFF_COLINFO: col range (%d, %d), XF index %d\n",
		(int) MS_OLE_GET_GUINT16(q->data + 0),
		(int) MS_OLE_GET_GUINT16(q->data + 2),
		(int) MS_OLE_GET_GUINT16(q->data + 6));
	break;

    case BIFF_XF: {
	unsigned short tp = MS_OLE_GET_GUINT16(q->data + 4);

	fprintf(stderr, "Got BIFF_XF: format record index %d ",
		(int) MS_OLE_GET_GUINT16(q->data + 2));
	if (tp & 0x04) {
	    fprintf(stderr, "(style XF)\n");
	} else {
	    fprintf(stderr, "(cell XF)\n");
	}
	break;
    }
#endif

    default: 
	break;
    }

    return 0;
}

static int handled_record (BiffQuery *q)
{ 
    if (q->opcode == BIFF_SST || 
	q->opcode == BIFF_CONTINUE ||  
	q->opcode == BIFF_LABELSST || 
	q->opcode == BIFF_MULRK || 
	q->opcode == BIFF_FORMULA) {
	return 1;
    }

    if (q->ms_op == 0x02) {
	if (q->ls_op == BIFF_LABEL ||
	    q->ls_op == BIFF_NUMBER ||
	    q->ls_op == BIFF_RK ||
	    q->ls_op == BIFF_STRING) {
	    return 1;
	}
    }

#ifdef FORMAT_INFO
    if (q->opcode == BIFF_COLINFO ||
	q->opcode == BIFF_XF) {
	return 1;
    }
    if (q->ms_op == 0x04) {
	if (q->ls_op == BIFF_FORMAT) {
	    return 1;
	}
    } 
#endif   

    if (q->ms_op == 0x08) {
	if (q->ls_op == BIFF_BOF) return 1;
    }

    return 0;
}

static int process_sheet (const char *filename, wbook *book, PRN *prn) 
{    
    int err = 0, gotbof = 0, eofcount = 0;
    long offset = book->byte_offsets[book->selected];

    MsOleStream *stream;
    MsOleErr result;
    BiffQuery *q;
    MsOle *file;

    if (ms_ole_open(&file, filename)) {
	return 1;
    }

    result = ms_ole_stream_open(&stream, file, "/", "workbook", 'r');

    if (result != MS_OLE_ERR_OK) {
	ms_ole_stream_close (&stream);

	result = ms_ole_stream_open(&stream, file, "/", "book", 'r');
	if (result != MS_OLE_ERR_OK) {
	    ms_ole_stream_close(&stream);
	    fprintf(stderr, _("No book or workbook streams found."));
	    return 1;
	}
    }

    fprintf(stderr, _("Reading file...\n"));
    q = ms_biff_query_new(stream);

    while (!gotbof && ms_biff_query_next(q)) {
	if (q->ls_op == BIFF_BOF) {
	    gotbof = 1;
	    break;
	}
    }

    if (!gotbof) {
	pprintf(prn, _("%s: No BOF record found"), filename);
	return 1;
    } 

    while (!err && ms_biff_query_next(q)) {
	dprintf("At %lu: q->opcode=0x%02x\n", (unsigned long) q->streamPos, q->opcode);
	if (q->opcode == BIFF_EOF) {
	    dprintf("got MSEOF at %lu\n", (unsigned long) stream->position);
	    eofcount++;

	    if (eofcount == 1) {
		if (stream->position < offset) {
		    /* skip to the worksheet we want? */
		    while (q->streamPos < offset && ms_biff_query_next(q)) ;
		    fprintf(stderr, "skipped forward to %lu\n", 
			    (unsigned long) q->streamPos);
		} else {
		    fprintf(stderr, "reading worksheet at %lu\n", 
			    (unsigned long) stream->position);
		}
	    }

	    if (eofcount == 2) {
		break;
	    } else {
		continue;
	    }
	} 

	if (handled_record(q)) {
	    err = process_item(q, book, prn);
	} else if (q->ms_op == 0x02 && q->ls_op == BIFF_ROW) {
	    dprintf("Got BIFF_ROW\n");
	} else if (q->opcode == BIFF_DBCELL) {
	    dprintf("Got BIFF_DBCELL\n");
	} else {
	    dprintf("skipping unhandled opcode 0x%02x\n", q->opcode);
	}
    }

    ms_biff_query_destroy(q);
    ms_ole_stream_close(&stream);
    ms_ole_destroy(&file);

    return err;    
}

static char *make_string (char *str) 
{
    char *ret = NULL;

    if (str != NULL) {
	ret = g_strdup_printf("\"%s", str);
	free(str);
    } else {
	ret = g_strdup("\"");
    }

    return ret;
}    

static char *
copy_unicode_string (unsigned char *src, int remlen, 
		     int *skip, int *slop) 
{
    int count = MS_OLE_GET_GUINT16(src);
    unsigned char flags = *(src + 2);
    int this_skip = 3, skip_to_next = 3;
    int csize = (flags & 0x01)? 2 : 1;

    dprintf("copy_unicode_string: count = %d, csize = %d\n",
	    count, csize);

    if (flags & 0x08) {
	dprintf(" contains Rich-Text info\n");
    }
    if (flags & 0x04) {
	dprintf(" contains Far-East info\n");
    }    

    skip_to_next += count * csize;

    if (flags & 0x08) {
	guint16 rich_text_info_len = 0;

	rich_text_info_len = 4 * MS_OLE_GET_GUINT16(src + 3);
	this_skip += 2;
	skip_to_next += 2 + rich_text_info_len;
    } 

    if (flags & 0x04) {
	guint32 far_east_info_len = 0;
	int far_east_offset = 3;

	if (flags & 0x08) {
	    far_east_offset = 5;
	}
	far_east_info_len = MS_OLE_GET_GUINT32(src + far_east_offset);
	this_skip += 4;
	skip_to_next += 4 + far_east_info_len;
    }

    /* skip for the next read */
    if (skip != NULL) {
	*skip = skip_to_next;
    }

    /* size check */
    if (slop != NULL) {
	if (remlen > 0 && this_skip + count > remlen) {
	    *slop = this_skip + count - remlen;
	} else {
	    *slop = 0;
	}
    }

    if (count > 64) {
	/* let's not mess with excessive strings */
	return g_strdup("bigstr");
    } else if (csize == 1) {
	return convert8to7(src + this_skip, count);
    } else { 
	return convert16to7(src + this_skip, count);
    }
}

static char *convert8to7 (const unsigned char *s, int count) 
{
    char *dest;

    if (count > USER_VLEN - 1) {
	count = USER_VLEN - 1;
    }

    dest = malloc(USER_VLEN);
    *dest = '\0';
    strncat(dest, s, count);
    iso_to_ascii(dest);

    if (*dest == '\0') {
	strcpy(dest, "varname");
    }

    dprintf("convert8to7: returning '%s'\n", dest);

    return dest;
}

static char *convert16to7 (const unsigned char *s, int count) 
{
    char *dest;
    int i, u, j = 0;

    dest = malloc(USER_VLEN);
    if (dest == NULL) {
	return NULL;
    }

    memset(dest, 0, USER_VLEN);

#if 1
    for (i=0; i<count && j<USER_VLEN-1; i++) {
	u = MS_OLE_GET_GUINT16(s);
	s += 2;
	if ((isalnum(u) || ispunct(u)) && u < 128) {
	    dest[j++] = u;
	}
    }
#else
    printf("convert16to7: count = %d\n", count);
    for (i=0; i<count; i++) {
	u = MS_OLE_GET_GUINT16(s);
	s += 2;
	if ((isalnum(u) || ispunct(u) || u == ' ') && u < 128) {
	    putchar(u);
	    if (j<USER_VLEN-1) {
		dest[j++] = u;
	    }
	    if (i % 56 == 0) {
		putchar('\n');
	    }	
	}
    }
    printf("\nconvert16to7: returning '%s'\n", dest);
#endif

    if (*dest == '\0') {
	strcpy(dest, "varname");
    }

    dprintf("convert16to7: returning '%s'\n", dest);

    return dest;    
}

static void row_init (struct sheetrow *row)
{
    row->last = 0;
    row->end = 0;
    row->cells = NULL;
}

static int allocate_row_col (int i, int j, wbook *book) 
{
    struct sheetrow *myrows = NULL;
    char **cells = NULL;
    int new_nrows, newcol, k;
    static int started;

    if (!started && i > book->row_offset) {
	book->row_offset = i;
	fprintf(stderr, "Missing rows: trying an offset of %d\n", i);
    }

    started = 1;

    dprintf("allocate: row=%d, col=%d, nrows=%d\n", i, j, nrows);

    if (i >= nrows) {
	new_nrows = (i / 16 + 1) * 16;

	myrows = realloc(rows, new_nrows * sizeof *myrows);

	if (myrows == NULL) {
	    return 1;
	}

	rows = myrows;

	for (k=nrows; k<new_nrows; k++) {
	    dprintf("allocate: initing rows[%d]\n", k);
	    row_init(&rows[k]);
	    dprintf("rows[%d].end=%d\n", i, rows[k].end);
	}
	nrows = new_nrows;
    }

    dprintf("allocate: col=%d and rows[%d].end = %d\n", j, i, rows[i].end);

    if (j >= rows[i].end) {
	newcol = (j / 16 + 1) * 16;
	dprintf("allocate: reallocing rows[%d].cells to size %d\n", i, newcol);
	cells = realloc(rows[i].cells, newcol * sizeof *cells);

	if (cells == NULL) {
	    return 1;
	}

	rows[i].cells = cells;

	for (k=rows[i].end; k<newcol; k++) {
	    rows[i].cells[k] = NULL;
	}
	rows[i].end = newcol;
    } 
 
    if (j > rows[i].last) {
	rows[i].last = j;
    }

    return 0;
}

static void free_sheet (void) 
{
    int i, j;

    dprintf("free_sheet(), nrows=%d\n", nrows);

    /* free shared string table */
    if (sst != NULL) {
	for (i=0; i<sstsize; i++) {
	    if (sst[i] != NULL) free(sst[i]);
	}
	free(sst);
    }

    /* free cells */
    if (rows != NULL) {
	for (i=0; i<nrows; i++) {
	    if (rows[i].cells == NULL) {
		dprintf("rows[%d].cells = NULL, skipping free\n", i);
		continue;
	    }
	    for (j=0; j<rows[i].end; j++) {
		if (rows[i].cells[j] != NULL) {
		    dprintf("Freeing rows[%d].cells[%d] at %p\n",
			    i, j, (void *) rows[i].cells[j]);
		    free(rows[i].cells[j]); 
		}
	    }
	    dprintf("Freeing rows[%d].cells at %p\n", i, (void *) rows[i].cells);
	    free(rows[i].cells);
	}
	free(rows);
	rows = NULL;
    }

    nrows = 0;
}

#define IS_STRING(v) ((v[0] == '"'))

/* check for full set of strings in first column to be read (which may
   be at an offset into the worksheet)
 */

static int first_col_strings (wbook *book)
{
    int i, j = book->col_offset;
    int startrow = book->row_offset + 1;
    int ret = 1;
    
    for (i=startrow; i<nrows; i++) {
	dprintf("book->row_offset=%d, i=%d\n", book->row_offset, i);
	dprintf("rows = %p\n", (void *) rows);
	if (rows == NULL || rows[i].cells == NULL || 
	    rows[i].cells[j] == NULL ||
	    !IS_STRING(rows[i].cells[j])) {
	    ret = 0;
	    break;
	}
	dprintf("first_col_strings: rows[%d].cells[%d]: '%s'\n", i, j,
		rows[i].cells[j]);
    }

    return ret;
}

#define obs_string(s) (!strcmp(s, "obs") || !strcmp(s, "id"))

static int fix_varname (char *vname)
{
    int i, n = strlen(vname);
    int bad = 0;

    for (i=1; i<n; i++) {
        if (!(isalpha((unsigned char) vname[i]))  
            && !(isdigit((unsigned char) vname[i]))
            && vname[i] != '_') {
	    vname[i] = '_';
	    bad++;
	}
    }

    return (bad == n);
}   

static int 
check_all_varnames (wbook *book, int ncols, const char *blank_col, int skip)
{
    int j, i = book->row_offset;
    int startcol = skip + book->col_offset;
    char *test;

    for (j=startcol; j<ncols; j++) { 
	if (blank_col[j]) {
	    continue;
	}
	if (rows[i].cells[j] == NULL) {
	    dprintf("got_varnames: rows[%d].cells[%d] is NULL\n", i, j);
	    return VARNAMES_NULL;
	}

	dprintf("got_varnames: rows[%d].cells[%d] is '%s'\n", i, j, 
		rows[i].cells[j]);

	if (!IS_STRING(rows[i].cells[j])) {
	    return VARNAMES_NOTSTR;
	}

	/* skip beyond the quote */
	test = rows[i].cells[j] + 1;

	/* "obs" or "id" is OK in the first col of the selection, 
	   but not thereafter */
	if (j == startcol && obs_string(test)) {
	    /* pass along */
	    ;
	} else {
	    int verr = check_varname(test);

	    if (verr == VARNAME_BADCHAR) {
		verr = fix_varname(test);
	    }
	    
	    if (verr) {
		return VARNAMES_INVALID;
	    }
	}
    }

    return VARNAMES_OK;
}

static int missval_string (const char *s)
{
    if (s + 1 == 0) {
	return 1;
    } else {
	char test[6] = {0};

	strncat(test, s + 1, 4);
	tailstrip(test);
	lower(test);
	if (!strcmp(test, "na") || 
	    !strcmp(test, "n.a.") ||
	    !strcmp(test, "..") ||
	    !strcmp(test, "?")) {
	    return 1;
	}
    }

    return 0;
}

struct string_err {
    int row;
    int column;
    char *str;
};

/* check for invalid data in the selected data block */

static int 
check_data_block (wbook *book, int ncols, const char *blank_col,
		  int skip, struct string_err *err)
{
    int startcol = skip + book->col_offset;
    int startrow = book->row_offset + 1;
    int j, i, ret = 0;

    err->row = 0;
    err->column = 0;
    err->str = NULL;

    for (j=startcol; j<ncols; j++) {
	if (blank_col[j]) {
	    continue;
	}
	for (i=startrow; i<nrows; i++) {
	    dprintf("data_block: looking at rows[%d], end = %d\n", i, rows[i].end);
	    if (rows[i].cells  == NULL) {
		dprintf("data_block: rows[%d].cells = NULL\n", i);
		ret = -1;
	    } else if (j >= rows[i].end) {
		dprintf("data_block: short row, fell off the end\n");
		ret = -1;
	    } else if (rows[i].cells[j] == NULL) {
		dprintf("data_block: rows[%d].cells[%d] = NULL\n", i, j);
		rows[i].cells[j] = g_strdup("-999.0");
		ret = -1;
	    } else if (IS_STRING(rows[i].cells[j])) {
		if (missval_string(rows[i].cells[j])) {
		    free(rows[i].cells[j]);
		    rows[i].cells[j] = g_strdup("-999.0");
		    ret = -1;
		} else {
		    err->row = i + 1;
		    err->column = j + 1;
		    err->str = g_strdup(rows[i].cells[j]);
		    return 1;
		}
	    }
	}
    }

    return ret;
}

/* determine the number of actual data columns, starting from a
   given offset into the worksheet, and allowing for the
   possibility that the first selected column contains string
   labels
*/

static int 
n_vars_from_col (int totcols, char *blank_col, int offset, int first_special)
{
    int i, nv = 1;

    if (first_special) {
	offset++;
    }

    for (i=offset; i<totcols; i++) {
	if (!blank_col[i]) nv++;
    }

    dprintf("n_vars_from_col: totcols=%d, nv=%d\n", totcols, nv);

    return nv;
}

static int transcribe_data (double **Z, DATAINFO *pdinfo, 
			    int row_offset, int totcols, int startcol, 
			    char *blank_col)
{
    int i, t, j = 1;

    for (i=startcol; i<totcols; i++) {

	if (blank_col[i]) {
	    continue;
	}

	if (j >= pdinfo->v) {
	    break;
	}

	pdinfo->varname[j][0] = 0;
	strncat(pdinfo->varname[j], rows[row_offset].cells[i] + 1, 
		USER_VLEN - 1);
	dprintf("accessing rows[%d].cells[%d] at %p\n",
		row_offset, i, (void *) rows[row_offset].cells[i]);
	dprintf("set varname[%d] = '%s'\n", j, pdinfo->varname[j]);

	for (t=0; t<pdinfo->n; t++) {
	    int ts = t + 1 + row_offset;

	    if (rows[ts].cells == NULL || i >= rows[ts].end ||
		rows[ts].cells[i] == NULL) {
		continue;
	    }
	    dprintf("accessing rows[%d].cells[%d] at %p\n", ts, i,
		    (void *) rows[ts].cells[i]);
	    dprintf("setting Z[%d][%d] = rows[%d].cells[%d] "
		    "= '%s'\n", j, t, i, ts, rows[ts].cells[i]);
	    Z[j][t] = atof(rows[ts].cells[i]);
	    if (Z[j][t] == -999.0) {
		Z[j][t] = NADBL;
	    }
	}
	j++;
    }

    return 0;
}

static int 
get_sheet_dimensions (int *totcols, int *datacols, char **blank_col,
		      PRN *prn)
{
    char *blanks = NULL;
    int i, j;

    *totcols = 0;
    *datacols = 0;
    *blank_col = NULL;

    /* trim any trailing blank rows */
    for (i=nrows-1; i>=0; i--) {
	if (rows[i].cells == NULL) {
	    nrows--;
	} else {
	    break;
	}
    }
	
    for (i=0; i<nrows; i++) {
	if (rows[i].cells != NULL) {
	    if (rows[i].last + 1 > *totcols) {
		*totcols = rows[i].last + 1;
	    }
	}
    }

    if (*totcols <= 0 || nrows < 1) {
	pputs(prn, _("No data found.\n"));
	pputs(prn, _(adjust_rc));
	return 1;
    }

    blanks = malloc(*totcols);
    if (blanks == NULL) {
	return E_ALLOC;
    }

    memset(blanks, 1, *totcols);

    for (i=0; i<nrows; i++) {
	if (rows[i].cells == NULL) {
	    continue;
	}
	for (j=0; j<=rows[i].last; j++) {
	    if (rows[i].cells[j] != NULL) {
		if (blanks[j]) {
		    blanks[j] = 0;
		}
	    }
	}
    }

    for (i=0; i<*totcols; i++) {
	if (!blanks[i]) {
	    *datacols += 1;
	}
    }

    printf("rows=%d, data cols=%d total cols=%d\n", nrows, 
	   *datacols, *totcols);

    if (*datacols < 1) {
	pputs(prn, _("No data found.\n"));
	pputs(prn, _(adjust_rc));
	return 1;
    }

    *blank_col = blanks;

    return 0;
}

int excel_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		    PRN *prn)
{
    wbook book;
    double **newZ = NULL;
    DATAINFO *newinfo;
    int label_strings, time_series = 0;
    int datacols, totcols;
    struct string_err strerr;
    char *blank_col = NULL;
    int i, t;
    int err = 0;

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	return 1;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    wbook_init(&book);

    if (excel_book_get_info(fname, &book)) {
	pputs(prn, _("Failed to get workbook info"));
	err = 1;
    } else if (book.nsheets == 0) {
	pputs(prn, _("No worksheets found"));
	err = 1;
    } else {
	wbook_print_info(&book);
    }

    if (!err) {
	if (book.nsheets > 1) {
	    wsheet_menu(&book, 1);
	} else {
	    wsheet_menu(&book, 0);
	}
	debug_print = book.debug;
	if (debug_print) {
	    print_version();
	}
    }

    dprintf("sheet selected=%d; import offsets: col=%d, row=%d\n",
	    book.selected, book.col_offset, book.row_offset);

    if (book.selected == -1) {
	/* canceled */
	err = -1; 
    }

    if (err) goto getout;

    /* processing for specific worksheet */
    err = process_sheet(fname, &book, prn);

    if (err) {
	const char *buf = gretl_print_get_buffer(prn);

	if (*buf == 0) {
	    pputs(prn, _("Failed to process Excel file"));
	    buf = gretl_print_get_buffer(prn);
	}
	fprintf(stderr, "%s\n", buf);
	goto getout;
    }

    /* locate any blank columns */
    err = get_sheet_dimensions(&totcols, &datacols, &blank_col, prn);

    if (err) goto getout;

    /* check feasibility of offsets */
    if (book.row_offset > nrows) {
	pputs(prn, _("Starting row is out of bounds.\n"));
	err = 1;
    } else if (book.col_offset > totcols) {
	pputs(prn, _("Starting column is out of bounds.\n"));
	err = 1;
    }

    if (err) goto getout;

    label_strings = first_col_strings(&book);
    if (label_strings) {
	puts("found label strings in first imported column");
    } else {
	puts("check for label strings in first imported column: not found");
    }

    /* any bad or missing variable names? */
    err = check_all_varnames(&book, totcols, blank_col, label_strings);

    if (err == VARNAMES_NULL || err == VARNAMES_NOTSTR) {
	pputs(prn, _("One or more variable names are missing.\n"));
	pputs(prn, _(adjust_rc));
    } else if (err == VARNAMES_INVALID) {
	invalid_varname(prn);
    }

    if (err) goto getout; 

    /* any bad data? */
    err = check_data_block(&book, totcols, blank_col, label_strings, &strerr);

    if (err == 1) {
	pprintf(prn, _("Expected numeric data, found string:\n"
		       "%s\" at row %d, column %d\n"),
		strerr.str, strerr.row, strerr.column);
	g_free(strerr.str);
	pputs(prn, _(adjust_rc));
	goto getout; 
    } else if (err == -1) {
	pputs(prn, _("Warning: there were missing values\n"));
	err = 0;
    }	    

    /* do we have a first column containing dates? */
    if (obs_column_heading(rows[book.row_offset].cells[book.col_offset])) {
	int pd = consistent_date_labels(nrows, book.row_offset, book.col_offset, NULL, 
					book.d1904);

	if (pd) {
	    time_series_setup(rows[1 + book.row_offset].cells[book.col_offset],
			      newinfo, pd, NULL, &time_series, &label_strings);
	}
    }

    /* number of variables and observations for import dataset */
    newinfo->v = n_vars_from_col(totcols, blank_col, book.col_offset,
				 time_series || label_strings);
    newinfo->n = nrows - 1 - book.row_offset;
    fprintf(stderr, "newinfo->v = %d, newinfo->n = %d\n",
	    newinfo->v, newinfo->n);

    /* create import dataset */
    err = start_new_Z(&newZ, newinfo, 0);
    if (err) {
	goto getout;
    }

    set_all_missing(newZ, newinfo);

    if (time_series) {
	ntodate_full(newinfo->endobs, newinfo->n - 1, newinfo);
    } else {
	strcpy(newinfo->stobs, "1");
	sprintf(newinfo->endobs, "%d", newinfo->n);
	newinfo->sd0 = 1.0;
	newinfo->pd = 1;
	newinfo->structure = CROSS_SECTION;
    } 

    /* OK: actually populate the dataset */
    transcribe_data(newZ, newinfo, book.row_offset, totcols,
		    book.col_offset + (time_series || label_strings),
		    blank_col);

    if (fix_varname_duplicates(newinfo)) {
	pputs(prn, _("warning: some variable names were duplicated\n"));
    }

    if (label_strings) {
	dataset_allocate_obs_markers(newinfo);
	if (newinfo->S != NULL) {
	    i = book.col_offset;
	    for (t=0; t<newinfo->n; t++) {
		int ts = t + 1 + book.row_offset;

		strncat(newinfo->S[t], rows[ts].cells[i] + 1, OBSLEN - 1);
	    }
	}
    }

    if (*pZ == NULL) {
	*pZ = newZ;
	*pdinfo = *newinfo;
    } else {
	err = merge_data(pZ, pdinfo, newZ, newinfo, prn);
    }

 getout:
    
    free(blank_col);
    wbook_free(&book);
    free_sheet();

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#ifdef WIN32
    if (fdb != NULL) {
	fclose(fdb);
    }
#endif

    return err;
}  

/*
 *  Copyright (c) by Allin Cottrell 2002
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
  Mostly borrowed from Gnumeric's Excel importer, written by Michael Meeks
*/

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "libgretl.h"

#include "importer.h"
#include "biff.h"

typedef struct _BiffBoundsheetData BiffBoundsheetData;

struct _BiffBoundsheetData
{
    guint16 index;
    guint32 streamStartPos;
    MsBiffFileType type;
    char *name;
};

BiffQuery *ms_biff_query_new (MsOleStream *ptr)
{
    BiffQuery *bq;

    if (ptr == NULL) return NULL;

    bq = g_new0 (BiffQuery, 1);
    bq->opcode        = 0;
    bq->length        = 0;
    bq->data_malloced = 0;
    bq->pos           = ptr;

    return bq;
}

/**
 * Returns 0 if has hit end
 **/
int ms_biff_query_next (BiffQuery *bq)
{
    guint8 tmp[4];
    int ans = 1;

    if (!bq || bq->pos->position >= bq->pos->size) {
	return 0;
    }

    if (bq->data_malloced) {
	g_free (bq->data);
	bq->data = 0;
	bq->data_malloced = 0;
    }

    bq->streamPos = bq->pos->position;
    if (!bq->pos->read_copy (bq->pos, tmp, 4)) {
	return 0;
    }

    bq->opcode = MS_OLE_GET_GUINT16(tmp);
    bq->length = MS_OLE_GET_GUINT16(tmp+2);
    bq->ms_op  = (bq->opcode >> 8);
    bq->ls_op  = (bq->opcode & 0xff);

    if (bq->length > 0 &&
	!(bq->data = bq->pos->read_ptr(bq->pos, bq->length))) {
	bq->data = g_new0 (guint8, bq->length);
	if (!bq->pos->read_copy(bq->pos, bq->data, bq->length)) {
	    ans = 0;
	    g_free(bq->data);
	    bq->data = NULL;
	    bq->length = 0;
	} else {
	    bq->data_malloced = 1;
	}
    }

    if (!bq->length) {
	bq->data = 0;
	return 1;
    }

    return ans;
}

void ms_biff_query_destroy (BiffQuery *bq)
{
    if (bq != NULL) {
	if (bq->data_malloced) {
	    g_free(bq->data);
	    bq->data = 0;
	    bq->data_malloced = 0;
	}
	g_free(bq);
    }
}

static void
get_xtn_lens (guint32 *pre_len, guint32 *end_len, const guint8 *ptr, 
	      gboolean ext_str, gboolean rich_str)
{
    *end_len = 0;
    *pre_len = 0;

    if (rich_str) { /* The data for this appears after the string */
	guint16 formatting_runs = MS_OLE_GET_GUINT16 (ptr);
	static int warned = FALSE;

	(*end_len) += formatting_runs * 4; /* 4 bytes per */
	(*pre_len) += 2;
	ptr        += 2;

	if (!warned) {
	    printf ("FIXME: rich string support unimplemented:"
		    "discarding %d runs\n", formatting_runs);
	}
	warned = TRUE;
    }
    if (ext_str) { /* NB this data always comes after the rich_str data */
	guint32 len_ext_rst = MS_OLE_GET_GUINT32 (ptr); /* A byte length */
	static int warned = FALSE;

	(*end_len) += len_ext_rst;
	(*pre_len) += 4;

	if (!warned) {
	    printf ("FIXME: extended string support unimplemented:"
		    "ignoring %u bytes\n", len_ext_rst);
	}
	warned = TRUE;
    }
}

static char *
get_chars (char const *ptr, guint length, gboolean high_byte)
{
    char *ans;
    guint32 lp;

    if (high_byte) {
	wchar_t *wc = g_new(wchar_t, length + 2);
	size_t retlength;

	ans = g_new(char, (length + 2) * 8);

	for (lp = 0; lp < length; lp++) {
	    guint16 c = MS_OLE_GET_GUINT16(ptr);

	    ptr += 2;
	    wc[lp] = c;
	}

	retlength = length;
	g_free (wc);
	if (retlength == (size_t)-1) {
	    retlength = 0;
	}
	ans[retlength] = 0;
	ans = g_realloc(ans, retlength + 2);
    } else {
	size_t outbytes = (length + 2) * 8;

	ans = g_new(char, outbytes + 1);
	for (lp = 0; lp < length; lp++) {
	    unsigned u = *ptr++;

	    ans[lp] = (u < 128)? u : '_';
	}
	ans[lp] = 0;
    }

    return ans;
}

static gboolean
biff_string_get_flags (const guint8 *ptr, gboolean *word_chars,
		       gboolean *extended, gboolean *rich)
{
    guint8 header;

    header = MS_OLE_GET_GUINT8(ptr);

    if (((header & 0xf2) == 0)) {
	*word_chars = (header & 0x1) != 0;
	*extended   = (header & 0x4) != 0;
	*rich       = (header & 0x8) != 0;
	return TRUE;
    } else { 
	*word_chars = 0;
	*extended   = 0;
	*rich       = 0;
	return FALSE;
    }
}

static char *
biff_get_text (guint8 const *pos, guint32 length, guint32 *byte_length)
{
    char *ans;
    const guint8 *ptr;
    guint32 byte_len;
    gboolean header;
    gboolean high_byte;
    gboolean ext_str;
    gboolean rich_str;
    guint32 pre_len, end_len;

    if (byte_length == NULL) {
	byte_length = &byte_len;
    }

    *byte_length = 0;

    if (!length) {
	return 0;
    }

    header = biff_string_get_flags (pos, &high_byte, &ext_str, &rich_str);
    if (header) {
	ptr = pos + 1;
	(*byte_length)++;
    } else {
	ptr = pos;
    }

    get_xtn_lens (&pre_len, &end_len, ptr, ext_str, rich_str);
    ptr += pre_len;
    (*byte_length) += pre_len + end_len;

    if (!length) {
	ans = g_new(char, 2);
	g_warning("Warning unterminated string floating");
    } else {
	(*byte_length) += (high_byte ? 2 : 1) * length;
	ans = get_chars((char *) ptr, length, high_byte);
    }

    return ans;
}

static BiffBoundsheetData *
biff_boundsheet_data_new (BiffQuery *q, MsBiffVersion ver)
{
    BiffBoundsheetData *ans = NULL;
    guint32 startpos;

    if (ver != MS_BIFF_V5 && ver != MS_BIFF_V7 && ver != MS_BIFF_V8) {
	printf ("Unknown BIFF Boundsheet spec. Assuming same as Biff7\n");
	ver = MS_BIFF_V7;
    }

    startpos = MS_OLE_GET_GUINT32(q->data);
    if (!(startpos == MS_OLE_GET_GUINT32(q->data))) {
	return NULL;
    }

    if (MS_OLE_GET_GUINT8(q->data + 4) == 0 &&
	((MS_OLE_GET_GUINT8(q->data + 5)) & 0x3) == 0) { /* worksheet, visible */
	ans = g_malloc(sizeof *ans);
	ans->streamStartPos = startpos;
	ans->name = biff_get_text (q->data + 7, 
				   MS_OLE_GET_GUINT8(q->data + 6), NULL);
    }

    return ans; 
}

static MsBiffBofData *ms_biff_bof_data_new (BiffQuery *q)
{
    MsBiffBofData *ans = g_new (MsBiffBofData, 1);

    if ((q->opcode & 0xff) == BIFF_BOF && (q->length >= 4)) {

	switch (q->opcode >> 8) {
	case 0: ans->version = MS_BIFF_V2;
	    break;
	case 2: ans->version = MS_BIFF_V3;
	    break;
	case 4: ans->version = MS_BIFF_V4;
	    break;
	case 8: {
	    switch (MS_OLE_GET_GUINT16 (q->data)) {
	    case 0x0600: ans->version = MS_BIFF_V8;
		break;
	    case 0x500: 
		ans->version = MS_BIFF_V7;
		break;
	    default:
		printf ("Unknown BIFF sub-number in BOF %x\n", q->opcode);
		ans->version = MS_BIFF_V_UNKNOWN;
	    }
	    break;
	}

	default:
	    printf ("Unknown BIFF number in BOF %x\n", q->opcode);
	    ans->version = MS_BIFF_V_UNKNOWN;
	    printf ("Biff version %d\n", ans->version);
	}

	switch (MS_OLE_GET_GUINT16(q->data + 2)) {
	case 0x0005: ans->type = MS_BIFF_TYPE_Workbook; break;
	case 0x0006: ans->type = MS_BIFF_TYPE_VBModule; break;
	case 0x0010: ans->type = MS_BIFF_TYPE_Worksheet; break;
	case 0x0020: ans->type = MS_BIFF_TYPE_Chart; break;
	case 0x0040: ans->type = MS_BIFF_TYPE_Macrosheet; break;
	case 0x0100: ans->type = MS_BIFF_TYPE_Workspace; break;
	default:
	    ans->type = MS_BIFF_TYPE_Unknown;
	    printf ("Unknown BIFF type in BOF %x\n", 
		    MS_OLE_GET_GUINT16 (q->data + 2));
	    break;
	}
    } else {
	printf ("Not a BOF !\n");
	ans->version = MS_BIFF_V_UNKNOWN;
	ans->type = MS_BIFF_TYPE_Unknown;
    }

    return ans;
}

static void
ms_biff_bof_data_destroy (MsBiffBofData *data)
{
    g_free(data);
}

static void
ms_excel_read_bof (BiffQuery *q,
		   MsBiffBofData **version)
{
    /* The first BOF seems to be OK, the rest lie ? */
    MsBiffVersion vv = MS_BIFF_V_UNKNOWN;
    MsBiffBofData *ver = *version;

    if (ver) {
	vv = ver->version;
	ms_biff_bof_data_destroy(ver);
    }

    *version = ver = ms_biff_bof_data_new(q);
    if (vv != MS_BIFF_V_UNKNOWN) {
	ver->version = vv;
    }

    if (ver->type == MS_BIFF_TYPE_Workbook) {
	if (ver->version >= MS_BIFF_V8) {
	    guint32 ver = MS_OLE_GET_GUINT32 (q->data + 4);
	    if (ver == 0x4107cd18) {
		printf ("Excel 2000 ?\n");
	    } else {
		printf ("Excel 97 +\n");
	    }
	} else if (ver->version >= MS_BIFF_V7)
	    printf ("Excel 95\n");
	else if (ver->version >= MS_BIFF_V5)
	    printf ("Excel 5.x\n");
	else if (ver->version >= MS_BIFF_V4)
	    printf ("Excel 4.x\n");
	else if (ver->version >= MS_BIFF_V3)
	    printf ("Excel 3.x\n");
	else if (ver->version >= MS_BIFF_V2)
	    printf ("Excel 2.x\n");
    } else if (ver->type == MS_BIFF_TYPE_Worksheet) {
	printf ("Got worksheet\n");
    } else if (ver->type == MS_BIFF_TYPE_Chart) {
	printf ("Chart.\n");
    } else if (ver->type == MS_BIFF_TYPE_VBModule ||
	       ver->type == MS_BIFF_TYPE_Macrosheet) {
	printf ("VB Module or Macrosheet.\n");

	while (ms_biff_query_next (q) && q->opcode != BIFF_EOF)
	    ;
	if (q->opcode != BIFF_EOF)
	    g_warning ("EXCEL: file format error.  Missing BIFF_EOF");
    } else {
	printf ("Unknown BOF (%x)\n", ver->type);
    }

    return;
}

static void
ms_excel_read_workbook (MsOle *file, BiffBoundsheetData ***bounds,
			int *nsheets)
{
    MsOleStream *stream;
    MsOleErr result;
    BiffQuery *q;
    MsBiffBofData *ver = NULL;
    char *problem_loading = NULL;

    result = ms_ole_stream_open (&stream, file, "/", "workbook", 'r');
    if (result != MS_OLE_ERR_OK) {
	ms_ole_stream_close (&stream);

	result = ms_ole_stream_open (&stream, file, "/", "book", 'r');
	if (result != MS_OLE_ERR_OK) {
	    ms_ole_stream_close (&stream);
	    fprintf (stderr, _("No book or workbook streams found."));
	    return;
	}
    }

    q = ms_biff_query_new(stream);

    while (problem_loading == NULL && ms_biff_query_next(q)) {

	if (0x1 == q->ms_op) {
	    switch (q->opcode) {
	    case BIFF_DSF:
	    case BIFF_XL9FILE:
	    case BIFF_REFRESHALL:
	    case BIFF_USESELFS:
	    case BIFF_TABID:
	    case BIFF_PROT4REV:
	    case BIFF_PROT4REVPASS:
	    case BIFF_CODENAME:
	    case BIFF_SUPBOOK:
		break;
	    default:
		fprintf(stderr, "Got unexpected BIFF token 0x%x\n", q->opcode);
	    }
	    continue;
	}

	switch (q->ls_op) {
	case BIFF_BOF:
#ifdef EDEBUG
	    fprintf(stderr, "Got BIFF_BOF\n");
#endif
	    ms_excel_read_bof(q, &ver);
	    break;

	case BIFF_EOF: 
#ifdef EDEBUG
	    fprintf(stderr, "Got BIFF_EOF\n");
#endif
	    break;

	case BIFF_BOUNDSHEET: {
	    BiffBoundsheetData *ans;

#ifdef EDEBUG
	    fprintf(stderr, "Got BIFF_BOUNDSHEET\n");
#endif
	    ans = biff_boundsheet_data_new(q, ver->version);
	    if (ans != NULL) {
		*bounds = g_realloc(*bounds, (*nsheets + 1) * sizeof **bounds);
		(*bounds)[*nsheets] = ans;
		*nsheets += 1;
	    }
	    break;
	}
	case BIFF_PALETTE:
	    break;

	case BIFF_FONT:
	    break;

	case BIFF_XF_OLD:
	case BIFF_XF:
	    break;

	case BIFF_SST:
	case BIFF_EXTSST: 
	    break;

	case BIFF_EXTERNSHEET: /* See: S59D82.HTM */
	    /* fprintf(stderr, "Got BIFF_EXTERNSHEET\n"); */
	    break;

	case BIFF_FORMAT: 
	case BIFF_EXTERNCOUNT:
	case BIFF_CODEPAGE: 
	case BIFF_OBJPROTECT:
	case BIFF_PROTECT:
	case BIFF_PASSWORD:
	    break;

	case BIFF_FILEPASS:
	    /* All records after this are encrypted */
	    problem_loading = g_strdup (_("Password protected workbooks "
					  "are not supported yet."));
	    break;

	case BIFF_STYLE:
	case BIFF_WINDOWPROTECT:
	    break;

	case BIFF_EXTERNNAME:
#ifdef EDEBUG
	    fprintf(stderr, "Got BIFF_EXTERNNAME\n");
#endif
	    break;

	case BIFF_NAME:
#ifdef EDEBUG
	    fprintf(stderr, "Got BIFF_NAME\n");
#endif
	    break;

	case BIFF_WRITEACCESS:
	case BIFF_HIDEOBJ:
	case BIFF_FNGROUPCOUNT:
	case BIFF_MMS:
	case BIFF_OBPROJ:
	case BIFF_BOOKBOOL:
	case BIFF_COUNTRY:
	case BIFF_INTERFACEHDR:
	case BIFF_INTERFACEEND:
	case BIFF_1904: /* 0, NOT 1 */
	case BIFF_WINDOW1:
	case BIFF_SELECTION: /* 0, NOT 10 */
	    break;
	case BIFF_DIMENSIONS:	/* 2, NOT 1,10 */
#ifdef EDEBUG
	    fprintf(stderr, "Got BIFF_DIMENSIONS\n");
#endif
	    /* ms_excel_biff_dimensions (q, wb); */
	    break;
	case BIFF_OBJ: /* See: ms-obj.c and S59DAD.HTM */
	case BIFF_SCL:
	    break;

	case BIFF_MS_O_DRAWING:
	case BIFF_MS_O_DRAWING_GROUP:
	case BIFF_MS_O_DRAWING_SELECTION:
	case BIFF_ADDMENU:
	    break;

	default:
	    /* fprintf(stderr, "ms_excel_unexpected_biff\n"); */
	    break;
	}
    }

    ms_biff_query_destroy(q);
    if (ver) {
	ms_biff_bof_data_destroy(ver);
    }
    ms_ole_stream_close(&stream);
}

/* public interface */

int excel_book_get_info (const char *fname, wbook *book)
{
    MsOleErr ole_error;
    MsOle *f;
    BiffBoundsheetData **bounds = NULL;
    int i, nsheets = 0;

    ole_error = ms_ole_open(&f, fname);

    if (ole_error != MS_OLE_ERR_OK) {
	char const *msg = (ole_error == MS_OLE_ERR_INVALID ||
			   ole_error == MS_OLE_ERR_FORMAT) ?
	    _("This file is not an 'OLE' file -- it may be too "
	      "old for gretl to read\n")
	    : _("Unexpected error reading the file\n");
	ms_ole_destroy (&f);
	fprintf(stderr, msg);
	return 1;
    }

    ms_excel_read_workbook(f, &bounds, &nsheets);
    ms_ole_destroy(&f);

    if (nsheets == 0 || bounds == NULL)
	return 1;

    book->sheetnames = g_malloc(nsheets * sizeof (char *));
    if (book->sheetnames == NULL) return 1;

    book->byte_offsets = g_malloc(nsheets * sizeof (unsigned));
    if (book->byte_offsets == NULL) return 1;

    book->nsheets = nsheets;

    for (i=0; i<nsheets; i++) {
	book->sheetnames[i] = (bounds[i])->name;
	book->byte_offsets[i] = (bounds[i])->streamStartPos;
	g_free(bounds[i]);
    }

    g_free(bounds);

    return 0;    
} 



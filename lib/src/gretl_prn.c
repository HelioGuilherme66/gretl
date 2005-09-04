/*
 *  Copyright (c) 2005 by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#include "libgretl.h"
#include <stdarg.h>

struct PRN_ {
    FILE *fp;
    FILE *fpaux;    
    char *buf;
    size_t bufsize;
    PrnFormat format;
    int fixed;
};

#undef PRN_DEBUG

/**
 * gretl_print_destroy:
 * @prn: pointer to gretl printing struct.
 *
 * Close a gretl printing struct and free any associated resources.
 */

void gretl_print_destroy (PRN *prn)
{
    int fpdup;

    if (prn == NULL) {
	return;
    }

    fpdup = (prn->fp == prn->fpaux);

    if (prn->fp != NULL && prn->fp != stdout && prn->fp != stderr) {
	fclose(prn->fp);
    }

    if (!fpdup && prn->fpaux != NULL && 
	prn->fpaux != stdout && prn->fpaux != stderr) {
	fclose(prn->fpaux);
    }

    if (prn->buf != NULL) {
#if PRN_DEBUG
  	fprintf(stderr, "gretl_print_destroy: freeing buffer at %p\n", 
		(void *) prn->buf); 
#endif
	free(prn->buf);
    }

    free(prn);
}

static PRN *real_gretl_print_new (PrnType ptype, 
				  const char *fname,
				  char *buf)
{
    PRN *prn = malloc(sizeof *prn);

    if (prn == NULL) {
	fprintf(stderr, _("gretl_prn_new: out of memory\n"));
	return NULL;
    }

#if PRN_DEBUG
    fprintf(stderr, "real_gretl_print_new: type=%d, fname='%s', buf=%p\n",
	    ptype, fname, (void *) buf);
#endif    

    prn->fp = NULL;
    prn->fpaux = NULL;
    prn->buf = NULL;
    prn->bufsize = 0;
    prn->format = GRETL_FORMAT_TXT;
    prn->fixed = 0;

    if (ptype == GRETL_PRINT_FILE) {
	prn->fp = gretl_fopen(fname, "w");
	if (prn->fp == NULL) {
	    fprintf(stderr, _("gretl_prn_new: couldn't open %s\n"), fname);
	    free(prn);
	    prn = NULL;
	}
    } else if (ptype == GRETL_PRINT_STDOUT) {
	prn->fp = stdout;
    } else if (ptype == GRETL_PRINT_STDERR) {
	prn->fp = stderr;
    } else if (ptype == GRETL_PRINT_BUFFER) {
	if (buf != NULL) {
	    prn->buf = buf;
	    prn->fixed = 1;
#if PRN_DEBUG
	    fprintf(stderr, "prn with fixed buffer\n");
#endif
	} else {
	    int p = pprintf(prn, "@init");

	    if (p < 0) {
		fprintf(stderr, _("gretl_prn_new: out of memory\n"));
		free(prn);
		prn = NULL;
	    }
#if PRN_DEBUG
	    fprintf(stderr, "new prn buffer, p = %d\n", p);
#endif
	}
    }

    return prn;
}

/**
 * gretl_print_new:
 * @ptype: code indicating the desired printing mode.
 * 
 * Create and initialize a gretl printing struct. If @ptype
 * is %GRETL_PRINT_BUFFER, output will go to an automatically
 * resized buffer; if @ptype is %GRETL_PRINT_STDOUT or
 * %GRETL_PRINT_STDERR output goes to %stdout or %stderr
 * respectively.
 *
 * If you want a named file associated with the struct, use 
 * #gretl_print_new_with_filename instead; if you want to
 * attach a fixed, pre-allocated text buffer, use
 * #gretl_print_new_with_buffer.
 *
 * Returns: pointer to newly created struct, or %NULL on failure.
 */

PRN *gretl_print_new (PrnType ptype)
{
    if (ptype == GRETL_PRINT_FILE) {
	fprintf(stderr, "gretl_print_new: needs a filename\n");
	return NULL;
    }

#if PRN_DEBUG
    fprintf(stderr, "gretl_print_new() called, type = %d\n", ptype);
#endif

    return real_gretl_print_new(ptype, NULL, NULL);
}

/**
 * gretl_print_new_with_filename:
 * @fname: name of the file to be opened for writing.
 * 
 * Create and initialize a gretl printing struct, with
 * output directed to the named file.  
 *
 * Returns: pointer to newly created struct, or %NULL on failure.
 */

PRN *gretl_print_new_with_filename (const char *fname)
{
    if (fname == NULL) {
	fprintf(stderr, _("gretl_prn_new: Must supply a filename\n"));
	return NULL;
    }

    return real_gretl_print_new(GRETL_PRINT_FILE, fname, NULL);
}

/**
 * gretl_print_new_with_buffer:
 * @buf: pre-allocated text buffer.
 * 
 * Creates and initializes a gretl printing struct, with
 * fixed text buffer @buf.  Fails if @buf is %NULL.  This is a
 * convenience function: you can't use #pprintf, #pputs or
 * #pputc with a printing struct obtained in this way.
 *
 * Note that @buf will be freed if and when #gretl_print_destroy
 * is called on the #PRN pointer obtained.
 *
 * Returns: pointer to newly created struct, or %NULL on failure.
 */

PRN *gretl_print_new_with_buffer (char *buf)
{
    if (buf == NULL) {
	return NULL;
    } else {
	return real_gretl_print_new(GRETL_PRINT_BUFFER, NULL, buf);
    }
}

/**
 * gretl_print_reset_buffer:
 * @prn: printing struct to operate on.
 * 
 * Write a NUL byte to the start of the buffer associated
 * with @prn, if any.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_print_reset_buffer (PRN *prn)
{
    int err = 0;

    if (prn != NULL && prn->buf != NULL) {
	*prn->buf = '\0';
    } else {
	err = 1;
    }

    return err;
}

/**
 * gretl_print_get_buffer:
 * @prn: printing struct.
 * 
 * Obtain a pointer to the buffer associated
 * with @prn, if any.  This pointer must not be
 * modified in any way.
 *
 * Returns: the buffer, or %NULL on failure.
 */

const char *gretl_print_get_buffer (PRN *prn)
{
    if (prn != NULL) {
	return prn->buf;
    } else {
	return NULL;
    }
}

/**
 * gretl_print_set_format:
 * @prn: printing struct.
 * @format: desired format code.
 * 
 * Sets the printing format on @prn.
 */

void gretl_print_set_format (PRN *prn, PrnFormat format)
{
    if (prn != NULL) {
	if (format == GRETL_FORMAT_RTF) {
	    format = GRETL_FORMAT_RTF | GRETL_FORMAT_DOC;
	}
	prn->format = format;
    }
}

/**
 * gretl_print_unset_doc_format:
 * @prn: printing struct.
 * 
 * Toggles the %GRETL_FORMAT_DOC flag on @prn.
 */

void gretl_print_toggle_doc_flag (PRN *prn)
{
    if (prn != NULL) {
	prn->format ^= GRETL_FORMAT_DOC;
    }
}

static int realloc_prn_buffer (PRN *prn, size_t blen)
{
    char *tmp;
    int err = 0;

#if PRN_DEBUG
    fprintf(stderr, "%d bytes left\ndoing realloc(%p, %d)\n",
	    prn->bufsize - blen, prn->buf, 2 * prn->bufsize);
#endif
    prn->bufsize *= 2; 

    tmp = realloc(prn->buf, prn->bufsize); 

    if (tmp == NULL) {
	err = 1;
    } else {
	prn->buf = tmp;
#if PRN_DEBUG
	fprintf(stderr, "realloc: prn->buf is %d bytes at %p\n",
		prn->bufsize, (void *) prn->buf);
#endif
    }

    memset(prn->buf + blen, 0, 1);

    return err;
}

/**
 * pprintf:
 * @prn: gretl printing struct.
 * @template: as in printf().
 * @Varargs: arguments to be printed.
 *
 * Multi-purpose printing function: can output to stream, to buffer
 * or to nowhere (silently discarding the output), depending on
 * how @prn was initialized.  It's not advisable to use this function
 * for large chunks of text: use #pputs instead.
 * 
 * Returns: 0 on successful completion, 1 on memory allocation
 * failure.
 * 
 */

int pprintf (PRN *prn, const char *template, ...)
{
    va_list args;
    size_t blen;
    int plen = 0;

    if (prn == NULL || prn->fixed) return 0;

    if (prn->fp != NULL) {
	va_start(args, template);
	plen = vfprintf(prn->fp, template, args);
	va_end(args);
	return plen;
    }

    if (strncmp(template, "@init", 5) == 0) {
	prn->bufsize = 2048;
	prn->buf = malloc(prn->bufsize);
#if PRN_DEBUG
  	fprintf(stderr, "pprintf: malloc'd %d bytes at %p\n", prn->bufsize,  
		(void *) prn->buf); 
#endif
	if (prn->buf == NULL) {
	    return -1;
	}
	memset(prn->buf, 0, 1);
	return 0;
    }

    if (prn->buf == NULL) return 0;

    blen = strlen(prn->buf);

    if (prn->bufsize - blen < 1024) {
	if (realloc_prn_buffer(prn, blen)) {
	    return -1;
	}
    }

    va_start(args, template);
#if PRN_DEBUG
    fprintf(stderr, "printing at %p\n", (void *) (prn->buf + blen));
#endif
    plen = vsprintf(prn->buf + blen, template, args);
    va_end(args);
#if PRN_DEBUG
    fprintf(stderr, "printed %d byte(s)\n", strlen(prn->buf) - blen);
#endif

    return plen;
}

/**
 * pputs:
 * @prn: gretl printing struct.
 * @s: constant string to print
 * 
 * Returns: 0 on successful completion, 1 on memory allocation
 * failure.
 */

int pputs (PRN *prn, const char *s)
{
    size_t blen;
    int slen, bytesleft;

    if (prn == NULL || prn->fixed) return 0;

    slen = strlen(s);

    if (prn->fp != NULL) {
	fputs(s, prn->fp);
	return slen;
    }

    if (prn->buf == NULL) return 0;

    blen = strlen(prn->buf);

    bytesleft = prn->bufsize - blen;

    while (prn->bufsize - blen < 1024 || bytesleft <= slen) {
	if (realloc_prn_buffer(prn, blen)) {
	    return -1;
	}
	bytesleft = prn->bufsize - blen;
    }	

#if PRN_DEBUG
    fprintf(stderr, "pputs: bufsize=%d, blen=%d, bytesleft=%d, copying %d bytes\n", 
	    (int) prn->bufsize, (int) blen, bytesleft, slen);
#endif

    strcpy(prn->buf + blen, s);

    return slen;
}

/**
 * pputc:
 * @prn: gretl printing struct.
 * @c: character to print
 * 
 * Returns: 0 on successful completion, 1 on memory allocation
 * failure.
 */

int pputc (PRN *prn, int c)
{
    size_t blen;

    if (prn == NULL || prn->fixed) return 0;

    if (prn->fp != NULL) {
	fputc(c, prn->fp);
	return 1;
    }

    if (prn->buf == NULL) return 0;

    blen = strlen(prn->buf);

    if (prn->bufsize - blen < 1024) {
	if (realloc_prn_buffer(prn, blen)) {
	    return -1;
	}
    }

#if PRN_DEBUG
    fprintf(stderr, "pputc: adding char at %p\n", (void *) (prn->buf + blen));
#endif
    prn->buf[blen] = c;
    prn->buf[blen + 1] = '\0';

    return 1;
}

/**
 * gretl_prn_newline:
 * @prn: gretl printing struct.
 * 
 * Print a line break, in the mode appropriate to the
 * format of @prn (plain text, TeX or RTF).
 */

void gretl_prn_newline (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

/**
 * gretl_print_flush_stream:
 * @prn: gretl printing struct.
 * 
 * If the output of @prn is directed to a stream, flush
 * that stream.
 */

void gretl_print_flush_stream (PRN *prn)
{
    if (prn != NULL && prn->fp != NULL) {
	fflush(prn->fp);
    }
}

/**
 * printing_to_standard_stream:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the output of @prn is directed to one of the
 * standard streams, %stdout or %stderr, else 0.
 */

int printing_to_standard_stream (PRN *prn)
{
    int ret = 0;

    if (prn != NULL && (prn->fp == stdout || prn->fp == stderr)) {
	ret = 1;
    }

    return ret;
}

/**
 * printing_is_redirected:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the output of @prn has been redirected
 * relative to its original setting, else 0.
 */

int printing_is_redirected (PRN *prn)
{
    int ret = 0;

    if (prn != NULL) {
	if (prn->fpaux != NULL || 
	    (prn->fp != NULL && prn->buf != NULL)) {
	    ret = 1;
	}
    }

    return ret;
}

/**
 * print_start_redirection:
 * @prn: gretl printing struct.
 * @fp: stream to which output should be redirected.
 * 
 * Redirects output of @prn to @fp.
 *
 * Returns: 0 on success, 1 on error.
 */

int print_start_redirection (PRN *prn, FILE *fp)
{
    int err = 0;

    if (prn != NULL) {
	/* save the current stream */
	prn->fpaux = prn->fp;
	/* hook output to specified file */
	prn->fp = fp;
    } else {
	err = 1;
    }

    return err;
}

/**
 * print_end_redirection:
 * @prn: gretl printing struct.
 * 
 * If the output of @prn has been diverted relative to
 * its original setting, terminate the redirection.
 *
 * Returns: 0 on success, 1 on error.
 */

int print_end_redirection (PRN *prn)
{
    int err = 0;

    if (prn != NULL && prn->fp != NULL) {
	fclose(prn->fp);
	prn->fp = prn->fpaux;
	prn->fpaux = NULL;
    } else {
	err = 1;
    }
    
    return err;
}

/**
 * plain_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is plain text, else 0.
 */

int plain_format (PRN *prn)
{
    return (prn != NULL && prn->format == GRETL_FORMAT_TXT);
}

/**
 * rtf_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is RTF, else 0.
 */

int rtf_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_RTF));
}

/**
 * rtf_doc_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is RTF, and the
 * RTF preamble should be printed, else 0.
 */

int rtf_doc_format (PRN *prn)
{
    return (prn != NULL && 
	    (prn->format & GRETL_FORMAT_RTF) &&
	    (prn->format & GRETL_FORMAT_DOC));
}

/**
 * tex_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is TeX, else 0.
 */

int tex_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_TEX));
}

/**
 * tex_doc_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is TeX and a complete
 * TeX document is being produced, else 0.
 */

int tex_doc_format (PRN *prn)
{
    return (prn != NULL && 
	    (prn->format & GRETL_FORMAT_TEX) &&
	    (prn->format & GRETL_FORMAT_DOC));
}

/**
 * tex_eqn_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is TeX and a model
 * is being represented as an equation, rather than in
 * tabular form, else 0.
 */

int tex_eqn_format (PRN *prn)
{
    return (prn != NULL && 
	    (prn->format & GRETL_FORMAT_TEX) &&
	    (prn->format & GRETL_FORMAT_EQN));
}

/**
 * csv_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is CSV, else 0.
 */

int csv_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_CSV));
}

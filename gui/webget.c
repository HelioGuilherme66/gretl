/*
 *   Copyright (c) by Allin Cottrell
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

/* webget.c for gretl -- based on parts of GNU Wget */

#include "gretl.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>

#ifdef OS_WIN32
# include <winsock.h>
#else
# include <sys/socket.h>
# include <netdb.h>
# include <netinet/in.h>
# include <arpa/inet.h>
#endif /* OS_WIN32 */

#include "webget.h"

#ifndef errno
extern int errno;
#endif
#ifndef h_errno
extern int h_errno;
#endif

#define DEFAULT_HTTP_PORT 80
#define DEFAULT_FTP_PORT 21
#define MINVAL(x, y) ((x) < (y) ? (x) : (y))

#define TEXTHTML_S "text/html"
#define HTTP_ACCEPT "*/*"

/* Some status code validation macros: */
#define H_20X(x)        (((x) >= 200) && ((x) < 300))
#define H_PARTIAL(x)    ((x) == HTTP_STATUS_PARTIAL_CONTENTS)
#define H_REDIRECTED(x) (((x) == HTTP_STATUS_MOVED_PERMANENTLY)	\
			 || ((x) == HTTP_STATUS_MOVED_TEMPORARILY))

/* HTTP/1.0 status codes from RFC1945, provided for reference.  */
/* Successful 2xx.  */
#define HTTP_STATUS_OK			200
#define HTTP_STATUS_CREATED		201
#define HTTP_STATUS_ACCEPTED		202
#define HTTP_STATUS_NO_CONTENT		204
#define HTTP_STATUS_PARTIAL_CONTENTS	206

/* Redirection 3xx.  */
#define HTTP_STATUS_MULTIPLE_CHOICES	300
#define HTTP_STATUS_MOVED_PERMANENTLY	301
#define HTTP_STATUS_MOVED_TEMPORARILY	302
#define HTTP_STATUS_NOT_MODIFIED	304

/* Client error 4xx.  */
#define HTTP_STATUS_BAD_REQUEST		400
#define HTTP_STATUS_UNAUTHORIZED	401
#define HTTP_STATUS_FORBIDDEN		403
#define HTTP_STATUS_NOT_FOUND		404

/* Server errors 5xx.  */
#define HTTP_STATUS_INTERNAL		500
#define HTTP_STATUS_NOT_IMPLEMENTED	501
#define HTTP_STATUS_BAD_GATEWAY		502
#define HTTP_STATUS_UNAVAILABLE		503

enum uflags {
    URELATIVE     = 0x0001,      /* Is URL relative? */
    UNOPROTO      = 0x0002,      /* Is URL without a protocol? */
    UABS2REL      = 0x0004,      /* Convert absolute to relative? */
    UREL2ABS      = 0x0008       /* Convert relative to absolute? */
};

enum {
    HG_OK, 
    HG_ERROR, 
    HG_EOF
};

enum header_get_flags { 
    HG_NONE = 0,
    HG_NO_CONTINUATIONS = 0x2 
};

/* Flags for show_progress(). */
enum spflags { 
    SP_NONE, 
    SP_INIT, 
    SP_FINISH 
};

extern const char *version_string;

/* prototypes */
static char *time_str (time_t *tm);
static void rbuf_initialize (struct rbuf *rbuf, int fd);
static int rbuf_peek (struct rbuf *rbuf, char *store);
static size_t rbuf_flush (struct rbuf *rbuf, char *where, int maxsize);

static int header_get (struct rbuf *rbuf, char **hdr, 
		       enum header_get_flags flags);
static int header_strdup (const char *header, void *closure);
static int header_extract_number (const char *header, void *closure);
static int skip_lws (const char *string);
static int header_process (const char *header, const char *name,
			   int (*procfun) (const char *, void *),
			   void *arg);
static int parse_http_status_line (const char *line, 
				   const char **reason_phrase_ptr);
static int http_process_range (const char *hdr, void *arg);
static int http_process_none (const char *hdr, void *arg);
static int http_process_type (const char *hdr, void *arg);
static int numdigit (long a);
static char *herrmsg (int error);
static uerr_t gethttp (struct urlinfo *u, struct http_stat *hs, int *dt);
static uerr_t http_loop (struct urlinfo *u, int *dt);
static struct urlinfo *newurl (void);
static void freeurl (struct urlinfo *u, int complete);
static int get_contents (int fd, FILE *fp, char **getbuf, long *len, 
			 long expected, struct rbuf *rbuf);
static int store_hostaddress (unsigned char *where, const char *hostname);
static uerr_t make_connection (int *sock, char *hostname, 
			       unsigned short port);
static int iread (int fd, char *buf, int len);
static int iwrite (int fd, char *buf, int len);
static char *print_option (int opt);
static void destroy_progress (GtkWidget *widget, ProgressData *pdata);
static ProgressData *progress_window (void);

#ifdef OS_WIN32
static void ws_cleanup (void)
{
    WSACleanup();
}

/* ........................................................... */

int ws_startup (void)
{
    WORD requested;
    WSADATA data;

    requested = MAKEWORD(1, 1);

    if (WSAStartup(requested, &data)) {
	fprintf(stderr, "Couldn't find usable socket driver.\n");
	return 1;
    }

    if (LOBYTE (requested) < 1 || (LOBYTE (requested) == 1 &&
				   HIBYTE (requested) < 1)) {
	fprintf(stderr, "Couldn't find usable socket driver.\n");
	WSACleanup();
	return 1;
    }
    atexit(ws_cleanup);
    return 0;
}
#endif

/* ........................................................... */

static char *time_str (time_t *tm)
{
  static char tms[15];
  struct tm *ptm;
  time_t tim;

  *tms = '\0';
  tim = time (tm);
  if (tim == -1)
    return tms;
  ptm = localtime (&tim);
  sprintf (tms, "%02d:%02d:%02d", ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
  return tms;
}

/* http header functions -- based on Wget */

/* ........................................................... */

static int rbuf_peek (struct rbuf *rbuf, char *store)
{
    if (!rbuf->buffer_left) {
	int res;
	rbuf->buffer_pos = rbuf->buffer;
	rbuf->buffer_left = 0;
	res = iread (rbuf->fd, rbuf->buffer, sizeof rbuf->buffer);
	if (res <= 0)
	    return res;
	rbuf->buffer_left = res;
    }
    *store = *rbuf->buffer_pos;
    return 1;
}

/* ........................................................... */

static int header_get (struct rbuf *rbuf, char **hdr, 
		       enum header_get_flags flags)
{
    int i;
    int bufsize = 80;

    *hdr = mymalloc(bufsize);
    if (*hdr == NULL) return HG_ERROR;
    for (i = 0; 1; i++) {
	int res;
	if (i > bufsize - 1)
	    *hdr = g_realloc(*hdr, (bufsize <<= 1));
	res = RBUF_READCHAR(rbuf, *hdr + i);
	if (res == 1) {
	    if ((*hdr)[i] == '\n') {
		if (!((flags & HG_NO_CONTINUATIONS)
		      || i == 0
		      || (i == 1 && (*hdr)[0] == '\r'))) {
		    char next;
		    /* If the header is non-empty, we need to check if
		       it continues on to the other line.  We do that by
		       peeking at the next character.  */
		    res = rbuf_peek (rbuf, &next);
		    if (res == 0)
			return HG_EOF;
		    else if (res == -1)
			return HG_ERROR;
		    /*  If the next character is HT or SP, just continue.  */
		    if (next == '\t' || next == ' ')
			continue;
		}
		/* The header ends.  */
		(*hdr)[i] = '\0';
		/* Get rid of '\r'.  */
		if (i > 0 && (*hdr)[i - 1] == '\r')
		    (*hdr)[i - 1] = '\0';
		break;
	    }
	}
	else if (res == 0)
	    return HG_EOF;
	else
	    return HG_ERROR;
    }
    return HG_OK;
}

/* ........................................................... */

static int header_extract_number (const char *header, void *closure)
{
    const char *p = header;
    long result;

    for (result = 0; isdigit((unsigned char) *p); p++)
	result = 10 * result + (*p - '0');
    if (*p)
	return 0;

    *(long *)closure = result;
    return 1;
}

/* ........................................................... */

static int header_strdup (const char *header, void *closure)
{
    *(char **)closure = g_strdup(header);
    return 1;
}

/* ........................................................... */

static int skip_lws (const char *string)
{
    const char *p = string;

    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
	++p;
    return p - string;
}

/* ........................................................... */

static int header_process (const char *header, const char *name,
			   int (*procfun) (const char *, void *),
			   void *arg)
{
    while (*name && (tolower (*name) == tolower (*header)))
	++name, ++header;
    if (*name || *header++ != ':')
	return 0;

    header += skip_lws (header);

    return ((*procfun) (header, arg));
}

/* end Wget http header functions */

/* further functions from Wget's http.c */

/* ........................................................... */

static void rbuf_initialize (struct rbuf *rbuf, int fd)
{
    rbuf->fd = fd;
    rbuf->buffer_pos = rbuf->buffer;
    rbuf->buffer_left = 0;
}

/* ........................................................... */

static int parse_http_status_line (const char *line, 
				   const char **reason_phrase_ptr)
/* Parse the HTTP status line, which is of format:

   HTTP-Version SP Status-Code SP Reason-Phrase

   The function returns the status-code, or -1 if the status line is
   malformed.  The pointer to reason-phrase is returned in RP.  
*/
{
    int mjr, mnr, statcode;
    const char *p;

    *reason_phrase_ptr = NULL;

    if (strncmp (line, "HTTP/", 5) != 0)
	return -1;
    line += 5;

    /* Calculate major HTTP version.  */
    p = line;
    for (mjr = 0; isdigit((unsigned char) *line); line++)
	mjr = 10 * mjr + (*line - '0');
    if (*line != '.' || p == line)
	return -1;
    ++line;

    /* Calculate minor HTTP version.  */
    p = line;
    for (mnr = 0; isdigit((unsigned char) *line); line++)
	mnr = 10 * mnr + (*line - '0');
    if (*line != ' ' || p == line)
	return -1;
    /* Will accept only 1.0 and higher HTTP-versions.  The value of
       minor version can be safely ignored.  */
    if (mjr < 1)
	return -1;
    ++line;

    /* Calculate status code.  */
    if (!(isdigit((unsigned char) *line) && 
	  isdigit((unsigned char) line[1]) && 
	  isdigit((unsigned char) line[2])))
	return -1;
    statcode = 100 * (*line - '0') + 10 * (line[1] - '0') + (line[2] - '0');

    /* Set up the reason phrase pointer.  */
    line += 3;
    /* RFC2068 requires SPC here, but we allow the string to finish
     here, in case no reason-phrase is present.  */
    if (*line != ' ') {
	if (!*line)
	    *reason_phrase_ptr = line;
	else
	    return -1;
    }
    else
	*reason_phrase_ptr = line + 1;

    return statcode;
}

struct http_process_range_closure {
    long first_byte_pos;
    long last_byte_pos;
    long entity_length;
};

/* ........................................................... */

static int http_process_range (const char *hdr, void *arg)
/* Parse the `Content-Range' header and extract the information it
   contains.  Returns 1 if successful, -1 otherwise.  */
{
    struct http_process_range_closure *closure
	= (struct http_process_range_closure *)arg;
    long num;

    if (!strncasecmp (hdr, "bytes", 5)) {
	hdr += 5;
	hdr += skip_lws (hdr);
	if (!*hdr)
	    return 0;
    }
    if (!isdigit((unsigned char) *hdr))
	return 0;
    for (num = 0; isdigit((unsigned char) *hdr); hdr++)
	num = 10 * num + (*hdr - '0');
    if (*hdr != '-' || !isdigit((unsigned char) *(hdr + 1)))
	return 0;
    closure->first_byte_pos = num;
    ++hdr;
    for (num = 0; isdigit((unsigned char) *hdr); hdr++)
	num = 10 * num + (*hdr - '0');
    if (*hdr != '/' || !isdigit((unsigned char) *(hdr + 1)))
	return 0;
    closure->last_byte_pos = num;
    ++hdr;
    for (num = 0; isdigit((unsigned char) *hdr); hdr++)
	num = 10 * num + (*hdr - '0');
    closure->entity_length = num;
    return 1;
}

/* ........................................................... */

static int http_process_none (const char *hdr, void *arg)
/* Place 1 to ARG if the HDR contains the word "none", 0 otherwise.
   Used for `Accept-Ranges'.  */
{
    int *where = (int *)arg;

    if (strstr (hdr, "none"))
	*where = 1;
    else
	*where = 0;
    return 1;
}

/* ........................................................... */

static int http_process_type (const char *hdr, void *arg)
/* Place the malloc-ed copy of HDR hdr, to the first `;' to ARG.  */
{
    char **result = (char **)arg;
    char *p;

    *result = g_strdup(hdr);
    p = strrchr (hdr, ';');
    g_free(*result); /* added 03/25/01 */
    if (p) {
	int len = p - hdr;

	*result = mymalloc(len + 1);
	memcpy(*result, hdr, len);
	(*result)[len] = '\0';
    } else
	*result = g_strdup(hdr);
    return 1;
}

#define FREEHSTAT(x) do				\
{						\
  free((x).newloc);				\
  free((x).remote_time);			\
  free((x).error);				\
  (x).newloc = (x).remote_time = (x).error = NULL;	\
} while (0)

/* ........................................................... */

static int numdigit (long a)
{
    int res = 1;

    while ((a /= 10) != 0)
	++res;
    return res;
}

/* ........................................................... */

static char *herrmsg (int error)
{
    if (error == HOST_NOT_FOUND
	|| error == NO_RECOVERY
	|| error == NO_DATA
	|| error == NO_ADDRESS
	|| error == TRY_AGAIN)
	return "Host not found";
    else
	return "Unknown error";
}

/* ........................................................... */

static uerr_t gethttp (struct urlinfo *u, struct http_stat *hs, int *dt)
{
    char *request, *type, *command, *path;
    char *pragma_h, *useragent, *range, *remhost;
    char *all_headers = NULL;
    int sock, hcount, num_written, all_length, remport, statcode;
    long contlen, contrange;
    uerr_t err;
    FILE *fp;
    struct rbuf rbuf;

    hs->len = 0L;
    hs->contlen = -1;
    hs->res = -1;
    hs->newloc = NULL;
    hs->remote_time = NULL;
    hs->error = NULL;

    err = make_connection(&sock, u->host, u->port);

    switch (err) {
    case HOSTERR:
	sprintf(u->errbuf, "%s: %s.\n", u->host, herrmsg(h_errno));
	return HOSTERR;
	break;
    case CONSOCKERR:
	sprintf(u->errbuf, "socket: %s\n", strerror(errno));
	return CONSOCKERR;
	break;
    case CONREFUSED:
	sprintf(u->errbuf, "Connection to %s:%hu refused.\n", u->host, u->port);
	close(sock);
	return CONREFUSED;
    case CONERROR:
	sprintf(u->errbuf, "connect: %s\n", strerror(errno));
	close(sock);
	return CONERROR;
	break;
    case NOCONERROR:
	break;
    default:
	abort();
	break;
    } 

    if (u->filesave) { /* save output to file */
	fp = fopen(*(u->local), "wb");
	if (fp == NULL) {
	    close(sock);
	    free(all_headers);
	    return FOPENERR;
	}
    } else 
	fp = NULL; /* use local buffer instead */

    path = u->path;
    command = (*dt & HEAD_ONLY)? "HEAD" : "GET";
    if (*dt & SEND_NOCACHE)
	pragma_h = "Pragma: no-cache\r\n";
    else
	pragma_h = "";

    range = NULL;
    useragent = mymalloc(16); 
    sprintf(useragent, "gretl-%s", version_string);
#ifdef OS_WIN32
    strcat(useragent, "w");
#endif
    remhost = u->host;
    remport = u->port;

    request = mymalloc(strlen(command) + strlen(path)
		       + strlen(useragent)
		       + strlen(remhost) + numdigit(remport)
		       + strlen(HTTP_ACCEPT)
		       + strlen(pragma_h)
		       + 64);
    sprintf(request, "%s %s HTTP/1.0\r\n"
	    "User-Agent: %s\r\n"
	    "Host: %s:%d\r\n"
	    "Accept: %s\r\n"
	    "%s\r\n",
	    command, path, useragent, remhost, remport, HTTP_ACCEPT, 
	    pragma_h); 
    free(useragent);
    /* Send the request to server */
    num_written = iwrite(sock, request, strlen(request));
    free(request); /* moved from within following conditional, 03/25/01 */
    if (num_written < 0) {
	close(sock);
	return WRITEFAILED;
    }
    contlen = contrange = -1;
    type = NULL;
    statcode = -1;
    *dt &= ~RETROKF;

    /* Before reading anything, initialize rbuf */
    rbuf_initialize(&rbuf, sock);

    all_headers = NULL;
    all_length = 0;
    /* Header-fetching loop */
    hcount = 0;
    while (1) {
	char *hdr;
	int status;

	++hcount;
	/* Get the header.  */
	status = header_get(&rbuf, &hdr,
			    /* Disallow continuations for status line.  */
			    (hcount == 1 ? HG_NO_CONTINUATIONS : HG_NONE));
	/* Check for errors.  */
	if (status == HG_EOF && *hdr) {
	    free(hdr);
	    free(type);
	    free(hs->newloc);
	    free(all_headers);
	    close(sock);
	    return HEOF;
	} else if (status == HG_ERROR) {
	    free(hdr);
	    free(type);
	    free(hs->newloc);
	    free(all_headers);
	    close(sock);
	    return HERR;
	}

	/* Check for status line.  */
	if (hcount == 1) {
	    const char *error;

	    /* Parse the first line of server response.  */
	    statcode = parse_http_status_line (hdr, &error);
	    hs->statcode = statcode;
	    /* Store the descriptive response */
	    if (statcode == -1) { /* malformed response */
		/* A common reason for "malformed response" error is the
		   case when no data was actually received.  Handle this
		   special case.  */
		if (!*hdr)
		    hs->error = g_strdup("No data received");
		else
		    hs->error = g_strdup("Malformed status line");
		free(hdr);
		break;
	    }
	    else if (!*error)
		hs->error = g_strdup("(no description)");
	    else
		hs->error = g_strdup(error);
	    goto done_header;
	}

	/* Exit on empty header */
	if (!*hdr) {
	    free(hdr);
	    break;
	}
	/* Try getting content-length */
	if (contlen == -1)
	    if (header_process(hdr, "Content-Length", header_extract_number,
			       &contlen))
		goto done_header;
	/* Try getting content-type */
	if (!type)
	    if (header_process (hdr, "Content-Type", http_process_type, &type))
		goto done_header;
	/* Try getting location */
	if (!hs->newloc)
	    if (header_process (hdr, "Location", header_strdup, &hs->newloc))
		goto done_header;
	/* Try getting last-modified */
	if (!hs->remote_time)
	    if (header_process (hdr, "Last-Modified", header_strdup,
				&hs->remote_time))
		goto done_header;
	/* Check for accept-ranges header.  If it contains the word
	   `none', disable the ranges */
	if (*dt & ACCEPTRANGES) {
	    int nonep;
	    if (header_process (hdr, "Accept-Ranges", 
				http_process_none, &nonep)) {
		if (nonep)
		    *dt &= ~ACCEPTRANGES;
		goto done_header;
	    }
	}
	/* Try getting content-range */
	if (contrange == -1) {
	    struct http_process_range_closure closure;
	    if (header_process(hdr, "Content-Range", 
			       http_process_range, &closure)) {
		contrange = closure.first_byte_pos;
		goto done_header;
	    }
	}
    done_header:
	free(hdr);
    }

    /* 20x responses are counted among successful by default */
    if (H_20X (statcode))
	*dt |= RETROKF;

    if (type && !strncasecmp (type, TEXTHTML_S, strlen (TEXTHTML_S)))
	*dt |= TEXTHTML;
    else
	/* We don't assume text/html by default */
	*dt &= ~TEXTHTML;

    hs->contlen = contlen;

    /* Return if redirected */
    if (H_REDIRECTED (statcode) || statcode == HTTP_STATUS_MULTIPLE_CHOICES) {
	/* RFC2068 says that in case of the 300 (multiple choices)
	   response, the server can output a preferred URL through
	   `Location' header; otherwise, the request should be treated
	   like GET.  So, if the location is set, it will be a
	   redirection; otherwise, just proceed normally.  */
	if (statcode == HTTP_STATUS_MULTIPLE_CHOICES && !hs->newloc)
	    *dt |= RETROKF;
	else {
	    close(sock);
	    free(type);
	    free(all_headers);
	    return NEWLOCATION;
	}
    }

    free(type);
    type = NULL;	

    /* Return if we have no intention of further downloading */
    if (!(*dt & RETROKF) || (*dt & HEAD_ONLY)) {
	hs->len = 0L;
	hs->res = 0;
	free(type);
	free(all_headers);
	close(sock);
	return RETRFINISHED;
    }

    /* Get the contents of the document */
    hs->res = get_contents(sock, fp, u->local, &hs->len, 
			   (contlen != -1 ? contlen : 0), &rbuf);

    if (fp != NULL) fclose(fp);

    free(all_headers);
    close(sock);
    if (hs->res == -2)
	return FWRITEERR;
    return RETRFINISHED;
}

#define MAXTRY 5

/* ........................................................... */

static uerr_t http_loop (struct urlinfo *u, int *dt)
{
    int count;
    char *tms;
    uerr_t err;
    struct http_stat hstat;	/* HTTP status */

    count = 0;
    *dt = 0 | ACCEPTRANGES;
    /* THE loop */
    do {
	++count;
	tms = time_str(NULL);

	*dt &= ~HEAD_ONLY;
	*dt &= ~SEND_NOCACHE;

	/* Try fetching the document, or at least its head.  :-) */
	err = gethttp(u, &hstat, dt);
	/* Time?  */
	tms = time_str(NULL);

	switch (err) {
	case HERR: case HEOF: case CONSOCKERR: case CONCLOSED:
	case CONERROR: case READERR: case WRITEFAILED:
	case RANGEERR:
	    FREEHSTAT(hstat);
	    continue;
	    break;
	case HOSTERR: case CONREFUSED: case PROXERR: case AUTHFAILED:
	    FREEHSTAT(hstat);
	    return err;
	    break;
	case FWRITEERR: case FOPENERR:
	    FREEHSTAT(hstat);
	    return err;
	    break;
	case RETRFINISHED:
	    break;
	default:
	    abort();
	}
	if (!(*dt & RETROKF)) {
	    FREEHSTAT(hstat);
	    return WRONGCODE;
	}

	FREEHSTAT(hstat);

	if (hstat.len == hstat.contlen)
	    return RETROK;
	else if (hstat.res == 0) { /* No read error */
	    if (hstat.contlen == -1)  
		return RETROK;
	    else	
		continue;
	}
	else {		          /* now hstat.res can only be -1 */
	    if (hstat.contlen == -1)
		continue;
	}
	break;
    } while (count < MAXTRY);
    return TRYLIMEXC;
}

/* other utility functions from Wget */

/* ........................................................... */

static size_t rbuf_flush (struct rbuf *rbuf, char *where, int maxsize)
/* Flush RBUF's buffer to WHERE.  Flush MAXSIZE bytes at most.
   Returns the number of bytes actually copied.  If the buffer is
   empty, 0 is returned.  */
{
    if (!rbuf->buffer_left)
	return 0;
    else {
	size_t howmuch = MINVAL(rbuf->buffer_left, (unsigned)maxsize);

	if (where)
	    memcpy(where, rbuf->buffer_pos, howmuch);
	rbuf->buffer_left -= howmuch;
	rbuf->buffer_pos += howmuch;
	return howmuch;
    }
}

/* ........................................................... */

static struct urlinfo *newurl (void)
/* Allocate a new urlinfo structure, fill it with default values and
   return a pointer to it.  */
{
    struct urlinfo *u;

    u = mymalloc(sizeof *u);
    memset(u, 0, sizeof *u);
    
    return u;
}

/* ........................................................... */

static void freeurl (struct urlinfo *u, int complete)
/* Perform a "deep" free of the urlinfo structure.  The structure
   should have been created with newurl, but need not have been used.
   If complete is non-0, free the pointer itself.  */
{
    if (u == NULL) return;
    free(u->url);
    free(u->host);
    free(u->path);
    free(u->user);
    free(u->passwd);
    if (complete) {
	free(u);
	u = NULL;
    }
    return;
}

/* ........................................................... */

static int show_progress (long res, long expected, enum spflags flags)
{
    static long offs;
    static ProgressData *pdata;

    if (expected == 0) return 0;
    if (flags == SP_FINISH) {
	if (pdata != NULL)
	    gtk_widget_destroy(GTK_WIDGET(pdata->window)); 
	return 0;
    }
    if (flags == SP_INIT) {
	char bytestr[48];

	offs = 0L;
	if ((pdata = progress_window()) == NULL)
	    return 0;
	gtk_progress_bar_update(GTK_PROGRESS_BAR(pdata->pbar), (gfloat) 0);
	sprintf(bytestr, "Grabbing %ld bytes", expected);
	gtk_label_set_text(GTK_LABEL(pdata->label), bytestr);
	while(gtk_events_pending())
	    gtk_main_iteration();
    }
    offs += res;
    if (pdata != NULL) {
	gtk_progress_bar_update(GTK_PROGRESS_BAR(pdata->pbar), 
				(gfloat) ((double) offs / expected));
	while(gtk_events_pending())
	    gtk_main_iteration();
    } else
	return -1;
	
    return 0;
}

/* ........................................................... */

static int get_contents (int fd, FILE *fp, char **getbuf, long *len, 
			 long expected, struct rbuf *rbuf)
{
    int i, res = 0;
    static char c[8192];

    *len = 0L;
    show_progress(res, expected, SP_INIT);
    if (rbuf && RBUF_FD(rbuf) == fd) {
	while ((res = rbuf_flush(rbuf, c, sizeof c)) != 0) {
	    if (fp == NULL) {
		/*  strncat(*getbuf, c, res); */
		memcpy(*getbuf, c, res);
	    } else 
		if (fwrite(c, 1, res, fp) < (unsigned)res)
		    return -2;
	    *len += res;
	    show_progress(res, expected, SP_NONE);
	}
    }
    /* Read from fd while there is available data. */
    i = (res)? 2 : 1;
    do {
	res = iread(fd, c, sizeof c);
	if (res > 0) {
	    if (fp == NULL) {
		if (i > 1) {
		    *getbuf = realloc(*getbuf, i * 8192);
		    if (*getbuf == NULL)
			return -2;
		}
		/*  strncat(*getbuf, c, res); */  
		memcpy(*getbuf + *len, c, res);
	    } else
		if (fwrite(c, 1, res, fp) < (unsigned)res)
		    return -2;
	    *len += res;
	    if (show_progress(res, expected, SP_NONE) < 0)
		break;
	}
	i++;
    } while (res > 0);
    if (res < -1)
	res = -1;
    show_progress(0, expected, SP_FINISH);
    return res;
}

/* ........................................................... */

static int store_hostaddress (unsigned char *where, const char *hostname)
{
    unsigned long addr;

    addr = (unsigned long)inet_addr(hostname);
    if ((int)addr != -1) {
	memcpy(where, &addr, 4);
	return 1;
    } else
	return 0;
}

#ifdef OS_WIN32
#define ECONNREFUSED WSAECONNREFUSED
#endif

/* ........................................................... */

static uerr_t make_connection (int *sock, char *hostname, unsigned short port)
/* Create an internet connection to HOSTNAME on PORT.  The created
   socket will be stored to *SOCK.  */
{
    struct sockaddr_in sock_name;

    if (!store_hostaddress((unsigned char *)&sock_name.sin_addr, hostname))
	return HOSTERR;

    sock_name.sin_family = AF_INET;
    sock_name.sin_port = g_htons(port);

    if ((*sock = socket(AF_INET, SOCK_STREAM, 0)) == -1)
	return CONSOCKERR;

    if (connect (*sock, (struct sockaddr *) &sock_name, sizeof sock_name)) {
	if (errno == ECONNREFUSED)
	    return CONREFUSED;
	else
	    return CONERROR;
    }
    return NOCONERROR;
}

/* ........................................................... */

static int iread (int fd, char *buf, int len)
/* Read at most LEN bytes from FD, storing them to BUF. */
{
    int res;

    do {
	res = READ(fd, buf, len);
    } while (res == -1 && errno == EINTR);

    return res;
}

/* ........................................................... */

static int iwrite (int fd, char *buf, int len)
/* Write LEN bytes from BUF to FD.  This is similar to iread(), but
   doesn't bother with select().  Unlike iread(), it makes sure that
   all of BUF is actually written to FD, so callers needn't bother
   with checking that the return value equals to LEN.  Instead, you
   should simply check for -1.  */
{
    int res = 0;

    while (len > 0) {
	do {
	    res = WRITE(fd, buf, len);
	} while (res == -1 && errno == EINTR);
	if (res <= 0)
	    break;
	buf += res;
	len -= res;
    }
    return res;
}

/* ........................................................... */

static char *print_option (int opt)
{
    switch (opt) {
    case LIST_DBS:
	return "LIST_DBS";
    case GRAB_IDX:
	return "GRAB_IDX";
    case GRAB_DATA:
	return "GRAB_DATA";
    case GRAB_NBO_DATA:
	return "GRAB_NBO_DATA";
    default:
	break;
    }
    return NULL;
} 

/* ........................................................... */

static int get_update_info (char **saver, char *errbuf, time_t filedate)
{
    uerr_t result;
    struct urlinfo *u;
    int dt, err = 0;
    char datestr[32];
    const char *cgi = "/gretl/cgi-bin/gretl_update.cgi";

    u = newurl();
    u->proto = URLHTTP;
    u->port = DEFAULT_HTTP_PORT;
    u->host = mymalloc(16);
    strcpy(u->host, paths.dbhost_ip);
    u->path = mymalloc(strlen(cgi) + 64);
    sprintf(u->path, "%s?opt=QUERY", cgi);

    strcat(u->path, "&date=");
    sprintf(datestr, "%lu", filedate);
    strcat(u->path, datestr);

    u->filesave = 0;
    u->local = saver;

    result = http_loop(u, &dt);

    if (result == RETROK) {
        errbuf[0] = '\0';
    } else {
	strcpy(errbuf, u->errbuf);
	err = 1;
    }
    freeurl(u, 1);
    return err;
}

/* ........................................................... */

#ifdef G_OS_WIN32
static size_t get_size (char *buf)
{
    size_t i, newsize = 0L;
    int pos;
    char line[60];

    while ((pos = haschar('\n', buf)) > 0) {
	strncpy(line, buf, pos);
	line[pos] = 0;
	sscanf(line, "%*s %u", &i);
	newsize += i;
	buf += pos + 1;
    }

    return newsize;
}
#endif /* G_OS_WIN32 */

/* ........................................................... */

int update_query (void)
{
    int err = 0, admin = 0;
    char *getbuf = NULL;
    char errbuf[80];
    char testfile[MAXLEN];
#ifndef OS_WIN32
    char hometest[MAXLEN];
    FILE *fp;
#endif
    struct stat fbuf;
    long filedate = 0L;

    sprintf(testfile, "%sgretl.stamp", paths.gretldir);

    if (stat(testfile, &fbuf)) {
	fprintf(stderr, "update_query: couldn't stat testfile '%s'\n", 
		testfile);
	return 1;
    } else {
	filedate = fbuf.st_mtime;
#ifndef OS_WIN32
	hometest[0] = '\0';
	if (getuid() != fbuf.st_uid) { 
	    /* user is not owner of gretl.stamp */
	    sprintf(hometest, "%s.gretl.stamp", paths.userdir);
	    if (!stat(hometest, &fbuf)) {
		filedate = fbuf.st_mtime;
	    }
	} else 
	    admin = 1;
#endif
    }

    getbuf = malloc(2048); 
    if (getbuf == NULL)
	return E_ALLOC;
    clear(getbuf, 2048);
    err = get_update_info(&getbuf, errbuf, filedate);
    
    if (err) 
	return 1;

    if (getbuf && strncmp(getbuf, "No new files", 12)) {
	char infotxt[512];

#ifdef G_OS_WIN32 
	sprintf(infotxt, "New files are available from the gretl web site.\n"
		"These files have a combined size of %u bytes.\n\nIf you "
		"would like to update your installation, please quit gretl\n"
		"and run the program titled \"gretl updater\".\n\nOnce the "
		"updater has completed you may restart gretl.",
		get_size(getbuf));

#else
	if (admin) {
	    strcpy(infotxt, "New files are available from the gretl web site\n"
		   "http://ricardo.ecn.wfu.edu/gretl/");
	    fp = fopen(testfile, "w");
	} else {
	    strcpy(infotxt, "You may want to let the system administrator know\n"
		   "that new files are available from the gretl web site\n"
		   "http://ricardo.ecn.wfu.edu/gretl/");
	    fp = fopen(hometest, "w");
	}
	if (fp != NULL) {
	    fprintf(fp, "This file is part of the gretl update notification "
		    "system\n");
	    fclose(fp);
	}
#endif /* G_OS_WIN32 */
	infobox(infotxt);
    }

    free(getbuf);
    return err;
} 

/* ........................................................... */

int retrieve_url (int opt, const char *dbase, const char *series, 
		  int filesave, char **saver, char *errbuf)
/* grab data from URL.  If filesave = 1 then data is stored to
   a local file whose name is given by "saver".  If filesave = 0
   then "save" is presumed to be a char buffer to which the data
   should be printed
*/
{
    uerr_t result;
    struct urlinfo *u;
    int dt;
    const char *cgi = "/gretl/cgi-bin/gretldata.cgi";
    size_t dblen = 0L;

    if (dbase != NULL)
	dblen = strlen(dbase);

    u = newurl();
    u->proto = URLHTTP;
    u->port = DEFAULT_HTTP_PORT;
    u->host = mymalloc(16);
    strcpy(u->host, paths.dbhost_ip);
    u->path = mymalloc(strlen(cgi) + dblen + 64);
    sprintf(u->path, "%s?opt=%s", cgi, print_option(opt));

    if (dblen) {
	strcat(u->path, "&dbase=");
	strcat(u->path, dbase);
    }
    if (series != NULL) {
	strcat(u->path, "&series=");
	strcat(u->path, series);
    }

    if (filesave) {
	u->filesave = 1;
	u->local = mymalloc(sizeof(char *));
	*(u->local) = mymalloc(strlen(*saver) + 1);
	strcpy(*(u->local), *saver);
    } else {
	u->filesave = 0;
	u->local = saver;
    }

    result = http_loop(u, &dt);
    freeurl(u, 1);

    if (result == RETROK) {
	errbuf[0] = 0;
	return 0;
    } else {
	strcpy(errbuf, u->errbuf);
	return 1;
    }
}

/* progress bar stuff for data download */

/* ........................................................... */

static void destroy_progress (GtkWidget *widget, ProgressData *pdata)
{
    pdata->window = NULL;
    g_free(pdata);
    pdata = NULL;
}

/* ........................................................... */

static ProgressData *progress_window (void)
{
    ProgressData *pdata;
    GtkWidget *align;
    GtkWidget *separator;
    GtkWidget *button;
    GtkWidget *vbox;

    pdata = mymalloc(sizeof *pdata);
    if (pdata == NULL) return NULL;

    pdata->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_policy(GTK_WINDOW(pdata->window), FALSE, FALSE, TRUE);

    gtk_signal_connect(GTK_OBJECT(pdata->window), "destroy",
			GTK_SIGNAL_FUNC(destroy_progress),
			pdata);
    gtk_window_set_title(GTK_WINDOW(pdata->window), "gretl download");
    gtk_container_set_border_width(GTK_CONTAINER(pdata->window), 0);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_container_add(GTK_CONTAINER(pdata->window), vbox);
    gtk_widget_show(vbox);

    /* Add a label */
    pdata->label = gtk_label_new("");
    gtk_widget_show(pdata->label);
    gtk_box_pack_start(GTK_BOX(vbox), pdata->label, FALSE, FALSE, 0);
        
    /* Create a centering alignment object */
    align = gtk_alignment_new(0.5, 0.5, 0, 0);
    gtk_box_pack_start(GTK_BOX(vbox), align, FALSE, FALSE, 5);
    gtk_widget_show(align);

     /* Create the GtkProgressBar */
    pdata->pbar = gtk_progress_bar_new();

    gtk_progress_set_format_string(GTK_PROGRESS(pdata->pbar), "%p%%");
    gtk_container_add(GTK_CONTAINER(align), pdata->pbar);
    gtk_progress_set_show_text(GTK_PROGRESS(pdata->pbar), TRUE);
    gtk_widget_show(pdata->pbar);

    separator = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), separator, FALSE, FALSE, 0);
    gtk_widget_show(separator);

    /* Add button to close progress bar window */
    button = gtk_button_new_with_label("Cancel");
    gtk_signal_connect_object(GTK_OBJECT(button), "clicked",
			      (GtkSignalFunc) gtk_widget_destroy,
			       GTK_OBJECT(pdata->window));
    gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);

    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    gtk_widget_show(pdata->window);

    return pdata;
}

#ifdef OS_WIN32

#include <windows.h>
#include <shellapi.h>
#include <string.h>

long GetRegKey (HKEY key, char *subkey, char *retdata)
{
    long err;
    HKEY hkey;

    err = RegOpenKeyEx(key, subkey, 0, KEY_QUERY_VALUE, &hkey);

    if (err == ERROR_SUCCESS) {
	long datasize = MAX_PATH;
	char data[MAX_PATH];

	RegQueryValue(hkey, NULL, (LPSTR)data, &datasize);

	lstrcpy(retdata, data);
	RegCloseKey(hkey);
    }

    return err;
}

int goto_url (const char *url)
{
    char key[MAX_PATH + MAX_PATH];
    int err = 0;

    /* if the ShellExecute() fails */
    if ((long)ShellExecute(NULL, "open", url, NULL, NULL, SW_SHOW) <= 32) {
	/* get the .htm regkey and lookup the program */
	if (GetRegKey(HKEY_CLASSES_ROOT, ".htm", key) == ERROR_SUCCESS) {
	    lstrcat(key,"\\shell\\open\\command");
	    if (GetRegKey(HKEY_CLASSES_ROOT, key, key) == ERROR_SUCCESS) {
		char *pos;
		pos = strstr(key,"\"%1\"");
		if (pos == NULL) {    /* if no quotes */
		    /* now check for %1, without the quotes */
		    pos = strstr(key, "%1");
		    if(pos == NULL) /* if no parameter */
			pos = key + lstrlen(key) - 1;
		    else
			*pos = '\0';    /* remove the parameter */
		}
		else
		    *pos = '\0';        /* remove the parameter */

		lstrcat(pos, " ");
		lstrcat(pos, url);
		if (WinExec(key, SW_SHOW) < 32) err = 1;
	    }
	}
    }
    else
	err = 0;

    return err;
}

#endif

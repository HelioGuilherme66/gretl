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

/* interact.c for gretl */

#include "libgretl.h"
#include "monte_carlo.h"
#include "var.h"
#include "johansen.h"
#include "gretl_func.h"
#include "compat.h"
#include "system.h"
#include "forecast.h"
#include "cmd_private.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_panel.h"
#include "texprint.h"
#include "gretl_xml.h"
#include "gretl_string_table.h"
#include "gretl_typemap.h"
#include "gretl_midas.h"
#include "dbread.h"
#include "gretl_foreign.h"
#include "boxplots.h"
#include "gretl_plot.h"
#include "kalman.h"
#include "flow_control.h"
#include "libglue.h"
#include "csvdata.h"
#include "gretl_zip.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <unistd.h> /* for getcwd() */
#include <errno.h>

/* for the "shell" command */
#ifdef WIN32
# include "gretl_win32.h"
#else
# ifdef HAVE_PATHS_H
#  include <paths.h>
# endif
#endif

#define CMD_DEBUG 0
#define ECHO_DEBUG 0

#include "tokenize.c"

#define bare_quote(p,s)   (*p == '"' && (p-s==0 || *(p-1) != '\\'))
#define starts_comment(p) (*p == '/' && *(p+1) == '*')
#define ends_comment(p)   (*p == '*' && *(p+1) == '/')

static int strip_inline_comments (char *s)
{
    int ret = 0;

    if (*s == '#') {
	/* the entire line is a comment */
	ret = 1;
    } else if (strstr(s, "#")) {
	int quoted = 0;
	int braced = 0;
	char *p = s;

	while (*p) {
	    if (bare_quote(p, s)) {
		quoted = !quoted;
	    } else if (!quoted) {
		if (*p == '{') {
		    braced++;
		} else if (*p == '}') {
		    braced--;
		}
	    }
	    if (!quoted && !braced) {
		if (*p == '#') {
		    *p = '\0';
		    break;
		}
	    }
	    p++;
	}
    }

    return ret;
}

/* filter_comments: strip comments out of line; return non-zero if
   the whole line is a comment */

static int filter_comments (char *s, CMD *cmd)
{
    char tmp[MAXLINE];
    char *p = s;
    int quoted = 0;
    int ignore = (cmd->flags & CMD_IGNORE);
    int j = 0, filt = 0;

    if (strlen(s) >= MAXLINE) {
	cmd->err = E_TOOLONG;
	return 0;
    }

    while (*p) {
	if (!quoted && !ignore && *p == '#') {
	    break;
	}
	if (!ignore && bare_quote(p, s)) {
	    quoted = !quoted;
	}
	if (!quoted) {
	    if (starts_comment(p)) {
		ignore = 1;
		p += 2;
	    } else if (ends_comment(p)) {
		if (!ignore) {
		    cmd->err = E_PARSE;
		    return 0;
		}
		ignore = 0;
		p += 2;
		p += strspn(p, " ");
	    }
	}
	if (!ignore && *p != '\r') {
	    tmp[j++] = *p;
	}
	if (*p) {
	    p++;
	}
    }

    tmp[j] = '\0';
    strcpy(s, tmp);
    tailstrip(s);

    if (*s == '\0') { 
	filt = 1;
    } else if (!ignore) {
	/* '#' comments */
	filt = strip_inline_comments(s);
	tailstrip(s);
    }

    if (filt) {
	/* the whole line is a comment */
	cmd->ci = CMD_COMMENT;
    }

    if (ignore) {
	/* the line ends in multi-line comment mode */
	cmd->flags |= CMD_IGNORE;
    } else {
	cmd->flags &= ~CMD_IGNORE;
    }

    return filt;
}

#define MODIFIES_LIST(c) (c == DIFF ||		\
			  c == DUMMIFY ||	\
			  c == LDIFF ||		\
			  c == SDIFF ||		\
			  c == LAGS ||		\
			  c == LOGS ||		\
			  c == SQUARE ||	\
	                  c == ORTHDEV)

static int has_param (const CMD *cmd)
{
    return cmd->param != NULL && *cmd->param != '\0';
}

/* Look for a line with an "implicit genr", such as
   y = 3*x, x += 10, etc. This is used in nls.c to
   assess auxiliary genrs in nls, mle, gmm.
*/

int plausible_genr_start (const char *s, const DATASET *dset)
{
    int ret = 0;

    if (strchr(s, '=') || strstr(s, "++") || strstr(s, "--")) {
	const char *ok = ".+-*/%^~|=[";
	char word[VNAMELEN] = {0};
	char fmt[20];

	sprintf(fmt, "%%%d[^[ .+*/%%^~|=-]", VNAMELEN - 1);

	if (sscanf(s, fmt, word)) {
	    s += strlen(word);
	    while (*s == ' ') s++;
	    if (strspn(s, ok) > 0 && check_varname(word) == 0) {
		ret = 1;
	    }
	}
    } else if (gretl_type_from_name(s, dset) != 0) {
	ret = 1;
    }

    return ret;
}

static int ends_foreign_block (const char *s)
{
    s += strspn(s, " \t");

    if (!strncmp(s, "end ", 4)) {
	s += 3;
	s += strspn(s, " \t");
	if (!strncmp(s, "foreign", 7)) {
	    return 1;
	} else if (!strncmp(s, "mpi", 3)) {
	    return 1;
	}
    }

    return 0;
}

/**
 * parse_command_line:
 * @line: the command line.
 * @cmd: pointer to command struct.
 * @dset: dataset struct.
 * @ptr: pointer for use with "compilation" of
 * conditionals in loops.
 *
 * Parses @line and fills out @cmd accordingly. 
 *
 * Returns: 0 on success, non-zero code on error.
 */

int parse_command_line (char *line, CMD *cmd, DATASET *dset, void *ptr) 
{
    gretl_cmd_clear(cmd);
    gretl_error_clear();

#if CMD_DEBUG
    fprintf(stderr, "parse_command_line: '%s' (nosub=%d)\n",
	    line, cmd_nosub(cmd) ? 1 : 0);
#endif

    if (cmd_nosub(cmd)) {
	cmd->flags &= ~CMD_SUBST;
    } else {
	int subst = 0;

	cmd->err = substitute_named_strings(line, &subst);
	if (cmd->err) {
	    return cmd->err;
	} else if (subst) {
	    /* record the fact that substitution has been done */
	    cmd->flags |= CMD_SUBST;
	} else {
	    cmd->flags &= ~CMD_SUBST;
	}
    }

#if CMD_DEBUG
    if (cmd->flags & CMD_SUBST) {
	fprintf(stderr, "after substitution: '%s'\n", line);
    }
#endif

    if ((cmd->context == FOREIGN || cmd->context == MPI) &&
	!ends_foreign_block(line)) {
	cmd->opt = OPT_NONE;
	cmd->ci = cmd->context;
	return 0;
    }

    if ((cmd->flags & CMD_SUBST) || !gretl_looping_currently()) {
	/* normalize line spaces */
	compress_spaces(line);
	
	/* trap lines that are nothing but comments */
	if (filter_comments(line, cmd)) {
	    return 0;
	}

	/* catch errors associated with comment syntax */
	if (cmd->err) {
	    return cmd->err;
	}
    }

    cmd->err = real_parse_command(line, cmd, dset, 0, ptr);

    if (cmd->err) {
	gretl_cmd_destroy_context(cmd);
    }

    return cmd->err;
}

#ifndef WIN32

static int gretl_shell_async (const char *arg, PRN *prn)
{
    GError *gerr = NULL;
    int err = 0;

    g_spawn_command_line_async(arg, &gerr);

    if (gerr != NULL) {
	pprintf(prn, "%s\n", gerr->message);
	g_error_free(gerr);
	err = 1;
    }    

    return err;
}

static int gretl_shell_sync (const char *arg, gchar **psout,
			     PRN *prn)
{
    gchar *sout = NULL;
    gchar *serr = NULL;
    GError *gerr = NULL;
    int status;
    gchar *argv[5];
    const char *theshell = getenv("SHELL");
    const char *namep;
    char shellnam[40];
    int err = 0;

    if (theshell == NULL) {
#ifdef HAVE_PATHS_H
	theshell =_PATH_BSHELL;
#else
	theshell = "/bin/sh"; 
#endif
    }

    namep = strrchr(theshell, '/');
    if (namep == NULL) {
	namep = theshell;
    }

    strcpy(shellnam, "-");
    strcat(shellnam, ++namep);
    if (strcmp(namep, "sh") != 0) {
	shellnam[0] = '+';
    }

    argv[0] = g_strdup(theshell);
    argv[1] = shellnam;
    argv[2] = g_strdup("-c");
    argv[3] = g_strdup(arg);
    argv[4] = NULL;

    g_spawn_sync(get_shelldir(), argv, NULL, 0, NULL, NULL,
		 &sout, &serr, &status, &gerr); 

    g_free(argv[0]);
    g_free(argv[2]);
    g_free(argv[3]);

    if (gerr != NULL) {
	if (prn != NULL) {
	    pprintf(prn, "%s\n", gerr->message);
	} else {
	    gretl_errmsg_set(gerr->message);
	}
	g_error_free(gerr);
	err = 1;
    }

    if (psout != NULL) {
	*psout = sout;
    } else if (sout != NULL) {
	pputs(prn, sout);
	g_free(sout);
    }

    if (serr != NULL) {
	pputs(prn, serr);
	g_free(serr);
    }

    return err;
}

/**
 * gretl_shell_grab:
 * @arg: command line to be executed.
 * @sout: location to receive output from command.
 *
 * Calls the shell to execute @arg syncronously and captures the
 * standard output, if any, in @sout.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_shell_grab (const char *arg, char **sout)
{
    return gretl_shell_sync(arg, sout, NULL);
}

static int gretl_shell (const char *arg, gretlopt opt, PRN *prn)
{
    int err = 0;
    
    if (arg == NULL || *arg == '\0') {
	return 0;
    }

    if (!libset_get_bool(SHELL_OK)) {
	gretl_errmsg_set(_("The shell command is not activated."));
	return 1;
    }

    arg += strspn(arg, " \t");

    if (opt & OPT_A) {
	/* "launch" */
	err = gretl_shell_async(arg, prn);
    } else {
	err = gretl_shell_sync(arg, NULL, prn);
    }

    return err;
}

#endif /* ! WIN32 */

#define SAFELEN 78 /* ? */

static void trim_to_length (char *s)
{
    int i, n = strlen(s);

    if (n < SAFELEN - 1) return;

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	    break;
	}
    }
}

void safe_print_line (const char *line, int *plen, PRN *prn)
{
    char tmp[SAFELEN];
    const char *q, *p = line;
    int n, m, rem, out = 0;
    int len0 = *plen;

    rem = n = strlen(line);

    while (out < n) {
	*tmp = 0;
	q = p;
	strncat(tmp, p, SAFELEN - 1);
	len0 = 0;
	trim_to_length(tmp - len0);
	len0 = 0;
	m = strlen(tmp);
	out += m;
	rem = n - out;
	p = q + m;
	if (rem > 0) {
	    pprintf(prn, "%s \\\n ", tmp);
	    *plen = 1;
	} else {
	    pprintf(prn, "%s", tmp);
	    *plen += m;
	}
    }
}

static void new_trim_to_length (char *s, int len)
{
    int n = strlen(s);

    if (n > len) {
	int i, quoted = 0;
	int bp0 = 0, bp1 = 0;

	for (i=1; i<n-1; i++) {
	    if (s[i] == '"' && s[i-1] != '\\') {
		quoted = !quoted;
	    } 
	    if (!quoted && s[i] == ' ') {
		if (i < len) {
		    bp0 = i;
		} else {
		    bp1 = i;
		    break;
		}
	    }
	}
	if (bp0 > 0) {
	    s[bp0] = '\0';
	} else if (bp1 > 0) {
	    s[bp1] = '\0';
	}
    }
}

#define TESTLEN 256
#define LINELEN 70

static void reflow_line (const char *line, const CMD *cmd,
			 const char *leader, PRN *prn)
{
    int maxline = LINELEN;

    if (leader != NULL) {
	maxline -= 2;
	pputs(prn, leader);
    }

    if (cmd != NULL && (cmd->ciflags & CI_EXPR)) {
	/* "genr"-type lines: be more generous? */
	maxline += 10;
    } else if (gretl_in_gui_mode()) {
	/* we can handle a little more width */
	maxline += 4;
    }

    if (strlen(line) < maxline) {
	pputs(prn, line);
    } else {
	const char *p = line;
	char buf[TESTLEN];
	int linenum = 0;

	while (*p) {
	    *buf = '\0';
	    strncat(buf, p, TESTLEN - 1);
	    if (linenum > 0 && leader == NULL) {
		new_trim_to_length(buf, maxline - 2);
	    } else {
		new_trim_to_length(buf, maxline);
	    }
	    p += strlen(buf);
	    if (!string_is_blank(buf)) {
		if (linenum > 0) {
		    pputs(prn, "  ");
		}
		pputs(prn, (*buf == ' ')? buf + 1 : buf);
		if (*p) {
		    pputs(prn, " \\\n");
		}
	    }
	    linenum++;
	}
    }
}

static int command_is_silent (const CMD *cmd, const char *line)
{
    if (cmd == NULL) {
	return 0;
    }
    
    if (cmd->ci == FUNCERR || cmd->ci == PRINTF ||
	(cmd->ci == PRINT && strchr(line, '"'))) {
	return 1;
    }

    if (!strcmp(line, "set echo off") || !strcmp(line, "flush")) {
	return 1;
    }

    if (!strncmp(line, "quit", 4) && string_is_blank(line + 4)) {
	return 1;
    }

    if (cmd->ci == SET && cmd->param != NULL &&
	!strcmp(cmd->param, "echo") &&
	gretl_function_depth() > 0) {
	return 1;
    }

    if (cmd->ci == OUTFILE && cmd->opt == OPT_C) {
	return 1;
    }

    if (*line == '!') {
	return 1;
    }

    return 0;
}

/*
 * real_echo_command:
 * @cmd: pointer to #CMD struct.
 * @line: "raw" command line associated with @cmd.
 * @recording: echo is going to command log (0/1).
 * @prn: pointer to gretl printing struct.
 *
 * Echoes the user command represented by @cmd and @line to
 * @prn.  This is used for two distinct purposes: to give 
 * visual feedback on the command supplied, and (in some 
 * contexts) to record a command that was executed interactively.
 */

static void real_echo_command (CMD *cmd, const char *line, 
			       int recording, PRN *prn)
{
    const char *leader = NULL;
    int compiling = 0;

    if (line == NULL || *line == '\0' || prn == NULL) {
	return;
    }

    if (cmd != NULL && cmd->ci >= NC) {
	return;
    }

    if (gretl_compiling_function() || gretl_compiling_loop()) {
	compiling = 1;
    }

#if ECHO_DEBUG
    if (cmd != NULL) {
	fprintf(stderr, "echo_cmd:\n*** line='%s'\n param='%s' parm2='%s'\n", 
		line, cmd->param, cmd->parm2);
	fprintf(stderr, " cmd->opt=%d, recording=%d, compiling=%d\n",
		cmd->opt, recording, compiling);
	fprintf(stderr, " cmd->ci = %d (%s), context = %d\n", cmd->ci, 
		gretl_command_word(cmd->ci), cmd->context);
	fprintf(stderr, " cmd->savename = '%s'\n", cmd->savename);
	if (cmd->list != NULL) {
	    printlist(cmd->list, "cmd->list");
	}
    }
#endif

    /* certain things don't get echoed at all, if not recording or
       compiling a function or loop */
    if (!recording && !compiling && command_is_silent(cmd, line)) {
	return;
    }

    /* print leading string before echo? */
    if (recording) {
	if (cmd != NULL && cmd->ci == STORE) {
	    leader = "# ";
	}
    } else if (compiling) {
	leader = "> ";
    } else {
	leader = "? ";
    }

    if (cmd != NULL && (cmd->context == FOREIGN || cmd->context == MPI)) {
	if (leader != NULL) {
	    pputs(prn, leader);
	}
	pputs(prn, line);
    } else {
	reflow_line(line, cmd, leader, prn);
    }

    pputc(prn, '\n');

    gretl_print_flush_stream(prn);
}

void gretl_echo_command (CMD *cmd, const char *line, PRN *prn)
{
    real_echo_command(cmd, line, 0, prn);
}

void gretl_record_command (CMD *cmd, const char *line, PRN *prn)
{
    real_echo_command(cmd, line, 1, prn);
}

static int set_var_info (const int *list,
			 const char *parm1,
			 const char *parm2,
			 gretlopt opt,
			 DATASET *dset)
{
    int v = list[1];
    int i, err = 0;

    if (dset == NULL || dset->varinfo == NULL) {
	return E_NODATA;
    } else if (v <= 0 || v >= dset->v) {
	return E_DATA;
    }

    if (opt & OPT_M) {
	err = gretl_list_set_midas(list, dset);
	if (err) {
	    return err;
	}
    }

    for (i=1; i<=list[0]; i++) {
	if (opt & OPT_D) {
	    series_set_discrete(dset, list[i], 1);
	} else if (opt & OPT_C) {
	    series_set_discrete(dset, list[i], 0);
	}
    }

    /* below: we'll accept multi-series lists, but the
       string-setting facility will apply to just the
       first member, as "representative" of the list
    */

    if (opt & OPT_I) {
	const char *s = get_optval_string(SETINFO, OPT_I);

	if (s == NULL) {
	    err = E_ARGS;
	} else {
	    series_record_label(dset, v, s);
	}
    } else if (parm1 != NULL) {
	/* backward compatibility */
	series_record_label(dset, v, parm1);	
    }

    if (opt & OPT_G) {
	const char *s = get_optval_string(SETINFO, OPT_G);

	if (s == NULL) {
	    err = E_ARGS;
	} else {
	    series_record_display_name(dset, v, s);
	}
    } else if (parm2 != NULL) {
	/* backward compatibility */
	series_record_display_name(dset, v, parm2);
    }    

    return err;
}

static void showlabels (const int *list, const DATASET *dset, PRN *prn)
{
    const char *label;
    int i, v, vmax, nl = 0;

    if (dset == NULL || dset->v == 0) {
	pprintf(prn, _("No series are defined\n"));
	return;
    }

    vmax = list == NULL ? dset->v - 1 : list[0];

    for (i=1; i<=vmax; i++) {
	v = list == NULL ? i : list[i];
	if (v >= 0 && v < dset->v) {
	    label = series_get_label(dset, v);
	    if (*label != '\0') {
		nl++;
	    }
	}
    } 

    if (nl == 0) {
	pprintf(prn, "No labels\n");
	return;
    } 

    pprintf(prn, _("Listing labels for variables:\n"));

    for (i=1; i<=vmax; i++) {
	v = list == NULL ? i : list[i];
	if (v >= 0 && v < dset->v) {
	    label = series_get_label(dset, v);
	    if (*label != '\0') {
		pprintf(prn, " %s: %s\n", dset->varname[v], label);
	    }
	}
    }

    pputc(prn, '\n');
}

static int outfile_redirect (PRN *prn, FILE *fp, gretlopt opt,
			     int *parms)
{
    int err;

    err = print_start_redirection(prn, fp);
    if (err) {
	return err;
    }
    
    if (opt & OPT_Q) {
	parms[0] = gretl_echo_on();
	parms[1] = gretl_messages_on();
	set_gretl_echo(0);
	set_gretl_messages(0);
    } else {
	parms[0] = parms[1] = -1;
    }

    return 0;
}

static void maybe_restore_vparms (int *parms)
{
    if (parms[0] == 1) {
	set_gretl_echo(1);
    }
    if (parms[1] == 1) {
	set_gretl_messages(1);
    }    
    parms[0] = parms[1] = -1;
}

static int cwd_is_workdir (void)
{
    char thisdir[MAXLEN];

    if (getcwd(thisdir, MAXLEN - 1) != NULL) {
	int n = strlen(thisdir);

	return strncmp(thisdir, gretl_workdir(), n) == 0;
    }

    return 0;
}

static int redirection_ok (PRN *prn)
{
    int fd = gretl_function_depth();
    
    if (fd == 0) {
	return 0;
    } else if (print_redirected_at_level(prn, fd)) {
	/* we may want to lift this ban in future? */
	return 0;
    } else {
	return 1;
    }
}

static int 
do_outfile_command (gretlopt opt, const char *fname, PRN *prn)
{
    static char outname[MAXLEN];
    static int vparms[2];
    int diverted = 0;
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (!(opt & (OPT_W | OPT_A | OPT_C))) {
	return E_ARGS;
    }

    diverted = print_redirection_level(prn) > 0;

    if (opt & OPT_C) {
	/* command to close outfile */
	if (!diverted) {
	    pputs(prn, _("Output is not currently diverted to file\n"));
	    err = 1;
	} else {
	    print_end_redirection(prn);
	    maybe_restore_vparms(vparms);
	    if (gretl_messages_on() && *outname != '\0') {
		pprintf(prn, _("Closed output file '%s'\n"), outname);
	    }
	}
	return err;
    }

    /* command to divert output to file */
    if (diverted && !redirection_ok(prn)) {
	gretl_errmsg_sprintf(_("Output is already diverted to '%s'"),
			     outname);
	return 1;
    } else if (fname == NULL || *fname == '\0') {
	return E_ARGS;
    } else if (!strcmp(fname, "null")) {
	if (gretl_messages_on()) {
	    pputs(prn, _("Now discarding output\n")); 
	}
	err = outfile_redirect(prn, NULL, opt, vparms);
	*outname = '\0';
    } else if (!strcmp(fname, "stderr")) {
	err = outfile_redirect(prn, stderr, opt, vparms);
	*outname = '\0';
    } else if (!strcmp(fname, "stdout")) {
	err = outfile_redirect(prn, stdout, opt, vparms);
	*outname = '\0';
    } else {
	/* should the stream be opened in binary mode on Windows? */
	char tmp[FILENAME_MAX];
	const char *name = tmp;
	FILE *fp;

	/* switches to workdir if needed */
	strcpy(tmp, fname);
	gretl_maybe_prepend_dir(tmp);

	if (opt & OPT_A) {
	    fp = gretl_fopen(tmp, "a");
	} else {
	    fp = gretl_fopen(tmp, "w");
	}

	if (fp == NULL) {
	    pprintf(prn, _("Couldn't open %s for writing\n"), tmp);
	    return 1;
	}

	if (gretl_messages_on()) {
	    if (cwd_is_workdir()) {
		name = fname;
	    }
	    if (opt & OPT_A) {
		pprintf(prn, _("Now appending output to '%s'\n"), name);
	    } else {
		pprintf(prn, _("Now writing output to '%s'\n"), name);
	    }
	}

	err = outfile_redirect(prn, fp, opt, vparms);
	if (err) {
	    fclose(fp);
	    remove(tmp);
	} else {
	    strcpy(outname, name);
	}
    }

    return err;
}

int call_pca_plugin (VMatrix *cmat, DATASET *dset, 
		     gretlopt opt, PRN *prn)
{
    int (*pca_from_cmatrix) (VMatrix *, DATASET *,
			     gretlopt, PRN *);

    gretl_error_clear();
    
    pca_from_cmatrix = get_plugin_function("pca_from_cmatrix");
    if (pca_from_cmatrix == NULL) {
        return 1;
    }
        
    return (*pca_from_cmatrix) (cmat, dset, opt, prn);
}

static int do_pca (int *list, DATASET *dset,
		   gretlopt opt, PRN *prn)
{
    int freelist = 0;
    int err = 0;

    if (list != NULL && list[0] == 0) {
	return 0;
    }

    if (list == NULL) {
	list = full_var_list(dset, NULL);
	freelist = 1;
    }

    if (list != NULL) {
	VMatrix *cmat = NULL;

	/* adding OPT_U ensures a uniform sample for the
	   correlation or covariance matrix */
	cmat = corrlist(PCA, list, dset, opt, &err);
	if (!err) {
	    err = call_pca_plugin(cmat, dset, opt, prn);
	    if (!err && (opt & (OPT_O | OPT_A))) {
		/* results saved as series */
		if (gretl_messages_on() && !gretl_looping_quietly()) {
		    pputs(prn, "Generated principal component series\n");
		}
	    }
	    free_vmatrix(cmat);
	}
	if (freelist) {
	    free(list);
	}
    }

    return err;
}

static void print_info (gretlopt opt, DATASET *dset, PRN *prn)
{
    if (dset != NULL && dset->descrip != NULL) {
	pprintf(prn, "%s\n", dset->descrip);
    } else {
	pputs(prn, _("No data information is available.\n"));
    }
}

/* After estimating a model, check its errcode member to see
   if anything went wrong, and reset gretl_errno to zero.

   If we're looping (that is, if a loop is in progress at the
   current level of function execution) and @loop_force is 0,
   that's all, but if not then:

   (a) print the model (this may require special handling inside
   loops); 

   (b) if the user has employed the "name <- command" mechanism,
   attach the supplied name to the model; 

   (c) conditionally add the model to the stack in objstack.c,
   and if this is done, signal the fact by setting the 'pmod'
   member of @ExecState.

   (d) if we're called by the GUI program and the model has
   been assigned a name, activate the callback that adds the 
   model to the GUI session.
*/

static int print_save_model (MODEL *pmod, DATASET *dset,
			     gretlopt opt, int loop_force,
			     PRN *prn, ExecState *s)
{
    int err = pmod->errcode;

    if (!err) {
	set_gretl_errno(0);
	if (!gretl_looping_currently() || loop_force) {
	    int havename = *s->cmd->savename != '\0';
	    int window = (opt & OPT_W) != 0;

	    if (havename) {
		gretl_model_set_name(pmod, s->cmd->savename);
	    }
	    printmodel(pmod, dset, opt, prn);
	    attach_subsample_to_model(pmod, dset);
	    s->pmod = maybe_stack_model(pmod, s->cmd, prn, &err);
	    if (!err && gretl_in_gui_mode() && s->callback != NULL && 
		(havename || window)) {
		s->callback(s, s->pmod, GRETL_OBJ_EQN);
	    }
	}
    } 

    return err;
}

static void save_var_vecm (ExecState *s)
{
    maybe_stack_var(s->var, s->cmd);

    if (s->callback != NULL && *s->cmd->savename != '\0' &&
	gretl_in_gui_mode()) {
	s->callback(s, s->var, GRETL_OBJ_VAR);
    }    
}

static void gui_save_system (ExecState *s)
{
    /* note: with GRETL_OBJ_SYS, the business of calling
       "maybe_stack" is handled within system.c, so here
       all we have to do is invoke the GUI callback, if
       appropriate
    */
    if (gretl_in_gui_mode() && s->callback != NULL && 
	*s->cmd->savename != '\0') {
	s->callback(s, s->sys, GRETL_OBJ_SYS);
    }    
}

static int model_test_check (CMD *cmd, DATASET *dset, PRN *prn)
{
    int err = last_model_test_ok(cmd->ci, cmd->opt, dset, prn);

    if (err == E_DATA && cmd->ci == RESTRICT && has_param(cmd)) {
	/* try for a not-yet estimated anonymous system */
	if (get_anonymous_equation_system() != NULL) {
	    gretl_error_clear();
	    err = 0;
	}
    }

    return err;
}

static int get_line_continuation (char *line, FILE *fp, PRN *prn)
{
    char tmp[MAXLINE];
    int err = 0;

    if (!strncmp(line, "quit", 4)) {
	return 0;
    }

    while (top_n_tail(line, MAXLINE, &err)) {
	if (err) {
	    break;
	}
	*tmp = '\0';
	if (fgets(tmp, sizeof tmp, fp) && *tmp != '\0') {
	    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
		pprintf(prn, _("Maximum length of command line "
			       "(%d bytes) exceeded\n"), MAXLINE);
		err = E_TOOLONG;
		break;
	    } else {
		strcat(line, tmp);
		compress_spaces(line);
	    }
	}
    }

    return err;
}

static int run_script (const char *fname, ExecState *s, 
		       DATASET *dset, gretlopt opt, 
		       PRN *prn)
{
    int indent = gretl_if_state_record();
    int echo = gretl_echo_on();
    int messages = gretl_messages_on();
    FILE *fp;
    int iferr, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), fname);
	return E_FOPEN;
    }

    strcpy(s->runfile, fname);

    if (opt & OPT_Q) {
	set_gretl_echo(0);
	set_gretl_messages(0);
    }

    if (gretl_echo_on()) {
	pprintf(prn, "run \"%s\"\n", fname);
    }

    while (fgets(s->line, MAXLINE - 1, fp) && !err) {
	err = get_line_continuation(s->line, fp, prn);
	if (!err) {
	    err = maybe_exec_line(s, dset, NULL);
	}
    }

    fclose(fp);

    if (opt & OPT_Q) {
	set_gretl_echo(echo);
	set_gretl_messages(messages);
    }    

    iferr = gretl_if_state_check(indent);
    if (iferr && !err) {
	err = iferr;
    }

    return err;
}

static int lib_try_http (const char *s, char *fname, int *http)
{
    int err = 0;

    if (strncmp(s, "http://", 7) == 0 ||
	strncmp(s, "https://", 8) == 0) {
#ifdef USE_CURL
	err = retrieve_public_file(s, fname);
	if (!err) {
	    *http = 1;
	}
#else
	gretl_errmsg_set(_("Internet access not supported"));
	err = E_DATA;
#endif
    }

    return err;
}

static int lib_clear_data (ExecState *s, DATASET *dset)
{
    int err = 0;

    if (dset->Z != NULL) {
	err = restore_full_sample(dset, NULL); 
	free_Z(dset); 
    }

    clear_model(s->model);
    clear_datainfo(dset, CLEAR_FULL);
    libgretl_session_cleanup(SESSION_CLEAR_DATASET);
    set_model_count(0);
    gretl_cmd_destroy_context(s->cmd);

    return err;
}

static int join_aggregation_method (const char *s, int *seqval,
				    char **auxname, int *err)
{
    int ret = -1;

    if (!strncmp(s, "seq:", 4)) {
	char *endptr;

	*seqval = (int) strtol(s + 4, &endptr, 10);
	if (*endptr == '\0' && *seqval != 0) {
	    ret = AGGR_SEQ;
	} else {
	    gretl_errmsg_sprintf(_("%s: invalid input '%s'\n"), "--seq", s + 4);
	    *err = E_DATA;
	}
    } else if (!strcmp(s, "count")) {
	ret = AGGR_COUNT;
    } else if (!strcmp(s, "avg")) {
	ret = AGGR_AVG;
    } else if (!strcmp(s, "sum")) {
	ret = AGGR_SUM;
    } else if (!strcmp(s, "min")) {
	ret = AGGR_MIN;
    } else if (!strcmp(s, "max")) {
	ret = AGGR_MAX;
    } else if (!strcmp(s, "none")) {
	ret = AGGR_NONE;
    } else if (!strncmp(s, "min(", 4) ||
	       !strncmp(s, "max(", 4)) {
	const char *p = strchr(s + 4, ')');

	if (p != NULL && strlen(p) == 1) {
	    int len = p - (s + 4);

	    if (len > 0) {
		*auxname = gretl_strndup(s + 4, len);
		if (*auxname == NULL) {
		    *err = E_ALLOC;
		} else {
		    ret = (s[1] == 'a')? AGGR_MAX : AGGR_MIN;
		}
	    }
	} else {
	    *err = E_PARSE;
	}
    } else {
	*err = E_PARSE;
    }

    return ret;
}

static int get_inner_key_id (const char *s, int n,
			     const DATASET *dset,
			     int *err)
{
    char vname[VNAMELEN];
    int id = -1;

    if (n == 0 || n >= VNAMELEN) {
	*err = E_PARSE;
    } else {
	*vname = '\0';
	strncat(vname, s, n);
	if (gretl_namechar_spn(vname) != n) {
	    gretl_errmsg_sprintf(_("field '%s' in command is invalid"), vname);
	    *err = E_PARSE;
	} else {
	    id = current_series_index(dset, vname);
	    if (id < 0) {
		*err = E_UNKVAR;
	    }
	}
    }

    return id;
}

static int *get_inner_keys (const char *s, DATASET *dset,
			    int *err)
{
    int *klist = NULL;
    int ikey1 = -1, ikey2 = -1;
    int nkeys = 0;

    if (strchr(s, ',') == NULL) {
	/* just one key, fine */
	ikey1 = current_series_index(dset, s);
	if (ikey1 < 0) {
	    *err = E_UNKVAR;
	} else {
	    nkeys = 1;
	}
    } else {
	/* we should have a double key */
	int n = strcspn(s, ",");

	ikey1 = get_inner_key_id(s, n, dset, err);

	if (!*err) {
	    s += n + 1;
	    n = strlen(s);
	    ikey2 = get_inner_key_id(s, n, dset, err);
	}

	if (!*err) {
	    nkeys = 2;
	}
    }

    if (!*err) {
	klist = gretl_list_new(nkeys);
	if (klist == NULL) {
	    *err = E_ALLOC;
	} else {
	    klist[1] = ikey1;
	    if (nkeys == 2) {
		klist[2] = ikey2;
	    }
	}
    }

    return klist;
}

static int check_join_import_names (char **S, int ns, 
				    DATASET *dset)
{
    int i, err = 0;
    
    for (i=0; i<ns && !err; i++) {
	if (S[i] == NULL || S[i][0] == '\0') {
	    err = E_DATA;
	} else if (current_series_index(dset, S[i]) < 0) {
	    err = check_varname(S[i]);
	    if (!err && gretl_type_from_name(S[i], NULL)) {
		err = E_TYPES;
	    }
	}
    }

    return err;
}

static char **names_from_array_arg (gretl_array *A, 
				    int *ns,
				    int *err)
{
    char **S = NULL;

    if (gretl_array_get_type(A) != GRETL_TYPE_STRINGS) {
	*err = E_TYPES;
    } else {
	S = gretl_array_get_strings(A, ns);
	if (S == NULL) {
	    *err = E_DATA;
	}
    }

    return S;
}

static char **strings_array_singleton (const char *s,
				       int *err)
{
    int len = strlen(s) + 1;
    char **S = strings_array_new_with_length(1, len);

    if (S == NULL) {
	*err = E_ALLOC;
    } else {
	strcat(S[0], s);
    }

    return S;
}

static int get_join_import_names (const char *s, 
				  DATASET *dset,
				  char ***pvnames,
				  int *pnvars)
{
    char **S = NULL;
    gretl_array *A = NULL;
    int ns = 0;
    int err = 0;

    if (s == NULL) {
	return E_PARSE;
    }

    if (strchr(s, ' ') != NULL) {
	/* @s should hold two or more names */
	S = gretl_string_split(s, &ns, NULL);
	if (S == NULL) {
	    err = E_DATA;
	}
    } else if ((A = get_array_by_name(s)) != NULL) {
	/* @s should be the name of an array of strings */
	S = names_from_array_arg(A, &ns, &err);
    } else {
	/* @s should be a legit series name */
	S = strings_array_singleton(s, &err);
	ns = 1;
    }

    if (S != NULL) {
	err = check_join_import_names(S, ns, dset);
    }

    if (!err) {
	if (A != NULL) {
	    /* copy strings "borrowed" from array */
	    *pvnames = strings_array_dup(S, ns);
	    if (*pvnames == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    /* grab strings allocated here */
	    *pvnames = S;
	}
	*pnvars = ns;
    }

    return err;
}

static int lib_join_data (ExecState *s,
			  char *newfile,
			  DATASET *dset,
			  gretlopt opt,
			  PRN *prn)
{
    gretlopt opts[] = { 
	OPT_I, /* ikey: inner key(s) */
	OPT_O, /* okey: outer key(s) */
	OPT_F, /* filter: filter expression */
	OPT_A, /* aggr: aggregation method */
	OPT_D, /* data: "payload" spec */
	OPT_K, /* tkey: outer time-key name,format */
	OPT_X, /* tconvert: date columns for conversion */
	OPT_T, /* tconv-fmt: format for "tconvert" */
	0 
    };
    char *okey = NULL, *filter = NULL;
    const char *param;
    char **vnames = NULL;
    char *dataname = NULL;
    char *auxname = NULL;
    char *tconvstr = NULL;
    char *tconvfmt = NULL;
    int *ikeyvars = NULL;
    int aggr = 0, seqval = 0;
    int tseries = 0;
    int nvars = 1;
    int i, err = 0;

    if (opt & OPT_K) {
	/* --tkey implies special handling of keys */
	if (opt & (OPT_I | OPT_O)) {
	    return E_BADOPT;
	} else if (!dataset_is_time_series(dset)) {
	    return E_PDWRONG;
	}
    }

    tseries = dataset_is_time_series(dset);
    param = s->cmd->parm2;

    err = get_join_import_names(param, dset, &vnames, &nvars);

    if (!err && nvars > 1) {
	/* multiple series: we can't handle the --data option */
	if (opt & OPT_D) {
	    err = E_BADOPT;
	}
    }

    for (i=0; opts[i] && !err; i++) {
	gretlopt jopt = opts[i];
	const char *param;

	if (opt & jopt) {
	    param = get_optval_string(JOIN, jopt);
	    if (param == NULL) {
		err = E_DATA;
	    } else if (jopt == OPT_I) {		
		/* --ikey: the inner key(s) string */
		ikeyvars = get_inner_keys(param, dset, &err);
	    } else if (jopt == OPT_O) {
		/* --okey: the outer key(s) string */
		okey = gretl_strdup(param);
	    } else if (jopt == OPT_F) {
		/* --filter: string specifying a row filter */
		filter = gretl_strdup(param);
	    } else if (jopt == OPT_A) {
		/* --aggr: aggregation */
		aggr = join_aggregation_method(param, &seqval, 
					       &auxname, &err);
	    } else if (jopt == OPT_D) {
		/* --data: string specifying the outer data series */
		dataname = gretl_strdup(param);
	    } else if (jopt == OPT_K) {
		/* --tkey: string specifying outer time key */
		okey = gretl_strdup(param);
	    } else if (jopt == OPT_X) {
		/* --tconvert: list of time/date cols */
		tconvstr = gretl_strdup(param);
	    } else if (jopt == OPT_T) {
		/* --tconv-fmt: format for tconvert columns */
		tconvfmt = gretl_strdup(param);
	    }
	}
    }

    if (!err && okey != NULL && ikeyvars == NULL && !(opt & OPT_K)) {
	/* We can't have an outer key but no inner one, unless
	   we're matching by the time-series structure of the
	   left-hand dataset (implied by OPT_K)
	*/
	gretl_errmsg_set(_("Inner key is missing"));
	err = E_PARSE;
    }

    if (!err && aggr != 0 && ikeyvars == NULL && !tseries) {
	/* aggregation requires ikeyvars, unless there's
	   an implicit time-series inner key
	*/
	gretl_errmsg_set(_("Inner key is missing"));
	err = E_ARGS;
    }

    if (!err) {
	err = gretl_join_data(newfile, 
			      (const char **) vnames, 
			      nvars, dset, 
			      ikeyvars, okey, filter,
			      dataname, aggr, seqval, 
			      auxname, tconvstr,
			      tconvfmt, opt, prn);
    }

    strings_array_free(vnames, nvars);
    free(ikeyvars);
    free(okey);
    free(filter);
    free(dataname);
    free(auxname);
    free(tconvstr);
    free(tconvfmt);

    return err;
}

#define ALLOW_GUI_OPEN 1

static int lib_open_append (ExecState *s, 
			    DATASET *dset, 
			    char *newfile,
			    PRN *prn)
{
    CMD *cmd = s->cmd;
    gretlopt opt = cmd->opt;
    PRN *vprn = prn;
    int quiet = (opt & OPT_Q);
    int http = 0;
    int dbdata = 0;
    int odbc = 0;
    int ftype;
    int err = 0;

    if (cmd->ci == JOIN && (dset == NULL || dset->v == 0)) {
	return E_NODATA;
    }

#if ALLOW_GUI_OPEN
    if (cmd->ci == OPEN && gretl_function_depth() > 0) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(cmd->ci));
	return E_DATA;
    }
#else
    if (cmd->ci == OPEN && (gretl_in_gui_mode() || gretl_function_depth() > 0)) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(cmd->ci));
	return E_DATA;
    }
#endif

    if (cmd->ci != JOIN && (opt & OPT_O)) {
	odbc = 1;
    }

    err = lib_try_http(cmd->param, newfile, &http);
    if (err) {
	errmsg(err, prn);
	return err;
    }

    if (!http && !odbc) {
	/* not using http or ODBC */
	err = get_full_filename(cmd->param, newfile, (opt & OPT_W)?
				OPT_W : OPT_NONE);
	if (err) {
	    errmsg(err, prn);
	    return err;
	}
    }

    if (opt & OPT_W) {
	ftype = GRETL_NATIVE_DB_WWW;
    } else if (odbc) {
	ftype = GRETL_ODBC;
    } else {
	ftype = detect_filetype(newfile, OPT_P);
    }

    if (cmd->ci == JOIN) {
	if (ftype == GRETL_CSV || ftype == GRETL_XML_DATA ||
	    ftype == GRETL_BINARY_DATA) {
	    dset->modflag = 0;
	    if (ftype != GRETL_CSV) {
		opt |= OPT_G;
	    }
	    err = lib_join_data(s, newfile, dset, opt, prn);
	} else {
	    gretl_errmsg_set("Only CSV and gdt[b] files are supported for now");
	    err = E_DATA;
	}
	if (err) {
	    errmsg(err, prn);
	}
	return err;
    }

    dbdata = (ftype == GRETL_NATIVE_DB || ftype == GRETL_NATIVE_DB_WWW ||
	      ftype == GRETL_RATS_DB || ftype == GRETL_PCGIVE_DB ||
	      ftype == GRETL_ODBC);

    if (cmd->ci == OPEN && !dbdata) {
	lib_clear_data(s, dset);
    } 

    if (quiet) {
	/* in case we hit any problems below... */
	vprn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    } 

    if (ftype == GRETL_CSV) {
	err = import_csv(newfile, dset, opt, vprn);
    } else if (SPREADSHEET_IMPORT(ftype)) {
	err = import_spreadsheet(newfile, ftype, cmd->list, cmd->parm2,
				 dset, opt, vprn);
    } else if (OTHER_IMPORT(ftype)) {
	err = import_other(newfile, ftype, dset, opt, vprn);
    } else if (ftype == GRETL_XML_DATA || ftype == GRETL_BINARY_DATA) {
	err = gretl_read_gdt(newfile, dset, opt, vprn);
    } else if (ftype == GRETL_ODBC) {
	err = set_odbc_dsn(cmd->param, vprn);
    } else if (dbdata) {
	err = set_db_name(newfile, ftype, vprn);
    } else {
	err = gretl_get_data(newfile, dset, opt, vprn);
    }

    if (vprn != prn) {
	if (err) {
	    /* The user asked for --quiet operation, but something
	       went wrong so let's print any info we got on
	       vprn.
	    */
	    const char *buf = gretl_print_get_buffer(vprn);

	    if (buf != NULL && *buf != '\0') {
		pputs(prn, buf);
	    }
	} else {
	    /* print minimal success message */
	    pprintf(prn, _("Read datafile %s\n"), newfile);
	}
	gretl_print_destroy(vprn);
    }

    if (err) {
	errmsg(err, prn);
	return err;
    }

    if (dset->v > 0 && !dbdata && !quiet) {
	list_series(dset, prn);
    }

    if (http) {
	remove(newfile);
    }

    if (dbdata || http || cmd->ci == JOIN) {
	/* signal to the gretlcli callback that we didn't do
	   a regular datafile open 
	*/
	*newfile = '\0';
    }

    return err;
}

static int check_clear_data (void)
{
#if ALLOW_GUI_OPEN
    if (gretl_function_depth() > 0) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(CLEAR));
	return E_DATA;
    }
#else
    if (gretl_in_gui_mode() || gretl_function_depth() > 0) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(CLEAR));
	return E_DATA;
    }
#endif

    return 0;
}

static void schedule_callback (ExecState *s)
{
    if (s->callback != NULL) {
	s->flags |= CALLBACK_EXEC;
    } 
}

static int callback_scheduled (ExecState *s)
{
    return (s->flags & CALLBACK_EXEC) ? 1 : 0;
}

static void callback_exec (ExecState *s, char *fname, int err)
{
    if (!err && s->callback != NULL) {
	if (s->cmd->ci == OPEN) {
	    s->callback(s, fname, 0);
	} else {
	    s->callback(s, NULL, 0);
	}
    }

    s->flags &= ~CALLBACK_EXEC;
    *s->cmd->savename = '\0';
}

static int do_end_restrict (ExecState *s, DATASET *dset)
{
    GretlObjType otype = gretl_restriction_get_type(s->rset);
    gretlopt ropt = gretl_restriction_get_options(s->rset);
    gretlopt opt = s->cmd->opt | ropt;
    int err = 0;

    if (opt & OPT_F) {
	/* restrict --full */
	if (otype == GRETL_OBJ_VAR) {
	    s->var = gretl_restricted_vecm(s->rset, dset, 
					   opt, s->prn, &err);
	    if (!err && s->var != NULL) {
		save_var_vecm(s);
	    }
	} else if (otype == GRETL_OBJ_EQN) {
	    err = gretl_restriction_finalize_full(s, s->rset, dset, 
						  opt, s->prn);
	    if (!err) {
		gretlopt printopt = OPT_NONE;

		if (opt & (OPT_Q | OPT_S)) {
		    printopt = OPT_Q;
		}
		print_save_model(s->pmod, dset, printopt, 1, 
				 s->prn, s);
	    }
	}
    } else {
	err = gretl_restriction_finalize(s->rset, dset, 
					 opt, s->prn);
    }

    s->rset = NULL;

    return err;
}

static int do_debug_command (ExecState *state, const char *param, 
			     gretlopt opt)
{
    int err = incompatible_options(opt, OPT_C | OPT_N | OPT_Q);

    if (err) {
	return err;
    }

    if (opt & (OPT_C | OPT_N)) {
	/* continue, next */
	if (!(state->flags & DEBUG_EXEC)) {
	    gretl_errmsg_set("Debugging is not in progress");
	    return E_DATA;
	} else {
	    /* handled in debug_command_loop */
	    return 0;
	}
    } else {
	/* OPT_Q quits debugging of the given function */
	return user_function_set_debug(param, !(opt & OPT_Q));
    } 
}

/* Given the name of a discrete variable, perform a command for each
   value of the discrete variable. Note that at present the only
   command supported in this way is SUMMARY.  
*/

static int do_command_by (CMD *cmd, DATASET *dset, PRN *prn)
{
    const char *byvar = get_optval_string(cmd->ci, OPT_B);
    series_table *st = NULL;
    gretl_matrix *xvals = NULL;
    const double *x;
    int i, v, nvals = 0;
    int single, err = 0;

    if (dset == NULL || byvar == NULL) {
	return E_DATA;
    }

    /* FIXME accept "unit" and "time"/"period" in place of actual
       variables for panel data? */

    v = current_series_index(dset, byvar);
    if (v < 0) {
	return E_UNKVAR;
    }

    x = (const double *) dset->Z[v];

    if (!series_is_discrete(dset, v) && !gretl_isdiscrete(dset->t1, dset->t2, x)) {
	gretl_errmsg_sprintf(_("The variable '%s' is not discrete"), byvar);
	return E_DATA;
    }

    single = cmd->list[0] == 1;

    xvals = gretl_matrix_values(x + dset->t1, dset->t2 - dset->t1 + 1, 
				OPT_S, &err);

    if (!err) {
	nvals = gretl_vector_get_length(xvals);
	if (nvals == 0) {
	    err = E_DATA;
	} else {
	    st = series_get_string_table(dset, v);
	}
    }

    if (!err && single) {
	pputc(prn, '\n');
	pprintf(prn, _("Summary statistics for %s, by value of %s"),
		dset->varname[cmd->list[1]], dset->varname[v]);
	pputc(prn, '\n');
    }

    for (i=0; i<nvals && !err; i++) {
	Summary *summ = NULL;
	char genline[64];
	double xi = gretl_vector_get(xvals, i);
	double *rv = NULL;

	gretl_push_c_numeric_locale();
	sprintf(genline, "%s == %g", byvar, xi);
	gretl_pop_c_numeric_locale();
	rv = generate_series(genline, dset, prn, &err);

	if (!err) {
	    summ = get_summary_restricted(cmd->list, dset, rv,
					  cmd->opt, prn, &err);
	}

	if (!err) {
	    if (i == 0) {
		pputc(prn, '\n');
	    }
	    if (single) {
		bufspace(2, prn);
	    }
	    if (st != NULL) {
		const char *s = series_table_get_string(st, xi);
		
		pprintf(prn, "%s = %s (n = %d):\n", byvar, s, summ->n);
	    } else {
		pprintf(prn, "%s = %g (n = %d):\n", byvar, xi, summ->n);
	    }
	    print_summary(summ, dset, prn);
	    free_summary(summ);
	}

	free(rv);
    }

    gretl_matrix_free(xvals);

    return err;
}

static void exec_state_prep (ExecState *s)
{
    s->flags &= ~CALLBACK_EXEC;
    s->pmod = NULL;
}

int gretl_delete_variables (int *list,
			    const char *param,
			    gretlopt opt,
			    DATASET *dset,
			    int *renumber,
			    PRN *prn)
{
    int err;

    err = incompatible_options(opt, OPT_T | OPT_D | OPT_F);

    if (!err) {
	if (opt & OPT_T) {
	    /* delete all vars of given type */
	    if (list != NULL || param != NULL) {
		err = E_BADOPT;
	    }
	} else if (opt & OPT_D) {
	    /* delete named vars from database */
	    if (list != NULL || param == NULL) {
		err = E_BADOPT;
	    }
	}
    }

    if (err) {
	return err;
    }

    if (opt & OPT_T) {
	const char *s = get_optval_string(DELEET, OPT_T);

	if (s == NULL) {
	    err = E_ARGS;
	} else {
	    GretlType type = gretl_type_from_string(s);

	    err = delete_user_vars_of_type(type, prn);
	}
    } else if (opt & OPT_D) {
	err = db_delete_series_by_name(param, prn);
    } else if (param != NULL) {
	err = gretl_delete_var_by_name(param, prn);
    } else if (list != NULL) {
	/* delete listed series from dataset */
	if (renumber == NULL && !(opt & OPT_F)) {
	    /* lacking the --force option */
	    pputs(prn, _("You cannot delete series in this context\n"));
	    err = E_DATA;
	} else {
	    err = dataset_drop_listed_variables(list, dset, 
						renumber, prn);
	}	    
    } else {
	err = E_DATA;
    }

    return err;
}

/* OMIT and ADD: if we're estimating a revised model, should
   we be saving it as the "last model", or are we just treating 
   the command as a stand-alone test? 
*/

static int add_omit_save (CMD *cmd)
{
    if (cmd->ci == ADD) {
	/* not saving if given the --lm option */
	return !(cmd->opt & OPT_L);
    } else {
	/* omit: not saving if given the --test-only option */
	return !(cmd->opt & OPT_W);
    }
}

static int VAR_omit_driver (CMD *cmd, DATASET *dset, PRN *prn)
{
    GRETL_VAR *var = get_last_model(NULL);
    int err = 0;

    if (cmd->opt & OPT_W) {
	/* Wald test using VCV */
	err = gretl_VAR_wald_omit_test(var, cmd->list, dset, 
				       cmd->opt, prn);
    } else {
	/* the full deal: estimate reduced system */
	GRETL_VAR *vnew;

	vnew = gretl_VAR_omit_test(var, cmd->list, dset, cmd->opt,
				   prn, &err);
	if (!err) {
	    err = maybe_stack_var(vnew, cmd);
	}
    }

    return err;
}

static int model_print_driver (MODEL *pmod, DATASET *dset,
			       int ci, const char *param,
			       gretlopt opt, PRN *prn)
{
    int err = incompatible_options(opt, OPT_R | OPT_C);

    if (!err) {
	char fname[FILENAME_MAX];

	*fname = '\0';

	if (param != NULL) {
	    /* the legacy mechanism */
	    strcpy(fname, param);
	} else if (opt & OPT_U) {
	    /* try for --output=filename, and if found let
	       the suffix determine the output type
	    */
	    const char *s = get_optval_string(ci, OPT_U);

	    if (s != NULL && *s != '\0') {
		strcpy(fname, s);
		if (has_suffix(fname, ".rtf")) {
		    opt |= OPT_R;
		} else if (has_suffix(fname, ".csv")) {
		    opt |= OPT_C;
		}
	    }
	}

	if (*fname == '\0') {
	    /* fallback */
	    const char *sfx = (opt & OPT_R)? "rtf" :
		(opt & OPT_C)? "csv" : "tex";
	    
	    sprintf(fname, "model_%d.%s", pmod->ID, sfx);
	}

	if (opt & OPT_R) {
	    err = rtfprint(pmod, dset, fname, opt);
	} else if (opt & OPT_C) {
	    err = csvprint(pmod, dset, fname, opt);
	} else {
	    gretlopt texopt = opt;

	    if (ci == EQNPRINT) {
		texopt |= OPT_E;
	    }		
	    err = texprint(pmod, dset, fname, texopt);
	}
	if (!err) {
	    pprintf(prn, _("Model printed to %s\n"), fname);
	}
    }

    return err;
}

#ifdef USE_CURL

static int install_function_package (const char *pkgname,
				     gretlopt opt,
				     PRN *prn)
{
    char *fname = NULL;
    int filetype = 0;
    int local = (opt & OPT_L);
    int http = 0;
    int err = 0;

    if (!strncmp(pkgname, "http://", 7) ||
	!strncmp(pkgname, "https://", 8)) {
	http = 1;
    }

    if (strstr(pkgname, ".gfn")) {
	filetype = 1;
    } else if (strstr(pkgname, ".zip")) {
	filetype = 2;
    } else if (local || http) {
	/* must have suitable suffix */
	err = E_DATA;
    } else {
	/* from gretl server: determine the correct suffix */
	fname = retrieve_remote_pkg_filename(pkgname, &err);
	if (!err) {
	    filetype = strstr(fname, ".zip") ? 2 : 1;
	}
    }

    if (!err) {
	if (http) {
	    /* get @fname as last portion of URL */
	    const char *p = strrchr(pkgname, '/');

	    if (p == NULL) {
		err = E_DATA;
	    } else {
		fname = gretl_strdup(p + 1);
	    }
	} else if (local) {
	    const char *p;
	    
	    gretl_maybe_switch_dir(pkgname);
	    p = strrchr(pkgname, SLASH);
	    if (p != NULL) {
		fname = gretl_strdup(p + 1);
	    }
	    
	}
    }

    if (!err && filetype) {
	const char *basename = fname != NULL ? fname : pkgname;
	const char *instpath = gretl_function_package_path();
	gchar *fullname;

	fullname = g_strdup_printf("%s%s", instpath, basename);

	if (local) {
	    err = gretl_copy_file(pkgname, fullname);
	} else if (http) {
	    /* get file from a specified server */
	    err = retrieve_public_file(pkgname, fullname);
	} else {
	    /* get file from default gretl server */
	    err = retrieve_remote_function_package(basename, fullname);
	}
	
	if (!err && filetype == 2) {
	    err = gretl_unzip_into(fullname, instpath);
	    if (!err) {
		/* delete the zipfile */
		gretl_remove(fullname);
	    }
	}
	
	g_free(fullname);

	if (!err && gretl_messages_on()) {
	    pprintf(prn, "Installed %s\n", basename);
	}
    }

    free(fname);
    
    return err;
}

#else /* !USE_CURL */

/* in this case we can only install from a local file */

static int install_function_package (const char *pkgname,
				     gretlopt opt,
				     PRN *prn)
{
    char *fname = NULL;
    int filetype = 0;
    int err = 0;

    if (!strncmp(pkgname, "http://", 7) ||
	!strncmp(pkgname, "https://", 8)) {
	gretl_errmsg_set(_("Internet access not supported"));
	return E_DATA;
    }

    if (strstr(pkgname, ".gfn")) {
	filetype = 1;
    } else if (strstr(pkgname, ".zip")) {
	filetype = 2;
    } else {
	/* must have suitable suffix */
	err = E_DATA;
    }

    if (!err) {
	/* get last portion of local filename */
	const char *p;

	gretl_maybe_switch_dir(pkgname);
	p = strrchr(pkgname, SLASH);
	if (p != NULL) {
	    fname = gretl_strdup(p + 1);
	}
    }

    if (!err && filetype) {
	const char *basename = fname != NULL ? fname : pkgname;
	const char *instpath = gretl_function_package_path();
	gchar *fullname;

	fullname = g_strdup_printf("%s%s", instpath, basename);

	/* copy file into place */
	err = gretl_copy_file(pkgname, fullname);
	
	if (!err && filetype == 2) {
	    err = gretl_unzip_into(fullname, instpath);
	    if (!err) {
		/* delete the zipfile */
		gretl_remove(fullname);
	    }
	}
	
	g_free(fullname);

	if (!err && gretl_messages_on()) {
	    pprintf(prn, "Installed %s\n", basename);
	}
    }

    free(fname);
    
    return err;
}

#endif

static void abort_execution (ExecState *s)
{
    *s->cmd->savename = '\0';
    gretl_cmd_destroy_context(s->cmd);
    errmsg(E_STOP, s->prn);
}

static int plot_ok;

void set_plot_produced (void)
{
    plot_ok = 1;
}

int is_plotting_command (CMD *cmd)
{
    if (GRAPHING_COMMAND(cmd->ci)) {
	return 1;
    } else if (cmd->ci == END &&
	       cmd->param != NULL &&
	       !strcmp(cmd->param, "plot")) {
	return 1;
    } else {
	return 0;
    }
}

static void maybe_schedule_graph_callback (ExecState *s)
{
    int gui_mode = gretl_in_gui_mode();

    if (graph_written_to_file()) {
	if (gui_mode && *s->cmd->savename != '\0') {
	    /* FIXME? */
	    pprintf(s->prn, "Warning: ignoring \"%s <-\"\n", s->cmd->savename);
	}
	report_plot_written(s->prn);
    } else if (gui_mode) {
	schedule_callback(s);
    }
}

static int execute_plot_call (CMD *cmd, DATASET *dset,
			      char *line, PRN *prn)
{
    gretlopt opt = cmd->opt;
    int err = 0;

    if (gretl_in_gui_mode() && *cmd->savename != '\0') {
	/* saving plot "as icon": add internal option to 
	   override production of a "gpttmp" file 
	*/
	opt |= OPT_G;
    }

#if 1
    if (!gretl_in_gui_mode() && getenv("CLI_NO_PLOTS")
	&& cmd->ci != END) {
	return 0;
    }
#endif
    
    if (cmd->ci == END) {
	/* end of a "plot" block */
	err = gretl_plot_finalize(line, dset, opt);
    } else if (opt & OPT_X) {
	err = matrix_command_driver(cmd->ci, cmd->list, cmd->param, 
				    dset, opt, prn);
    } else if (cmd->ci == GNUPLOT) {
	if (opt & OPT_I) {
	    err = gnuplot_process_file(opt, prn);
	} else if (opt & OPT_C) {
	    err = xy_plot_with_control(cmd->list, cmd->param, 
				       dset, opt);
	} else {
	    err = gnuplot(cmd->list, cmd->param, dset, opt);
	}
    } else if (cmd->ci == SCATTERS) {
	err = multi_scatters(cmd->list, dset, opt);
    } else if (cmd->ci == BXPLOT) {
	err = boxplots(cmd->list, cmd->param, dset, opt);
    } else if (cmd->ci == HFPLOT) {
	err = hf_plot(cmd->list, cmd->param, dset, opt);
    }

    return err;
}

static void maybe_print_error_message (CMD *cmd, int err, PRN *prn)
{
    if (gretl_function_depth() > 0) {
	; /* defer printing */
    } else if (cmd->flags & CMD_CATCH) {
	/* print only if messages on */
	if (gretl_messages_on()) {
	    errmsg(err, prn);
	}
    } else {
	/* otherwise go ahead and print */
	errmsg(err, prn);
    }
}

int gretl_cmd_exec (ExecState *s, DATASET *dset)
{
    CMD *cmd = s->cmd;
    char *line = s->line;
    MODEL *model = s->model;
    PRN *prn = s->prn;
    char readfile[MAXLEN];
    int *listcpy = NULL;
    int err = 0;

    exec_state_prep(s);
    plot_ok = 0;

    if (gretl_in_gui_mode() && check_for_stop()) {
	/* the GUI user clicked the "Stop" button */
	abort_execution(s);
	return E_STOP;
    }

    if (NEEDS_MODEL_CHECK(cmd->ci)) {
	err = model_test_check(cmd, dset, prn);
    } else if (MODIFIES_LIST(cmd->ci)) {
	if (cmd->list[0] == 0) {
	    /* no-op */
	    return 0;
	} else {
	    /* list is potentially modified -> make a copy */
	    listcpy = gretl_list_copy(cmd->list);
	    if (listcpy == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    *readfile = '\0';

    if (cmd->ci == OLS && dataset_is_panel(dset)) {
	cmd->ci = PANEL;
	cmd->opt |= OPT_P; /* panel pooled OLS flag */
    }

#if 0
    fprintf(stderr, "gretl_cmd_exec: '%s' (ci %d) \n", line, cmd->ci);
#endif

    switch (cmd->ci) {

    case APPEND:
    case JOIN:
    case OPEN:
	err = lib_open_append(s, dset, readfile, prn);
	if (!err) {
	    schedule_callback(s);
	}
	break;

    case CLEAR:
	err = check_clear_data();
	if (!err) {
	    if (gretl_in_gui_mode()) {
		schedule_callback(s);
	    } else {
		lib_clear_data(s, dset);
	    }
	}
	break;

    case FLUSH:
	if (gretl_in_gui_mode()) {
	    schedule_callback(s);
	}
	break;

    case ANOVA:
	err = anova(cmd->list, dset, cmd->opt, prn);
	break;

    case ADF:
	err = adf_test(cmd->order, cmd->list, dset, cmd->opt, prn);
	break;

    case KPSS:
	err = kpss_test(cmd->order, cmd->list, dset, cmd->opt, prn);
	break;

    case LEVINLIN:
	err = llc_test_driver(cmd->param, cmd->list, dset, 
			      cmd->opt, prn);
	break;

    case COINT:
	err = engle_granger_test(cmd->order, cmd->list, dset, 
				 cmd->opt, prn);
	break;

    case COINT2:
	err = johansen_test_simple(cmd->order, cmd->list, 
				   dset, cmd->opt, prn);
	break;

    case CORR:
	err = incompatible_options(cmd->opt, OPT_U | OPT_S | OPT_K);
	if (err) {
	    break;
	}
	if (cmd->opt & OPT_K) {
	    err = kendall_tau(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->opt & OPT_S) {
	    err = spearman_rho(cmd->list, dset, cmd->opt, prn);
	} else {
	    err = gretl_corrmx(cmd->list, dset, cmd->opt, prn);
	}
	break;

    case CORRGM:
	err = corrgram(cmd->list[1], cmd->order, 0, dset, 
		       cmd->opt, prn);
	break;

    case XCORRGM:
	err = xcorrgram(cmd->list, cmd->order, dset, 
			cmd->opt, prn);
	break;

    case PERGM:
	err = periodogram(cmd->list[1], cmd->order, dset, 
			  cmd->opt, prn);
	break;

    case FRACTINT:
	err = fractint(cmd->list[1], cmd->order, dset, 
		       cmd->opt, prn);
	break;	

    case FUNDEBUG:
	err = do_debug_command(s, cmd->param, cmd->opt);
	break;

    case BREAK:
    case ENDLOOP:
	pprintf(prn, _("You can't end a loop here, "
		       "you haven't started one\n"));
	err = 1;
	break;

    case FCAST:
	err = do_forecast(cmd->param, dset, cmd->opt, prn);
	break;

    case FREQ:
	if (cmd->opt & OPT_X) {
	    err = matrix_freq_driver(cmd->list, cmd->opt, prn);
	} else {
	    err = freqdist(cmd->list[1], dset, cmd->opt, prn);
	}
	break;

    case DISCRETE:
	err = list_makediscrete(cmd->list, dset, cmd->opt);
	break;

    case ESTIMATE:
	err = estimate_named_system(cmd->param, cmd->parm2, dset, 
				    cmd->opt, prn);
	break;

    case FUNC:
	err = gretl_start_compiling_function(cmd->param, prn);
	break;

    case GENR:
    case EVAL:
	if (cmd->flags & CMD_CATCH) {
	    err = generate(cmd->vstart, dset, cmd->gtype,
			   cmd->opt | OPT_C, prn);
	    if (err == E_BADCATCH) {
		cmd->flags ^= CMD_CATCH;
	    }
	} else {
	    err = generate(cmd->vstart, dset, cmd->gtype,
			   cmd->opt, prn);
	}
	break;

    case PCA:
	err = do_pca(cmd->list, dset, cmd->opt, prn);
	break;

    case DATA:
	err = db_get_series(cmd->param, dset, cmd->opt, prn);
	break;

    case DATAMOD:
	err = modify_dataset(dset, cmd->auxint, cmd->list, 
			     cmd->parm2, prn);
	if (!err) { 
	    schedule_callback(s);
	} 
	break;

    case DIFF:
    case LDIFF:
    case SDIFF:
	err = list_diffgenr(listcpy, cmd->ci, dset);
	if (!err) {
	    maybe_list_series(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case ORTHDEV:
	err = list_orthdev(listcpy, dset);
	if (!err) {
	    maybe_list_series(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case DUMMIFY:
	err = list_dumgenr(&listcpy, dset, cmd->opt);
	if (!err) {
	    maybe_list_series(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case LAGS:
	err = list_laggenr(&listcpy, 1, cmd->order, NULL,
			   dset, 0, cmd->opt); 
	if (!err) {
	    maybe_list_series(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case LOGS:
	err = list_loggenr(listcpy, dset);
	if (!err) {
	    maybe_list_series(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case SQUARE:
	err = list_xpxgenr(&listcpy, dset, cmd->opt);
	if (!err) {
	    maybe_list_series(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case TEXTPLOT:
	err = textplot(cmd->list, dset, cmd->opt, prn);
	break;

    case RMPLOT:
    case HURST:
	if (cmd->ci == RMPLOT) {
	    err = rmplot(cmd->list, dset, cmd->opt, prn);
	} else {
	    err = hurstplot(cmd->list, dset, cmd->opt, prn);
	}
	break;

    case QQPLOT:
	err = qq_plot(cmd->list, dset, cmd->opt);
	break;	

    case INFO:
	print_info(cmd->opt, dset, prn);
	break;

    case RENAME:
	err = dataset_rename_series(dset, cmd->auxint, cmd->parm2);
	if (!err) {
	    maybe_list_series(dset, prn);
	}
	break;

    case SET:
	err = execute_set(cmd->param, cmd->parm2, dset, cmd->opt, prn);
	break;

    case SETINFO:
	err = set_var_info(cmd->list, cmd->param, cmd->parm2, 
			   cmd->opt, dset);
	break;

    case SETMISS:
	err = set_miss(cmd->list, cmd->param, dset, prn);
        break;

    case LABELS:
	if (cmd->opt) {
	    err = read_or_write_var_labels(cmd->opt, dset, prn);
	    if (!err && (cmd->opt & (OPT_D | OPT_F))) {
		schedule_callback(s);
	    }
	} else {
	    showlabels(cmd->list, dset, prn);
	}
	break;

    case MARKERS:
	err = read_or_write_obs_markers(cmd->opt, dset, prn);
	if (!err && (cmd->opt & (OPT_D | OPT_F))) {
	    schedule_callback(s);
	}
	break;

    case VARLIST:
	if (cmd->opt & OPT_T) {
	    list_user_vars_of_type(dset, prn);
	} else if (cmd->opt & OPT_S) {
	    print_scalars(prn);
	} else if (cmd->opt & OPT_A) {
	    list_ok_dollar_vars(dset, prn);
	} else {
	    list_series(dset, prn);
	}
	break;

    case PRINT:
	if (cmd->opt & OPT_L) {
	    err = printdata(NULL, cmd->param, dset, OPT_NONE, prn);
	} else if (cmd->param != NULL) {
	    /* directly printing a string literal */
	    pputs(prn, cmd->param);
	    pputc(prn, '\n');
	} else {
	    err = printdata(cmd->list, cmd->parm2, dset, cmd->opt, prn);
	}
	break;

    case PRINTF:
    case SPRINTF:
    case SSCANF:
	err = do_printscan_command(cmd->ci, cmd->param, cmd->parm2,
				   cmd->vstart, dset, prn); 	 
	break;

    case PVAL:
	err = batch_pvalue(cmd->param, dset, prn);
	break;

    case SUMMARY:
	err = incompatible_options(cmd->opt, OPT_B | OPT_W);
	if (err) {
	    break;
	}
	if (cmd->opt & OPT_B) {
	    err = do_command_by(cmd, dset, prn);
	} else if (cmd->opt & OPT_X) {
	    err = matrix_command_driver(cmd->ci, cmd->list, cmd->param,
					dset, cmd->opt, prn);
	} else {
	    err = list_summary_driver(cmd->list, dset, cmd->opt, prn);
	}
	break; 

    case XTAB:
	if (cmd->opt & OPT_X) {
	    err = crosstab_from_matrix(cmd->opt, prn);
	} else {
	    err = crosstab(cmd->list, dset, cmd->opt, prn);
	}
	break;

    case MAHAL:
	err = mahalanobis_distance(cmd->list, dset, cmd->opt, prn);
	break;

    case MEANTEST:
	err = means_test(cmd->list, dset, cmd->opt, prn);
	break;	

    case VARTEST:
	err = vars_test(cmd->list, dset, prn);
	break;

    case RUNS:
	err = runs_test(cmd->list[1], dset, cmd->opt, prn);
	break;

    case SPEARMAN:
	err = spearman_rho(cmd->list, dset, cmd->opt, prn);
	break;

    case DIFFTEST:
	err = diff_test(cmd->list, dset, cmd->opt, prn);
	break;

    case OUTFILE:
	err = do_outfile_command(cmd->opt, cmd->param, prn);
	break;

    case SETOBS:
	err = set_obs(cmd->param, cmd->parm2, dset, cmd->opt);
	if (!err) {
	    if (dset->n > 0) {
		if (!(cmd->opt & (OPT_I | OPT_G))) {
		    print_smpl(dset, 0, prn);
		}
		schedule_callback(s);
	    } else {
		pprintf(prn, _("data frequency = %d\n"), dset->pd);
	    }
	}
	break;

    case SETOPT:
	err = set_options_for_command(cmd->param, cmd->parm2, cmd->opt);
	if (!err && gretl_messages_on()) {
	    pprintf(prn, "Set option(s) for command \"%s\"\n", cmd->param);
	}
	break;

    case SMPL:
	if (cmd->opt == OPT_F) {
	    err = restore_full_sample(dset, s);
	} else if ((cmd->opt & OPT_T) && (cmd->opt & OPT_U)) {
	    err = perma_sample(dset, cmd->opt, prn, NULL);
	} else if (cmd->opt) {
	    err = restrict_sample(cmd->param, cmd->list, dset, 
				  s, cmd->opt, prn, NULL);
	} else { 
	    err = set_sample(cmd->param, cmd->parm2, dset);
	}
	if (!err) {
	    print_smpl(dset, get_full_length_n(), prn);
	}	
	break;

    case INSTALL:
	if (cmd->opt & (OPT_R | OPT_P)) {
	    err = uninstall_function_package(cmd->param, cmd->opt, prn);
	} else {
	    err = install_function_package(cmd->param, cmd->opt, prn);
	}
	break;

    case MAKEPKG:
	err = create_and_write_function_package(cmd->param, cmd->opt, prn);
	break;

    case STORE:
	if (dset == NULL || dset->Z == NULL) {
	    err = E_NODATA;
	} else if (!has_param(cmd)) {
	    pputs(prn, _("store: no filename given\n"));
	    err = E_PARSE;
	}
	if (!err) {
	    err = write_data(cmd->param, cmd->list, dset, cmd->opt, prn);
	}
	break;

    case SHELL:
	err = gretl_shell(cmd->vstart, cmd->opt, prn);
	break;

    case OLS:
    case WLS:
	clear_model(model);
	*model = lsq(cmd->list, dset, cmd->ci, cmd->opt);
	err = print_save_model(model, dset, cmd->opt, 0, prn, s);
	break;
	
    case MPOLS:
	clear_model(model);
	*model = mp_ols(cmd->list, dset);
	err = print_save_model(model, dset, cmd->opt, 0, prn, s);
	break;

    case AR:
    case AR1:
    case ARMA:
    case ARCH:
	clear_model(model);
	if (cmd->ci == AR) {
	    *model = ar_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == AR1) {
	    *model = ar1_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == ARMA) {
	    *model = arma(cmd->list, cmd->auxlist, dset, 
			  cmd->opt, prn);
	} else {
	    *model = arch_model(cmd->list, cmd->order, dset,
				cmd->opt);
	}
	err = print_save_model(model, dset, cmd->opt, 0, prn, s);
	break;

    case ARBOND:
    case PANEL:	
    case DPANEL:
	if (!dataset_is_panel(dset)) {
	    gretl_errmsg_set(_("This estimator requires panel data"));
	    err = E_DATA;
	    break;
	}
    case GARCH:
    case HECKIT:
    case HSK:
    case INTREG:
    case IVREG:
    case LAD:
    case LOGISTIC:
    case LOGIT:
    case POISSON:
    case NEGBIN:
    case PROBIT:
    case QUANTREG:
    case TOBIT:
    case DURATION:
    case BIPROBIT:
    case MIDASREG:
	clear_model(model);
	if (cmd->ci == LOGIT || cmd->ci == PROBIT) {
	    *model = logit_probit(cmd->list, dset, cmd->ci, cmd->opt, prn);
	} else if (cmd->ci == HSK) {
	    *model = hsk_model(cmd->list, dset, cmd->opt);
	} else if (cmd->ci == LOGISTIC) {
	    *model = logistic_driver(cmd->list, dset, cmd->opt);
	} else if (cmd->ci == TOBIT) {
	    *model = tobit_driver(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == POISSON || cmd->ci == NEGBIN) {
	    *model = count_model(cmd->list, cmd->ci, dset, cmd->opt, prn);
	} else if (cmd->ci == HECKIT) {
	    *model = heckit_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == IVREG) {
	    *model = ivreg(cmd->list, dset, cmd->opt);
	} else if (cmd->ci == LAD) {
	    *model = lad(cmd->list, dset);
	} else if (cmd->ci == QUANTREG) {
	    *model = quantreg_driver(cmd->param, cmd->list, dset,
				     cmd->opt, prn);
	} else if (cmd->ci == DURATION) {
	    *model = duration_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == GARCH) {
	    *model = garch(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == PANEL) {
	    *model = panel_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == ARBOND) {
	    *model = arbond_model(cmd->list, cmd->param, dset, 
				  cmd->opt, prn);
	} else if (cmd->ci == DPANEL) {
	    *model = dpd_model(cmd->list, cmd->auxlist, cmd->param, 
			       dset, cmd->opt, prn);
	} else if (cmd->ci == INTREG) {
	    *model = interval_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == BIPROBIT) {
	    *model = biprobit_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == MIDASREG) {
	    *model = midas_model(cmd->list, cmd->param, dset,
				 cmd->opt, prn);
	} else {
	    /* can't happen */
	    err = 1;
	    break;
	}
	err = print_save_model(model, dset, cmd->opt, 0, prn, s);
	break;

    case GMM:
    case MLE:
    case NLS:
	err = nl_parse_line(cmd->ci, cmd->vstart, dset, prn);
	if (!err) {
	    gretl_cmd_set_context(cmd, cmd->ci);
	} 
	break;

    case FOREIGN:
    case MPI:
	if (cmd->context == FOREIGN || cmd->context == MPI) {
	    err = foreign_append(line, cmd->context);
	} else {
	    err = foreign_start(cmd->ci, cmd->param, cmd->opt, prn);
	    if (!err) {
		gretl_cmd_set_context(cmd, cmd->ci);
	    }
	}
	break;

    case KALMAN:
	/* tokenizer: @line arg OK for now */
	err = kalman_parse_line(line, dset, cmd->opt, prn);
	if (!err && (cmd->opt == OPT_NONE)) {
	    gretl_cmd_set_context(cmd, cmd->ci);
	}
	break;

    case PLOT:
	if (!cmd->context) {
	    err = gretl_plot_start(cmd->param, dset);
	} else {
	    err = gretl_plot_append_line(line, dset);
	}
	if (!err && !cmd->context) {
	    gretl_cmd_set_context(cmd, cmd->ci);
	}
	break;	

    case ADD:
    case OMIT:
	if (get_last_model_type() == GRETL_OBJ_VAR) {
	    err = VAR_omit_driver(cmd, dset, prn);
	} else if (add_omit_save(cmd)) {
	    MODEL mymod;
	    
	    gretl_model_init(&mymod, dset);
	    if (cmd->ci == ADD) {
		err = add_test_full(model, &mymod, cmd->list, 
				    dset, cmd->opt, prn);
	    } else {
		err = omit_test_full(model, &mymod, cmd->list, 
				     dset, cmd->opt, prn);
	    }
	    if (!err) {
		gretlopt popt = OPT_NONE;

		if (cmd->opt & (OPT_I | OPT_Q)) {
		    popt = OPT_Q;
		} else if (cmd->opt & OPT_O) {
		    popt = OPT_O; /* --vcv printing option */
		}
		clear_model(model);
		*model = mymod;
		print_save_model(model, dset, popt, 1, prn, s);
	    } 
	} else if (cmd->ci == ADD) {
	    err = add_test(model, cmd->list, dset, cmd->opt, prn);
	} else {
	    err = omit_test(model, cmd->list, dset, cmd->opt, prn);
	}
	if (err == E_NOOMIT) {
	    /* auto-omit was a no-op */
	    err = 0;
	}	
	break;	

    case COEFFSUM:
    case CUSUM:
    case RESET:
    case CHOW:
    case QLRTEST:
    case VIF:
	if (cmd->ci == COEFFSUM) {
	    err = gretl_sum_test(cmd->list, model, dset, prn);
	} else if (cmd->ci == CUSUM) {
	    err = cusum_test(model, dset, cmd->opt, prn);
	} else if (cmd->ci == RESET) {
	    err = reset_test(model, dset, cmd->opt, prn);
	} else if (cmd->ci == CHOW) {
	    err = chow_test_driver(cmd->param, model, dset, cmd->opt, prn);
	} else if (cmd->ci == QLRTEST) {
	    err = QLR_test(model, dset, cmd->opt, prn);
	} else if (cmd->ci == VIF) {
	    err = vif_test(model, dset, prn);
	} 
	break;

    case NORMTEST:
	err = gretl_normality_test(cmd->list[1], dset, cmd->opt, prn);
	break;

    case HAUSMAN:
	err = panel_hausman_test(model, dset, cmd->opt, prn);
	break;

    case MODTEST:
	err = model_test_driver(cmd->order, dset, cmd->opt, prn);
	break;

    case LEVERAGE:
	err = leverage_test(model, dset, cmd->opt, prn);
	if (!err && (cmd->opt & OPT_S) && !(cmd->opt & OPT_Q)) {
	    maybe_list_series(dset, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	if (model->errcode == E_NAN) {
	    pprintf(prn, _("Couldn't format model\n"));
	} else {
	    err = model_print_driver(model, dset, cmd->ci,
				     cmd->param, cmd->opt, 
				     prn);
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model */
	if (s->rset == NULL) {
	    if (!has_param(cmd)) {
		/* if param is non-blank, we're restricting a named system */
		err = model_test_check(cmd, dset, prn);
		if (err) break;
	    }
	    s->rset = restriction_set_start(cmd->param, cmd->opt, &err);
	    if (!err) {
		/* FIXME redundant? */
		gretl_cmd_set_context(cmd, RESTRICT);
	    }
	} else {
	    err = restriction_set_parse_line(s->rset, line, dset);
	    if (err) {
		s->rset = NULL;
	    }	
	}
	break;

    case SYSTEM:
	if (s->sys == NULL) {
	    /* no equation system is defined currently */
	    s->sys = equation_system_start(cmd->param, cmd->savename, 
					   cmd->opt, &err);
	    if (!err) {
		gretl_cmd_set_context(cmd, SYSTEM);
	    }
	} else {
	    /* tokenize: use of @line OK here? */
	    err = system_parse_line(s->sys, line, dset);
	    if (err) {
		s->sys = NULL;
	    } 
	}
	break;

    case EQUATION:
	if (cmd->opt & OPT_M) {
	    err = equation_system_append_multi(s->sys, cmd->param, dset);
	} else {
	    err = equation_system_append(s->sys, cmd->list);
	}
	if (err) {
	    s->sys = NULL;
	}
	break;

    case END:
	if (!strcmp(cmd->param, "system")) {
	    err = equation_system_finalize(s->sys, dset, cmd->opt, prn);
	    if (!err) {
		gui_save_system(s);
	    }
	    /* clear for next use */
	    s->sys = NULL;
	} else if (!strcmp(cmd->param, "mle") || 
		   !strcmp(cmd->param, "nls") ||
		   !strcmp(cmd->param, "gmm")) {
	    clear_model(model);
	    *model = nl_model(dset, cmd->opt, prn);
	    err = print_save_model(model, dset, cmd->opt, 0, prn, s);
	} else if (!strcmp(cmd->param, "restrict")) {
	    err = do_end_restrict(s, dset);
	} else if (!strcmp(cmd->param, "foreign")) {
	    err = foreign_execute(dset, cmd->opt, prn);
	} else if (!strcmp(cmd->param, "kalman")) {
	    err = kalman_parse_line(line, dset, cmd->opt, prn);
	} else if (!strcmp(cmd->param, "mpi")) {
	    err = foreign_execute(dset, cmd->opt, prn);
	} else if (!strcmp(cmd->param, "plot")) {
	    err = execute_plot_call(cmd, dset, line, prn);
	} else {
	    err = 1;
	}
	break;

    case VAR:
    case VECM:
	if (cmd->ci == VAR) {
	    s->var = gretl_VAR(cmd->order, cmd->auxlist, cmd->list, 
			       dset, cmd->opt, prn, &err);
	} else {
	    s->var = gretl_VECM(cmd->order, cmd->auxint, cmd->list, 
				dset, cmd->opt, prn, &err);
	}
	if (!err && s->var != NULL) {
	    save_var_vecm(s);
	}
	break;

    case RUN:
    case INCLUDE:
	if (cmd->ci == RUN) {
	    err = get_full_filename(cmd->param, readfile, OPT_S);
	} else {
	    err = get_full_filename(cmd->param, readfile, OPT_I);
	    cmd->opt |= OPT_Q;
	}
	if (err) { 
	    break;
	} 
	if (gretl_messages_on()) {
	    pprintf(prn, " %s\n", readfile);
	}
	if (cmd->ci == INCLUDE && gretl_is_xml_file(readfile)) {
	    err = load_user_XML_file(readfile, prn);
	    break;
	} else if (cmd->ci == INCLUDE && gfn_is_loaded(readfile)) {
	    break;
	}
	if (!strcmp(readfile, s->runfile)) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    err = 1;
	    break;
	}
	err = run_script(readfile, s, dset, cmd->opt, prn);
	break;

    case FUNCERR:
    case FUNCRET:
	if (gretl_function_depth() == 0) {
	    gretl_errmsg_sprintf("'%s': can only be used within a function",
				 gretl_command_word(cmd->ci));
	    err = 1;
	} else if (cmd->ci == FUNCERR) {
	    err = E_FUNCERR;
	} 
	break;

    case DELEET:
	err = gretl_delete_variables(cmd->list, cmd->param, cmd->opt,
				     dset, NULL, prn);
	break;

    case MODPRINT:
	err = do_modprint(cmd->param, cmd->parm2, cmd->opt, prn);
	break;

    case GNUPLOT:
    case BXPLOT:
    case SCATTERS:
    case HFPLOT:
	err = execute_plot_call(cmd, dset, NULL, prn);
	break;

    case MODELTAB:
    case GRAPHPG:
	if (gretl_in_gui_mode()) {
	    schedule_callback(s);
	} else {
	    pprintf(prn, _("%s: command not available\n"), 
		    gretl_command_word(cmd->ci));
	}
	break;

    default:
	{
	    const char *word = gretl_command_word(cmd->ci);

	    if (*word != '\0') {
		pprintf(prn, _("Sorry, the %s command is not yet implemented "
			       "in libgretl\n"), word);
	    } else {
		pprintf(prn, "What?\n");
	    }
	}
	err = 1;
	break;
    }

    if (listcpy != NULL) {
	free(listcpy);
    }

    if (err == E_OK) {
	err = 0;
    }

    if (!err && plot_ok) {
	maybe_schedule_graph_callback(s);
    }

    if (callback_scheduled(s)) {
	callback_exec(s, readfile, err);
    } 

 bailout:

    if (err) {
	maybe_print_error_message(cmd, err, prn);
	err = process_command_error(s, err);
    }

    if (err) {
	gretl_cmd_destroy_context(cmd);
    } else {
	/* this is a no-op if there's no warning */
	warnmsg(prn);
    }

    return err;
}

/* called by functions, and by scripts executed from within
   functions */

int maybe_exec_line (ExecState *s, DATASET *dset, int *loopstart)
{
    int err = 0;

    if (string_is_blank(s->line)) {
	return 0;
    }

    if (gretl_compiling_loop()) {
	err = get_command_index(s->line, LOOP, s->cmd);
    } else {
	/* FIXME last arg to parse_command_line() ? */
	err = parse_command_line(s->line, s->cmd, dset, NULL);
	if (loopstart != NULL && s->cmd->ci == LOOP) {
	    *loopstart = 1;
	}
    }

    if (err) {
        errmsg(err, s->prn);
        return err;
    }

    gretl_exec_state_transcribe_flags(s, s->cmd);

    if (s->cmd->ci < 0) {
	return 0; /* nothing there, or a comment */
    }

    if (s->cmd->ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	err = gretl_loop_append_line(s, dset);
	if (err) {
	    errmsg(err, s->prn);
	    return err;
	}
	return 0;
    } 

    s->pmod = NULL; /* be on the safe side */

    if (s->cmd->ci == FUNCERR) {
	err = E_FUNCERR;
    } else {
	/* note: error messages may be printed to s->prn */
	err = gretl_cmd_exec(s, dset);
    }

    return err;
}

/**
 * get_command_index:
 * @line: command line.
 * @cmode: compilation mode: LOOP or FUNC
 * @cmd: pointer to gretl command struct.
 *
 * Parse @line and assign to the %ci field of @cmd the index number of
 * the command embedded in @line.  Note: this is a "lite" version of
 * parse_command_line().  It is used when commands are being stacked
 * for execution within a loop.  Command options are not parsed out of
 * @line.
 *
 * Returns: 1 on error, otherwise 0.
 */

int get_command_index (char *line, int cmode, CMD *cmd)
{
    int err = 0;

    gretl_cmd_clear(cmd);

#if CMD_DEBUG
    fprintf(stderr, "get_command_index: line='%s'\n", line);
#endif

    if ((cmd->context == FOREIGN || cmd->context == MPI) &&
	!ends_foreign_block(line)) {
	cmd->opt = OPT_NONE;
	cmd->ci = cmd->context;
	return 0;
    }

    if (filter_comments(line, cmd)) {
	return 0;
    }

    err = real_parse_command(line, cmd, NULL, cmode, NULL);

    if (!err && cmd->ci == 0) {
	/* maybe genr via series name? */
	const char *s = cmd->toks[0].s;

	if (s != NULL) {
	    if (*s == '$' || *s == '@') s++;
	    if (strlen(s) == gretl_namechar_spn(s)) {
		cmd->ci = GENR;
	    }
	}  
    }

    if (!err && cmd->ci == 0) {
	/* FIXME watch out for fallout! (2012-03-01) */
	cmd->ci = CMD_NULL;
	err = E_PARSE;
    }

    if (err) {
	return err;
    }

    if (cmd->ci == END) {
	cmd->context = 0;
    } else if (cmd->context) {
	cmd->ci = cmd->context;
    }	

    if (cmd->ci == NLS || cmd->ci == MLE ||
	cmd->ci == GMM || cmd->ci == FOREIGN ||
	cmd->ci == KALMAN || cmd->ci == PLOT ||
	cmd->ci == MPI) {
	cmd->context = cmd->ci;
    }

#if CMD_DEBUG
    fprintf(stderr, " get_command_index: cmd->ci set to %d\n", cmd->ci);
#endif

    return 0;
}

void gretl_cmd_set_context (CMD *cmd, int ci)
{
    cmd->context = ci;
}

void gretl_cmd_destroy_context (CMD *cmd)
{
    if (cmd->context == FOREIGN || cmd->context == MPI) {
	foreign_destroy();
    }
    cmd->context = 0;
    *cmd->savename = '\0';
}

gretlopt gretl_cmd_get_opt (const CMD *cmd)
{
    return cmd->opt;
}

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt)
{
    cmd->opt = opt;
}

const char *gretl_cmd_get_savename (CMD *cmd)
{
    return cmd->savename;
}

void gretl_exec_state_init (ExecState *s,
			    ExecFlags flags,
			    char *line,
			    CMD *cmd,
			    MODEL *model, 
			    PRN *prn)
{
    s->flags = flags;

    s->line = line;
    if (s->line != NULL) {
	*s->line = '\0';
    }    

    s->cmd = cmd;
    if (s->cmd != NULL) {
	s->cmd->ci = 0;
    }    

    *s->runfile = '\0';

    s->model = model;
    s->prn = prn;

    s->pmod = NULL;
    s->sys = NULL;
    s->rset = NULL;
    s->var = NULL;
    s->in_comment = 0;
    s->padded = 0;

    if (flags == FUNCTION_EXEC) {
	/* On entry to function execution we check if there's
	   a 'last model' in place. If so, we want to make
	   this invisible within the function, but set things
	   up so that we can restore it as last model on
	   exit from the function -- the idea being that
	   excuting a function should not change the 'last
	   model' state at caller level. To achieve this we
	   need to take out a 'private' reference to the
	   model, stored in the ExecState, and then remove
	   it from last model position for the present.
	*/
	s->prev_model = get_last_model(&s->prev_type);
	if (s->prev_model != NULL) {
	    gretl_object_ref(s->prev_model, s->prev_type);
	    set_as_last_model(NULL, GRETL_OBJ_NULL);
	}
	s->prev_model_count = get_model_count();
    } else {
	s->prev_model = NULL;
	s->prev_type = GRETL_OBJ_NULL;
	s->prev_model_count = -1;
    }

    s->submask = NULL;
    s->callback = NULL;
}

void function_state_init (CMD *cmd, ExecState *state, int *indent0)
{
    cmd->list = NULL;
    cmd->auxlist = NULL;
    cmd->param = NULL;
    cmd->parm2 = NULL;
    /* FIXME tokenize more needed? */

    state->cmd = NULL;
    state->model = NULL;
    state->submask = NULL;

    state->padded = 0;

    *indent0 = gretl_if_state_record();
}

static EXEC_CALLBACK gui_callback;

void gretl_exec_state_set_callback (ExecState *s, EXEC_CALLBACK callback,
				    gretlopt opt)
{
    s->callback = callback;
    s->pmod = NULL;

    if (opt & OPT_G) {
	gui_callback = callback;
    }
}

EXEC_CALLBACK get_gui_callback (void)
{
    return gui_callback;
}

void gretl_exec_state_clear (ExecState *s)
{
    gretl_cmd_free(s->cmd);

    if (s->flags & FUNCTION_EXEC) {
	/* Restore whatever was the 'last model' before 
	   function execution. Note that this includes
	   the case where there was no 'last model', in
	   which case we restore the null state. Drop
	   the extra refcount for the model we put into
	   last model position (if any), so we don't end 
	   up leaking memory.
	*/
	set_as_last_model(s->prev_model, s->prev_type);
	if (s->prev_model != NULL) {
	    gretl_object_unref(s->prev_model, s->prev_type);
	}
	/* restore the previous model count */
	if (s->prev_model_count >= 0) {
	    set_model_count(s->prev_model_count);
	}
    }

    destroy_working_model(s->model);

    s->prev_model = NULL;
    s->prev_type = GRETL_OBJ_NULL;
    s->prev_model_count = -1;

    free_subsample_mask(s->submask);
}

void gretl_exec_state_uncomment (ExecState *s)
{
    s->in_comment = 0;
    s->cmd->flags &= ~CMD_IGNORE;
}

void gretl_exec_state_transcribe_flags (ExecState *s, CMD *cmd)
{
    s->in_comment = (cmd_ignore(cmd))? 1 : 0;
}

void gretl_exec_state_set_model (ExecState *s, MODEL *pmod)
{
    s->pmod = pmod;
}

int process_command_error (ExecState *s, int err)
{
    int ret = err;

    if (err) {
	if (gretl_compiling_function() ||
	    gretl_compiling_loop()) {
	    ; /* pass the error through */
	} else if (s->cmd->flags & CMD_CATCH) {
	    /* local "continue on error" */
	    set_gretl_errno(err);
	    s->cmd->flags ^= CMD_CATCH;
	    ret = 0;
	}
    }

    if (ret && print_redirection_level(s->prn) > 0) {
	print_end_redirection(s->prn);
	pputs(s->prn, _("An error occurred when 'outfile' was active\n"));
    }

    return ret;
}

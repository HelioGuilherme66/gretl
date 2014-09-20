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

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_help.h"

static int maybe_need_recode (void)
{
    const gchar *cset = NULL;
    int utf = g_get_charset(&cset);

    return !utf;
}

static int recode_print_line (const char *s, PRN *prn)
{
    gchar *trs;
    gsize bytes;
    GError *err = NULL;

    trs = g_locale_from_utf8(s, -1, NULL, &bytes, &err);  

    if (err != NULL) {
	pprintf(prn, "%s\n", err->message);
	g_error_free(err);
    } else {
	pputs(prn, trs);
    }

    if (trs != NULL) {
	g_free(trs);
    }

    return 0;
}

/* list the topics available in the functions help file */

static int func_help_topics (const char *helpfile, PRN *prn)
{
    char line[128], word[12];
    FILE *fp;
    int j, k;

    if ((fp = gretl_fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return E_FOPEN;
    } 

    j = 1;
    k = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
	if (!strncmp(line, "## ", 3)) {
	    /* sub-heading */
	    tailstrip(line);
	    if (k++ > 0) {
		pputc(prn, '\n');
	    }
	    pprintf(prn, "\n%s:\n", line + 3);
	    j = 1;
	} else if (*line == '#') {
	    /* actual entry */
	    sscanf(line + 2, "%10s", word);
	    pprintf(prn, "%-10s", word);
	    if (j % 7 == 0) {
		pputc(prn, '\n');
	    } else {
		pputc(prn, ' ');
	    }
	    j++;
	}
    } 

    pputs(prn, _("\n\nFor help on a specific function, type: help funname"));
    pputs(prn, _(" (e.g. help qrdecomp)\n"));

    fclose(fp);
    
    return 0;
}

static void output_help_line (const char *line, PRN *prn, int recode)
{
    if (recode > 0) {
	recode_print_line(line, prn);
    } else {
	pputs(prn, line);
    }
}

/* check in the CLI helpfile for a line of the form "  @s: ...", 
   where @s has been recognized as a libset variable */

static int got_setvar_line (const char *s, int n, char *line)
{
    if (!strncmp(line, "  ", 2) &&
	!strncmp(line + 2, s, n) &&
	line[2 + n] == ':') {
	return 1;
    } else {
	return 0;
    }
}

#define HELPLEN 128

/* special case: the user has done "help set @setvar" */	

static int do_set_help (const char *setvar, FILE *fp, 
			int recode, PRN *prn)
{
    char s[9], line[HELPLEN];
    int n = strlen(setvar);
    int ok = 0, count = 0;

    while (count < 2 && fgets(line, HELPLEN, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%8s", s);
	if (!strcmp(s, "set")) {
	    while (count < 2 && fgets(line, HELPLEN, fp)) {
		if (got_setvar_line(setvar, n, line)) {
		    pputc(prn, '\n');
		    output_help_line(line + 2, prn, recode);
		    count++;
		} else if (count > 0) {
		    if (string_is_blank(line)) {
			/* reached the end of the item */
			pputc(prn, '\n');
			count++;
		    } else {
			output_help_line(line + 2, prn, recode);
		    }
		} else if (*line == '#') {
		    /* overshot: not found */
		    count = 999;
		}
	    }
	}
    }

    ok = count > 0 && count < 999;

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), setvar);
    }

    return ok;
}

static int output_help_text (const char *cmd, FILE *fp,
			     int recode, PRN *prn)
{
    char word[16], line[HELPLEN];
    int ok = 0;

    while (fgets(line, sizeof line, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%15s", word);
	if (!strcmp(cmd, word)) {
	    ok = 1;
	    pprintf(prn, "\n%s\n", word);
	    while (fgets(line, sizeof line, fp)) {
		if (*line == '#') {
		    break;
		}
		output_help_line(line, prn, recode);
	    }
	    break;
	}
    }

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), cmd);
    }

    return ok;
}

static void do_help_on_help (PRN *prn)
{
    int i, j;

    pputs(prn, _("\nValid gretl commands are:\n"));
    j = 1;
    for (i=1; i<NC; i++) {
	if (HIDDEN_COMMAND(i)) {
	    continue;
	}
	pprintf(prn, "%-9s", gretl_command_word(i));
	if (j % 8 == 0) {
	    pputc(prn, '\n');
	} else {
	    pputc(prn, ' ');
	}
	j++;
    }

    pputs(prn, _("\n\nFor help on a specific command, type: help cmdname"));
    pputs(prn, _(" (e.g. help smpl)\n"));
    pputs(prn, _("You can also do 'help functions' for a list of functions\n"));
}

/**
 * cli_help:
 * @hlpword: the word (usually a command) on which help is wanted.
 * @param: parameter usable when @hlpword is "set".
 * @opt: may include %OPT_F to give priority to functions
 * rather than commands.
 * @prn: pointer to gretl printing struct.
 *
 * Searches in the gretl helpfile for help on @hlpword and, 
 * if help is found, prints it to @prn.  If @hlpword is %NULL, 
 * lists the valid commands.
 *
 * Returns: 0 on success, 1 if the helpfile was not found or the
 * requested topic was not found.
 */

int cli_help (const char *hlpword, const char *param,
	      gretlopt opt, PRN *prn)
{
    static int recode = -1;
    char helpfile[FILENAME_MAX];
    FILE *fp;
    int noword, funhelp = (opt & OPT_F);
    int ok = 0;

    noword = (hlpword == NULL || *hlpword == '\0');

    /* @param should be NULL unless hlpword is "set" */
    if ((noword || strcmp(hlpword, "set")) && param != NULL) {
	return 1;
    }

    if (noword && !funhelp) {
	do_help_on_help(prn);
	return 0; /* handled */
    }

    if ((noword && funhelp) || !strcmp(hlpword, "functions")) {
	sprintf(helpfile, "%s%s", gretl_home(), _("genrcli.hlp"));
	return func_help_topics(helpfile, prn); /* handled */
    }

    if (!funhelp) {
	int ci = gretl_command_number(hlpword);
	
	if (ci == SET && param != NULL) {
	    ok = is_libset_var(param);
	    if (!ok) {
		pprintf(prn, _("%s: sorry, no help available.\n"), param);
		return 1;
	    }
	} else {
	    ok = ci > 0;
	}
    } 

    if (ok) {
	strcpy(helpfile, helpfile_path(GRETL_CLI_HELPFILE));
    } else if (genr_function_word(hlpword)) {
	sprintf(helpfile, "%sgenrcli.hlp", gretl_home());
    } else if (gretl_is_public_user_function(hlpword)) {
	return user_function_help(hlpword, OPT_NONE, prn);
    } else {
	pprintf(prn, _("\"%s\" is not a gretl command.\n"), hlpword);
	return 1;
    }

    if ((fp = gretl_fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return 1;
    } 

    if (!gretl_in_gui_mode() && recode < 0) {
	recode = maybe_need_recode();
    }

    /* actually output the relevant help text */

    if (param != NULL) {
	ok = do_set_help(param, fp, recode, prn);
    } else {
	ok = output_help_text(hlpword, fp, recode, prn);
    }

    fclose(fp);

    return 0;
}

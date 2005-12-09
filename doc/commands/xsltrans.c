/*
 *  Copyright (c) 2004 by Allin Cottrell
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

/* Formatter for gretl commands info stored as "generic" XML: takes
   a purely content-based XML representation of the info relating to
   the gretl commands (conforming to gretl_commands.dtd) and uses
   XSL transformation to turn this into output suitable for:

   * the "cmdref" chapter of the gretl manual (docbook XML); and 
   * the "online" help files.
   
   Uses the XSL stylesheets gretlman.xsl and gretltxt.xsl.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

enum {
    OUTPUT_ALL,
    OUTPUT_DOCBOOK,
    OUTPUT_DOCBOOK_STANDALONE,
    OUTPUT_HLP,
    OUTPUT_CLI_HLP,
    OUTPUT_CMD_HLP,
    OUTPUT_GUI_HLP,
    OUTPUT_TEST
};

#define ROOTNODE "commandlist"
#define UTF const xmlChar *

static void full_fname (const char *fname, const char *dir,
			char *targ)
{
    if (dir == NULL || *dir == '\0') {
	strcpy(targ, fname);
    } else {
	sprintf(targ, "%s/%s", dir, fname);
    }
}

static void build_params (char const **params, int output, 
			  int xrefs, const char *lang)
{
    int i = 0;

    if (output == OUTPUT_CMD_HLP && !xrefs) {
	return;
    }

    if (strcmp(lang, "en")) {
	params[0] = "lang";
	params[1] = lang;
	i = 2;
    }    
    
    if (output == OUTPUT_DOCBOOK_STANDALONE) {
	params[i++] = "standalone";
	params[i++] = "\"true\"";
    } else if (output == OUTPUT_GUI_HLP) {
	params[i++] = "hlp";
	params[i++] = "\"gui\"";
    } 

    if (xrefs) {
	params[i++] = "xrefs";
	params[i++] = "\"true\"";
    }
	
}

int apply_xslt (xmlDocPtr doc, int output, int xrefs, const char *lang, 
		const char *docdir)
{
    xsltStylesheetPtr style;
    xmlDocPtr result;
    char styname[FILENAME_MAX];
    char const *xsl_params[12] = {0};
    FILE *fp;
    int err = 0;

    xmlIndentTreeOutput = 1;

    if (output == OUTPUT_TEST) {
	/* not ready yet: may some day produce nicely formatted help
	   text for use with GtkTextView */
	full_fname("gretlhlp.xsl", docdir, styname);
	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_DOCBOOK, 0, lang);
	    result = xsltApplyStylesheet(style, doc, xsl_params);
	    if (result == NULL) {
		err = 1;
	    } else {
		fp = fopen("cmdlist2.txt", "w");
		if (fp == NULL) {
		    err = 1;
		} else {
		    xsltSaveResultToFile(fp, result, style);
		    fclose(fp);
		}
		xsltFreeStylesheet(style);
		xmlFreeDoc(result);
	    }	    
	}
    }	

    /* make "full" DocBook XML output */
    if (output == OUTPUT_ALL || output == OUTPUT_DOCBOOK) {
	full_fname("gretlman.xsl", docdir, styname);
	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_DOCBOOK, 0, lang);
	    result = xsltApplyStylesheet(style, doc, xsl_params);
	    if (result == NULL) {
		err = 1;
	    } else {
		fp = fopen("cmdlist.xml", "w");
		if (fp == NULL) {
		    err = 1;
		} else {
		    xsltSaveResultToFile(fp, result, style);
		    fclose(fp);
		}
		xsltFreeStylesheet(style);
		xmlFreeDoc(result);
	    }	    
	}
    }

    /* make plain text "hlp" output */
    if (output == OUTPUT_ALL || output == OUTPUT_HLP) {
	full_fname("gretltxt.xsl", docdir, styname);
	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    /* cli version */
	    build_params(xsl_params, OUTPUT_CLI_HLP, 0, lang);
	    result = xsltApplyStylesheet(style, doc, xsl_params);
	    if (result == NULL) {
		err = 1;
	    } else {
		fp = fopen("clilist.txt", "w");
		if (fp == NULL) {
		    err = 1;
		} else {
		    xsltSaveResultToFile(fp, result, style);
		    fclose(fp);
		}
		xmlFreeDoc(result);
	    }
	    /* script help, gui version (gtk2 only) */
	    build_params(xsl_params, OUTPUT_CMD_HLP, xrefs, lang);
	    result = xsltApplyStylesheet(style, doc, xsl_params);
	    if (result == NULL) {
		err = 1;
	    } else {
		fp = fopen("cmdlist.txt", "w");
		if (fp == NULL) {
		    err = 1;
		} else {
		    xsltSaveResultToFile(fp, result, style);
		    fclose(fp);
		}
		xmlFreeDoc(result);
	    }
	    /* gui help */
	    build_params(xsl_params, OUTPUT_GUI_HLP, xrefs, lang);
	    result = xsltApplyStylesheet(style, doc, xsl_params);
	    if (result == NULL) {
		err = 1;
	    } else {
		fp = fopen("guilist.txt", "w");
		if (fp == NULL) {
		    err = 1;
		} else {
		    xsltSaveResultToFile(fp, result, style);
		    fclose(fp);
		}
		xmlFreeDoc(result);
	    }
	    xsltFreeStylesheet(style);
	}
    }

    return err;
}

char *get_abbreviated_lang (char *lang, const char *full_lang)
{
    if (!strcmp(full_lang, "italian")) {
	strcpy(lang, "'it'");
    }
    else if (!strcmp(full_lang, "spanish")) {
	strcpy(lang, "'es'");
    }
    else if (!strcmp(full_lang, "french")) {
	strcpy(lang, "'fr'");
    }

    return lang;
}

int parse_commands_data (const char *fname, int output, 
			 int xrefs, const char *docdir) 
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    char *tmp = NULL;
    char lang[8] = "en";
    int err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	err = 1;
	goto bailout;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
	err = 1;
	goto bailout;
    }

    if (xmlStrcmp(cur->name, (UTF) ROOTNODE)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", ROOTNODE);
	err = 1;
	goto bailout;
    }

    tmp = (char *) xmlGetProp(cur, (UTF) "language");
    if (tmp != NULL) {
	get_abbreviated_lang(lang, tmp);
	free(tmp);
    }

    apply_xslt(doc, output, xrefs, lang, docdir);

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

static char *get_docdir (char *ddir, const char *fname)
{
    char *p;

    p = strrchr(fname, '/');
    if (p == NULL) return NULL;

    strcpy(ddir, fname);
    p = strrchr(ddir, '/');
    *p = 0;

    return ddir;
}

static void usage (void)
{
    fputs("Please give the name of an XML file to parse\n", stderr);
    exit(EXIT_FAILURE);
}

int main (int argc, char **argv)
{
    const char *fname = NULL;
    char docdir[FILENAME_MAX];
    int output = OUTPUT_ALL;
    int xrefs = 0;
    int i, err;

    if (argc < 2) {
	usage();
    }

    *docdir = '\0';

    for (i=1; i<argc; i++) {
	if (!strcmp(argv[i], "--docbook")) {
	    output = OUTPUT_DOCBOOK;
	} else if (!strcmp(argv[i], "--docbook-standalone")) {
	    output = OUTPUT_DOCBOOK_STANDALONE;
	} else if (!strcmp(argv[i], "--hlp")) {
	    output = OUTPUT_HLP;
	} else if (!strcmp(argv[i], "--all")) {
	    output = OUTPUT_ALL;
	} else if (!strcmp(argv[i], "--xrefs")) {
	    xrefs = 1;
	} else if (!strcmp(argv[i], "--test")) {
	    output = OUTPUT_TEST;
	} else if (!strncmp(argv[i], "--docdir=", 9)) {
	    strcpy(docdir, argv[i] + 9);
	} else {
	    fname = argv[i];
	}
    }

    if (fname == NULL) {
	usage();
    }

    if (*docdir == '\0') {
	get_docdir(docdir, fname);
    }

    fprintf(stderr, "%s: input file '%s', docdir '%s'\n", 
	    argv[0], fname, docdir);

    err = parse_commands_data(fname, output, xrefs, docdir);

    return err;
}

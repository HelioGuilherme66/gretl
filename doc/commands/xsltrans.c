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
    OUTPUT_GUI_HLP
};

#define ROOTNODE "commandlist"
#define UTF const xmlChar *

static void full_fname (const char *fname, const char *dir,
			char *targ)
{
    if (dir == NULL) {
	strcpy(targ, fname);
    } else {
	sprintf(targ, "%s%s", dir, fname);
    }
}

static void build_params (char const **params, int output, 
			  const char *lang)
{
    int i = 0;

    if (strcmp(lang, "en")) {
	params[0] = "lang";
	params[1] = lang;
	i = 2;
    }    
    
    if (output == OUTPUT_DOCBOOK_STANDALONE) {
	params[i++] = "standalone";
	params[i++] = "\"true\"";
    }

    else if (output == OUTPUT_GUI_HLP) {
	params[i++] = "hlp";
	params[i++] = "\"gui\"";
    }
}

int apply_xslt (xmlDocPtr doc, int output, const char *lang, 
		const char *docdir)
{
    xsltStylesheetPtr style;
    xmlDocPtr result;
    char styname[FILENAME_MAX];
    char const *xsl_params[5] = {0};
    FILE *fp;
    int err = 0;

    xmlIndentTreeOutput = 1;

    /* make "full" DocBook XML output */
    if (output == OUTPUT_ALL || output == OUTPUT_DOCBOOK) {
	full_fname("gretlman.xsl", docdir, styname);
	style = xsltParseStylesheetFile(styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_DOCBOOK, lang);
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

    /* make "standalone" DocBook XML output */
    if (output == OUTPUT_ALL || output == OUTPUT_DOCBOOK_STANDALONE) {
	full_fname("gretlman.xsl", docdir, styname);
	style = xsltParseStylesheetFile(styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_DOCBOOK_STANDALONE, lang);
	    result = xsltApplyStylesheet(style, doc, xsl_params);
	    if (result == NULL) {
		err = 1;
	    } else {
		fp = fopen("cmdlist_standalone.xml", "w");
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
	style = xsltParseStylesheetFile(styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    /* cli version */
	    build_params(xsl_params, OUTPUT_CLI_HLP, lang);
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
	    /* gui version */
	    build_params(xsl_params, OUTPUT_GUI_HLP, lang);
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
			 const char *docdir) 
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

    tmp = xmlGetProp(cur, (UTF) "language");
    if (tmp != NULL) {
	get_abbreviated_lang(lang, tmp);
	free(tmp);
    }

    apply_xslt(doc, output, lang, docdir);

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

static char *get_docdir (const char *fname)
{
    char *docdir, *p;

    p = strrchr(fname, '/');
    if (p == NULL) return NULL;

    docdir = malloc(strlen(fname) + 1);
    if (docdir == NULL) return NULL;

    strcpy(docdir, fname);
    p = strrchr(docdir, '/');
    *(p + 1) = 0;

    return docdir;
}

int main (int argc, char **argv)
{
    const char *fname;
    char *docdir;
    int output = OUTPUT_ALL;
    int err;

    if (argc < 2) {
	fputs("Please give the name of an XML file to parse\n", stderr);
	exit(EXIT_FAILURE);
    }

    if (argc == 3) {
	fname = argv[2];
	if (!strcmp(argv[1], "--docbook")) {
	    output = OUTPUT_DOCBOOK;
	}
	if (!strcmp(argv[1], "--docbook-standalone")) {
	    output = OUTPUT_DOCBOOK_STANDALONE;
	}
	else if (!strcmp(argv[1], "--hlp")) {
	    output = OUTPUT_HLP;
	}
    } else {
	fname = argv[1];
    }

    docdir = get_docdir(fname);

    err = parse_commands_data(fname, output, docdir);

    return err;
}

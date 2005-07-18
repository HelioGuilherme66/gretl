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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "gretl.h"
#include "graph_page.h"

#ifdef G_OS_WIN32 
# include <io.h>
#else
# include <unistd.h>
# include <sys/stat.h>
#endif

#define GRAPHS_MAX 8

enum {
    PS_OUTPUT,
    PDF_OUTPUT
} output_types;

typedef struct _graphpage graphpage;

struct _graphpage {
    int output;
    int color;
    int ngraphs;
    char **fnames;
};

static graphpage gpage;
static char gpage_base[FILENAME_MAX];
static char gpage_tex_base[FILENAME_MAX];

static void gpage_filenames_init (const char *base)
{
    if (base == NULL) {
	strcpy(gpage_base, paths.userdir);
	strcat(gpage_base, "gretl_graphpage");
	strcpy(gpage_tex_base, "gretl_graphpage");
    } else {
	const char *p;

	strcpy(gpage_base, base);
	if (has_suffix(gpage_base, ".tex")) {
	    gpage_base[strlen(gpage_base) - 4] = '\0';
	}
	p = strrchr(base, SLASH);
	if (p != NULL) {
	    strcpy(gpage_tex_base, p + 1);
	} else {
	    strcpy(gpage_tex_base, base);
	}
    }
}

static char *gpage_fname (const char *ext, int i)
{
    static char fname[MAXLEN];

    strcpy(fname, gpage_base);

    if (i > 0) {
	char num[6];

	sprintf(num, "_%d", i);
	strcat(fname, num);
    }

    if (ext != NULL) {
	strcat(fname, ext);
    }

    return fname;
}

static void gpage_errmsg (char *msg, int gui)
{
    if (gui) {
	errbox(msg);
    } else {
	gretl_errmsg_set(msg);
    }
}

static void graph_page_init (void)
{
    gpage.color = 0;
    gpage.ngraphs = 0;
    gpage.fnames = NULL;
}

static void doctop (FILE *fp)
{
    int letter = in_usa();

    fprintf(fp, "\\documentclass%s{article}\n", (letter)? "" : "[a4paper]");
    fprintf(fp, "\\usepackage[%s]{graphicx}\n", 
	    (gpage.output == PS_OUTPUT)? "dvips" : "pdftex");
}

/* A4 is 210mm * 297mm */

static int geomline (int ng, FILE *fp)
{
    double width = 7.0, height = 10.0;
    double tmarg = 0.5, lmarg = 0.75;
    char unit[3];

    if (in_usa()) {
	strcpy(unit, "in");
    } else {
	strcpy(unit, "mm");
	width = 190;
	height = 277.0;
	tmarg = lmarg = 10.0;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    fprintf(fp, "\\usepackage[body={%g%s,%g%s},"
	    "top=%g%s,left=%g%s,nohead]{geometry}\n\n",
	    width, unit, height, unit, 
	    tmarg, unit, lmarg, unit);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    return 0;
}

static void common_setup (FILE *fp)
{
    fputs("\\begin{document}\n\n"
	  "\\thispagestyle{empty}\n\n"
	  "\\vspace*{\\stretch{1}}\n\n"
	  "\\begin{center}\n", fp);
}

static int oddgraph (int ng, int i)
{
    return (ng % 2) && (i == ng - 1);
}

static int tex_graph_setup (int ng, FILE *fp)
{
    char fname[FILENAME_MAX];
    double scale = 1.0;
    double vspace = 1.0;
    int i;

    if (ng > GRAPHS_MAX) {
	fprintf(stderr, "ng (%d) out of range\n", ng);
	return 1;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    if (ng == 1) {
	sprintf(fname, "%s_1", gpage_tex_base);
	fprintf(fp, "\\includegraphics[scale=1.2]{%s}\n", fname);
    } else if (ng == 2) {
	sprintf(fname, "%s_%d", gpage_tex_base, 1);
	fprintf(fp, "\\includegraphics{%s}\n\n", fname);
	fprintf(fp, "\\vspace{%gin}\n\n", vspace);
	sprintf(fname, "%s_%d", gpage_tex_base, 2);
	fprintf(fp, "\\includegraphics{%s}\n\n", fname);
    } else if (ng == 3) {
	scale = 0.9;
	vspace = 0.25;
	for (i=0; i<3; i++) {
	    sprintf(fname, "%s_%d", gpage_tex_base, i + 1);
	    fprintf(fp, "\\includegraphics[scale=%g]{%s}\n\n",
		    scale, fname);
	    fprintf(fp, "\\vspace{%gin}\n", vspace);
	}
    } else {
	if (ng > 6) {
	    scale = 0.85;
	    vspace = 0.20;
	} else {
	    scale = 0.9;
	    vspace = 0.25;
	}	    
	fputs("\\begin{tabular}{cc}\n", fp);
	for (i=0; i<ng; i++) {
	    sprintf(fname, "%s_%d", gpage_tex_base, i + 1);
	    if (oddgraph(ng, i)) {
		fprintf(fp, "\\multicolumn{2}{c}{\\includegraphics[scale=%g]{%s}}",
			scale, fname);
	    } else {
		fprintf(fp, "\\includegraphics[scale=%g]{%s}",
			scale, fname);
		if (i % 2 == 0) {
		    fputs(" &\n  ", fp);
		} else if (i < ng - 1) {
		    fprintf(fp, " \\\\ [%gin]\n", vspace);
		}
	    }
	}
	fputs("\\end{tabular}\n", fp);
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif 

    return 0;
}

static void common_end (FILE *fp)
{
    fputs("\\end{center}\n\n"
	  "\\vspace*{\\stretch{2}}\n\n"
	  "\\end{document}\n", fp);
}

static int make_graphpage_tex (void)
{
    char *fname;
    FILE *fp;
    int err = 0;

    fname = gpage_fname(".tex", 0);

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    doctop(fp);

    err = geomline(gpage.ngraphs, fp);

    if (!err) {
	common_setup(fp);
	err = tex_graph_setup(gpage.ngraphs, fp);
    }

    if (!err) {
	common_end(fp);
    }

    fclose(fp);

    return err;
}

int graph_page_add_file (const char *fname)
{
    char **fnames;
    int ng = gpage.ngraphs + 1;

    if (ng > GRAPHS_MAX) {
	gpage_errmsg(_("The graph page is full"), 1);
	return 1;
    }

    fnames = myrealloc(gpage.fnames, ng * sizeof *fnames);
    if (fnames == NULL) return 1;

    fnames[ng - 1] = g_strdup(fname);
    gpage.fnames = fnames;
    gpage.ngraphs = ng;

    return 0;
}

static int gnuplot_compile (const char *fname)
{
    char plotcmd[MAXLEN];
    int err = 0;

# ifdef WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", paths.gnuplot, fname);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
# else
    sprintf(plotcmd, "%s \"%s\"", paths.gnuplot, fname);
    err = gretl_spawn(plotcmd);  
# endif /* WIN32 */

    return err;
}

static int gp_make_outfile (const char *gfname, int i, double scale)
{
    char line[128];
    char *fname;
    FILE *fp, *fq;
    int err = 0;

    fp = gretl_fopen(gfname, "r");
    if (fp == NULL) return 1;

    fname = gpage_fname(".plt", 0);

    fq = gretl_fopen(fname, "w");
    if (fq == NULL) {
	fclose(fp);
	return 1;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    
    if (gpage.output == PDF_OUTPUT) {
	fputs("set term pdf\n", fq);
	if (scale != 1.0) {
	    fprintf(fq, "set size %g,%g\n", scale, scale);
	}	
	fname = gpage_fname(".pdf", i);
    } else {
#ifdef ENABLE_NLS
	fprint_gnuplot_encoding("postscript", fq);
#endif
	fprintf(fq, "set term postscript eps%s\n", (gpage.color)? " color" : "");
	if (scale != 1.0) {
	    fprintf(fq, "set size %g,%g\n", scale, scale);
	}
	fname = gpage_fname(".ps", i);
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fprintf(fq, "set output '%s'\n", fname);

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "set out", 7)) continue;
	if (!strncmp(line, "set term", 8)) continue;
	if (!strncmp(line, "set size", 8)) continue;
	fputs(line, fq);
    }

    fclose(fp);
    fclose(fq);

    fname = gpage_fname(".plt", 0);
    err = gnuplot_compile(fname);

    return err;
}

#if defined(G_OS_WIN32)

static int get_dvips_path (char *path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, "dvips.exe", NULL, MAXLEN, path, &p);

    return (ret == 0);
}

#elif !defined(OLD_GTK)

#include <signal.h>

static int spawn_dvips (char *texsrc)
{
    GError *error = NULL;
    gchar *sout = NULL;
    gchar *argv[5];
    char outfile[32];
    int ok, status;
    int ret = 0;

    sprintf(outfile, "%s.ps", texsrc);

    argv[0] = "dvips";
    argv[1] = "-o";
    argv[2] = outfile;
    argv[3] = texsrc;
    argv[4] = NULL;

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (paths.userdir, /* working dir */
		       argv,
		       NULL,    /* envp */
		       G_SPAWN_SEARCH_PATH,
		       NULL,    /* child_setup */
		       NULL,    /* user_data */
		       &sout,   /* standard output */
		       NULL,    /* standard error */
		       &status, /* exit status */
		       &error);

    if (!ok) {
	errbox(error->message);
	g_error_free(error);
	ret = LATEX_EXEC_FAILED;
    } else if (status != 0) {
	gchar *errmsg;

	errmsg = g_strdup_printf("%s\n%s", 
				 _("Failed to process TeX file"),
				 sout);
	errbox(errmsg);
	g_free(errmsg);
	ret = 1;
    }

    if (sout != NULL) g_free(sout);

    return ret;
}

#endif

int dvips_compile (char *texshort)
{
#ifdef G_OS_WIN32
    static char dvips_path[MAXLEN];
#endif
#if defined(G_OS_WIN32) || defined(OLD_GTK)
    char tmp[MAXLEN];
#endif
    int err = 0;

#if defined(G_OS_WIN32)
    if (*dvips_path == 0 && get_dvips_path(dvips_path)) {
	DWORD dw = GetLastError();
	win_show_error(dw);
	return 1;
    }

    sprintf(tmp, "\"%s\" -o %s.ps %s", dvips_path, texshort, texshort);
    if (winfork(tmp, paths.userdir, SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE)) {
	return 1;
    }
#elif defined(OLD_GTK)
    sprintf(tmp, "cd \"%s\" && dvips -o %s.ps %s", paths.userdir, texshort, texshort);
    err = system(tmp);
#else
    err = spawn_dvips(texshort);
#endif 

    return err;
}

static int latex_compile_graph_page (void)
{
    char *fname;
    int err;

    /* ensure we don't get stale output */
    if (gpage.output == PS_OUTPUT) {
	fname = gpage_fname(".ps", 0);
    } else {
	fname = gpage_fname(".pdf", 0);
    }
    remove(fname);

    err = latex_compile(gpage_base);

    if (err == LATEX_ERROR) {
	char *fname = gpage_fname(".log", 0);

	view_file(fname, 0, 1, 78, 350, VIEW_FILE);
    }

    if (gpage.output == PS_OUTPUT && !err) {
	err = dvips_compile(gpage_base);
    }    

    return err;
}

static int make_gp_output (void)
{
    char *fname;
    double scale = 1.0;
    int i;
    int err = 0;

    if (gpage.ngraphs == 3) {
	scale = 0.8;
    } else if (gpage.ngraphs > 3) {
	scale = 0.75;
    }

    for (i=0; i<gpage.ngraphs && !err; i++) {
	err = gp_make_outfile(gpage.fnames[i], i + 1, scale);
    }

    fname = gpage_fname(".plt", 0);
    remove(fname);

    return err;
}

static int real_display_gpage (void)
{
#ifndef G_OS_WIN32
    char *viewer;
#endif
    char *fname;
    int err = 0;

    if (gpage.output == PDF_OUTPUT) {
	fname = gpage_fname(".pdf", 0);
    } else {
	fname = gpage_fname(".ps", 0);
    }

#ifdef G_OS_WIN32
    if ((int) ShellExecute(NULL, "open", fname, NULL, NULL, SW_SHOW) <= 32) {
	DWORD dw = GetLastError();
	win_show_error(dw);
	err = 1;
    }
#else
    viewer = (gpage.output == PDF_OUTPUT)? viewpdf : viewps;
    err = gretl_fork(viewer, fname);
#endif

    return err;
}

static void gpage_cleanup (void)
{
    char *fname;
    int i;

    for (i=0; i<gpage.ngraphs; i++) {
	if (gpage.output == PDF_OUTPUT) {
	    fname = gpage_fname(".pdf", i + 1);
	} else {
	    fname = gpage_fname(".ps", i + 1);
	}
	remove(fname);
    }

    fname = gpage_fname(".tex", 0);
    remove(fname);
    fname = gpage_fname(".dvi", 0);
    remove(fname);
    fname = gpage_fname(".log", 0);
    remove(fname);
    fname = gpage_fname(".aux", 0);
    remove(fname);
}

int display_graph_page (void)
{
    char *latex_orig = NULL;
    int err = 0;

    if (gpage.ngraphs == 0) {
	gpage_errmsg(_("The graph page is empty"), 1);
	return 1;
    }

    gpage_filenames_init(NULL);

    if (!strncmp(latex, "pdf", 3)) {
	if (gnuplot_has_pdf()) {
	    gpage.output = PDF_OUTPUT;
	} else {
	    latex_orig = g_strdup(latex);
	    strcpy(latex, latex_orig + 3);
	    gpage.output = PS_OUTPUT;
	}
    } else {
	gpage.output = PS_OUTPUT;
    }

    /* write the LaTeX driver file */
    err = make_graphpage_tex();

    if (!err) {
	/* transform individual plot files and compile 
	   using gnuplot */
	err = make_gp_output();
    }

    if (!err) {
	err = latex_compile_graph_page();
    }

    if (!err) {
	/* compile LaTeX and display output */
	err = real_display_gpage();
    }

    if (latex_orig != NULL) {
	strcpy(latex, latex_orig);
	g_free(latex_orig);
    }

    gpage_cleanup();

    return err;
}

void clear_graph_page (void)
{
    int i;

    for (i=0; i<gpage.ngraphs; i++) {
	free(gpage.fnames[i]);
    }
    free(gpage.fnames);

    graph_page_init();
}

int graph_page_get_n_graphs (void)
{
    return gpage.ngraphs;
}

int save_graph_page (const char *fname)
{
    char *latex_orig = NULL;
    int err = 0;

    gpage_filenames_init(fname);

    if (!strncmp(latex, "pdf", 3)) {
	if (gnuplot_has_pdf()) {
	    gpage.output = PDF_OUTPUT;
	} else {
	    latex_orig = g_strdup(latex);
	    strcpy(latex, latex_orig + 3);
	    gpage.output = PS_OUTPUT;
	}
    } else {
	gpage.output = PS_OUTPUT;
    }

    /* write the LaTeX driver file */
    err = make_graphpage_tex();

    if (!err) {
	/* transform individual plot files and compile 
	   using gnuplot */
	err = make_gp_output();
    }

    if (latex_orig != NULL) {
	strcpy(latex, latex_orig);
	g_free(latex_orig);
    }

    return err;
}




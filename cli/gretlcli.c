/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* interactive client program for libgretl - 
   uses the GNU readline library if available */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>

#include "libgretl.h"
#include "var.h"
#include "gretl_restrict.h"
#include "gretl_func.h"
#include "modelspec.h"
#include "libset.h"

#ifdef WIN32
# include <windows.h>
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
#endif 

#ifdef HAVE_READLINE
# include <readline/readline.h>
/* readline functions from complete.c */
extern char *rl_gets (char **line_read, const char *prompt);
extern void initialize_readline (void);
#endif /* HAVE_READLINE */

char prefix[MAXLEN];
char runfile[MAXLEN];
char cmdfile[MAXLEN];
char datfile[MAXLEN];
char outfile[MAXLEN];
char hdrfile[MAXLEN];
char syscmd[MAXLEN];
double **Z;                   /* data set */
MODEL **models;               /* holds ptrs to model structs */
DATAINFO *datainfo;           /* info on data set */
FREQDIST *freq;               /* struct for freq distributions */
CMD cmd;                      /* struct for command characteristics */
PATHS paths;                  /* useful paths */
PRN *cmdprn;
MODELSPEC *modelspec;
MODEL tmpmod;
FILE *dat, *fb;
int i, j, dot, opt, err, errfatal, batch;
int runit, loopstack, looprun;
int data_status, runfile_open;
int echo_off;               /* suppress command echoing */
int plot_count;             /* graphs via gnuplot */
int ignore;                 /* trap for comments */
int order;                  /* VAR lag order */
int lines[1];               /* for gnuplot command */
char *line;                 /* non-Readline command line */
char texfile[MAXLEN];
char response[3];
char linebak[MAXLEN];      /* for storing comments */
char *line_read;

gretl_equation_system *sys;
gretl_restriction_set *rset;

static void exec_line (char *line, LOOPSET **ploop, PRN *prn); 

void usage(void)
{
    logo();
    printf(_("\nYou may supply the name of a data file on the command line.\n"
	     "Options:\n"
	     " -b or --batch     Process a command script and exit.\n"
	     " -r or --run       Run a script then hand control to command line.\n"
	     " -p or --pvalue    Determine p-values interactively.\n"
	     " -h or --help      Print this info and exit.\n"
	     " -v or --version   Print version info and exit.\n"
	     " -e or --english   Force use of English rather than translation.\n"
	     "Example of batch mode usage:\n"
	     " gretlcli -b myfile.inp >myfile.out\n"
	     "Example of run mode usage:\n"
	     " gretlcli -r myfile.inp\n"));
    exit(EXIT_SUCCESS);
}

#ifndef WIN32

int make_userdir (PATHS *ppaths) 
{
    DIR *dir = NULL;
    int err = 0;
    
    if ((dir = opendir(ppaths->userdir)) == NULL) {
	err = mkdir(ppaths->userdir, 0755);
	if (err) {
	    fprintf(stderr, _("Couldn't create user directory %s\n"), 
		    ppaths->userdir);
	} else {
	    fprintf(stderr, _("Created user directory %s\n"), ppaths->userdir);
	}
    } else {
	closedir(dir);
    }

    return err;
}

#endif /* WIN32 */

void gretl_abort (char *line)
{
    fprintf(stderr, _("\ngretlcli: error executing script: halting\n"));
    fprintf(stderr, "> %s\n", line);
    exit(EXIT_FAILURE);
}

void noalloc (const char *str)
{
    fprintf(stderr, _("Couldn't allocate memory for %s\n"), str);
    exit(EXIT_FAILURE);
}

int model_test_start (int test_ci, int model_id, PRN *prn)
{
    int m, err = 0;

    if (model_id != 0) {
	m = modelspec_index_from_model_id(modelspec, model_id);
    } else {
	m = modelspec_last_index(modelspec);
    }

#ifdef MSPEC_DEBUG
    fprintf(stderr, "model_test_start: test_ci=%d, model_id=%d, m=%d\n",
	    test_ci, model_id, m);
#endif

    if (m < 0) { 
	if (model_id == 0) {
	    pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	    err = 1;
	} else {
	    pprintf(prn, _("Can't do this: there is no model %d\n"), model_id);
	    err = 1;
	}
    }    
    else if (!command_ok_for_model(test_ci, 
				   model_ci_from_modelspec(modelspec, m))) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
	err = 1;
    }
    else if (model_sample_issue(NULL, modelspec, m, datainfo)) {
	pputs(prn, _("Can't do: the current data set is different from "
		     "the one on which\nthe reference model was estimated\n"));
	err = 1;
    }

    return err;
}

void file_get_line (void)
{
    clear(line, MAXLINE);
    fgets(line, MAXLINE - 1, fb);

    if (*line == '\0') {
	strcpy(line, "quit");
    } else {
	*linebak = 0;
	strncat(linebak, line, MAXLEN-1);
    }

    if (!strncmp(line, "noecho", 6)) {
	echo_off = 1;
    }
    if (!echo_off && cmd.ci == RUN && batch && *line == '(') {
	printf("%s", line);
	*linebak = 0;
    }
}

void fn_get_line (void)
{
    clear(line, MAXLINE);
    
    gretl_function_get_line(line, MAXLINE, &Z, datainfo);

    if (*line == '\0') {
	return;
    }

    *linebak = 0;
    strncat(linebak, line, MAXLEN-1);

    if (!strncmp(line, "noecho", 6)) {
	echo_off = 1;
    }
    if (!echo_off && cmd.ci == RUN && batch && *line == '(') {
	printf("%s", line);
	*linebak = 0;
    }
}

unsigned char gp_flags (int batch, gretlopt opt)
{
    unsigned char flags = 0;

    if (batch) flags |= GP_BATCH;

    if (opt & OPT_M) flags |= GP_IMPULSES;
    else if (opt & OPT_Z) flags |= GP_DUMMY;
    else if (opt & OPT_S) flags |= GP_OLS_OMIT;

    return flags;
}

#ifdef ENABLE_NLS

void nls_init (void)
{
# ifdef WIN32
    char gretldir[MAXLEN], LOCALEDIR[MAXLEN];

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretl", "gretldir", gretldir)) {
        return;
    }
    sprintf(LOCALEDIR, "%s\\locale", gretldir);
# endif /* WIN32 */

    setlocale (LC_ALL, "");
    bindtextdomain (PACKAGE, LOCALEDIR);
    textdomain (PACKAGE); 
    iso_gettext("@CLI_INIT");

    putenv("LC_NUMERIC=");
    setlocale(LC_NUMERIC, "");
    reset_local_decpoint();
}

static void force_english (void)
{
    setlocale (LC_ALL, "C");

# ifdef WIN32
    SetEnvironmentVariable("LC_ALL", "C");
    putenv("LC_ALL=C");
# endif
}

#endif /* ENABLE_NLS */

static int clear_data (void)
{
    int err = 0;

    *paths.datfile = 0;

    err = restore_full_sample(&Z, &datainfo, OPT_C); 

    if (Z != NULL) {
	free_Z(Z, datainfo); 
    }
    clear_datainfo(datainfo, CLEAR_FULL);

    Z = NULL;
    data_status = 0;

    clear_model(models[0]);
    clear_model(models[1]);

    free_modelspec(modelspec);
    modelspec = NULL;

    reset_model_count();

    return err;
}

static void get_an_input_line (void)
{
    if (gretl_executing_function()) {
	/* reading from compiled function */
	fn_get_line();
    } else if (runit || batch) {
	/* reading from script file */
	file_get_line();
    } else {
	/* normal interactive use */
#ifdef HAVE_READLINE
	rl_gets(&line_read, (loopstack || gretl_compiling_function())? 
		"> " : "? ");
	if (line_read == NULL) {
	    strcpy(line, "quit");
	} else {
	    strcpy(line, line_read);
	}
#else
	printf("%s", (loopstack || gretl_compiling_function())? 
	       "> " : "? ");
	fflush(stdout);
	file_get_line(); /* "file" = stdin here */
#endif
    }
}

static int maybe_get_input_line_continuation (char *tmp)
{
    int err = 0;

    if (!strncmp(line, "quit", 4)) {
	return 0;
    }

    /* allow for backslash continuation of lines */
    while (top_n_tail(line)) {
	tmp[0] = '\0';

	if (gretl_executing_function()) {
	    gretl_function_get_line(tmp, MAXLEN - 1, &Z, datainfo);
	} else if (batch || runit) {
	    fgets(tmp, MAXLEN - 1, fb);
	} else {
#ifdef HAVE_READLINE
	    rl_gets(&line_read, "> ");
	    strcpy(tmp, line_read);
#else
	    fgets(tmp, MAXLEN - 1, stdin); 
#endif /* HAVE_READLINE */
	}

	if (*tmp != '\0') {
	    if (strlen(line) + strlen(tmp) > MAXLEN - 1) {
		err = 1;
		break;
	    } else {
		strcat(line, tmp);
		compress_spaces(line);
	    }
	}
    }

    return err;
}

int main (int argc, char *argv[])
{
    int cli_get_data = 0;
    int cmd_overflow = 0;
    char filearg[MAXLEN];
    char tmp[MAXLINE];
    LOOPSET *loop = NULL;
    PRN *prn;

#ifdef WIN32
    strcpy(tmp, argv[0]);
#endif

#ifdef ENABLE_NLS
    nls_init();
#endif

    datainfo = datainfo_new();
    if (datainfo == NULL) {
	noalloc(_("data information"));
    }
    
    if (argc > 1) {
	int english = 0;
	int opt = parseopt((const char **) argv, argc, filearg, &english);

	switch (opt) {
	case OPT_BATCH:
	    batch = 1;
	    if (*filearg == '\0') usage();
	    strcpy(runfile, filearg);
	    cli_get_data = 1;
	    break;
	case OPT_HELP:
	case OPT_DBOPEN:
	case OPT_WEBDB:
	    usage();
	    break;
	case OPT_PVALS:
	    interact_pvalue();
	    exit(EXIT_SUCCESS);
	    break;
	case OPT_VERSION:
	    logo();
	    exit(EXIT_SUCCESS);
	    break;
	case OPT_RUNIT:
	    runit = 1;
	    if (*filearg == '\0') usage();
	    strcpy(runfile, filearg); 
	    cli_get_data = 1;
	    break;
	default:
	    break;
	}

#ifdef ENABLE_NLS
	if (english) {
	    force_english();
	    if (argc == 2) {
		cli_get_data = 1;
	    }
	}
#endif
    } else {
	cli_get_data = 1;
    }

#ifdef WIN32
    if (!batch) {
	bind_textdomain_codeset (PACKAGE, "CP850");
    }
#endif

    logo();     /* print version info */
    session_time(stdout);
    fb = stdin; /* may be reset later with "run" command */

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    line = malloc(MAXLINE);
    if (line == NULL) {
	noalloc(_("command line"));
    } 

    set_paths(&paths, 1, 0); /* 1 = defaults, 0 = not gui */
#ifdef WIN32
    cli_read_registry(tmp, &paths);
    set_paths(&paths, 0, 0); /* not defaults; use registry info */
#else
    make_userdir(&paths);
#endif /* WIN32 */

    if (!batch) {
	strcpy(cmdfile, paths.userdir);
	strcat(cmdfile, "session.inp");
	cmdprn = gretl_print_new(GRETL_PRINT_FILE, cmdfile);
	if (cmdprn == NULL) {
	    printf(_("Can't open file to save commands\n"));
	    return EXIT_FAILURE;
	}
    }

    if (!cli_get_data) {
	char given_file[MAXLEN];

	strcpy(given_file, filearg);
	strcpy(paths.datfile, filearg);

	err = detect_filetype(paths.datfile, &paths, prn);

	if (err == GRETL_UNRECOGNIZED || err == GRETL_NATIVE_DB ||
	    err == GRETL_RATS_DB) { 
	    exit(EXIT_FAILURE);
	}

	if (err == GRETL_NATIVE_DATA) {
	    err = gretl_get_data(&Z, &datainfo, paths.datfile, &paths, 
				 DATA_NONE, prn);
	} 
	else if (err == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, &datainfo, paths.datfile, &paths, 
			      DATA_NONE, prn, 0);
	} 
	else if (err == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, paths.datfile, &paths, prn);
	} 
	else if (err == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, paths.datfile, prn);
	} 
	else if (err == GRETL_SCRIPT) { /* maybe it's a script file? */
	    runit = 1;
	    strcpy(runfile, paths.datfile); 
	    clear(paths.datfile, MAXLEN);
	    cli_get_data = 1;
	}

	if (!cli_get_data) {
	    if (err) {
		errmsg(err, prn);
		if (err == E_FOPEN) show_paths(&paths);
		return EXIT_FAILURE;
	    }
	    data_status = 1;
	    if (!batch) { 
		pprintf(cmdprn, "open %s\n", given_file);
	    }
	}
    }

    /* allocate memory for models */
    models = malloc(2 * sizeof *models);
    if (models == NULL) noalloc("models"); 

    models[0] = gretl_model_new();
    models[1] = gretl_model_new();

    if (models[0] == NULL || models[1] == NULL) {
	noalloc("models"); 
    }

    libgretl_init(&cmd);

    /* print list of variables */
    if (data_status) {
	varlist(datainfo, prn);
    }

    /* check for help file */
    if (!batch) {
	dat = fopen(paths.helpfile, "r");
	if (dat != NULL) { 
	    printf(_("\n\"help\" gives a list of commands\n"));
	    fclose(dat);
	} else {
	    printf(_("help file %s is not accessible\n"), 
		   paths.helpfile);
	    show_paths(&paths);
	}
    } 

    if (!batch && !runit && !data_status) {
	fprintf(stderr, _("Type \"open filename\" to open a data set\n"));
    }

#ifdef HAVE_READLINE
    if (!batch) initialize_readline();
#endif

    if (batch || runit) {
	sprintf(line, "run %s\n", runfile);
	exec_line(line, &loop, prn);
    }

    /* should we stop immediately on error, in batch mode? */
    if (getenv("GRETL_KEEP_GOING") == NULL) {
	errfatal = batch;  /* exit on error in batch mode */
    }

    /* main command loop */
    while (strcmp(cmd.cmd, "quit")) {
	char linecopy[MAXLEN];

	if (err && batch && errfatal) {
	    gretl_abort(linecopy);
	}

	if (looprun) { 
	    if (loop_exec(loop, line, &Z, &datainfo,
			  models, &paths, 
			  &echo_off, prn)) {
		return 1;
	    }
	    gretl_loop_destroy(loop);
	    loop = NULL;
	    looprun = errfatal = 0;
	} else {
	    get_an_input_line();
	}

	cmd_overflow = maybe_get_input_line_continuation(tmp);
	if (cmd_overflow) {
	    fprintf(stderr, _("Maximum length of command line "
			      "(%d bytes) exceeded\n"), MAXLEN);
	    break;
	} else {
	    strcpy(linecopy, line);
	    exec_line(line, &loop, prn);
	}
    } /* end of get commands loop */

    /* leak check -- try explicitly freeing all memory allocated */

    free_Z(Z, datainfo);

    free_model(models[0]);
    free_model(models[1]);
    free(models);

    if (data_status) {
	free_datainfo(datainfo);
    }

    if (runfile_open && fb != NULL) {
	fclose(fb);
    }

    free(line);

    if (modelspec != NULL) {
	free_modelspec(modelspec);
    }

    if (!batch) {
	remove(paths.plotfile);
    }

    gretl_print_destroy(prn);

    libgretl_cleanup(&cmd);

    return 0;
}

static void printf_strip (char *s, int loopstack)
{
    int i, n;

    while (isspace((unsigned char) *s)) s++;

    n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (isspace(s[i]) || s[i] == '\r') s[i] = '\0';
	else break;
    }

    if (loopstack) {
	printf("> %s\n", s);
    } else {
	printf("%s\n", s);
    }
}

static int handle_user_defined_function (char *line, int *fncall)
{
    int ufunc = gretl_is_user_function(line);
    int err = 0;

    /* allow for nested function calls */
    if (ufunc && gretl_compiling_function()) {
	return 0;
    }

    /* an actual function call */
    else if (ufunc) {
	err = gretl_function_start_exec(line);
	*fncall = 1;
    } 

    return err;
}

static void exec_line (char *line, LOOPSET **ploop, PRN *prn) 
{
    LOOPSET *loop = *ploop;
    int chk, nulldata_n, renumber;
    int dbdata = 0, do_arch = 0, do_nls = 0;
    char s1[12], s2[12];
    int fncall = 0;
    double rho;

    if (*line == '\0') return;

    /* catch any user-defined functions */
    err = handle_user_defined_function(line, &fncall);
    if (err) {
	errmsg(err, prn);
	return;
    } else if (fncall) {
	return;
    }

    if (gretl_compiling_function()) {
	err = gretl_function_append_line(line);
	if (err) errmsg(err, prn);
	return;
    }  
    
    /* are we ready for this? */
    if (!data_status && !ignore && 
	!ready_for_command(line)) {
	fprintf(stderr, _("You must open a data file first\n"));
	err = 1;
	return;
    }

    err = catchflags(line, &cmd.opt);
    if (err) {
	errmsg(err, prn);
	return;
    }

    compress_spaces(line);

    /* if we're stacking commands for a loop, parse "lightly" */
    if (loopstack) {
	get_cmd_ci(line, &cmd);
    } else {
	getcmd(line, datainfo, &cmd, &ignore, &Z, 
	       (runit)? NULL : cmdprn);
    }

    /* tell libgretl if we should pause between screens */
    gretl_set_text_pause(!batch);

    /* if in batch mode, echo comments from input */
    if (batch && cmd.ci == CMD_COMMENT && !echo_off) {
	printf_strip(linebak, loopstack);
    }

    if ((err = cmd.errcode)) {
	errmsg(err, prn);
	return;
    }

    if (cmd.ci < 0) return; /* there's nothing there */ 

    if (sys != NULL && cmd.ci != END && cmd.ci != EQUATION &&
	cmd.ci != SYSTEM) {
	printf(_("Command '%s' ignored; not valid within equation system\n"), 
	       line);
	gretl_equation_system_destroy(sys);
	sys = NULL;
	return;
    }

    if (cmd.ci == LOOP && !batch && !runit) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }
   
    if (loopstack || cmd.ci == LOOP) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd.ci, loop)) {
	    printf(_("Command '%s' ignored; not available in loop mode\n"), line);
	} else {
	    if (!echo_off) {
		echo_cmd(&cmd, datainfo, line, (batch || runit)? 1 : 0, 
			 0, 1, cmdprn);
	    }
	    loop = add_to_loop(line, cmd.ci, cmd.opt,
			       datainfo, &Z, loop, 
			       &loopstack, &looprun);
	    if (loop == NULL) {
		print_gretl_errmsg(prn);
		err = 1;
	    } 
	    *ploop = loop;
	}
	return;
    }

    if (!echo_off) {
	echo_cmd(&cmd, datainfo, line, (batch || runit)? 1 : 0, 0, 0,
		 cmdprn);
    }

    check_for_loop_only_options(cmd.ci, cmd.opt, prn);

    switch (cmd.ci) {

    case ADF: case COINT: case COINT2:
    case CORR:
    case CRITERIA: case CRITICAL: case DATA:
    case DIFF: case LDIFF: case LAGS: case LOGS:
    case MULTIPLY:
    case GRAPH: case PLOT: case LABEL:
    case INFO: case LABELS: case VARLIST:
    case PRINT: 
    case SUMMARY:
    case MEANTEST: case VARTEST: case STORE:
    case RUNS: case SPEARMAN: case OUTFILE: case PCA:
	err = simple_commands(&cmd, line, &Z, datainfo, &paths,
			      prn);
	if (err) errmsg(err, prn);
	break;

    case ADD:
    case OMIT:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
    plain_add_omit:
	clear_model(models[1]);
	if (cmd.ci == ADD || cmd.ci == ADDTO) {
	    err = auxreg(cmd.list, models[0], models[1], 
			 &Z, datainfo, AUX_ADD, NULL, cmd.opt, prn);
	} else {
	    err = omit_test(cmd.list, models[0], models[1],
			    &Z, datainfo, cmd.opt, prn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	} else {
	    /* for command-line use, we keep a "stack" of 
	       two models, and recycle the places */
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(&models[0], &models[1]);
	    } 
	    clear_model(models[1]);
	}
	break;	

    case ADDTO:
    case OMITFROM:
	i = atoi(cmd.param);
	if ((err = model_test_start(cmd.ci, i, prn))) break;
	if (i == (models[0])->ID) goto plain_add_omit;
	err = re_estimate(modelspec_get_command_by_id(modelspec, i), 
			  &tmpmod, &Z, datainfo);
	if (err) {
	    pprintf(prn, _("Failed to reconstruct model %d\n"), i);
	    break;
	} 
	clear_model(models[1]);
	tmpmod.ID = i;
	if (cmd.ci == ADDTO) {
	    err = auxreg(cmd.list, &tmpmod, models[1], 
			 &Z, datainfo, AUX_ADD, NULL, cmd.opt, prn);
	} else {
	    err = omit_test(cmd.list, &tmpmod, models[1],
			    &Z, datainfo, cmd.opt, prn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	    break;
	} else {
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(&models[0], &models[1]);
	    }
	    clear_model(models[1]);
	}
	clear_model(&tmpmod);
	break;

    case AR:
	clear_model(models[0]);
	*models[0] = ar_func(cmd.list, atoi(cmd.param), &Z, 
			     datainfo, cmd.opt, prn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	    break;
	}
	break;

    case ARCH:
	order = atoi(cmd.param);
	clear_model(models[1]);
	*models[1] = arch(order, cmd.list, &Z, datainfo, 
			  NULL, cmd.opt, prn);
	if ((err = (models[1])->errcode)) 
	    errmsg(err, prn);
	if ((models[1])->ci == ARCH) {
	    do_arch = 1;
	    swap_models(&models[0], &models[1]); 
	}
	clear_model(models[1]);
	break;

    case ARMA:
	clear_model(models[0]);
#ifdef HAVE_X12A
	if (cmd.opt & OPT_X) {
	    *models[0] = arma_x12(cmd.list, (const double **) Z, datainfo,
				  ((cmd.opt & OPT_V) ? prn : NULL), &paths); 
	} else {
	    *models[0] = arma(cmd.list, (const double **) Z, datainfo, 
			      (cmd.opt & OPT_V)? prn : NULL);
	}
#else
	*models[0] = arma(cmd.list, (const double **) Z, datainfo, 
			  (cmd.opt & OPT_V)? prn : NULL);
#endif
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	} else {	
	    printmodel(models[0], datainfo, cmd.opt, prn);
	}	
	break;

    case GARCH:
	clear_model(models[0]);
	*models[0] = garch(cmd.list, &Z, datainfo, cmd.opt, prn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	} else {	
	    printmodel(models[0], datainfo, cmd.opt, prn);
	}	
	break;

    case CHOW:
        if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = chow_test(line, models[0], &Z, datainfo, prn, NULL);
	if (err) errmsg(err, prn);
	break;

    case COEFFSUM:
        if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = sum_test(cmd.list, models[0], &Z, datainfo, prn);
	if (err) errmsg(err, prn);
	break;

    case CUSUM:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = cusum_test(models[0], &Z, datainfo, prn, &paths, NULL);
	if (err) errmsg(err, prn);
	break;

    case RESET:
        if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = reset_test(models[0], &Z, datainfo, prn, NULL);
	if (err) errmsg(err, prn);
	break;
	
    case CORC:
    case HILU:
    case PWE:
	err = hilu_corc(&rho, cmd.list, &Z, datainfo, 
			NULL, 1, cmd.ci, prn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	clear_model(models[0]);
	*models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, rho);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, prn); 
	break;

    case LAD:
	clear_model(models[0]);
	*models[0] = lad(cmd.list, &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, prn);
	break;

    case CORRGM:
	order = atoi(cmd.param);
	err = corrgram(cmd.list[1], order, &Z, datainfo, &paths,
		       batch, prn);
	if (err) pputs(prn, _("Failed to generate correlogram\n"));
	break;

    case DELEET:
	if (complex_subsampled()) {
	    pputs(prn, _("Can't delete a variable when in sub-sample"
			 " mode\n"));
	    break;
	}	
	if (cmd.list[0]) {
	    err = dataset_drop_listed_vars(cmd.list, &Z, datainfo, 
					   &renumber);
	} else {
	    err = dataset_drop_vars(1, &Z, datainfo);
	    renumber = 0;
	}
	if (err) {
	    pputs(prn, _("Failed to shrink the data set"));
	} else {
	    if (renumber) {
		pputs(prn, _("Take note: variables have been renumbered"));
		pputc(prn, '\n');
	    }
	    varlist(datainfo, prn);
	}
	break;

    case RENAME:
	err = rename_var_by_id(cmd.str, cmd.param, datainfo);
	if (err) {
	    errmsg(err, prn);
	} else {
	    varlist(datainfo, prn);
	}
	break;

    case END:
	if (!strcmp(cmd.param, "system")) {
	    err = gretl_equation_system_finalize(sys, &Z, datainfo, prn);
	    if (err) errmsg(err, prn);
	    sys = NULL;
	} 
	else if (!strcmp(cmd.param, "nls")) {
	    clear_model(models[0]);
	    *models[0] = nls(&Z, datainfo, prn);
	    if ((err = (models[0])->errcode)) {
		errmsg(err, prn);
		break;
	    }
	    do_nls = 1;
	    printmodel(models[0], datainfo, cmd.opt, prn);
	}
	else if (!strcmp(cmd.param, "restrict")) {
	    err = gretl_restriction_set_finalize(rset, prn);
	    if (err) errmsg(err, prn);
	    rset = NULL;
	}  
	else {
	    err = 1;
	}
	break;

    case ENDLOOP:
	pputs(prn, _("You can't end a loop here, "
		     "you haven't started one\n"));
	err = 1;
	break;

    case EQUATION:
	err = gretl_equation_system_append(sys, cmd.list);
	if (err) {
	    gretl_equation_system_destroy(sys);
	    sys = NULL;
	    errmsg(err, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	strcpy(texfile, cmd.param);
	if ((err = model_test_start(cmd.ci, 0, prn))) 
	    break;
	if (cmd.ci == EQNPRINT) {
	    err = eqnprint(models[0], datainfo, &paths, 
			   texfile, cmd.opt);
	} else {
	    err = tabprint(models[0], datainfo, &paths, 
			   texfile, cmd.opt);
	}
	if (err) {
	    pputs(prn, _("Couldn't open tex file for writing\n"));
	} else {
	    pprintf(prn, _("Model printed to %s\n"), texfile);
	}
	break;

    case FCAST:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = fcast(line, models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    pputs(prn, _("Error retrieving fitted values\n"));
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	varlist(datainfo, prn);
	break;

    case FCASTERR:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = fcast_with_errs(line, models[0], &Z, datainfo, prn,
			      &paths, cmd.opt); 
	if (err) errmsg(err, prn);
	break;

    case FIT:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = fcast("fcast autofit", models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    pputs(prn, _("Error retrieving fitted values\n"));
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	pputs(prn, _("Retrieved fitted values as \"autofit\"\n"));
	varlist(datainfo, prn);
	if (dataset_is_time_series(datainfo)) {
	    plotvar(&Z, datainfo, "time");
	    cmd.list = realloc(cmd.list, 4 * sizeof(int));
	    cmd.list[0] = 3; 
	    if ((models[0])->ci == ARMA) {
		cmd.list[1] = (models[0])->list[4];
	    } else {
		cmd.list[1] = (models[0])->list[1];
	    }
	    cmd.list[2] = varindex(datainfo, "autofit");
	    cmd.list[3] = varindex(datainfo, "time");
	    lines[0] = 1;
	    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
			  &paths, &plot_count, gp_flags(batch, 0));
	    if (err) {
		pputs(prn, _("gnuplot command failed\n"));
	    }
	}
	break;
		
    case FREQ:
	freq = freqdist(&Z, datainfo, cmd.list[1], 1);
	if (freq == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if ((err = get_gretl_errno())) 
	    errmsg(err, prn);
	else {
	    printfreq(freq, prn); 
	    if (!batch) {
		if (plot_freq(freq, &paths, NORMAL))
		    pputs(prn, _("gnuplot command failed\n"));
	    }
	    free_freq(freq);
	}
	break;

    case FUNC:
	err = gretl_start_compiling_function(line);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case GENR:
	err = generate(&Z, datainfo, line, models[0]);
	if (err) { 
	    errmsg(err, prn);
	} else {
	    print_gretl_msg(prn);
	}	
	break;

    case GNUPLOT:
	if ((cmd.opt & OPT_Z) && 
	    (cmd.list[0] != 3 || 
	     !isdummy(Z[cmd.list[3]], datainfo->t1, datainfo->t2))) { 
	    pputs(prn, _("You must supply three variables, the last of "
			 "which is a dummy variable\n(with values 1 or 0)\n"));
	    break;
	}
	if ((cmd.opt & OPT_M) || (cmd.opt & OPT_Z) || (cmd.opt & OPT_S)) { 
	    err = gnuplot(cmd.list, NULL, NULL, &Z, datainfo,
			  &paths, &plot_count, gp_flags(batch, cmd.opt));
	} else {
	    lines[0] = (cmd.opt != 0);
	    err = gnuplot(cmd.list, lines, cmd.param, 
			  &Z, datainfo, &paths, &plot_count, 
			  gp_flags(batch, 0));
	}
	if (err < 0) {
	    pputs(prn, _("gnuplot command failed\n"));
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), paths.plotfile);
	}
	break;

    case SCATTERS:
	err = multi_scatters(cmd.list, atoi(cmd.param), &Z, 
			     datainfo, &paths, &plot_count, 
			     gp_flags(batch, cmd.opt));
	if (err) {
	    pputs(prn, _("scatters command failed\n"));
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), paths.plotfile);
	}
	break;

    case HAUSMAN:
	err = model_test_start(cmd.ci, 0, prn);
	if (!err) {
	    err = hausman_test(models[0], &Z, datainfo, prn);
	}
	break;

    case HCCM:
    case HSK:
	clear_model(models[0]);
	if (cmd.ci == HCCM) {
	    *models[0] = hccm_func(cmd.list, &Z, datainfo);
	} else {
	    *models[0] = hsk_func(cmd.list, &Z, datainfo);
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, prn);
	break;

    case HELP:
	help(cmd.param, paths.helpfile, prn);
	break;

    case IMPORT:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pputs(prn, _("import command is malformed\n"));
	    break;
	}
	if (data_status) {
	    clear_data();
	}
	if (cmd.opt) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else {
	    err = import_csv(&Z, &datainfo, datfile, &paths, prn);
	}
	if (!err) { 
	    data_status = 1;
	    print_smpl(datainfo, 0, prn);
	    varlist(datainfo, prn);
	    pputs(prn, _("You should now use the \"print\" command "
			 "to verify the data\n"));
	    pputs(prn, _("If they are OK, use the \"store\" command "
			 "to save them in gretl format\n"));
	}
	break;

    case OPEN:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pputs(prn, _("'open' command is malformed\n"));
	    break;
	}

	chk = detect_filetype(datfile, &paths, prn);
	dbdata = (chk == GRETL_NATIVE_DB || chk == GRETL_RATS_DB);

	if (data_status && !batch && !dbdata &&
	    strcmp(datfile, paths.datfile)) {
	    fprintf(stderr, _("Opening a new data file closes the "
			      "present one.  Proceed? (y/n) "));
	    fgets(response, 2, stdin);
	    if (*response != 'y' && *response != 'Y') {
		fprintf(stderr, 
			_("OK, staying with current data set\n"));
		break;
	    }
	}

	if (data_status && !dbdata) {
	    clear_data();
	}

	if (chk == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, datfile, &paths, prn);
	} else if (chk == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else if (chk == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, &datainfo, datfile, &paths, 
			      data_status, prn, 0);
	} else if (dbdata) {
	    err = set_db_name(datfile, chk, &paths, prn);
	} else {
	    err = gretl_get_data(&Z, &datainfo, datfile, &paths, 
				 data_status, prn);
	}

	if (err) {
	    errmsg(err, prn);
	    break;
	}
	strncpy(paths.datfile, datfile, MAXLEN-1);
	data_status = 1;
	if (datainfo->v > 0 && !dbdata) {
	    varlist(datainfo, prn);
	}
	paths.currdir[0] = '\0';
	break;

    case LEVERAGE:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;	
	err = leverage_test(models[0], &Z, datainfo, prn, NULL, cmd.opt);
	if (err > 1) {
	    errmsg(err, prn);
	} else if (cmd.opt) {
	    varlist(datainfo, prn);
	}
	break;

    case VIF:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	err = vif_test(models[0], &Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case LMTEST:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	/* non-linearity (squares) */
	if ((cmd.opt & OPT_S) || (cmd.opt & OPT_O) || !cmd.opt) {
	    clear_model(models[1]);
	    err = auxreg(NULL, models[0], models[1], 
			 &Z, datainfo, AUX_SQ, NULL, OPT_NONE, prn);
	    clear_model(models[1]);
	    if (err) errmsg(err, prn);
	    if (cmd.opt == OPT_S) break;
	    if (!err && !batch && page_pause()) break; 
	}
	/* non-linearity (logs) */
	if ((cmd.opt & OPT_L) || (cmd.opt & OPT_O) || !cmd.opt) {
	    err = auxreg(NULL, models[0], models[1], 
			 &Z, datainfo, AUX_LOG, NULL, OPT_NONE, prn);
	    clear_model(models[1]); 
	    if (err) errmsg(err, prn);
	    if (cmd.opt == OPT_L) break;
	    if (!err && !batch && page_pause()) break;
	}
	/* autocorrelation */
	if ((cmd.opt & OPT_M) || (cmd.opt & OPT_O)) {
	    int order = atoi(cmd.param);

	    err = autocorr_test(models[0], order, &Z, datainfo, prn, NULL);
	    if (err) errmsg(err, prn);
	}
	/* heteroskedasticity */
	if ((cmd.opt & OPT_W) || !cmd.opt) {
	    err = whites_test(models[0], &Z, datainfo, prn, NULL);
	    if (err) errmsg(err, prn);
	    /* need to take more action in case of err? */
	}
	break;

    case LOGISTIC:
    case LOGIT:
    case PROBIT:
    case TOBIT:
	clear_model(models[0]);
	if (cmd.ci == LOGIT || cmd.ci == PROBIT) {
	    *models[0] = logit_probit(cmd.list, &Z, datainfo, cmd.ci);
	} else if (cmd.ci == LOGISTIC) {
	    *models[0] = logistic_model(cmd.list, &Z, datainfo, cmd.param);
	} else {
	    *models[0] = tobit_model(cmd.list, &Z, datainfo,
				     (cmd.opt & OPT_V)? prn : NULL);
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, prn);
	break;

    case NLS:
	err = nls_parse_line(line, (const double **) Z, datainfo);
	if (err) errmsg(err, prn);
	break;

    case NULLDATA:
	nulldata_n = atoi(cmd.param);
	if (nulldata_n < 2) {
	    pputs(prn, _("Data series length count missing or invalid\n"));
	    err = 1;
	    break;
	}
	if (nulldata_n > 1000000) {
	    pputs(prn, _("Data series too long\n"));
	    err = 1;
	    break;
	}
	if (data_status) {
	    clear_data();
	}	
	err = open_nulldata(&Z, datainfo, data_status, 
			    nulldata_n, prn);
	if (err) { 
	    pputs(prn, _("Failed to create empty data set\n"));
	} else {
	    data_status = 1;
	}
	break;

    case OLS:
    case WLS:
    case POOLED:
	clear_model(models[0]);
	*models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, 0.0);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    clear_model(models[0]);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, prn);
	break;

#ifdef ENABLE_GMP
    case MPOLS:
	err = mp_ols(cmd.list, cmd.param, &Z, datainfo, prn);
	if (err) {
	    pputs(prn, _("mpols command failed\n"));
	    errmsg(err, prn);
	}
	break;
#endif

    case PANEL:	
	err = set_panel_structure(cmd.opt, datainfo, prn);
	break;

    case PERGM:
	err = periodogram(cmd.list[1], &Z, datainfo, &paths,
			  batch, cmd.opt, prn);
	if (err) pputs(prn, _("Failed to generate periodogram\n"));
	break;

    case PRINTF:
	err = do_printf(line, &Z, datainfo, models[0], prn);
	break;

    case PVALUE:
	if (batch || runit || (sscanf(line, "%s %s", s1, s2) == 2))
	    batch_pvalue(line, (const double **) Z, datainfo, prn);
	else interact_pvalue();
	break;

    case QUIT:
	if (batch) {
	    pputs(prn, _("Done\n"));
	    break;
	}
	if (runit) {
	    runit = 0;
	    if (fb != NULL) {
		fclose(fb);
	    }
	    fb = stdin;
	    runfile_open = 0;
	    strcpy(cmd.cmd, "endrun"); /* overwrite "quit" */
	    break;
	}
	printf(_("commands saved as %s\n"), cmdfile);
	gretl_print_destroy(cmdprn);

	if (cmd.param[0] == 'x') break;

	printf(_("type a filename to store output (enter to quit): "));
	*outfile = '\0';
	fgets(outfile, MAXLEN-1, stdin); 
	top_n_tail(outfile);

	if (*outfile != '\n' && *outfile != '\r' && strcmp(outfile, "q")) {
	    printf(_("writing session output to %s%s\n"), 
		   paths.userdir, outfile);
#ifdef WIN32
	    sprintf(syscmd, "\"%s\\gretlcli\" -b \"%s\" > \"%s%s\"", 
		    paths.gretldir, cmdfile, paths.userdir, outfile);
	    /* WinExec(syscmd, SW_SHOWMINIMIZED); */
	    system(syscmd);
#else
	    sprintf(syscmd, "gretlcli -b \"%s\" > \"%s%s\"", 
		    cmdfile, paths.userdir, outfile);
	    gretl_spawn(syscmd);
#endif
	    printf("%s\n", syscmd);
	} 
	break;

    case RHODIFF:
	if (!cmd.list[0]) {
	    pputs(prn, _("This command requires a list of variables\n"));
	    break;
	}
	err = rhodiff(cmd.param, cmd.list, &Z, datainfo);
	if (err) errmsg(err, prn);
	else varlist(datainfo, prn);
	break;

    case RUN:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pputs(prn, _("Command is malformed\n"));
	    break;
	}
	if ((fb = fopen(runfile, "r")) == NULL) {
	    fprintf(stderr, _("Couldn't open script \"%s\"\n"), runfile);
	    if (runit) {
		fb = stdin;
		runit = 0;
	    } else {
		if (batch) exit(EXIT_FAILURE);
	    }
	} else {
	    fprintf(stderr, _("%s opened OK\n"), runfile);
	    pprintf(cmdprn, "run \"%s\"\n", runfile);
	    runfile_open = 1;
	    if (!batch) runit = 1;
	}
	break;

    case SET:
	err = parse_set_line(line, &echo_off, prn);
	if (err) errmsg(err, prn);
	break;

    case SETOBS:
	err = set_obs(line, datainfo, cmd.opt);
	if (err) errmsg(err, prn);
	else {
	    if (datainfo->n > 0) {
		print_smpl(datainfo, 0, prn);
	    } else {
		pprintf(prn, _("setting data frequency = %d\n"), datainfo->pd);
	    }
	}
	break;

    case SETMISS:
	set_miss(cmd.list, cmd.param, Z, datainfo, prn);
	break;

    case SHELL:
#ifdef WIN32
	WinExec(line + 1, SW_SHOWNORMAL);
#else		
	shell(line + 1);
#endif
	break;

    case SIM:
	err = simulate(line, &Z, datainfo);
	if (err) { 
	    errmsg(err, prn);
	} else {
	    print_gretl_msg(prn);
	}
	break;

    case SMPL:
	if (cmd.opt) {
	    err = restore_full_sample(&Z, &datainfo, cmd.opt);
	    if (err) {
		errmsg(err, prn);
		break;
	    } else {
		err = restrict_sample(line, &Z, &datainfo, 
				      cmd.list, cmd.opt);
	    }
	} else if (!strcmp(line, "smpl full") ||
		   !strcmp(line, "smpl --full")) {
	    err = restore_full_sample(&Z, &datainfo, OPT_C);
	} else { 
	    err = set_sample(line, datainfo);
	}

	if (err) {
	    errmsg(err, prn);
	} else {
	    print_smpl(datainfo, get_full_length_n(), prn);
	}
	break;

    case SQUARE:
	if (cmd.opt) {
	    chk = xpxgenr(cmd.list, &Z, datainfo, 1, 1);
	} else {
	    chk = xpxgenr(cmd.list, &Z, datainfo, 0, 1);
	}
	if (chk < 0) {
	    pputs(prn, _("Failed to generate squares\n"));
	    err = 1;
	} else {
	    pputs(prn, _("Squares generated OK\n"));
	    varlist(datainfo, prn);
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model */
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	if (rset == NULL) {
	    rset = restriction_set_start(line, models[0], datainfo);
	    if (rset == NULL) {
		err = 1;
		errmsg(err, prn);
	    }
	} else {
	    err = restriction_set_parse_line(rset, line);
	    if (err) {
		errmsg(err, prn);
		rset = NULL;
	    }	
	}
	break;

    case SYSTEM:
	/* system of equations */
	if (sys == NULL) {
	    sys = system_start(line);
	    if (sys == NULL) {
		err = 1;
		errmsg(err, prn);
	    }
	} else {
	    err = system_parse_line(sys, line, datainfo);
	    if (err) {
		errmsg(err, prn);
		sys = NULL;
	    }	
	}
	break;

    case TESTUHAT:
	if ((err = model_test_start(cmd.ci, 0, prn))) break;
	if (genr_fit_resid(models[0], &Z, datainfo, GENR_RESID, 1)) {
	    pputs(prn, _("Out of memory attempting to add variable\n"));
	    err = 1;
	    break;
	}
	freq = freqdist(&Z, datainfo, datainfo->v - 1, (models[0])->ncoeff);	
	dataset_drop_vars(1, &Z, datainfo);
	if (freq == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if ((err = get_gretl_errno())) { 
	    errmsg(err, prn);
	} else {
	    printfreq(freq, prn); 
	    free_freq(freq);
	}
	break;

    case TSLS:
	clear_model(models[0]);
	*models[0] = tsls_func(cmd.list, atoi(cmd.param), 
			       &Z, datainfo, cmd.opt);
	if ((err = (models[0])->errcode)) {
	    errmsg((models[0])->errcode, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, prn);
	break;

    case VAR:
	order = atoi(cmd.param);
	err = simple_var(order, cmd.list, &Z, datainfo, !batch, 
			 cmd.opt, prn);
	break;

    case VARDUP:
	err = 1;
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		       "in gretlcli\n"), cmd.cmd);
	err = 1;
	break;
    }

    if (err && gretl_executing_function()) {
	gretl_function_error();
    }

    if (!err && (is_model_cmd(cmd.cmd) || do_nls || do_arch)
	&& !is_quiet_model_test(cmd.ci, cmd.opt)) { 

	attach_subsample_to_model(models[0], datainfo);

#ifdef MSPEC_DEBUG
	fprintf(stderr, "\ngretlcli: saving spec: model.ID = %d, model_count = %d\n",
		(models[0])->ID, get_model_count());
#endif

	err = modelspec_save(models[0], &modelspec);
    }

}







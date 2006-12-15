/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
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

/* interact.h for gretl */

#ifndef INTERACT_H
#define INTERACT_H

#define MAXSAVENAME 32
#define CMD_NULL    -1
#define CMD_COMMENT -2

typedef struct CMD_ CMD;
typedef struct ExecState_ ExecState;

typedef enum {
    CMD_BATCH_MODE     = 1 << 0,
    CMD_ECHO_TO_STDOUT = 1 << 1,
    CMD_STACKING       = 1 << 2
} CmdEchoFlags;

typedef enum {
    OPT_BATCH = 1,
    OPT_HELP,
    OPT_VERSION,
    OPT_RUNIT,
    OPT_DBOPEN,
    OPT_WEBDB,
    OPT_DUMP
} ProgramOptions;

typedef enum {
    ENGLISH = 1,
    BASQUE
} ForcedLangs;

typedef enum {
    CONSOLE_EXEC      = 1 << 0,
    SCRIPT_EXEC       = 1 << 1,
    SESSION_EXEC      = 1 << 2,
    INCLUDE_EXEC      = 1 << 3,
    FUNCTION_EXEC     = 1 << 4
} ExecFlags;

#define HIDDEN_COMMAND(c) (c == ADDTO || \
                           c == FUNCERR || \
                           c == OMITFROM || \
                           c == REMEMBER)
    
/* functions follow */

int gretl_cmd_init (CMD *cmd);

void gretl_cmd_free (CMD *cmd);

CMD *gretl_cmd_new (void);

void gretl_cmd_destroy (CMD *cmd);

void gretl_cmd_set_context (CMD *cmd, int ci);

void gretl_cmd_destroy_context (CMD *cmd);

char *gretl_cmd_get_savename (char *sname);

gretlopt gretl_cmd_get_opt (const CMD *cmd);

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt);

int parse_command_line (char *line, CMD *cmd, double ***pZ, 
			DATAINFO *pdinfo); 

int get_command_index (char *line, CMD *cmd, const DATAINFO *pdinfo);

int command_number (const char *cmd);

int help (const char *cmdword, const char *helpfile, PRN *prn);

int parseopt (const char **argv, int argc, char *fname, 
	      int *force_lang);

int gretl_shell (const char *arg);

void echo_cmd (const CMD *cmd, const DATAINFO *pdinfo, const char *line, 
	       unsigned char flags, PRN *prn);

void echo_function_call (const char *line, unsigned char flags, PRN *prn);

int gretl_cmd_exec (ExecState *s, double ***pZ, DATAINFO **ppdinfo,
		    PRN *prn);

int call_pca_plugin (VMatrix *corrmat, double ***pZ,
		     DATAINFO *pdinfo, gretlopt *pflag,
		     PRN *prn);

int ready_for_command (const char *line);

void safe_print_line (const char *line, PRN *prn);

#endif /* INTERACT_H */



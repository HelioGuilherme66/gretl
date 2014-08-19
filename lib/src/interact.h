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

/* interact.h for gretl */

#ifndef INTERACT_H
#define INTERACT_H

#define MAXSAVENAME 32
#define CMD_NULL    -1
#define CMD_COMMENT -2
#define CMD_MASKED  -3

typedef struct CMD_ CMD;
typedef struct ExecState_ ExecState;

typedef enum {
    OPT_BATCH   = 1 << 0,
    OPT_HELP    = 1 << 1,
    OPT_VERSION = 1 << 2,
    OPT_RUNIT   = 1 << 3,
    OPT_DBOPEN  = 1 << 4,
    OPT_WEBDB   = 1 << 5,
    OPT_FNPKG   = 1 << 6,
    OPT_DUMP    = 1 << 7,
    OPT_DEBUG   = 1 << 8,
    OPT_QUIET   = 1 << 9,
    OPT_ENGLISH = 1 << 10,
    OPT_BASQUE  = 1 << 11,
    OPT_MAKEPKG = 1 << 12
} ProgramOptions;

typedef enum {
    CONSOLE_EXEC      = 1 << 0,
    SCRIPT_EXEC       = 1 << 1,
    INCLUDE_EXEC      = 1 << 2,
    FUNCTION_EXEC     = 1 << 3,
    DEBUG_EXEC        = 1 << 4,
    CALLBACK_EXEC     = 1 << 5,
    INIT_EXEC         = 1 << 6
} ExecFlags;

#define HIDDEN_COMMAND(c) (c == FUNCERR || c == FUNCRET)
    
/* functions follow */

int gretl_cmd_init (CMD *cmd);

void gretl_cmd_free (CMD *cmd);

CMD *gretl_cmd_new (void);

void gretl_cmd_destroy (CMD *cmd);

void gretl_cmd_set_context (CMD *cmd, int ci);

void gretl_cmd_destroy_context (CMD *cmd);

const char *gretl_cmd_get_savename (CMD *cmd);

gretlopt gretl_cmd_get_opt (const CMD *cmd);

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt);

int filter_comments (char *s, CMD *cmd);

int parse_command_line (char *line, CMD *cmd, DATASET *dset, void *ptr); 

int get_command_index (char *line, CMD *cmd);

int command_number (const char *cmd);

void echo_command (CMD *cmd, const DATASET *dset, 
		   const char *line, PRN *prn);

void gretl_record_command (CMD *cmd, const DATASET *dset, 
			   const char *line, PRN *prn);

void safe_print_line (const char *line, int *plen, PRN *prn);

int gretl_cmd_exec (ExecState *s, DATASET *dset);

int call_pca_plugin (VMatrix *cmat, DATASET *dset, 
		     gretlopt opt, PRN *prn);

int gretl_shell_grab (const char *arg, char **sout);

void set_plot_produced (void);

#endif /* INTERACT_H */



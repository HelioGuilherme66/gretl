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

/* cmd_private.c for libgretl: functions used in a small set of
   libgretl C files concerned with command execution.
*/

#include "libgretl.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "cmd_private.h"

#define PMDEBUG 0

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
    s->more = NULL;

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
           executing a function should not change the 'last
           model' state at caller level. To achieve this we
           need to take out a 'private' reference to the
           model, stored in the ExecState, and then remove
           it from last model position for the present.
        */
        s->prev_model = get_last_model(&s->prev_type);
        if (s->prev_model != NULL) {
#if PMDEBUG
            fprintf(stderr, "ExecState %p: set prev_model %p\n",
                    (void *) s, s->prev_model);
#endif
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

static int try_compile_func_genr (ExecState *s,
                                  DATASET *dset,
                                  void *ptr,
                                  int *done)
{
    GENERATOR **pgen = ptr;
    GretlType gtype = s->cmd->gtype;
    const char *line = s->cmd->vstart;
    gretlopt gopt = OPT_NONE;
    int err = 0;

    if (s->cmd->opt & OPT_O) {
	gopt |= OPT_O;
    }

    *pgen = genr_compile(line, dset, gtype, gopt, s->prn, &err);
    if (!err && *pgen != NULL) {
        *done = 1;
    } else if (err == E_EQN) {
	/* may be a non-compilable special such as "genr time",
           or perhaps a bare declaration */
        gretl_error_clear();
	err = 0;
    }

    return err;
}

/* Called by functions, and by scripts executed from within
   functions. Augmented 2022-08-11 to support a 3rd argument,
   for use in a function that is being called from a loop or
   internal iteration: we'll attempt to "compile" GENR
   statements in this context, for reuse when the function is
   next called.
*/

int maybe_exec_line (ExecState *s, DATASET *dset, void *ptr)
{
    int done = 0;
    int err = 0;

    if (string_is_blank(s->line)) {
        return 0;
    }

    if (gretl_compiling_loop()) {
        err = get_command_index(s, LOOP);
    } else {
        err = parse_command_line(s, dset, ptr);
        if (s->cmd->ci == GENR && !err && ptr != NULL) {
            if (!(s->cmd->flags & CMD_SUBST)) {
                err = try_compile_func_genr(s, dset, ptr, &done);
            }
        }
    }

    if (err) {
	errmsg(err, s->prn);
	if (s->cmd->flags & CMD_CATCH) {
	    set_gretl_errno(err);
	    return 0;
	} else {
	    return err;
	}
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
    } else if (!done) {
        /* note: error messages may be printed to s->prn */
        err = gretl_cmd_exec(s, dset);
    }

    return err;
}

#define STRUCTURE_DEBUG 0

/* Takes an array of "stmt" structs (statements or lines) from a
   function or loop and examines them for structure, in terms of
   conditionality (if/elif/else/endif). The object of the exercise is
   to figure out which line we should skip to in case of a false
   if-condition, and on reaching the end of a block enabled by a true
   if-condition. We thereby construct a set of efficiency-enhancing
   "gotos", which are recorded in the "next" member of each relevant
   line.

   In the case of functions we also determine and record embedded loop
   structure. (Note the asymmetry: a function may contain one or more
   loops but a loop cannot contain a function definition.)

   The @context argument should be FUNC or LOOP to distinguish the
   cases. The @name argument is used only in debugging, to identify
   the function we're working on; in the LOOP case we just pass
   "loop".
*/

int statements_get_structure (stmt *lines,
			      int n_lines,
			      int context,
			      const char *name)
{
    stmt *line;
    int *match_start = NULL;
    int *match_end = NULL;
    int *next_from_loop = NULL;
    int d_max = 0, ld_max = 0;
    int d = 0, ld = 0;
    int last_endif = 0;
    int first_if = -1;
    int target, src;
    int i, j, ci;
    int err = 0;

    /* first pass: determine the maximum depth of
       conditionality and looping
    */
    for (i=0; i<n_lines; i++) {
	ci = lines[i].ci;
	if (context == FUNC && ci == LOOP) {
            if (++ld > ld_max) {
                ld_max = ld;
            }
        } else if (context == FUNC && ci == ENDLOOP) {
            ld--;
        } else if (ld > 0) {
            continue;
	} else if (ci == IF) {
            if (first_if < 0) {
                first_if = i;
            }
            if (++d > d_max) {
                d_max = d;
            }
	} else if (ci == ENDIF) {
            last_endif = i;
            d--;
	}
    }

#if STRUCTURE_DEBUG
    fprintf(stderr, "\n%s: max if-depth %d, max loop-depth %d\n",
            name, d_max, ld_max);
#endif

    if (d != 0 || ld != 0) {
	fprintf(stderr, "broken structure: d = %d, ld = %d at end\n", d, ld);
	return E_PARSE;
    }

    if (d_max > 0) {
        match_start = gretl_list_new(d_max);
	fprintf(stderr, "match_start allocated for d_max = %d\n", d_max);
        match_end = gretl_list_new(d_max);
    }
    if (ld_max > 0) {
        next_from_loop = gretl_list_new(ld_max);
    }

    /* second pass: analysis */
    for (i=0; i<n_lines && !err; i++) {
	line = &lines[i];
	if (context == FUNC && line->ci == LOOP) {
            ld++;
            next_from_loop[ld] = i;
	} else if (context == FUNC && line->ci == ENDLOOP) {
            if (ld == 0) {
                err = 1;
            } else {
                j = next_from_loop[ld];
                lines[j].next = i;
                ld--;
            }
        } else if (ld > 0) {
            continue;
	} else if (line->ci == IF) {
	    d++;
	    match_start[d] = i;
	} else if (line->ci == ENDIF) {
	    if (d == 0) {
		err = 1;
	    } else {
                line->next = -d;
		j = match_start[d];
		lines[j].next = i;
                d--;
	    }
	} else if (line->ci == ELIF) {
	    if (d == 0) {
		err = 1;
	    } else {
                if (lines[i-1].ci != ENDIF) {
                    lines[i-1].next = -d;
                }
		j = match_start[d];
		lines[j].next = i;
		match_start[d] = i;
	    }
	} else if (line->ci == ELSE) {
	    if (d == 0) {
		err = 1;
	    } else {
		if (lines[i-1].ci != ENDIF) {
                    lines[i-1].next = -d;
		}
		j = match_start[d];
		lines[j].next = i;
		match_start[d] = i;
	    }
	}
    }

#if 0 /* not yet: apparently this can be dodgy (e.g. dbnomics_sample.inp */
    /* third pass: fill in goto's for true-block terminators */
    d = target = src = 0;
    for (i=last_endif; i>=first_if; i--) {
        line = &lines[i];
        if (line->ci == ENDIF) {
            d++;
	    fprintf(stderr, "match_start: writing to element %d\n", d);
            target = match_start[d] = line->next;
            src = match_end[d] = i;
        } else if (line->ci == IF) {
            target = (d == 0)? 0 : match_start[d];
            src = (d == 0)? 0 : match_end[d];
	    d--;
        } else if (target < 0 && line->next == target) {
            line->next = src;
        }
    }
#endif

    free(match_start);
    free(match_end);
    free(next_from_loop);

#if STRUCTURE_DEBUG
    /* display what we figured out */
    fputc('\n', stderr);
    for (i=0; i<n_lines; i++) {
	line = &lines[i];
        j = line->next;
	if (j <= 0) {
	    continue;
	}
	if (line->ci == IF || line->ci == ELIF || line->ci == ELSE) {
	    fprintf(stderr, "L%d ('%s'): next-on-false = %d ('%s')\n",
		    i, line->s, j, lines[j].s);
        } else if (context == FUNC && line->ci == LOOP) {
	    fprintf(stderr, "L%d ('%s'): end-of-loop = %d\n", i, line->s, j);
	} else {
	    fprintf(stderr, "L%d ('%s'): next-on-true = %d ('%s')\n",
		    i, line->s, j, lines[j].s);
	}
    }
#endif

    return err;
}

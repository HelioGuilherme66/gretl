/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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

/* errors.c - error messages for gretl */

#include "libgretl.h"
#include "libset.h"

int gretl_errno;
char gretl_errmsg[ERRLEN];
char gretl_msg[ERRLEN];

const char *gretl_error_messages[] = {
    NULL,
    NULL,
    N_("Data error"),                                            /* E_DATA = 2 */
    N_("Exact or near collinearity encountered"),                /* E_SINGULAR */
    N_("Insufficient degrees of freedom for regression"),        /* E_DF */
    N_("Y-prime * Y equals zero"),                               /* E_YPY */
    N_("Dependent variable is all zeros, aborting regression"),  /* E_ZERO */
    N_("Total sum of squares was not positive"),                 /* E_TSS */
    N_("Sum of squared residuals negative!"),                    /* E_ESS */
    N_("Unbalanced parentheses in genr command"),                /* E_UNBAL */
    N_("Sorry, command not available for this estimator"),       /* E_NOTIMP */
    N_("Unspecified error -- FIXME"),                            /* E_UNSPEC */
    N_("Syntax error in genr formula"),                          /* E_SYNTAX */
    N_("This command won't work with the current periodicity"),  /* E_PDWRONG */
    N_("Error attempting to open file"),                         /* E_FOPEN */
    N_("Out of memory error"),                                   /* E_ALLOC */
    N_("No formula supplied in genr"),                           /* E_EQN */
    N_("Unknown variable name in command"),                      /* E_UNKVAR */
    N_("The observations specified for the regression "
       "exceed those in the data set"),                          /* E_NODATA */
    N_("Command has insufficient arguments"),                    /* E_ARGS */
    N_("This command is implemented only for OLS models"),       /* E_OLSONLY */
    N_("Invalid argument for function"),                         /* E_INVARG */
    N_("Invalid sample split for Chow test"),                    /* E_SPLIT */
    N_("Syntax error in command line"),                          /* E_PARSE */
    N_("No independent variables left after omissions"),         /* E_NOVARS */
    N_("No independent variables were omitted"),                 /* E_NOOMIT */
    N_("Can't do this: some vars in original model "
       "have been redefined"),                                   /* E_VARCHANGE */
    N_("No new independent variables were added"),               /* E_NOADD */
    N_("One or more \"added\" vars were already present"),       /* E_ADDDUP */
    N_("Error generating logarithms"),                           /* E_LOGS */
    N_("Error generating squares"),                              /* E_SQUARES */
    N_("Error generating lagged variables"),                     /* E_LAGS */
    N_("Attempting to take square root of negative number"),     /* E_SQRT */
    N_("Excessive exponent in genr formula"),                    /* E_HIGH */
    N_("Weight variable is all zeros, aborting regression"),     /* E_WTZERO */
    N_("Weight variable contains negative values"),              /* E_WTNEG */
    N_("Need valid starting and ending observations"),           /* E_OBS */
    N_("You must include a constant in this sort of model"),     /* E_NOCONST */
    N_("There were missing observations for the added "
       "variable(s).\nReset the sample and rerun the original "
       "regression first"),                                      /* E_MISS */
    N_("The statistic you requested is not available"),          /* E_BADSTAT */
    N_("Missing sub-sample information; can't merge data"),      /* E_NOMERGE */
    N_("The convergence criterion was not met"),                 /* E_NOCONV */
    N_("The operation was canceled"),                            /* E_CANCEL */
    N_("Missing values encountered"),                            /* E_MISSDATA */
    N_("Not a Number in calculation"),                           /* E_NAN */
    N_("Matrices not conformable for operation"),                /* E_NONCONF */
    N_("Data types not conformable for operation"),              /* E_TYPES */
    N_("Wrong data type"),                                       /* E_DATATYPE */
    N_("Incompatible options"),                                  /* E_BADOPT */
    NULL,                                                        /* E_DB_DUP */
    NULL,                                                        /* E_OK */
    NULL                                                         /* E_MAX */
};

/**
 * get_errmsg:
 * @errcode: gretl error code (see #error_codes).
 * @errtext: pre-allocated string or NULL.
 * @prn: gretl printing struct.
 *
 * Print an error message, given an error code number.  The message
 * is printed to the string variable errtext, if it is non-NULL,
 * or otherwise to the printing struct @prn.
 * 
 * Returns: the error text string, or NULL if @errtext is NULL.
 */

char *get_errmsg (const int errcode, char *errtext, PRN *prn)
{
    char *msg = NULL;

    if (errcode > 0 && errcode < E_MAX) {
	if (gretl_error_messages[errcode] != NULL) {
	    if (errtext != NULL) {
		strcpy(errtext, _(gretl_error_messages[errcode]));
		msg = errtext;
	    } else {
		pprintf(prn, "%s\n", _(gretl_error_messages[errcode]));
	    }
	}
    } else {
	fprintf(stderr, "get_errmsg: out of bounds errcode %d\n", 
		errcode);
    }

    return msg;
}

/**
 * errmsg:
 * @errcode: gretl error code (see #error_codes).
 * @prn: gretl printing struct.
 *
 * Print an error message looked up from a given an error code number, 
 * or a more specific error message if available.  
 * 
 */

void errmsg (const int errcode, PRN *prn)
{
    if (*gretl_errmsg == '\0') {
	get_errmsg(errcode, NULL, prn);
    } else {
	pprintf(prn, "%s\n", gretl_errmsg);
    }
}

int get_gretl_errno (void)
{
    return gretl_errno;
}

const char *get_gretl_errmsg (void)
{
    return gretl_errmsg;
}

int print_gretl_errmsg (PRN *prn)
{
    int ret = 0;

    if (*gretl_errmsg != '\0') {
	pprintf(prn, "%s\n", gretl_errmsg);
	ret = 1;
    } else if (get_errmsg(gretl_errno, NULL, prn)) {
	ret = 1;
    }

    return ret;
}

/**
 * gretl_errmsg_set:
 * @str: an error message.
 *
 * If gretl_errmsg is currently blank, copy the given string into
 * the message space.
 * 
 */

void gretl_errmsg_set (const char *str)
{
    if (*gretl_errmsg == '\0') {
	strncat(gretl_errmsg, str, ERRLEN - 1);
    }
}

void gretl_errmsg_clear (void)
{
    *gretl_errmsg = '\0';
}


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

enum tx_objects {
    TX_SA,    /* save seasonally adjusted series */
    TX_TR,    /* save trend/cycle */
    TX_IR,    /* save irregular component */
    TX_LN,    /* save linearized series */
    TRIGRAPH, /* graph showing some/all of the above */
    TEXTOUT,  /* for full text output */
    TX_MAXOPT /* sentinel */
};

typedef struct _common_opt_info common_opt_info;
typedef struct _x13a_opts x13a_opts;
typedef struct _tx_request tx_request;

struct _common_opt_info {
    GtkWidget *check;
    char save;
    unsigned short v;
    char savename[VNAMELEN];
};

struct _x13a_opts {
    int logtrans; /* log transformation */
    int outliers; /* outlier detection */
    int trdays;   /* trading days */
    int wdays;    /* working days */
    int easter;   /* Easter effect */
    int seats;    /* Use SEATS with x13as (vs X-11) */
    int airline;  /* use the airline model (vs auto model selection) */
    int output;   /* which output series is wanted? */
    int verbose;  /* verbosity level */
    int save_spc; /* option to save spc file content */
    double critical; /* critical value for outliers (vs auto) */
    int *savelist;   /* spec for series to save */
    guint8 *aspec;   /* arima spec, as array of int */
};

struct _tx_request {
    int prog;          /* tramo vs x13as */
    GtkWidget *dialog;
    void (*helpfunc);
    common_opt_info opts[TX_MAXOPT];
    char yname[VNAMELEN];
    void *gui;
    gretlopt *popt;
    int savevars;
    int pd;
    int seasonal_ok;
    x13a_opts xopt;
};

int add_tramo_options (tx_request *request, GtkWidget *vbox);

int print_tramo_options (tx_request *request, FILE *fp);

const char *get_tramo_save_string (int i);

void sensitize_tx_entry (GtkToggleButton *b, GtkWidget *w);

void update_tx_savename (GtkEntry *entry, char *name);

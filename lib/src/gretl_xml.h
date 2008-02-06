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

#ifndef GRETL_XML_H
#define GRETL_XML_H

#ifndef GRETLCLI

#define XUC const xmlChar *

void gretl_xml_header (FILE *fp);

int gretl_xml_open_doc_root (const char *fname,
			     const char *rootname,
			     xmlDocPtr *pdoc, 
			     xmlNodePtr *pnode);

void gretl_xml_put_int (const char *tag, int i, FILE *fp);

void gretl_xml_put_double (const char *tag, double x, FILE *fp);

void gretl_xml_put_double_array (const char *tag, double *x, int n,
				 FILE *fp);

void gretl_xml_put_strings_array (const char *tag, const char **strs, int n,
				  FILE *fp);

void gretl_xml_put_named_list (const char *name, const int *list, FILE *fp);

void gretl_xml_put_tagged_list (const char *tag, const int *list, FILE *fp);

int gretl_xml_put_tagged_string (const char *tag, const char *str, 
				 FILE *fp);

int gretl_xml_put_raw_string (const char *str, FILE *fp);

void gretl_xml_put_matrix (const gretl_matrix *m, const char *name, 
			   FILE *fp);

int gretl_xml_get_prop_as_int (xmlNodePtr node, const char *tag,
			       int *i);

int gretl_xml_get_prop_as_char (xmlNodePtr node, const char *tag,
				char *c);

int gretl_xml_get_prop_as_uchar (xmlNodePtr node, const char *tag,
				 unsigned char *u);

int gretl_xml_get_prop_as_double (xmlNodePtr node, const char *tag,
				  double *x);

int gretl_xml_get_prop_as_string (xmlNodePtr node, const char *tag,
				  char **pstr);

int gretl_xml_get_prop_as_bool (xmlNodePtr node, const char *tag);

int gretl_xml_node_get_int (xmlNodePtr node, xmlDocPtr doc, int *i);

int gretl_xml_node_get_double (xmlNodePtr node, xmlDocPtr doc, 
			       double *x);

int *gretl_xml_node_get_list (xmlNodePtr node, xmlDocPtr doc, int *err);

int gretl_xml_node_get_string (xmlNodePtr node, xmlDocPtr doc, 
			       char **pstr);

int gretl_xml_node_get_trimmed_string (xmlNodePtr node, xmlDocPtr doc, 
				       char **pstr);

int *gretl_xml_get_int_array (xmlNodePtr node, xmlDocPtr doc,
			      int *nelem, int *err);

double *gretl_xml_get_double_array (xmlNodePtr node, xmlDocPtr doc,
				    int *nelem, int *err);

cmplx *gretl_xml_get_cmplx_array (xmlNodePtr node, xmlDocPtr doc,
				  int *nelem, int *err);

char **gretl_xml_get_strings_array (xmlNodePtr node, xmlDocPtr doc,
				    int *nelem, int *err);

gretl_matrix *gretl_xml_get_matrix (xmlNodePtr node, xmlDocPtr doc, 
				    int *err);

int gretl_xml_get_submask (xmlNodePtr node, xmlDocPtr doc, char **pmask);

#endif /* !GRETLCLI */

int gretl_write_matrix_as_gdt (const char *fname, 
			       const gretl_matrix *X,
			       const char **varnames, 
			       const char **labels);

int gretl_write_gdt (const char *fname, const int *list, 
		     const double **Z, const DATAINFO *pdinfo, 
		     GretlDataFormat fmt, PATHS *ppaths);

int gretl_read_gdt (double ***pZ, DATAINFO *pdinfo, char *fname,
		    PATHS *ppaths, gretlopt opt, PRN *prn);

char *gretl_get_gdt_description (const char *fname);

int load_user_matrix_file (const char *fname);

int load_user_XML_file (const char *fname);

#endif /* GRETL_XML_H */

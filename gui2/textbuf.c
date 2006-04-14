/*
 *  Copyright (c) by Allin Cottrell
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

#include "gretl.h"
#include "textbuf.h"

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
# include <gtksourceview/gtksourcelanguage.h>
# include <gtksourceview/gtksourcelanguagesmanager.h>
#endif

#define GUIDE_PAGE 999

enum {
    PLAIN_TEXT,
    BLUE_TEXT,
    RED_TEXT
};

#define gui_help(r) (r == GUI_HELP || r == GUI_HELP_EN)

void text_set_cursor (GtkWidget *w, GdkCursorType cspec)
{
    GdkWindow *win = gtk_text_view_get_window(GTK_TEXT_VIEW(w),
                                              GTK_TEXT_WINDOW_TEXT);

    if (cspec == 0) {
	gdk_window_set_cursor(win, NULL);
    } else {
	GdkCursor *cursor = gdk_cursor_new(cspec);

	gdk_window_set_cursor(win, cursor);
	gdk_cursor_destroy(cursor);
    } 
}

void cursor_to_top (windata_t *vwin)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w)); 
    GtkTextIter start;
    GtkTextMark *mark;

    gtk_text_buffer_get_start_iter(buf, &start);
    gtk_text_buffer_place_cursor(buf, &start);
    mark = gtk_text_buffer_create_mark(buf, NULL, &start, FALSE);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				 mark, 0.0, FALSE, 0, 0);
}

gint get_char_width (GtkWidget *widget)
{
    PangoLayout *pl;
    PangoContext *pc;
    GtkRcStyle *style;
    int width;

    pc = gtk_widget_get_pango_context(widget);
    style = gtk_widget_get_modifier_style(widget);
    pango_context_set_font_description(pc, style->font_desc);

    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, "X", 1);
    pango_layout_get_pixel_size(pl, &width, NULL);

    g_object_unref(G_OBJECT(pl));

    return width;
}

gchar *textview_get_text (GtkTextView *view)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    tbuf = gtk_text_view_get_buffer(view);
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);

    return gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
}

int viewer_char_count (windata_t *vwin)
{
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    return gtk_text_buffer_get_char_count(tbuf);
}

void text_paste (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old;
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gchar *undo_buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));

    old = g_object_get_data(G_OBJECT(vwin->w), "undo");
    g_free(old);

    g_object_set_data(G_OBJECT(vwin->w), "undo", undo_buf);

    gtk_text_buffer_paste_clipboard(buf, gtk_clipboard_get(GDK_NONE),
				    NULL, TRUE);
}

void text_undo (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old = NULL;

#ifdef USE_GTKSOURCEVIEW
    if (vwin->sbuf != NULL) {
	if (gtk_source_buffer_can_undo(vwin->sbuf)) {
	    gtk_source_buffer_undo(vwin->sbuf);
	} else {
	    errbox(_("No undo information available"));
	}
	return;
    }
#endif
    
    old = g_object_steal_data(G_OBJECT(vwin->w), "undo");

    if (old == NULL) {
	errbox(_("No undo information available"));
    } else {
	GtkTextBuffer *buf;
	GtkTextIter start, end;
	GtkTextMark *ins;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

	ins = gtk_text_buffer_get_insert(buf);

	gtk_text_buffer_get_start_iter(buf, &start);
	gtk_text_buffer_get_end_iter(buf, &end);
	gtk_text_buffer_delete(buf, &start, &end);

	gtk_text_buffer_insert(buf, &start, old, strlen(old));
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				     ins, 0.0, TRUE, 0.1, 0.0);
	g_free(old);
    }
}

#ifdef USE_GTKSOURCEVIEW

static int 
gtk_source_buffer_load_file (GtkSourceBuffer *sbuf, 
			     const char *fname)
{
    FILE *fp;
    GtkTextIter iter;    
    char readbuf[MAXSTR], *chunk = NULL;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) return 1;

    gtk_source_buffer_begin_not_undoable_action (sbuf);

    gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sbuf), "", 0);
    gtk_text_buffer_get_iter_at_offset(GTK_TEXT_BUFFER(sbuf), &iter, 0);

    memset(readbuf, 0, sizeof readbuf);

    while (fgets(readbuf, sizeof readbuf, fp)) {
	int len;

# ifdef ENABLE_NLS
	if (!g_utf8_validate(readbuf, -1, NULL)) {
	    chunk = my_locale_to_utf8(readbuf);
	    if (chunk == NULL) {
		continue;
	    }
	} else {
	    chunk = readbuf;
	}
# else
	chunk = readbuf;
# endif /* ENABLE_NLS */

	/* check that this works */
	len = strlen(chunk);
	if (chunk[len - 2] == '\r') {
	    chunk[len - 2] = '\n';
	    chunk[len - 1] = '\0';
	}

	gtk_text_buffer_insert(GTK_TEXT_BUFFER(sbuf), &iter, chunk, -1);
	memset(readbuf, 0, sizeof readbuf);
	if (chunk != readbuf) {
	    g_free(chunk);
	    chunk = NULL;
	}
    }

    fclose(fp);
	
    gtk_source_buffer_end_not_undoable_action(sbuf);

    gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sbuf), FALSE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(GTK_TEXT_BUFFER(sbuf), &iter);
    gtk_text_buffer_place_cursor(GTK_TEXT_BUFFER(sbuf), &iter);

    return 0;
}

void source_buffer_insert_file (GtkSourceBuffer *sbuf, const char *filename,
				int role)
{
    GtkSourceLanguagesManager *manager;    
    GtkSourceLanguage *language = NULL;
		
    manager = g_object_get_data(G_OBJECT (sbuf), "languages-manager");

    if (role == GR_PLOT) {
	language = 
	    gtk_source_languages_manager_get_language_from_mime_type 
	    (manager, "application/x-gnuplot");
    } else {
	language = 
	    gtk_source_languages_manager_get_language_from_mime_type 
	    (manager, "application/x-gretlsession");
    }

    if (language == NULL) {
	g_object_set(G_OBJECT(sbuf), "highlight", FALSE, NULL);
    } else {
	g_object_set(G_OBJECT(sbuf), "highlight", TRUE, NULL);
	gtk_source_buffer_set_language(sbuf, language);
    }

    gtk_source_buffer_load_file(sbuf, filename);
}

void create_source (windata_t *vwin, GtkSourceBuffer **buf, 
		    int hsize, int vsize, gboolean editable)
{
    GtkSourceLanguagesManager *lm;
    GtkSourceBuffer *sbuf;
    GtkSourceTagStyle *tagstyle;
    GdkColormap *cmap;
    GdkColor blue;

    /* set up paren-matching in blue */

    cmap = gdk_colormap_get_system ();
    gdk_color_parse ("blue", &blue);
    gdk_colormap_alloc_color (cmap, &blue, FALSE, TRUE);

    lm = gtk_source_languages_manager_new ();
    tagstyle = gtk_source_tag_style_new ();
    
    sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new(NULL));
    g_object_ref (lm);
    g_object_set_data_full (G_OBJECT (sbuf), "languages-manager",
			    lm, (GDestroyNotify) g_object_unref); 
    g_object_unref (lm); 

    tagstyle->mask = GTK_SOURCE_TAG_STYLE_USE_FOREGROUND;
    tagstyle->foreground = blue;
    g_object_set_data_full (G_OBJECT (sbuf), "tag-style",
			    tagstyle, 
			    (GDestroyNotify) gtk_source_tag_style_free); 
    gtk_source_buffer_set_bracket_match_style(sbuf, tagstyle);
    gtk_source_buffer_set_check_brackets(sbuf, TRUE);

    vwin->w = gtk_source_view_new_with_buffer(sbuf);
    *buf = sbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vwin->w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vwin->w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vwin->w), 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->w), fixed_font);

    hsize *= get_char_width(vwin->w);
    hsize += 48;

    gtk_window_set_default_size (GTK_WINDOW(vwin->dialog), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), editable);

    g_object_unref(cmap);
}

#endif /* USE_GTKSOURCEVIEW */

static gchar *my_utf_string (char *t)
{
    static gchar *s = NULL;
    GError *error = NULL;
    gsize r_bytes, w_bytes;
    unsigned char *c;
    const char *fc;
    const char *smb;
    gchar *from_codeset = NULL;
    
    if (t == NULL || *t == '\0') return t;

    if (g_utf8_validate(t, -1, NULL)) return t;   
    
    /* so we got a non-UTF-8 */

    smb = getenv("SMB_CODESET");
    if (smb != NULL && *smb != '\0') {
	from_codeset = g_strdup(smb);
    } else {
    	g_get_charset(&fc);
    	if (fc) {
	    from_codeset = g_strdup(fc);
    	} else {
	    from_codeset = g_strdup("ISO-8859-1");
	}
    }
    
    if (!strcmp(from_codeset, "ISO-")) {
	g_free(from_codeset);
	from_codeset = g_strdup("ISO-8859-1");
    }  
  
    if (s) g_free(s);

    for (c = (unsigned char *) t; *c != 0; c++) {
	if (*c < 32 && *c != '\n') {
	    *c = ' ';
	}
    }

    s = g_convert(t, strlen(t), "UTF-8", from_codeset, &r_bytes, &w_bytes,
		  &error);

    g_free(from_codeset);

    if (s == NULL) {
	s = g_strdup(t);
	for (c = (unsigned char *) s; *c != 0; c++) {
	    if (*c > 128) {
		*c = '?';
	    }
	}
    }

    if (error) {
        printf("DBG: %s. Codeset for system is: %s\n",
	       error->message, from_codeset);
        printf("DBG: You should set the environment variable "
	       "SMB_CODESET to ISO-8859-1\n");
	g_error_free(error);
    }

    return s;
}

static GtkTextTagTable *gretl_tags_new (void)
{
    GtkTextTagTable *table;
    GtkTextTag *tag;

    table = gtk_text_tag_table_new(); 

    tag = gtk_text_tag_new("bluetext");
    g_object_set(tag, "foreground", "blue", NULL);
    gtk_text_tag_table_add(table, tag);
    
    tag = gtk_text_tag_new("redtext");
    g_object_set(tag, "foreground", "red", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("title");
    g_object_set(tag, "justification", GTK_JUSTIFY_CENTER,
		 "pixels_above_lines", 20,
		 "family", "sans",
		 "size", 15 * PANGO_SCALE, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("sansbold");
    g_object_set(tag, "family", "sans",
		 "weight", PANGO_WEIGHT_BOLD, 
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("italic");
    g_object_set(tag, "family", "sans",
		 "style", PANGO_STYLE_ITALIC, 
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("replaceable");
    g_object_set(tag, "family", "sans",
		 "style", PANGO_STYLE_ITALIC, 
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("literal");
    g_object_set(tag, "family", "monospace", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("text");
    g_object_set(tag, "family", "sans", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("indented");
    g_object_set(tag, "left_margin", 30, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("code");
    g_object_set(tag, "family", "monospace", NULL);
    gtk_text_tag_table_add(table, tag);

    return table;
}

static GtkTextBuffer *gretl_text_buf_new (void)
{
    static GtkTextTagTable *tags = NULL;
    GtkTextBuffer *tbuf; 

    if (tags == NULL) {
	tags = gretl_tags_new();
    }

    tbuf = gtk_text_buffer_new(tags);

    return tbuf;
}

void text_buffer_insert_colorized_buffer (GtkTextBuffer *tbuf, PRN *prn)
{
    GtkTextIter iter;    
    int thiscolor, nextcolor;
    const char *pbuf;
    char readbuf[MAXSTR];

    pbuf = gretl_print_get_buffer(prn);

    thiscolor = PLAIN_TEXT;
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    bufgets_init(pbuf);

    while (bufgets(readbuf, sizeof readbuf, pbuf)) {

	if (ends_with_backslash(readbuf)) {
	    nextcolor = BLUE_TEXT;
	} else {
	    nextcolor = PLAIN_TEXT;
	}

	if (*readbuf == '#' || *readbuf == '?' || *readbuf == '>') {
	    thiscolor = BLUE_TEXT;
	} 

	if (thiscolor == BLUE_TEXT) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
						     readbuf, -1,
						     "bluetext", NULL);
	} else {
	    gtk_text_buffer_insert(tbuf, &iter, readbuf, -1);
	}

	/* bufgets strips newlines */
	gtk_text_buffer_insert(tbuf, &iter, "\n", 1);

	thiscolor = nextcolor;
    }
}

void text_buffer_insert_file (GtkTextBuffer *tbuf, const char *fname, 
			      int role)
{
    FILE *fp;
    GtkTextIter iter;    
    int thiscolor, nextcolor;
    char readbuf[MAXSTR], *chunk;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) return;

    thiscolor = nextcolor = PLAIN_TEXT;

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    memset(readbuf, 0, sizeof readbuf);

    while (fgets(readbuf, sizeof readbuf, fp)) {
#ifdef ENABLE_NLS
	chunk = my_utf_string(readbuf);
#else
	chunk = readbuf;
#endif

	nextcolor = PLAIN_TEXT;
	
	if (role == SCRIPT_OUT && ends_with_backslash(chunk)) {
	    nextcolor = BLUE_TEXT;
	}

	if (*chunk == '?') {
	    thiscolor = (role == CONSOLE)? RED_TEXT : BLUE_TEXT;
	} else if (*chunk == '#') {
	    thiscolor = BLUE_TEXT;
	} 

	switch (thiscolor) {
	case PLAIN_TEXT:
	    gtk_text_buffer_insert(tbuf, &iter, chunk, -1);
	    break;
	case BLUE_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
						     chunk, -1,
						     "bluetext", NULL);
	    break;
	case RED_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
						     chunk, -1,
						     "redtext", NULL);
	    break;
	}

	thiscolor = nextcolor;
	memset(readbuf, 0, sizeof readbuf);
    }

    fclose(fp);
}

static void insert_link (GtkTextBuffer *tbuf, GtkTextIter *iter, 
			 const gchar *text, gint page,
			 const char *indent)
{
    GtkTextTag *ltag;

    if (page == GUIDE_PAGE) {
	ltag = gtk_text_buffer_create_tag(tbuf, NULL, "foreground", "blue", 
					  "family", "sans", NULL);
    } else if (indent != NULL) {
	ltag = gtk_text_buffer_create_tag(tbuf, NULL, "foreground", "blue", 
					  "left_margin", 30, NULL);
    } else {
	ltag = gtk_text_buffer_create_tag(tbuf, NULL, "foreground", "blue", NULL);
    }

    g_object_set_data(G_OBJECT(ltag), "page", GINT_TO_POINTER(page));
    gtk_text_buffer_insert_with_tags(tbuf, iter, text, -1, ltag, NULL);
}

static void follow_if_link (GtkWidget *tview, GtkTextIter *iter, gpointer p)
{
    GSList *tags = NULL, *tagp = NULL;

    tags = gtk_text_iter_get_tags(iter);

    for (tagp = tags; tagp != NULL; tagp = tagp->next) {
	GtkTextTag *tag = tagp->data;
	gint page = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(tag), "page"));

	if (page != 0) {
	    if (page == GUIDE_PAGE) {
		display_pdf_help(NULL, 1, NULL);
	    } else {
		plain_text_cmdref(p, page, NULL);
	    }
	    break;
	}
    }

    if (tags) {
	g_slist_free (tags);
    }
}

/* Help links can be activated by pressing Enter */

static gboolean cmdref_key_press (GtkWidget *tview, GdkEventKey *ev,
				  gpointer p)
{
    GtkTextIter iter;
    GtkTextBuffer *tbuf;

    switch (ev->keyval) {
    case GDK_Return: 
    case GDK_KP_Enter:
	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));
	gtk_text_buffer_get_iter_at_mark(tbuf, &iter, 
					 gtk_text_buffer_get_insert(tbuf));
	follow_if_link(tview, &iter, p);
	break;
    default:
	break;
    }

    return FALSE;
}

/* Help links can be activated by clicking */

static gboolean cmdref_event_after (GtkWidget *tview, GdkEvent *ev,
				    gpointer p)
{
    GtkTextIter start, end, iter;
    GtkTextBuffer *buffer;
    GdkEventButton *event;
    gint x, y;

    if (ev->type != GDK_BUTTON_RELEASE) {
	return FALSE;
    }

    event = (GdkEventButton *) ev;

    if (event->button != 1)
	return FALSE;

    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));

    /* don't follow a link if the user has selected something */
    gtk_text_buffer_get_selection_bounds(buffer, &start, &end);
    if (gtk_text_iter_get_offset(&start) != gtk_text_iter_get_offset(&end))
	return FALSE;

    gtk_text_view_window_to_buffer_coords(GTK_TEXT_VIEW(tview), 
					  GTK_TEXT_WINDOW_WIDGET,
					  event->x, event->y, &x, &y);

    gtk_text_view_get_iter_at_location(GTK_TEXT_VIEW(tview), &iter, x, y);

    follow_if_link(tview, &iter, p);

    return FALSE;
}

static GdkCursor *hand_cursor = NULL;
static GdkCursor *regular_cursor = NULL;

static void
set_cursor_if_appropriate (GtkTextView *tview, gint x, gint y)
{
    static gboolean hovering_over_link = FALSE;
    GSList *tags = NULL, *tagp = NULL;
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    gboolean hovering = FALSE;

    tbuf = gtk_text_view_get_buffer(tview);
    gtk_text_view_get_iter_at_location(tview, &iter, x, y);
  
    tags = gtk_text_iter_get_tags(&iter);

    for (tagp = tags; tagp != NULL; tagp = tagp->next) {
	GtkTextTag *tag = tagp->data;
	gint page = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(tag), "page"));

	if (page != 0) { 
	    hovering = TRUE;
	    break;
        }
    }

    if (hovering != hovering_over_link) {
	hovering_over_link = hovering;

	if (hovering_over_link) {
	    gdk_window_set_cursor(gtk_text_view_get_window(tview, GTK_TEXT_WINDOW_TEXT), 
				  hand_cursor);
	} else {
	    gdk_window_set_cursor(gtk_text_view_get_window(tview, GTK_TEXT_WINDOW_TEXT), 
				  regular_cursor);
	}
    }

    if (tags) {
	g_slist_free(tags);
    }
}

static gboolean 
cmdref_motion_notify (GtkWidget *tview, GdkEventMotion *event)
{
    gint x, y;

    gtk_text_view_window_to_buffer_coords(GTK_TEXT_VIEW(tview), 
					  GTK_TEXT_WINDOW_WIDGET,
					  event->x, event->y, &x, &y);
    set_cursor_if_appropriate(GTK_TEXT_VIEW(tview), x, y);
    gdk_window_get_pointer(tview->window, NULL, NULL, NULL);

    return FALSE;
}

static gboolean
cmdref_visibility_notify (GtkWidget *tview,  GdkEventVisibility *e)
{
    gint wx, wy, bx, by;
  
    gdk_window_get_pointer(tview->window, &wx, &wy, NULL);
    gtk_text_view_window_to_buffer_coords(GTK_TEXT_VIEW(tview), 
					  GTK_TEXT_WINDOW_WIDGET,
					  wx, wy, &bx, &by);
    set_cursor_if_appropriate(GTK_TEXT_VIEW(tview), bx, by);

    return FALSE;
}

static void maybe_connect_help_signals (windata_t *hwin, int en)
{
    int done = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(hwin->w), 
						 "sigs_connected"));

    if (hand_cursor == NULL) {
	hand_cursor = gdk_cursor_new(GDK_HAND2);
    }

    if (regular_cursor == NULL) {
	regular_cursor = gdk_cursor_new(GDK_XTERM);
    }    

    if (!done) {
	gpointer en_ptr = GINT_TO_POINTER(en);

	g_signal_connect(hwin->w, "key-press-event", 
			 G_CALLBACK(cmdref_key_press), en_ptr);
	g_signal_connect(hwin->w, "event-after", 
			 G_CALLBACK(cmdref_event_after), en_ptr);
	g_signal_connect(hwin->w, "motion-notify-event", 
			 G_CALLBACK(cmdref_motion_notify), NULL);
	g_signal_connect(hwin->w, "visibility-notify-event", 
			 G_CALLBACK(cmdref_visibility_notify), NULL);
	g_object_set_data(G_OBJECT(hwin->w), "sigs_connected", 
			  GINT_TO_POINTER(1));
    }
}

static void maybe_set_help_tabs (windata_t *hwin)
{
    int done = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(hwin->w), 
						 "tabs_set"));

    if (!done) {
	PangoTabArray *tabs;
	
	tabs = pango_tab_array_new(1, TRUE);
	pango_tab_array_set_tab(tabs, 0, PANGO_TAB_LEFT, 50);
	gtk_text_view_set_tabs(GTK_TEXT_VIEW(hwin->w), tabs);
	pango_tab_array_free(tabs);
	g_object_set_data(G_OBJECT(hwin->w), "tabs_set", GINT_TO_POINTER(1));
    }
}

static void cmdref_title_page (windata_t *hwin, GtkTextBuffer *tbuf, int en)
{
    const char *header = N_("Gretl Command Reference");
    GtkTextIter iter;
    int i, j, k;

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
					     (en)? header : _(header), -1,
					     "title", NULL);
    gtk_text_buffer_insert(tbuf, &iter, "\n\n\n", -1);

    j = 1;

    for (i=1; i<NC; i++) {
	const char *word;

	if (HIDDEN_COMMAND(i)) {
	    continue;
	}

	word = gretl_command_word(i);
	insert_link(tbuf, &iter, gretl_command_word(i), i, NULL);
	if (j++ % 8 == 0) {
	    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
	} else {
	    int n = 10 - strlen(word);

	    for (k=0; k<n; k++) {
		gtk_text_buffer_insert(tbuf, &iter, " ", -1);
	    }
	}
    }

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->w), tbuf);

    maybe_connect_help_signals(hwin, en);
    maybe_set_help_tabs(hwin);
}

static gint help_popup_click (GtkWidget *w, gpointer p)
{
    windata_t *hwin = (windata_t *) p;
    int action = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    gpointer enp = NULL;
    int page = 0;

    if (hwin->role == CLI_HELP_EN) {
	enp = GINT_TO_POINTER(1);
    }

    if (action == 2) {
	page = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(hwin->w), 
						 "backpage"));
    }

    plain_text_cmdref(enp, page, NULL);

    return FALSE;
}

GtkWidget *build_help_popup (windata_t *hwin)
{
    const char *items[] = {
	N_("Index"),
	N_("Back")
    };
    GtkWidget *pmenu = gtk_menu_new();
    GtkWidget *item;
    int i;

    for (i=0; i<2; i++) {
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_object_set_data(G_OBJECT(item), "action", GINT_TO_POINTER(i+1));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(help_popup_click),
			 hwin);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }

    return pmenu;
}

gboolean 
help_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p)
{
    GdkModifierType mods;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods);

    if (mods & GDK_BUTTON3_MASK) {
	windata_t *hwin = (windata_t *) p;

	if (hwin->active_var == 0) {
	    return TRUE;
	}

	if (hwin->popup) {
	    gtk_widget_destroy(hwin->popup);
	    hwin->popup = NULL;
	}

	hwin->popup = build_help_popup(hwin);

	if (hwin->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(hwin->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    gtk_signal_connect(GTK_OBJECT(hwin->popup), "destroy",
			       GTK_SIGNAL_FUNC(gtk_widget_destroyed), 
			       &hwin->popup);
	}

	return TRUE;
    }

    return FALSE;
}

enum {
    INSERT_NONE,
    INSERT_REF,
    INSERT_FIG,
    INSERT_REPL,
    INSERT_LIT,
    INSERT_ITAL,
    INSERT_TEXT,
    INSERT_PDFLINK
};

static void insert_help_figure (GtkTextBuffer *tbuf, GtkTextIter *iter,
				const char *fig)
{
    char figfile[FILENAME_MAX];
    GdkPixbuf *pixbuf;

    sprintf(figfile, "%shelpfigs%c%s.png", paths.gretldir,
	    SLASH, fig);

    pixbuf = gdk_pixbuf_new_from_file(figfile, NULL);

    if (pixbuf != NULL) {
	gtk_text_buffer_insert_pixbuf(tbuf, iter, pixbuf);
	g_object_unref(G_OBJECT(pixbuf));
    }
}

static void insert_tagged_text (GtkTextBuffer *tbuf, GtkTextIter *iter,
				const char *s, int ins, const char *indent)
{
    const char *ftag = NULL;

    switch (ins) {
    case INSERT_ITAL:
	ftag = "italic";
	break;
    case INSERT_REPL:
	ftag = "replaceable";
	break;
    case INSERT_LIT:
	ftag = "literal";
	break;
    default:
	break;
    }

    if (ftag != NULL) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 ftag, indent, NULL);
    }
}

static int get_instruction_and_string (const char *p, char *str)
{
    int ins = INSERT_NONE;
    *str = '\0';

    if (!strncmp(p, "ref", 3)) {
	ins = INSERT_REF;
    } else if (!strncmp(p, "fig", 3)) {
	ins = INSERT_FIG;
    } else if (!strncmp(p, "itl", 3)) {
	ins = INSERT_ITAL;
    } else if (!strncmp(p, "var", 3)) {
	ins = INSERT_REPL;
    } else if (!strncmp(p, "lit", 3)) {
	ins = INSERT_LIT;
    } else if (!strncmp(p, "pdf", 3)) {
	ins = INSERT_PDFLINK;
    }

    if (ins != INSERT_NONE) {
	int i = 0;

	p += 5;
	while (*p) {
	    if (*p == '"' && *(p+1) == '>') {
		str[i] = 0;
		break;
	    } else {
		str[i++] = *p++;
	    }
	}
    }

    return ins;
}

static int get_code_skip (const char *s)
{
    int skip = 5;

    while (*s) {
	if (*s == '\n') {
	    skip++;
	    break;
	} else if (isspace(*s)) {
	    skip++;
	}
	s++;
    } 

    return skip;
}

static void
insert_text_with_markup (GtkTextBuffer *tbuf, GtkTextIter *iter,
			 const char *s)
{
    char targ[128];
    const char *p;
    int ins;

    const char *indent = NULL;
    const char *code = NULL;

    while ((p = strstr(s, "<"))) {
	int skip = 0;

	if (code != NULL) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "code", indent, NULL);
	} else {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "text", indent, NULL);
	}

	p++;

	if (*p == '@') {
	    /* "atomic" markup */
	    ins = get_instruction_and_string(p + 1, targ);
	    if (ins == INSERT_REF) {
		insert_link(tbuf, iter, targ, gretl_command_number(targ), indent);
	    } else if (ins == INSERT_PDFLINK) {
		insert_link(tbuf, iter, targ, GUIDE_PAGE, indent);
	    } else if (ins == INSERT_FIG) {
		insert_help_figure(tbuf, iter, targ);
	    } else if (ins != INSERT_NONE) {
		insert_tagged_text(tbuf, iter, targ, ins, indent);
	    }
	    skip = 8 + strlen(targ);
	} else if (!strncmp(p, "indent", 6)) {
	    indent = "indented";
	    skip = 7 + (*(p+7) == '\n');
	} else if (!strncmp(p, "/indent", 7)) {
	    indent = NULL;
	    skip = 8 + (*(p+8) == '\n');
	} else if (!strncmp(p, "code", 4)) {
	    code = "code";
	    skip = get_code_skip(p + 5);
	} else if (!strncmp(p, "/code", 5)) {
	    code = NULL;
	    skip = 6 + (*(p+6) == '\n');
	} else {
	    /* literal "<" */
	    gtk_text_buffer_insert(tbuf, iter, "<", 1);
	}

	s = p + skip;
    }

    if (code != NULL) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "code", indent, NULL);
    } else {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "text", indent, NULL);
    }
}

static char *grab_topic_buffer (const char *s)
{
    const char *p = strstr(s, "\n# ");
    char *buf;

    if (p != NULL) {
	buf = g_strndup(s, p - s);
    } else {
	buf = g_strdup(s);
    }

    return buf;
}

/* pull the appropriate chunk of help text out of the buffer attached
   to the help viewer and display it */

void set_help_topic_buffer (windata_t *hwin, int hcode, int pos, int en)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    char line[256];
    gchar *hbuf;
    char *topicbuf;

    tbuf = gretl_text_buf_new();

    if (pos == 0) {
	/* cli help with no topic selected */
	cmdref_title_page(hwin, tbuf, en);
	cursor_to_top(hwin);
	hwin->active_var = 0;
	return;
    }
    
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    hbuf = (gchar *) hwin->data + pos;
    bufgets_init(hbuf);

    if (bufgets(line, sizeof line, hbuf) == NULL) {
	return;
    }

    if (gui_help(hwin->role)) {
	/* topic heading: descriptive string */
	gchar *p = quoted_help_string(line);

	gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
						 p, -1,
						 "sansbold", NULL);
	free(p);
    } else {
	/* topic heading: plain command word */
	char hword[9];

	sscanf(line + 2, "%8s", hword);
	gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
						 hword, -1,
						 "redtext", NULL);
    }

    gtk_text_buffer_insert(tbuf, &iter, "\n", 1);

    topicbuf = grab_topic_buffer(hbuf + strlen(line) + 1);
    if (topicbuf == NULL) {
	return;
    }

    insert_text_with_markup(tbuf, &iter, topicbuf);
    free(topicbuf);

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->w), tbuf);
    g_object_set_data(G_OBJECT(hwin->w), "backpage", 
		      GINT_TO_POINTER(hwin->active_var));
    maybe_connect_help_signals(hwin, en);
    cursor_to_top(hwin);
    hwin->active_var = hcode;
}

#ifndef USE_GTKSOURCEVIEW

void correct_line_color (windata_t *vwin)
{
    GtkTextBuffer *buf;
    GtkTextIter start, end;
    gint linelen;
    gchar *txt;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gtk_text_buffer_get_iter_at_mark(buf, &end, 
				     gtk_text_buffer_get_insert(buf));
    linelen = gtk_text_iter_get_chars_in_line(&end);
    start = end;
    gtk_text_iter_backward_chars(&start, linelen);

    txt = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

    if (*txt == '#') {
	gtk_text_buffer_apply_tag_by_name (buf, "bluetext",
					   &start, &end);
    }
    g_free(txt);
}

#endif /* not USE_GTKSOURCEVIEW */

GtkWidget *create_text (GtkWidget *dlg, GtkTextBuffer **buf, 
			int hsize, int vsize, gboolean editable)
{
    GtkTextBuffer *tbuf = gretl_text_buf_new();
    GtkWidget *w = gtk_text_view_new_with_buffer(tbuf);

    *buf = tbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(w), 4);

    gtk_widget_modify_font(GTK_WIDGET(w), fixed_font);

    hsize *= get_char_width(w);
    hsize += 48;

    gtk_window_set_default_size(GTK_WINDOW(dlg), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), editable);

    return w;
}

void text_table_setup (GtkWidget *vbox, GtkWidget *w)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), w); 
    gtk_widget_show(w);
    gtk_widget_show(sw);
}

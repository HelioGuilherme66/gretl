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

#include "libgretl.h"
#include "texprint.h"
#include "libset.h"
#include "gretl_string_table.h"

#include <glib.h>

#ifdef ENABLE_NLS

static int numeric_c_locale_depth = 0;
static char *numeric_locale = NULL;
static int native_dot = -1;

/**
 * gretl_push_c_numeric_locale:
 *
 * Description: Saves the current %LC_NUMERIC locale and sets it to "C".
 * This way you can safely read/write floating point numbers all in the
 * same format, using '.' as the decimal character.  You should make sure 
 * that code between gretl_push_c_numeric_locale() and gretl_pop_c_numeric_locale()
 * doesn't do any setlocale calls, or locale may end up in a strange setting.
 * Also make sure to always pop the C numeric locale after you've pushed it.
 * The calls can be nested.
 **/

void gretl_push_c_numeric_locale (void)
{
    if (native_dot == -1) {
	struct lconv *lc = localeconv();

	native_dot = (*lc->decimal_point == '.');
    }

    if (native_dot == 1) {
	return;
    }
    
    if (numeric_c_locale_depth == 0) {
	free(numeric_locale);
	numeric_locale = gretl_strdup(setlocale(LC_NUMERIC, NULL));
	setlocale(LC_NUMERIC, "C");
    }
    
    numeric_c_locale_depth++;
}

/**
 * gretl_pop_c_numeric_locale:
 *
 * Description:  Restores the LC_NUMERIC locale to what it was 
 * before the matching gretl_push_c_numeric_locale(). If these calls
 * were nested, then this is a no-op until we get to the most outermost
 * layer. Code in between these should not do any setlocale calls
 * to change the %LC_NUMERIC locale or things may come out very strange.
 **/

void gretl_pop_c_numeric_locale (void)
{
    if (numeric_c_locale_depth == 0) {
	return;
    }

    numeric_c_locale_depth--;

    if (numeric_c_locale_depth == 0) {
	setlocale(LC_NUMERIC, numeric_locale);
	free(numeric_locale);
	numeric_locale = NULL;
    }
}

#else

void gretl_push_c_numeric_locale (void)
{
    return;
}

void gretl_pop_c_numeric_locale (void)
{
    return;
}

#endif /* ENABLE_NLS or not */

/**
 * doing_nls:
 *
 * Returns: 1 if NLS translation is in effect, 0 otherwise.
 */

int doing_nls (void)
{
    static int called, nls;

    if (!called) {
	nls = (strcmp("_Open data", _("_Open data")) ||
	       strcmp("Test statistic", _("Test statistic")) ||
	       strcmp("annual", _("annual")));
	called = 1;
    }

    return nls;
}

#ifdef ENABLE_NLS

static int decpoint;

/**
 * reset_local_decpoint:
 *
 * Uses localeconv() to determine the representation of the decimal
 * point in the current locale.
 *
 * Returns: the decimal character for the current locale.
 */

int reset_local_decpoint (void)
{
    struct lconv *lc;

    lc = localeconv();
    decpoint = *lc->decimal_point;
    set_atof_point(decpoint);
    return decpoint;
}

/**
 * get_local_decpoint:
 *
 * Returns: the decimal character for the current locale.
 */

int get_local_decpoint (void)
{
    if (decpoint == 0) {
	decpoint = reset_local_decpoint();
    }
    return decpoint;
}

#else

int reset_local_decpoint (void)
{
    return '.';
}

int get_local_decpoint (void)
{
    return '.';
}

#endif /* ENABLE_NLS or not */

/* fudges for strings that should not be in utf-8 under some
   conditions: under GTK translations always come out in utf-8 in the
   GUI, but when we're sending stuff to stderr we may have to put it
   into ISO-8859-N or some MS Code Page.
*/

static int gretl_cset_maj;
static int gretl_cset_min;
#ifdef WIN32
static int gretl_cpage;
#endif
static int native_utf8;

void set_native_utf8 (int s)
{
    native_utf8 = s;
    
    if (native_utf8) {
	set_stdio_use_utf8();
    }
}


/* Use g_get_charset() to determine the current local character set,
   and record this information.  If we get an "ISO-XXXX-Y" locale,
   record the numerical elements as gretl_cset_maj and gretl_cset_min
   respectively.  If we get a Windows "CPXXXX" reading, record the
   codepage as gretl_cpage.
*/

void set_gretl_charset (void)
{
    const char *charset = NULL;
    char gretl_charset[32];

#ifdef MAC_NATIVE
    /* FIXME - why is this necessary? */
    native_utf8 = 1;
#else
    native_utf8 = g_get_charset(&charset);
#endif

    if (native_utf8) {
	set_stdio_use_utf8();
    }

    *gretl_charset = '\0';

    if (!native_utf8 && charset != NULL && *charset != '\0') {
	char *p;

	strncat(gretl_charset, charset, 31);
	gretl_lower(gretl_charset);
	p = strstr(gretl_charset, "iso");
	if (p != NULL) {
	    char numstr[6] = {0};

	    while (*p && !isdigit((unsigned char) *p)) p++;
	    strncat(numstr, p, 4);
	    gretl_cset_maj = atoi(numstr);
	    if (strlen(p) > 4) {
		p += 4;
		while (*p && !isdigit((unsigned char) *p)) p++;
		gretl_cset_min = atoi(p);
	    }
	    
	    if (gretl_cset_maj < 0 || gretl_cset_maj > 9000) {
		gretl_cset_maj = gretl_cset_min = 0;
	    } else if (gretl_cset_min < 0 || gretl_cset_min > 30) {
		gretl_cset_maj = gretl_cset_min = 0;
	    }
	} 
#ifdef WIN32
	if (p == NULL) {
	    sscanf(gretl_charset, "cp%d", &gretl_cpage);
	}
#endif
    }

#ifdef WIN32
    fprintf(stderr, "codepage = %d\n", gretl_cpage);
    if (gretl_cpage != 1250) {
	char *e = getenv("GRETL_CPAGE");

	if (e != NULL && !strcmp(e, "1250")) {
	    gretl_cpage = 1250;
	    fprintf(stderr, "revised codepage to 1250\n");
	}
    }
#endif
}

#ifdef WIN32

static void set_cp_from_locale (const char *loc)
{
    const char *p = strrchr(loc, '.');

    if (p != NULL && strlen(p) > 3 && isdigit(p[1])) {
        gretl_cpage = atoi(p + 1);
        fprintf(stderr, "set_cp_from_locale: CP = %d\n", gretl_cpage);
    }
}

#endif /* WIN32 */

#ifdef ENABLE_NLS

static const char *get_gretl_charset (void)
{
    static char cset[12];

# ifdef WIN32
    if (gretl_cpage > 0) {
	sprintf(cset, "CP%d", gretl_cpage);
	return cset;
    }
# endif

    if (gretl_cset_maj > 0 && gretl_cset_min > 0) {
	sprintf(cset, "ISO-%d-%d", gretl_cset_maj, gretl_cset_min);
	return cset;
    } 

    return NULL;
}

#endif /* ENABLE_NLS */

int iso_latin_version (void)
{
    char *lang = NULL;

    if (gretl_cset_maj == 8859 &&
	(gretl_cset_min == 1 || 
	 gretl_cset_min == 2 ||
	 gretl_cset_min == 5 ||
	 gretl_cset_min == 9 ||
	 gretl_cset_min == 15)) {
	return gretl_cset_min;
    }

#ifdef WIN32
    if (gretl_cpage == 1252) {
	return 1;
    } else if (gretl_cpage == 1250) {
	return 2;
    } else if (gretl_cpage == 1251) {
	return 5;
    } else if (gretl_cpage == 1254) {
	return 9;
    }
#endif

    /* Polish, Russian, Turkish, Ukrainian: UTF-8 locale? */
    lang = getenv("LANG");
    if (lang != NULL) {
	if (!strncmp(lang, "pl", 2)) {
	    return 2;
	} else if (!strncmp(lang, "ru", 2)) {
	    return 5;
	} else if (!strncmp(lang, "tr", 2)) {
	    return 9;
	} else if (!strncmp(lang, "uk", 2)) {
	    return 5;
	}
    } 

    return 1;
}

int chinese_locale (void)
{
#ifdef WIN32
    return (gretl_cpage == 950);
#elif defined(ENABLE_NLS)
    char *loc = setlocale(LC_ALL, NULL);

    return (loc != NULL && !strncmp(loc, "zh", 2));
#else
    return 0;
#endif
}

int japanese_locale (void)
{
#ifdef WIN32
    return (gretl_cpage == 932);
#elif defined(ENABLE_NLS)
    char *loc = setlocale(LC_ALL, NULL);

    return (loc != NULL && !strncmp(loc, "ja", 2));
#else
    return 0;
#endif
}

int east_asian_locale (void)
{
#ifdef WIN32
    return (gretl_cpage == 950 || gretl_cpage == 932);
#elif defined(ENABLE_NLS)
    char *loc = setlocale(LC_ALL, NULL);

    return (loc != NULL && (!strncmp(loc, "zh", 2) ||
			    !strncmp(loc, "ja", 2)));
#else
    return 0;
#endif    
}

#ifdef ENABLE_NLS

char *iso_gettext (const char *msgid)
{
    static int iso_switch = -1;
    static const char *cset;
    static int cli;
    char *ret;

    /* command line program: switching of codesets should not be 
       required, since unlike the GUI program there's no need
       to force UTF-8 as the default.  
    */

    if (!strcmp(msgid, "@CLI_INIT")) {
	cli = 1;
	return NULL;
    }

    if (cli) { 
	return gettext(msgid);
    }

    /* iso_switch: we'll reckon that if the system character set is
       not UTF-8, and is an ISO-8859-N or Windows CP12NN 8-bit set,
       then we should probably recode when printing translated strings
       in the context of writing CSV files.  

       If iso_switch is non-zero (once it's determinate) this makes the
       I_() gettext macro use the system encoding, otherwise I_() is
       equivalent to plain _(), which always spits out UTF-8.
    */

    if (iso_switch < 0) {
	/* not yet determined */
	cset = get_gretl_charset();
	if (cset == NULL) {
	    fprintf(stderr, "get_gretl_charset: using UTF-8\n");
	} else {
	    fprintf(stderr, "get_gretl_charset gave %s\n", cset);
	}
	iso_switch = (cset != NULL);
    }

    if (iso_switch) {
	bind_textdomain_codeset(PACKAGE, cset);
	ret = gettext(msgid);
	bind_textdomain_codeset(PACKAGE, "UTF-8");
    } else {
	ret = gettext(msgid);
    }

    return ret;
} 

/* return translated @msgid in UTF-8, when this is not
   the default operation for plain gettext()
*/

static char *force_utf8_gettext (const char *msgid)
{
    static const char *cset;
    char *ret;

    if (cset == NULL) {
	cset = get_gretl_charset();
    }

    bind_textdomain_codeset(PACKAGE, "UTF-8");
    ret = gettext(msgid);
    bind_textdomain_codeset(PACKAGE, cset);

    return ret;
} 

/* Cases to consider below:

   (1) The system-native encoding is UTF-8: in that case all output
   should be in UTF-8, with the possible exception of RTF; and all
   output _will_ be in UTF-8 if we just call plain gettext().

   (2) The system-native encoding is not UTF-8 (this means
   MS Windows -- or possibly ancient Linux?). Then we have two
   sub-cases.

   (2a) We're in the GUI program: plain gettext() will have been
   coerced to produce UTF-8 for the sake of GTK. But we should
   probably arrange for RTF and CSV output to be in the locale
   encoding, and also any output going to the console.

   (2b) We're at the command line (gretlcli): gettext() output will
   be in the system locale, but under the "new scheme" (December
   2011) TeX output should be coerced to UTF-8. All other output
   should use plain gettext().
*/

enum {
    GETTEXT_DEFAULT,
    GETTEXT_FORCE_UTF8,
    GETTEXT_FORCE_LOCALE
};

static int gettext_mode;

void set_alt_gettext_mode (PRN *prn)
{
    gettext_mode = GETTEXT_DEFAULT;

    /* As of 2014-09-07, we'll handle RTF by (a) getting
       gettext to write UTF-8 (forcing it if necessary)
       then (b) using utf8_to_rtf() to convert to ASCII
       plus \uXXXX codes, as per the spec for RTF >= 1.5.
       (It would be nice if gettext were able to generate
       the latter directly, but it can't.)

       Note: this means that RTF should not be written to
       file directly: a bufferized PRN should be used
       first so that the buffer can be sent through
       utf8_to_rtf(); then it can be written to file.
    */

    if (prn != NULL && !native_utf8) {
	if (gretl_in_gui_mode()) {
	    if (printing_to_standard_stream(prn)) {
		gettext_mode = GETTEXT_FORCE_LOCALE;
	    }
	} else if (tex_format(prn) || rtf_format(prn)) {
	    /* CLI mode, writing TeX or RTF */
	    gettext_mode = GETTEXT_FORCE_UTF8;
	}
    }
}

char *alt_gettext (const char *msgid)
{
    if (gettext_mode == GETTEXT_DEFAULT) {
	return gettext(msgid);
    } else if (gettext_mode == GETTEXT_FORCE_UTF8) {
	return force_utf8_gettext(msgid);
    } else if (gettext_mode == GETTEXT_FORCE_LOCALE) {
	return iso_gettext(msgid);
    } else {
	return gettext(msgid);
    }
}

#else /* !ENABLE_NLS */

void set_alt_gettext_mode (PRN *prn)
{
    return;
}

#endif /* ENABLE_NLS */

#ifdef WIN32

struct localeinfo {
    int id;
    const char *code;
};

/* the following are the strings that are accepted by setlocale()
   on win32 */

static struct localeinfo locales[] = {
    { LANG_AUTO,  NULL },
    { LANG_C,     "english" },
    { LANG_SQ,    "albanian" },
    { LANG_EU,    "basque" },
    { LANG_BG,    "bulgarian" },    
    { LANG_CA,    "catalan" },
    { LANG_ZH_TW, "chinese-traditional" },
    { LANG_CS,    "czech" },
    { LANG_FR,    "french" },
    { LANG_GL,    "galician" },
    { LANG_DE,    "german" },
    { LANG_EL,    "greek" },
    { LANG_IT,    "italian" },
    { LANG_JA,    "japanese" },
    { LANG_PL,    "polish" },
    { LANG_PT,    "portuguese" },
    { LANG_PT_BR, "portuguese-brazilian" },
    { LANG_RO,    "romanian" },
    { LANG_RU,    "russian" },
    { LANG_ES,    "spanish" },
    { LANG_TR,    "turkish" },
    { LANG_UK,    "ukrainian" },
    { LANG_MAX,    NULL }
};

const char *locale_code_from_id (int langid)
{
    int i;

    for (i=0; i<LANG_MAX; i++) {
	if (langid == locales[i].id) {
	    return locales[i].code;
	}
    }

    return NULL;
}

#endif /* WIN32 */

struct langinfo {
    int id;
    const char *name;
    const char *code;
};

static struct langinfo langs[] = {
    { LANG_AUTO,  "Automatic",            NULL    },
    { LANG_C,     "English",              "C"     },
    { LANG_SQ,    "Albanian",             "sq_AL" },
    { LANG_EU,    "Basque",               "eu_ES" },
    { LANG_BG,    "Bulgarian",            "bg_BG" },    
    { LANG_CA,    "Catalan",              "ca_ES" },
    { LANG_ZH_TW, "Chinese (Taiwan)",     "zh_TW" },
    { LANG_CS,    "Czech",                "cs_CZ" },
    { LANG_FR,    "French",               "fr_FR" },
    { LANG_GL,    "Galician",             "gl_ES" },
    { LANG_DE,    "German",               "de_DE" },
    { LANG_EL,    "Greek",                "el_GR" },
    { LANG_IT,    "Italian",              "it_IT" },
    { LANG_JA,    "Japanese",             "ja_JP" },
    { LANG_PL,    "Polish",               "pl_PL" },
    { LANG_PT,    "Portuguese",           "pt_PT" },
    { LANG_PT_BR, "Portuguese (Brazil)",  "pt_BR" },
    { LANG_RO,    "Romanian",             "ro_RO" },
    { LANG_RU,    "Russian",              "ru_RU" },
    { LANG_ES,    "Spanish",              "es_ES" },
    { LANG_TR,    "Turkish",              "tr_TR" },
    { LANG_UK,    "Ukrainian",            "uk_UA" },
    { LANG_MAX,    NULL,                   NULL   }
};

const char *lang_string_from_id (int langid)
{
    int i;

    for (i=0; i<LANG_MAX; i++) {
	if (langid == langs[i].id) {
	    return langs[i].name;
	}
    }

    return NULL;
}

int lang_id_from_name (const char *s)
{
    if (s != NULL || *s != '\0') {
	int i;

	for (i=0; i<LANG_MAX; i++) {
	    if (!strcmp(s, langs[i].name)) {
		return langs[i].id;
	    }
	}
    }

    return 0;
}

const char *lang_code_from_id (int langid)
{
    int i;

    for (i=0; i<LANG_MAX; i++) {
	if (langid == langs[i].id) {
	    return langs[i].code;
	}
    }

    return NULL;
}

#ifdef WIN32

static const char *lang_code_from_windows_name (char *s)
{
    const char *ret = NULL;
    int i, n;

    if (strstr(s, "razil")) {
	return "pr_BR";
    } else if (!strncmp(s + 1, "hinese", 6)) {
	return "zh_TW";
    }

    for (i=1; i<LANG_MAX; i++) {
	n = strlen(langs[i].name);
	if (!strncmp(s, langs[i].name, n)) {
	    ret = langs[i].code;
	    break;
	}
    }

    if (ret == NULL && strchr(s, '-')) {
	/* Windows : "en-US" in place of "en_US"? */
	char *p = strchr(s, '-');

	*p = '_';
	for (i=1; i<LANG_MAX; i++) {
	    n = strlen(langs[i].code);
	    if (!strncmp(s, langs[i].code, n)) {
		ret = langs[i].code;
		break;
	    }
	}
    }

    return ret;
}

static char *win32_set_numeric (const char *lang)
{
    char *set = NULL;
    int i;

    for (i=LANG_SQ; i<LANG_MAX; i++) {
	if (!strcmp(lang, langs[i].code) ||
	    !strncmp(lang, langs[i].code, 2)) {
	    set = setlocale(LC_NUMERIC, locales[i].code);
	    if (set != NULL) {
		break;
	    }
	}
    }

    return set;
}

#else /* !WIN32 */

# ifdef ENABLE_NLS

static char *other_set_numeric (const char *lang)
{
    char *set = setlocale(LC_NUMERIC, lang);

    if (set == NULL) {
	char lfix[32];

	sprintf(lfix, "%s.UTF-8", lang);
	set = setlocale(LC_NUMERIC, lfix);
    }

    return set;
}

# endif /* ENABLE_NLS */

#endif /* WIN32 or not */

void set_lcnumeric (int langid, int lcnumeric)
{
#ifndef ENABLE_NLS
    return;
#else
    if (lcnumeric && langid != LANG_C) {
	const char *lang;
	char *set = NULL;

	if (langid == LANG_AUTO) {
	    /* respect the system LANG setting */
	    lang = getenv("LANG");
	} else {
	    /* fake it from user preference */
	    lang = lang_code_from_id(langid);
	}

	if (lang != NULL) {
# ifdef WIN32
	    set = win32_set_numeric(lang);
# else
	    set = other_set_numeric(lang);
# endif
	} 
	if (set == NULL) {
	    setlocale(LC_NUMERIC, "");
	    gretl_setenv("LC_NUMERIC", "");
	}
    } else {
	/* either lcnumeric is not chosen, or we're in LANG_C */
	setlocale(LC_NUMERIC, "C");
	gretl_setenv("LC_NUMERIC", "C");
    }

    reset_local_decpoint();
#endif
}

#ifdef ENABLE_NLS

static int 
set_locale_with_workaround (int langid, const char *lcode,
			    char **locp)
{
    char *test = setlocale(LC_ALL, lcode);

# ifndef WIN32    
    if (test == NULL) {
	char lfix[32];

	sprintf(lfix, "%s.UTF-8", lcode);
	test = setlocale(LC_ALL, lfix);
    }
# endif

    if (test != NULL) {
	fprintf(stderr, "setlocale: '%s' -> '%s'\n", lcode, test);
	if (strcmp("_File", _("_File")) == 0) {
	    const char *langstr;
	    char tmp[64];
	    
	    langstr = lang_string_from_id(langid);
	    sscanf(langstr, "%s", tmp);
	    gretl_lower(tmp);
	    gretl_setenv("LANGUAGE", tmp);
	}
    }

    if (locp != NULL && test != NULL) {
	*locp = gretl_strdup(test);
    }

    return test == NULL;
}

# ifdef WIN32
# define get_setlocale_string(i) (locale_code_from_id(i))
# else
# define get_setlocale_string(i) (lang_code_from_id(i))
# endif

#endif /* ENABLE_NLS */

/* @langstr should be the English name of the selected language
   as displayed in the GUI (e.g. "German", "French") 
*/

int test_locale (const char *langstr)
{
#ifndef ENABLE_NLS
    return 1;
#else
    const char *lcode;
    char *orig, ocpy[64];
    int langid, err = 0;

    langid = lang_id_from_name(langstr);
    lcode = get_setlocale_string(langid);
    orig = setlocale(LC_ALL, NULL);

    gretl_error_clear();

    *ocpy = '\0';
    strncat(ocpy, orig, 63);

    err = set_locale_with_workaround(langid, lcode, NULL);

    if (err) {
	gretl_errmsg_sprintf(_("%s: locale is not supported "
			       "on this system"), lcode);
    } else {
	setlocale(LC_ALL, ocpy); /* restore the original locale */
    } 

    return err;
#endif
}

#ifdef ENABLE_NLS

static void record_locale (char *locale)
{
    int done = 0;
    
# ifdef WIN32
    /* LANG probably not present, use setlocale output */
    if (locale != NULL) {
	const char *s = lang_code_from_windows_name(locale);

	if (s != NULL) {
	    gretl_insert_builtin_string("lang", s);
	    done = 1;
	}
    }
# else
    char *lang = getenv("LANG");

    if (lang != NULL) {
	/* prefer using LANG */
	if (strrchr(lang, '.') == NULL) {
	    gretl_insert_builtin_string("lang", lang);
	} else {
	    char *tmp = gretl_strdup(lang);
	    char *p = strrchr(tmp, '.');

	    *p = '\0';
	    gretl_insert_builtin_string("lang", tmp);
	    free(tmp);
	}
	done = 1;
    } else if (locale != NULL) {
	/* use locale as fallback */
	if (strrchr(locale, '.') == NULL) {
	    gretl_insert_builtin_string("lang", locale);
	} else {
	    char *p = strrchr(locale, '.');

	    *p = '\0';
	}
	gretl_insert_builtin_string("lang", locale);
	done = 1;
    }
# endif

    if (!done) {
	gretl_insert_builtin_string("lang", "unknown");
    }
}

#endif

int force_language (int langid)
{
#ifndef ENABLE_NLS
    return 1;
#else
    const char *lcode = NULL;
    char *locale = NULL;
    int err = 0;

    if (langid == LANG_AUTO) {
	/* note: avoid getting long spew from Windows */
	locale = gretl_strdup(setlocale(LC_COLLATE, NULL));
	goto record;
    }

    if (langid == LANG_C) {
	gretl_setenv("LANGUAGE", "english");
	gretl_setenv("LANG", "C");
	setlocale(LC_ALL, "C");
    } else {
	lcode = get_setlocale_string(langid);
	if (lcode != NULL) {  
# ifdef WIN32
	    locale = gretl_strdup(setlocale(LC_ALL, lcode));
            fprintf(stderr, "lcode='%s', newloc='%s'\n", lcode, locale);
	    if (locale != NULL) {
		set_cp_from_locale(locale);
	    } else {
		err = 1;
	    }
# else
	    err = set_locale_with_workaround(langid, lcode, &locale);
# endif
	}
    }

# if defined(WIN32)
    if (langid == LANG_C) {
	gretl_setenv("LC_ALL", "C");
	textdomain("none");
    } else if (lcode != NULL) {
        lcode = lang_code_from_id(langid);
	if (lcode != NULL) {
	    gretl_setenv("LC_ALL", lcode);
	    gretl_setenv("LANG", lcode);
	}
    }
# elif defined(OS_OSX)
    if (langid != LANG_C) {
	lcode = lang_code_from_id(langid);
	if (lcode != NULL) {
	    gretl_setenv("LANGUAGE", lcode);
	}
    }
# endif

 record:

    record_locale(locale);
    free(locale);

    return err;
#endif /* ENABLE_NLS */
}

static void 
iso_to_ascii_translate (char *targ, const char *src, int latin)
{
    char *p;
    const char *q;

    p = targ;
    q = src;

    if (latin == 1) {
	while (*q) {
	    unsigned char c = *q;

	    if (c == '\t' || c == '\n' || (c >= 32 && c <= 126)) {
		*p++ = c;
	    } else if (c >= 192 && c <= 198) {
		*p++ = 'A';
	    } else if (c == 199) {
		*p++ = 'C';
	    } else if (c >= 200 && c <= 203) {
		*p++ = 'E';
	    } else if (c >= 204 && c <= 207) {
		*p++ = 'I';
	    } else if (c == 208) {
		*p++ = 'D';
	    } else if (c == 209) {
		*p++ = 'N';
	    } else if (c >= 210 && c <= 214) {
		*p++ = 'O';
	    } else if (c == 216) {
		*p++ = 'O';
	    } else if (c >= 217 && c <= 220) {
		*p++ = 'U';
	    } else if (c == 221) {
		*p++ = 'Y';
	    } else if (c >= 224 && c <= 230) {
		*p++ = 'a';
	    } else if (c == 231) {
		*p++ = 'c';
	    } else if (c >= 232 && c <= 235) {
		*p++ = 'e';
	    } else if (c >= 236 && c <= 239) {
		*p++ = 'i';
	    } else if (c == 240) {
		*p++ = 'd';
	    } else if (c == 241) {
		*p++ = 'n';
	    } else if (c >= 242 && c <= 246) {
		*p++ = 'o';
	    } else if (c == 248) {
		*p++ = 'o';
	    } else if (c >= 249 && c <= 252) {
		*p++ = 'u';
	    } else if (c == 253) {
		*p++ = 'y';
	    }
	    q++;
	}
    } else if (latin == 2) {
	while (*q) {
	    unsigned char c = *q;

	    if (c == '\t' || c == '\n' || (c >= 32 && c <= 126)) {
		*p++ = c;
	    }

#ifndef WIN32
	    if (c==161 || c==193 || c==194 || c==195 || c==196) {
		*p++ = 'A';
	    }
#else
	    if (c==165 || c==193 || c==194 || c==195 || c==196) {
		*p++ = 'A';
	    }
#endif
	    else if (c==198 || c==199 || c==200) {
		*p++ = 'C';
	    }
	    else if (c==207 || c==208) {
		*p++ = 'D';
	    }
	    else if (c==201 || c==202 || c==203 || c==204) {
		*p++ = 'E';
	    }
	    else if (c==205 || c==206) {
		*p++ = 'I';
	    }
#ifndef WIN32
	    else if (c==163 || c==165 || c==197) {
		*p++ = 'L';
	    }
#else
	    else if (c==163 || c==188 || c==197) {
		*p++ = 'L';
	    }
#endif
	    else if (c==209 || c==210) {
		*p++ = 'N';
	    }
	    else if (c==211 || c==212 || c==213 || c==214) {
		*p++ = 'O';
	    }
	    else if (c==192 || c==216) {
		*p++ = 'R';
	    }
#ifndef WIN32
	    else if (c==166 || c==169 || c==170) {
		*p++ = 'S';
	    }
#else
	    else if (c==138 || c==140 || c==170) {
		*p++ = 'S';
	    }
#endif
#ifndef WIN32
	    else if (c==171 || c==222) {
		*p++ = 'T';
	    }
#else
	    else if (c==141 || c==222) {
		*p++ = 'T';
	    }
#endif
	    else if (c==217 || c==218 || c==219 || c==220) {
		*p++ = 'U';
	    }
	    else if (c==221) {
		*p++ = 'Y';
	    }
#ifndef WIN32
	    else if (c==172 || c==174 || c==175) {
		*p++ = 'Z';
	    }
#else
	    else if (c==142 || c==143 || c==175) {
		*p++ = 'Z';
	    }
#endif
#ifndef WIN32
	    else if (c==177 || c==225 || c==226 || c==227 || c==228) {
		*p++ = 'a';
	    }
#else
	    else if (c==185 || c==225 || c==226 || c==227 || c==228) {
		*p++ = 'a';
	    }
#endif
	    else if (c==230 || c==231 || c==232) {
		*p++ = 'c';
	    }
	    else if (c==239 || c==240) {
		*p++ = 'd';
	    }
	    else if (c==233 || c==234 || c==235 || c==236) {
		*p++ = 'e';
	    }
	    else if (c==237 || c==238) {
		*p++ = 'i';
	    }
#ifndef WIN32
	    else if (c==179 || c==181 || c==229) {
		*p++ = 'l';
	    }
#else
	    else if (c==179 || c==190 || c==229) {
		*p++ = 'l';
	    }
#endif
	    else if (c==241 || c==242) {
		*p++ = 'n';
	    }
	    else if (c==243 || c==244 || c==245 || c==246) {
		*p++ = 'o';
	    }
	    else if (c==224 || c==248) {
		*p++ = 'r';
	    }
#ifndef WIN32
	    else if (c==182 || c==185 || c==186 || c==223) {
		*p++ = 's';
	    }
#else
	    else if (c==154 || c==156 || c==186 || c==223) {
		*p++ = 's';
	    }
#endif
#ifndef WIN32
	    else if (c==187 || c==254) {
		*p++ = 't';
	    }
#else
	    else if (c==157 || c==254) {
		*p++ = 't';
	    }
#endif
	    else if (c==249 || c==250 || c==251 || c==252) {
		*p++ = 'u';
	    }
	    else if (c==253) {
		*p++ = 'y';
	    }
#ifndef WIN32
	    else if (c==188 || c==190 || c==191) {
		*p++ = 'z';
	    }
#else
	    else if (c==158 || c==159 || c==191) {
		*p++ = 'z';
	    }
#endif
	    q++;
	}
    }

    *p = '\0';
}

/* If @maxlen > 0 we limit the write to @targ to at most
   @maxlen bytes (excluding the terminating nul byte).
   If @sub > 0 we write this byte to @targ in place of
   UTF-8 characters that we can't represent in ASCII,
   provided they are lower than 0x0180.
*/

char *u8_to_ascii_convert (char *targ, const char *src,
			   int maxlen, char sub)
{
    int prevspace = 0;
    const char *q = src;
    char *p = targ;
    gunichar u;
    int c, skip;
    int len = 0;

    *p = '\0';

    /* If sub == 0 we assume we're doing varnames and
       so we skip all characters that are not valid in
       a gretl varname. But if sub > 0 we pass through
       all printable ASCII characters.
    */

    while (q && *q) {
	skip = 0;
	c = *q;
	if (sub > 0 && ((c >= 32 && c <= 126) || c == 9 || c == 10)) {
	    /* ASCII printables */
	    *p++ = c;
	    q++;
	} else if (c >= 0x0030 && c <= 0x0039) {
	    /* digits 0-9 */
	    *p++ = c;
	    q++;
	} else if (c >= 0x0041 && c <= 0x005A) {
	    /* upper-case ASCII letters */
	    *p++ = c;
	    q++;
	} else if (c >= 0x0061 && c <= 0x007A) {
	    /* lower-case ASCII letters */
	    *p++ = c;
	    q++;
	} else if (c == 0x005F) {
	    /* underscore */
	    *p++ = c;
	    q++;
	} else if (c == 0x0020) {
	    if (!prevspace) {
		prevspace = 1;
		*p++ = '_';
	    } else {
		skip = 1;
	    }
	    q++;
	} else {
	    /* handle Latin-1 and Latin-2, only */
	    u = g_utf8_get_char(q);
	    if (u >= 0x0180) {
		skip = 1; /* can't handle */
	    } else if ((u >= 0x00C0 && u <= 0x00C6) || u == 0x0102 || u == 0x0104) {
		*p++ = 'A';
	    } else if (u == 0x00C7 || u == 0x0106 || u == 0x010C) {
		*p++ = 'C';
	    } else if ((u >= 0x00C8 && u <= 0x00CB) || u == 0x0118 || u == 0x011A) {
		*p++ = 'E';
	    } else if (u >= 0x00CC && u <= 0x00CF) {
		*p++ = 'I';
	    } else if (u == 0x00D0 || u == 0x010E || u == 0x0110 || u == 0x010E) {
		*p++ = 'D';
	    } else if (u == 0x00D1 || u == 0x0143 || u == 0x0147) {
		*p++ = 'N';
	    } else if (u == 0x00D8 || (u >= 0x00D2 && u <= 0x00D6) || u == 0x0150) {
		*p++ = 'O';
	    } else if ((u >= 0x00D9 && u <= 0x00DC) || u == 0x016E || u == 0x0170) {
		*p++ = 'U';
	    } else if (u == 0x00DD) {
		*p++ = 'Y';
	    } else if (u == 0x00DE || u == 0x0164) {
		*p++ = 'T';
	    } else if (u == 0x00DF) {
		*p++ = 's';
	    } else if ((u >= 0x00E0 && u <= 0x00E6) || u == 0x0103) {
		*p++ = 'a';
	    } else if (u == 0x00E7 || u == 0x0107) {
		*p++ = 'c';
	    } else if ((u >= 0x00E8 && u <= 0x00EB) || u == 0x0119 || u == 0x011B) {
		*p++ = 'e';
	    } else if (u >= 0x00EC && u <= 0x00EF) {
		*p++ = 'i';
	    } else if (u == 0x00F0 || u == 0x0111 || u == 0x010F) {
		*p++ = 'd';
	    } else if (u == 0x00F1 || u == 0x0144 || u == 0x0148) {
		*p++ = 'n';
	    } else if (u == 0x00F8 || u == 0x0151 || (u >= 0x00F2 && u <= 0x00F6)) {
		*p++ = 'o';
	    } else if ((u >= 0x00F9 && u <= 0x00FC) || u == 0x016F || u == 0x0171) {
		*p++ = 'u';
	    } else if (u == 0x00FD || u == 0x00FF) {
		*p++ = 'y';
	    } else if (u == 0x00FE || u == 0x0163) {
		*p++ = 't';
	    } else if (u == 0x0141 || u == 0x013D || u == 0x0139) {
		*p++ = 'L';
	    } else if (u == 0x0142 || u == 0x013E || u == 0x013A) {
		*p++ = 'l';
	    } else if (u == 0x0154 || u == 0x0158) {
		*p++ = 'R';
	    } else if (u == 0x0155 || u == 0x0159) {
		*p++ = 'r';
	    } else if (u == 0x0160 || u == 0x015E) {
		*p++ = 'S';
	    } else if (u == 0x0161 || u == 0x015F) {
		*p = 's';
	    } else if (u == 0x0179 || u == 0x017D || u == 0x0178) {
		*p = 'Z';
	    } else if (u == 0x017A || u == 0x017E || u == 0x017C) {
		*p = 'z';
	    } else if (sub > 0) {
		*p = sub;
	    } else {
		skip = 1;
	    }
	    q = g_utf8_next_char(q);
	}
	if (c != 0x0020) {
	    prevspace = 0;
	}
	if (!skip) len++;
	if (maxlen > 0 && len == maxlen) {
	    break;
	}
    }

    *p = '\0';

    return targ;
}

static char *real_iso_to_ascii (char *s, int latin)
{
    char *tmp;

    tmp = malloc(strlen(s) + 1);
    if (tmp == NULL) {
	return NULL;
    }

    if (latin != 1 && latin != 2) {
	/* fallback?? */
	latin = 1;
    }

    iso_to_ascii_translate(tmp, s, latin);

    strcpy(s, tmp);
    free(tmp);

    return s;
}

char *iso_to_ascii (char *s) 
{
    return real_iso_to_ascii(s, 1);
}

char *sprint_l2_to_ascii (char *targ, const char *s, size_t len)
{
    iso_to_ascii_translate(targ, s, 2);

    return targ;
}

char *asciify_utf8_varname (char *s)
{
    char *tmp = malloc(32);

    if (tmp != NULL) {
	u8_to_ascii_convert(tmp, s, 31, 0);
	strcpy(s, tmp);
	free(tmp);
    }

    return s;
}

enum {
    ENC_ISO_LATIN,
    ENC_CODEPAGE
};

/* Construct a codeset string to pass to g_convert, to convert UTF-8
   encoded strings for gnuplot formats that can't handle UTF-8.  At
   this point we don't check if gnuplot can deal with the target
   codeset via its "set encoding" command, we're just trying to get
   the text localized: the graph may work even without an explicit
   "set encoding".
*/

static char *get_gp_encoding_set (char *s, int targ)
{
    int latin = iso_latin_version();

    if (targ == ENC_ISO_LATIN) {
	strcpy(s, "ISO-8859-");
	if (latin == 2) {
	    strcat(s, "2");
	} else if (latin == 5) {
	    strcat(s, "5");
	} else if (latin == 9) {
	    strcat(s, "9");
	} else if (latin == 15) {
	    strcat(s, "15");
	} else {
	    /* default is ISO-8859-1 */
	    strcat(s, "1");
	}
    } else if (chinese_locale()) {
	strcpy(s, "CP950");
    } else if (japanese_locale()) {
	strcpy(s, "CP932");
    } else {
	/* ENC_CODEPAGE */
	strcpy(s, "CP125");
	if (latin == 2) {
	    strcat(s, "0");
	} else if (latin == 5) {
	    strcat(s, "1");
	} else if (latin == 9) {
	    strcat(s, "4");
	} else {
	    /* default is CP1252 */
	    strcat(s, "2");
	} 
    }

    return s;
}

/* convert from UTF-8 to ISO latin, for, e.g., EPS graphs */

char *utf8_to_latin (const char *s)
{
    char to_set[12];
    gsize read, wrote;
    GError *err = NULL;
    char *ret = NULL;

    get_gp_encoding_set(to_set, ENC_ISO_LATIN);

    ret = g_convert(s, -1, to_set, "UTF-8",
		    &read, &wrote, &err);

    if (err != NULL) {
	gretl_errmsg_set(err->message);
	g_error_free(err);
    }

    return ret;
}

char *utf8_to_rtf (const char *s)
{
    const char *nextp, *p = s;
    short int k;
    PRN *prn;
    char *ret = NULL;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (prn == NULL) {
	return NULL;
    }
 
    while (*p) {
	nextp = g_utf8_next_char(p);
	if (nextp - p > 1) {
	    k = (short) g_utf8_get_char(p);
	    pprintf(prn, "\\u%d?", k);
	} else {
	    pputc(prn, *p);
	}
	p = nextp;
    }

    ret = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);

    return ret;
}

#define ascii_ctrl(a) (a == '\t' || a == '\n' || \
                       a == '\r' || a == CTRLZ)

int gretl_is_ascii (const char *buf)
{
    int a;

    while (*buf) {
	a = *buf;
	if (a > 126 || (a < 32 && !(ascii_ctrl(a)))) {
	    return 0;
	}
	buf++;
    }

    return 1;
}

/* We want to print @str in a field of @width (visible) characters,
   but @str may contain multi-byte characters. In that case, determine 
   the adjustment to @width that is needed to avoid underrun and
   return the adjusted value.
*/

int get_utf_width (const char *str, int width)
{
    /* the number of "invisible" bytes */
    int invis = strlen(str) - g_utf8_strlen(str, -1);

    return width + invis;
}

int get_translated_width (const char *str)
{
    int w = strlen(str);

    w += w - g_utf8_strlen(str, -1);

    return w;
}

/* utility functionality: recoding of an entire file:
   we start with a couple of static "helpers"
*/

static gchar *file_get_content (const char *fname,
				gsize *bytes,
				PRN *prn,
				int *err)
{
    GError *gerr = NULL;
    gchar *buf = NULL;
    int ok = 0;

#ifdef WIN32
    gchar *tmp = NULL;

    *err = maybe_recode_path(fname, &tmp, 1);
    if (!*err) {
	if (tmp != NULL) {
	    ok = g_file_get_contents(tmp, &buf, bytes, &gerr);
	    g_free(tmp);
	} else {
	    ok = g_file_get_contents(fname, &buf, bytes, &gerr);
	}
    }
#else
    ok = g_file_get_contents(fname, &buf, bytes, &gerr);
#endif

    if (ok) {
	pprintf(prn, "got content, %" G_GSIZE_FORMAT " bytes\n", *bytes);
    } else {
	*err = E_FOPEN;
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	}
    }

    return buf;
}

static int file_set_content (const char *fname, 
			     const gchar *buf,
			     gsize buflen)
{
    GError *gerr = NULL;
    int ok = 0;
    int err = 0;

#ifdef WIN32
    gchar *tmp = NULL;

    err = maybe_recode_path(fname, &tmp, 1);
    if (!err) {
	if (tmp != NULL) {
	    ok = g_file_set_contents(tmp, buf, buflen, &gerr);
	    g_free(tmp);
	} else {
	    ok = g_file_set_contents(fname, buf, buflen, &gerr);
	}
    }
#else
    ok = g_file_set_contents(fname, buf, buflen, &gerr);
#endif

    if (!ok) {
	err = E_FOPEN;
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	}
    }

    return err;
}

/**
 * gretl_recode_file:
 * @path1: path to original file.
 * @path2: path to file to be written.
 * @from_set: the codeset of the original file.
 * @to_set: the codeset for the recoded file.
 * @prn: gretl printer (for a few comments) or NULL.
 *
 * Returns: 0 on success or non-zero code on error.
 */

int gretl_recode_file (const char *path1, const char *path2,
		       const char *from_set, const char *to_set,
		       PRN *prn)
{
    gchar *buf = NULL;
    gsize bytes = 0;
    int err = 0;

    /* get entire content of original file */
    buf = file_get_content(path1, &bytes, prn, &err);

    if (!err) {
	GError *gerr = NULL;
	gchar *trbuf = NULL;
	gsize written = 0;

	/* recode the buffer */
	trbuf = g_convert(buf, bytes, to_set, from_set,
			  NULL, &written, &gerr);

	if (gerr != NULL) {
	    err = E_DATA;
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	} else {
	    /* write recoded text to file */
	    pprintf(prn, "recoded: %" G_GSIZE_FORMAT " bytes\n", written);
	    err = file_set_content(path2, trbuf, written);
	}

	g_free(trbuf);
    }

    g_free(buf);

    return err;
}


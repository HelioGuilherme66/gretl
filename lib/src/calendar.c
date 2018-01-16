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

/**
 * SECTION:calendar
 * @short_description: functions for working with dates
 * @title: Calendar
 * @include: libgretl.h
 *
 * Here we have various functions dealing with calendar dates;
 * for the most part these are designed for handling daily
 * time-series data.
 */

static int days_in_month[2][13] = {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
};

/* Jan 1, 0001, was a Monday on the proleptic Gregorian calendar */
#define DAY1 1

/* Note on the GLib API used below: where "julian" occurs in the names
   of GLib calendrical functions it refers to the "Julian day"; that is,
   the number of days since some fixed starting point, as used by
   astronomers. This is quite distinct from the Julian calendar.
   However, GLib takes Julian day 1 to be the first of January in
   AD 1, as opposed to the astronomical starting point in 4713 BC,
   so these are not strictly Julian days, and in our own functions
   which use the same concept we call them "epoch days".
*/

static int leap_year (int yr)
{
    return (!(yr % 4) && (yr % 100)) || !(yr % 400);
}

static int valid_ymd (int y, int m, int d, int julian)
{
    int ok = 1;
    
    if (!g_date_valid_dmy(d, m, y)) {
	ok = 0;
	if (julian && y > 0 && m == 2 && d == 29 && y%4 == 0) {
	    ok = 1;
	}
    }
	
    return ok;
}

static int day_of_week_from_ymd (int y, int m, int d, int julian)
{
    GDate date;
    int wd;

    if (!valid_ymd(y, m, d, julian)) {
	return -1;
    }

    g_date_clear(&date, 1);

    if (julian) {
	guint32 ed = epoch_day_from_julian_ymd(y, m, d);
	
	g_date_set_julian(&date, ed);
    } else {
	g_date_set_dmy(&date, d, m, y);
    }

    wd = g_date_get_weekday(&date);

    /* switch to Sunday == 0 */
    return wd == G_DATE_SUNDAY ? 0 : wd;
}

/**
 * epoch_day_from_ymd:
 * @y: year (1 <= y <= 9999).
 * @m: month (1 <= m <= 12).
 * @d: day of month (1 <= d <= 31).
 *
 * Returns: the epoch day number, which equals 1 for the first of
 * January in the year AD 1 on the proleptic Gregorian calendar,
 * or 0 on error.
 */

guint32 epoch_day_from_ymd (int y, int m, int d)
{
    GDate date;

    if (!g_date_valid_dmy(d, m, y)) {
	return 0;
    }

    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);

    return g_date_get_julian(&date);
}

/**
 * ymd_bits_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @y: location to receive year.
 * @m: location to receive month.
 * @m: location to receive day.
 *
 * Returns: 0 on success, non-zero on error.
 */

int ymd_bits_from_epoch_day (guint32 ed, int *y, int *m, int *d)
{
    GDate date;

    if (!g_date_valid_julian(ed)) {
	return E_INVARG;
    }

    g_date_clear(&date, 1);
    g_date_set_julian(&date, ed);

    *y = g_date_get_year(&date);
    *m = g_date_get_month(&date);
    *d = g_date_get_day(&date);

    return 0;
}

/**
 * julian_ymd_bits_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @y: location to receive year.
 * @m: location to receive month.
 * @m: location to receive day.
 *
 * Follows the algorithm of E.G. Richards (2013), "Calendars," In S.E.
 * Urban & P.K. Seidelmann, eds. Explanatory Supplement to the Astronomical
 * Almanac, 3rd ed. (pp. 585-624), Mill Valley, CA: University Science Books
 * (as set out on https://en.wikipedia.org/wiki/Julian_day).
 *
 * There are other algorithms for this purpose on the internet but they are
 * mostly wrong (at least, not right for all dates); many of them fail
 * the round-trip test (date -> epoch day -> date) for some dates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int julian_ymd_bits_from_epoch_day (guint32 ed, int *py,
				    int *pm, int *pd)
{
    int y = 4716;
    int p = 1461;
    int f = ed + 1721425 + 1401;
    int e = 4 * f + 3;
    int g = (e % p)/4;
    int h = 5 * g + 2;

    /* The addition of 1721425 above translates from our
       "epoch day" to Julian Day Number; the addition of 1401
       to the JDN is specified by Richards' algorithm.
    */

    *pd = (h % 153)/5 + 1;
    *pm = (h/153 + 2) % 12 + 1;
    *py = e/p - y + (14 - *pm)/12;

    return 0;
}

/**
 * epoch_day_from_julian_ymd:
 * @y: year (y >= 1).
 * @m: month (1 <= m <= 12).
 * @d: day of month (1 <= d <= 31).
 *
 * The @y, @m and @d arguments are assumed to refer to a date on
 * the Julian calendar. The conversion algorithm is taken from
 * https://en.wikipedia.org/wiki/Julian_day, where it appears to
 * be credited to the Department of Computer Science at UT, Austin.
 *
 * Returns: the epoch day number, which equals 1 for the first of
 * January in the year AD 1 on the proleptic Gregorian calendar,
 * or 0 on error.
 */

guint32 epoch_day_from_julian_ymd (int y, int m, int d)
{
    int a = (14 - m)/12;
    int jd;

    y = y + 4800 - a;
    m = m + 12*a - 3;

    jd = d + (153*m + 2)/5 + 365*y + y/4 - 32083;

    if (jd <= 1721425) {
	/* prior to AD 1 */
	return 0;
    } else {
	return (guint32) jd - 1721425;
    }
}

/**
 * ymd_extended_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @julian: non-zero to use Julian calendar, otherwise Gregorian.
 * @err: location to receive error code.
 *
 * Returns: a string on the pattern YYYY-MM-DD (ISO 8601 extended
 * date format) given the epoch day number, which equals 1 for the
 * first of January in the year 1 AD, or NULL on error.
 */

char *ymd_extended_from_epoch_day (guint32 ed, int julian, int *err)
{
    char *ret = NULL;
    int y, m, d;
    int myerr;

    if (julian) {
	myerr = julian_ymd_bits_from_epoch_day(ed, &y, &m, &d);
    } else {
	myerr = ymd_bits_from_epoch_day(ed, &y, &m, &d);
    }

    if (!myerr) {
	ret = calloc(12, 1);
	if (ret == NULL) {
	    myerr = E_ALLOC;
	} else {
	    sprintf(ret, "%04d-%02d-%02d", y, m, d);
	}
    }

    if (err != NULL) {
	*err = myerr;
    }

    return ret;
}

/**
 * ymd_basic_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @julian: non-zero to use Julian calendar, otherwise Gregorian.
 * @err: location to receive error code.
 *
 * Returns: an 8-digit number on the pattern YYYYMMDD (ISO 8601 basic
 * date format) given the epoch day number, which equals 1 for the
 * first of January in the year 1 AD, or #NADBL on error.
 */

double ymd_basic_from_epoch_day (guint32 ed, int julian, int *err)
{
    int y = 0, m = 0, d = 0;

    if (julian) {
	*err = julian_ymd_bits_from_epoch_day(ed, &y, &m, &d);
    } else {
	*err = ymd_bits_from_epoch_day(ed, &y, &m, &d);
    }

    if (*err) {
	return NADBL;
    } else {
	return 10000*y + 100*m + d;
    }
}

/**
 * weekday_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 *
 * Returns: the weekday (Sunday = 0) corrsponding to @ed,
 * or -1 on error.
 */

int weekday_from_epoch_day (guint32 ed)
{
    GDate date;
    int wd;

    if (!g_date_valid_julian(ed)) {
	return -1;
    }

    g_date_clear(&date, 1);
    g_date_set_julian(&date, ed);
    wd = g_date_get_weekday(&date);

    return wd == G_DATE_SUNDAY ? 0 : wd;
}

/**
 * get_epoch_day:
 * @datestr: string representation of calendar date, in form
 * YY[YY]-MM-DD.
 *
 * Returns: the epoch day number, or -1 on failure.
 */

guint32 get_epoch_day (const char *datestr)
{
    GDate date;
    int y, m, d, nf = 0;
    int ydigits = 0;

    if (strchr(datestr, '-')) {
	ydigits = strcspn(datestr, "-");
	nf = sscanf(datestr, YMD_READ_FMT, &y, &m, &d);
    } else if (strchr(datestr, '/')) {
	ydigits = strcspn(datestr, "/");
	nf = sscanf(datestr, "%d/%d/%d", &y, &m, &d);
    } else if (strlen(datestr) == 8) {
	ydigits = 4;
	nf = sscanf(datestr, "%4d%2d%2d", &y, &m, &d);
    }

    if (nf != 3 || (ydigits != 4 && ydigits != 2)) {
	return 0;
    }

    if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    if (!g_date_valid_dmy(d, m, y)) {
	return 0;
    }

    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);

    return g_date_get_julian(&date);
}

/**
 * calendar_obs_number:
 * @datestr: string representation of calendar date, in form
 * YY[YY]/MM/DD.
 * @dset: pointer to dataset information.
 *
 * Returns: The zero-based observation number for the given
 * date within the current data set.
 */

int calendar_obs_number (const char *datestr, const DATASET *dset)
{
    guint32 ed0 = (guint32) dset->sd0;
    guint32 t = get_epoch_day(datestr);

#ifdef CAL_DEBUG
    fprintf(stderr, "calendar_obs_number: '%s' gave epoch day = %u\n",
	    datestr, t);
#endif

    if (t <= 0 || t < ed0) {
	return -1;
    } else if (t == ed0) {
	return 0;
    }

    /* subtract starting day for dataset */
    t -= ed0;

    if (t <= 0) {
	return -1;
    }

    if (dset->pd == 52) {
	/* weekly data */
	t /= 7;
    } else if (dset->pd == 5 || dset->pd == 6) {
	/* daily, 5- or 6-day week: subtract number of irrelevant days */
	int startday = (ed0 - 1 + DAY1) % 7;
	int wkends = (t + startday - 1) / 7;

#ifdef CAL_DEBUG
	printf("calendar_obs_number: ed0=%d, date=%s, t=%d, startday=%d, wkends=%d\n",
	       (int) ed0, date, (int) t, startday, wkends);
#endif

	if (dset->pd == 5) {
	    t -= (2 * wkends);
	} else {
	    t -= wkends;
	}
    }

    return (int) t;
}

/* convert from $t in dataset to epoch day */

static int t_to_epoch_day (int t, guint32 start, int wkdays)
{
    int startday = (start - 1 + DAY1) % 7;
    int wkends = (t + startday - 1) / wkdays;

    if (wkdays == 5) {
	wkends *= 2;
    }

    return start + t + wkends;
}

/**
 * epoch_day_from_t:
 * @t: 0-based observation index.
 * @dset: pointer to dataset.
 *
 * Returns: the epoch day based on calendrical information
 * in @dset. In case of error 0 is returned.
 */

guint32 epoch_day_from_t (int t, const DATASET *dset)
{
    guint32 d0 = (guint32) dset->sd0;
    guint32 dt = 0;

    if (dset->pd == 52) {
	dt = d0 + 7 * t;
    } else if (dset->pd == 7) {
	dt = d0 + t;
    } else {
	dt = t_to_epoch_day(t, d0, dset->pd);
    }

    return dt;
}

/**
 * calendar_date_string:
 * @targ: string to be filled out.
 * @t: zero-based index of observation.
 * @dset: pointer to dataset.
 *
 * Writes to @targ the calendar representation of the date of
 * observation @t, in the form YY[YY]/MM/DD.
 *
 * Returns: 0 on success, non-zero on error.
 */

int calendar_date_string (char *targ, int t, const DATASET *dset)
{
    guint32 d0, dt = 0;
    int y, m, d;
    int err = 0;

    *targ = '\0';
    d0 = (guint32) dset->sd0;

    if (dset->pd == 52) {
	dt = d0 + 7 * t;
    } else if (dset->pd == 7) {
	dt = d0 + t;
    } else {
	/* 5- or 6-day data */
	if (t == 0 && dset->pd == 5) {
	    int wd = weekday_from_epoch_day(d0);

	    if (wd == 0 || wd == 6) {
		gretl_errmsg_sprintf(_("Invalid starting date for %d-day data"), dset->pd);
		return E_DATA;
	    }
	}
	dt = t_to_epoch_day(t, d0, dset->pd);
    }

    err = ymd_bits_from_epoch_day(dt, &y, &m, &d);

    if (!err) {
	if (strlen(dset->stobs) == 8) {
	    sprintf(targ, YMD_WRITE_Y2_FMT, y % 100, m, d);
	} else {
	    sprintf(targ, YMD_WRITE_Y4_FMT, y, m, d);
	}
    }

    return err;
}

/**
 * MS_excel_date_string:
 * @targ: date string to be filled out.
 * @mst: MS Excel-type date code: days since base.
 * @pd: periodicity of data (or 0 if unknown).
 * @d1904: set to 1 if the base is 1904/01/01; otherwise
 * the base is assumed to be 1899/12/31.
 *
 * Writes to @targ the calendar representation of the date of
 * observation @mst, in the form YYYY-MM-DD if @pd is 0, 5,
 * 6, 7 or 52 (unknown, daily, or weekly frequency), otherwise
 * in the appropriate format for annual, quarterly or monthly
 * according to @pd.
 *
 * Returns: 0.
 */

int MS_excel_date_string (char *targ, int mst, int pd, int d1904)
{
    int y = (d1904)? 1904 : 1900;
    int d = (d1904)? 2 : 1;
    int m = 1;
    int leap, drem;

    *targ = '\0';

    if (mst == 0) {
	/* date coincident with base */
	if (d1904) {
	    d = 1;
	} else {
	    y = 1899;
	    m = 12;
	    d = 31;
	}
    } else if (mst > 0) {
	/* date subsequent to base */
	drem = mst + d1904;

	while (1) {
	    int yd = 365 + leap_year(y);

	    /* MS nincompoopery */
	    if (y == 1900) yd++;

	    if (drem > yd) {
		drem -= yd;
		y++;
	    } else {
		break;
	    }
	}

	leap = leap_year(y) + (y == 1900);

	for (m=1; m<=12; m++) {
	    int md = days_in_month[leap][m];

	    if (drem > md) {
		drem -= md;
	    } else {
		d = drem;
		break;
	    }
	}
    } else {
	/* mst < 0, date prior to base */
	drem = -(mst + d1904);

	y = (d1904)? 1903 : 1899;

	while (1) {
	    int yd = 365 + leap_year(y);

	    if (drem > yd) {
		drem -= yd;
		y--;
	    } else {
		break;
	    }
	}

	leap = leap_year(y);

	for (m=12; m>0; m--) {
	    int md = days_in_month[leap][m];

	    if (drem >= md) {
		drem -= md;
	    } else {
		d = md - drem;
		break;
	    }
	}
    }

    if (pd == 1) {
	sprintf(targ, "%d", y);
    } else if (pd == 12) {
	sprintf(targ, "%d:%02d", y, m);
    } else if (pd == 4) {
	int q = 1 + m / 3.25;

	sprintf(targ, "%d:%d", y, q);
    } else {
	sprintf(targ, YMD_WRITE_Y4_FMT, y, m, d);
    }

    return 0;
}

/**
 * get_dec_date:
 * @datestr: calendar representation of date: YYYY-MM-DD.
 *
 * Returns: representation of date as year plus fraction of year.
 */

double get_dec_date (const char *datestr)
{
    GDate date;
    int y, m, d, nf;
    int ydigits = 0;
    double num, den;

    nf = sscanf(datestr, YMD_READ_FMT, &y, &m, &d);

    if (nf == 3) {
	ydigits = strcspn(datestr, "-");
    } else if (strchr(datestr, '/') != NULL) {
	/* backward compatibility */
	ydigits = strcspn(datestr, "/");
	nf = sscanf(datestr, "%d/%d/%d", &y, &m, &d);
    }

    if (nf != 3 || (ydigits != 4 && ydigits != 2)) {
	return NADBL;
    }

    if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    if (!g_date_valid_dmy(d, m, y)) {
	return NADBL;
    }

    den = 365 + g_date_is_leap_year(y);
    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);
    num = g_date_get_day_of_year(&date);

    return y + num / den;
}

/**
 * day_of_week:
 * @y: year.
 * @m: month, 1 to 12.
 * @d: day in month, 1 to 31.
 * @julian: non-zero to use Julian calendar, otherwise Gregorian.
 * @err: location to receive error code.
 *
 * Returns: the day of the week for the supplied date
 * (Sunday = 0, Monday = 1, ...) or %NADBL on failure
 * (the date is invalid).
 */

double day_of_week (int y, int m, int d, int julian, int *err)
{
    int wd = day_of_week_from_ymd(y, m, d, julian);

    return wd < 0 ? NADBL : wd;
}

#define day_in_calendar(w, d) (((w) == 6 && d != 0) || \
			       ((w) == 5 && d != 0 && d != 6))

/**
 * weekday_from_date:
 * @datestr: calendar representation of date, [YY]YY/MM/DD
 *
 * Returns: day of week as integer, Sunday = 0.
 */

int weekday_from_date (const char *datestr)
{
    int y, m, d;
    int ydigits;

    if (sscanf(datestr, YMD_READ_FMT, &y, &m, &d) != 3) {
	return -1;
    }

    ydigits = strcspn(datestr, "-");
    if (ydigits != 4 && ydigits != 2) {
	return -1;
    }

    if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    return day_of_week_from_ymd(y, m, d, 0);
}

/**
 * day_starts_month:
 * @d: day of month, 1-based
 * @m: month number, 1-based
 * @y: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 * @pad: location to receive 1 if the first day of the month
 * can reasonably be padded by a missing value (Jan 1), or NULL.
 *
 * Returns: 1 if the day is the "first day of the month",
 * allowance made for the possibility of a 5- or 6-day week,
 * else 0.
 */

int day_starts_month (int d, int m, int y, int wkdays, int *pad)
{
    int ret = 0;

    if (wkdays == 7) {
	if (d == 1) {
	    ret = 1;
	} else if (pad != NULL && m == 1 && d == 2) {
	    /* second of January */
	    *pad = 1;
	    ret = 1;
	}
    } else {
	/* 5- or 6-day week: check for first weekday or non-Sunday */
	int i, idx = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<6; i++) {
	   if (day_in_calendar(wkdays, idx % 7)) {
	       break;
	   }
	   idx++;
	}
	if (d == i) {
	    ret = 1;
	} else if (pad != NULL && m == 1 && d == i + 1) {
	    /* January again */
	    *pad = 1;
	    ret = 1;
	}
    }

    return ret;
}

/**
 * day_ends_month:
 * @d: day of month, 1-based
 * @m: month number, 1-based
 * @y: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: 1 if the day is the "last day of the month",
 * allowance made for the possibility of a 5- or 6-day week, else 0.
 */

int day_ends_month (int d, int m, int y, int wkdays)
{
    int ret = 0;
    int leap = (m == 2)? leap_year(y) : 0;
    int dm = days_in_month[leap][m];

    if (wkdays == 7) {
	ret = (d == dm);
    } else {
	/* 5- or 6-day week: check for last weekday or non-Sunday */
	int i, idx = day_of_week_from_ymd(y, m, dm, 0);

	for (i=dm; i>0; i--) {
	    if (day_in_calendar(wkdays, idx % 7)) {
		break;
	    }
	    idx--;
	}
	ret = (d == i);
    }

    return ret;
}

/**
 * get_days_in_month:
 * @m: month number, 1-based
 * @y: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 * @julian: non-zero for Julian calendar, otherwise Gregorian.
 *
 * Returns: the number of (relevant) days in the month, allowance
 * made for the possibility of a 5- or 6-day week.
 */

int get_days_in_month (int m, int y, int wkdays, int julian)
{
    int dm, leap = 0;
    int ret = 0;

    if (m == 2) {
	leap = julian ? (y%4 == 0) : leap_year(y);
    }
    dm = days_in_month[leap][m];

    if (wkdays == 7) {
	ret = dm;
    } else {
	int i, idx = day_of_week_from_ymd(y, m, 1, julian);

	for (i=0; i<dm; i++) {
	    if (day_in_calendar(wkdays, idx % 7)) {
		ret++;
	    }
	    idx++;
	}
    }

    return ret;
}

/**
 * days_in_month_before:
 * @y: 4-digit year
 * @m: month number, 1-based
 * @d: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the number of relevant days in the month prior to
 * the supplied date, allowing for the possibility of a 5- or
 * 6-day week.
 */

int days_in_month_before (int y, int m, int d, int wkdays)
{
    int ret = 0;

    if (wkdays == 7) {
	ret = d - 1;
    } else {
	int i, idx = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<d; i++) {
	    if (day_in_calendar(wkdays, idx % 7)) {
		ret++;
	    }
	    idx++;
	}
    }

    return ret;
}

/**
 * days_in_month_after:
 * @y: 4-digit year
 * @m: month number, 1-based
 * @d: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the number of relevant days in the month after
 * the supplied date, allowing for the possibility of a 5- or
 * 6-day week.
 */

int days_in_month_after (int y, int m, int d, int wkdays)
{
    int leap = (m == 2)? leap_year(y) : 0;
    int dm = days_in_month[leap][m];
    int ret = 0;

    if (wkdays == 7) {
	ret = dm - d;
    } else {
	int i, wd = day_of_week_from_ymd(y, m, dm, 0);

	for (i=dm; i>d; i--) {
	    if (day_in_calendar(wkdays, wd)) {
		ret++;
	    }
	    if (wd > 0) {
		wd--;
	    } else {
		wd = 6;
	    }
	}
    }

    return ret;
}

/**
 * date_to_daily_index:
 * @datestr: date in format YYYY-MM-DD.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the zero-based index of the specified day
 * within the specified month and year. In the case
 * of 5- or 6-day data index zero does not necessarily
 * correspond to the first day of the month but rather
 * to the first relevant day.
 */

int date_to_daily_index (const char *datestr, int wkdays)
{
    int y, m, d, seq = 0;

    if (sscanf(datestr, YMD_READ_FMT, &y, &m, &d) != 3) {
	return -1;
    }

    if (wkdays == 7) {
	seq = d - 1;
    } else {
	int leap = leap_year(y);
	int n = days_in_month[leap][m];
	int i, idx = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<=n; i++) {
	    if (d == i) {
		break;
	    }
	    if (day_in_calendar(wkdays, idx % 7)) {
		seq++;
	    }
	    idx++;
	}
    }

    return seq;
}

/**
 * daily_index_to_date:
 * @targ: location to receive the date (YYYY-MM-DD).
 * @y: year.
 * @m: month.
 * @idx: zero-based index of day within month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Fills out @targ with the calendar data implied by
 * the specification of @y, @m, @seq and @wkdays,
 * provided this specification corresponds to an actual
 * calendar date.

 * Returns: 0 on successful completion, non-zero if
 * there is no such calendar date.
 */

int daily_index_to_date (char *targ, int y, int m, int idx,
			 int wkdays)
{
    int day = 0;

    *targ = '\0';

    if (m < 1 || m > 12 || idx < 0 || idx > 30) {
	fprintf(stderr, "daily_index_to_date: y=%d, m=%d, seq=%d\n",
		y, m, idx);
	return E_INVARG;
    }

    if (wkdays == 7) {
	day = idx + 1;
    } else {
	int leap = leap_year(y);
	int n = days_in_month[leap][m];
	int wd = day_of_week_from_ymd(y, m, 1, 0);
	int i, seq = 0;

	for (i=1; i<=n; i++) {
	    if (day_in_calendar(wkdays, wd % 7)) {
		if (seq == idx) {
		    day = i;
		    break;
		}
		seq++;
	    }
	    wd++;
	}
    }

    if (day <= 0) {
	return E_DATA;
    } else {
	sprintf(targ, YMD_WRITE_FMT, y, m, day);
	return 0;
    }
}

/**
 * n_hidden_missing_obs:
 * @dset: dataset information.
 * @t1: first observation.
 * @t2: last observation.
 *
 * For daily data with user-supplied data strings,
 * determine the number of "hidden" missing observations
 * in the range @t1 to @t2 inclusive. This is
 * the difference between the actual number of
 * observations and the number that should be there,
 * according to the calendar. Allowance is made for
 * 5- or 6-day data, via the data frequency given
 * in @dset.
 *
 * Returns: number of hidden observations.
 */

int n_hidden_missing_obs (const DATASET *dset, int t1, int t2)
{
    int n_present = t2 - t1 + 1;
    int cal_n;

    if (!dated_daily_data(dset) || dset->S == NULL) {
	return 0;
    }

    t1 = calendar_obs_number(dset->S[t1], dset);
    t2 = calendar_obs_number(dset->S[t2], dset);

    cal_n = t2 - t1 + 1;

    return cal_n - n_present;
}

/**
 * guess_daily_pd:
 * @dset: dataset information.
 *
 * Based on user-supplied daily date strings recorded in
 * @dset, try to guess whether the number of observations
 * per week is 5, 6 or 7 (given that some observations
 * may be missing).
 *
 * Returns: best quess at data frequency.
 */

int guess_daily_pd (const DATASET *dset)
{
    int t, wd, pd = 5;
    int wdbak = -1;
    int havesat = 0;
    int gotsat = 0, gotsun = 0;
    int contig = 0;

    wd = weekday_from_date(dset->S[0]);
    if (6 - wd < dset->n) {
	havesat = 1;
    }

    for (t=0; t<dset->n && t<28; t++) {
	wd = weekday_from_date(dset->S[t]);
	if (wd == 0) {
	    gotsun = 1;
	} else if (wd == 6) {
	    gotsat = 1;
	}
	if ((wdbak + 1) % 7 == wd) {
	    contig++;
	}
	wdbak = wd;
    }

    if (gotsat && gotsun) {
	pd = 7;
    } else if (contig > 10) {
	if (gotsun) pd = 7;
	else if (gotsat) pd = 6;
    } else if (dset->n > 7) {
	if (!gotsun && !gotsat) {
	    pd = 5;
	} else if (!gotsun) {
	    pd = 6;
	}
    } else if (havesat && !gotsat) {
	pd = 5;
    } else {
	pd = 7;
    }

    return pd;
}

/**
 * iso_basic_to_extended:
 * @b: source array of YYYYMMDD values.
 * @y: array to hold year values.
 * @m: array to hold month values.
 * @d: array to hold day-of-week values.
 * @n: length of all the above arrays.
 *
 * Given the array @b of ISO 8601 "basic" daily dates (YYYYMMDD as
 * doubles), fill out the arrays @y, @m and @d with year, month
 * and day.
 *
 * Returns: 0.
 */

int iso_basic_to_extended (const double *b, double *y, double *m,
			   double *d, int n)
{
    int bi, yi, mi, di;
    int i, julian;

    for (i=0; i<n; i++) {
	julian = 0;
	if (na(b[i])) {
	    y[i] = m[i] = NADBL;
	    if (d != NULL) {
		d[i] = NADBL;
	    }
	} else {
	    bi = (int) b[i];
	    if (bi < 0) {
		julian = 1;
		bi = -bi;
	    }
	    yi = bi / 10000;
	    mi = (bi - 10000*yi) / 100;
	    di = bi - 10000*yi - 100*mi;
	    /* now check for legit date */
	    if (!valid_ymd(yi, mi, di, julian)) {
		y[i] = m[i] = NADBL;
		if (d != NULL) {
		    d[i] = NADBL;
		}
	    } else {
		y[i] = yi;
		m[i] = mi;
		if (d != NULL) {
		    d[i] = di;
		}
	    }
	}
    }

    return 0;
}

/**
 * easterdate:
 * @year: year for which we want Easter date (Gregorian).
 *
 * Algorithm taken from Wikipedia page
 * https://en.wikipedia.org/wiki/Computus
 * under the heading "Anonymous Gregorian algorithm".
 *
 * Returns: the date of Easter in the Gregorian calendar as
 * (month + day/100). Note that April the 10th is, under
 * this convention, 4.1; hence, 4.2 is April the 20th, not
 * April the 2nd (which would be 4.02).
 */

double easterdate (int year)
{
    int a = year % 19;
    int b = year / 100;
    int c = year % 100;
    int d = b / 4;
    int e = b % 4;
    int f = (b + 8) / 25;
    int g = (b - f + 1) / 3;
    int h = (19 * a + b - d - g + 15) % 30;
    int i = c / 4;
    int k = c % 4;
    int L = (32 + 2 * e + 2 * i - h - k) % 7;
    int m = (a + 11 * h + 22 * L) / 451 ;

    int month = (h + L - 7 * m + 114) / 31;
    int day = ((h + L - 7 * m + 114) % 31) + 1;

    return month + day * 0.01;
}

/**
 * dayspan:
 * @ed1: first epoch day.
 * @ed2: last epoch day.
 * @wkdays: relevant days per week (5, 6 or 7).
 * @err: location to receive error code.
 *
 * Returns: The number of days in the interval @ed1 to
 * @ed2, inclusive, taking account of the number of daily
 * observations per week, @wkdays. If @wkdays = 6 Sundays
 * are disregarded; if @wkdays = 5 both Saturdays and
 * Sundays are disregarded.
 */

int day_span (guint32 ed1, guint32 ed2, int wkdays, int *err)
{
    int n = 0;
    
    if (!g_date_valid_julian(ed1) ||
	!g_date_valid_julian(ed2) ||
	ed2 < ed1) {
	*err = E_INVARG;
    } else if (wkdays == 7) {
	/* simple! */
	n = ed2 - ed1 + 1;
    } else {
	GDate date;
	guint32 i;
	int idx;

	g_date_clear(&date, 1);
	g_date_set_julian(&date, ed1);
	idx = g_date_get_weekday(&date) - 1;
	
	for (i=ed1; i<=ed2; i++) {
	    if (day_in_calendar(wkdays, idx % 7)) {
		n++;
	    }
	    idx++;
	}
    }

    return n;
}

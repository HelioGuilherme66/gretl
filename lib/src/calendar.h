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

#ifndef CALENDAR_H
#define CALENDAR_H

#define FOUR_DIGIT_YEAR(y) ((y < 50)? y + 2000 : y + 1900)

long get_epoch_day (const char *date);

int get_day_of_week (const char *date);

int day_starts_month (int d, int m, int y, int wkdays, int *pad);

int day_ends_month (int d, int m, int y, int wkdays);

int get_days_in_month (int mon, int yr, int wkdays);

int days_in_month_before (int yr, int mon, int day, int wkdays);

int days_in_month_after (int yr, int mon, int day, int wkdays);

int daily_obs_number (const char *date, const DATAINFO *pdinfo);

void daily_date_string (char *str, int t, const DATAINFO *pdinfo);

double get_dec_date (const char *date);

int n_hidden_missing_obs (const DATAINFO *pdinfo);

int guess_daily_pd (const DATAINFO *pdinfo);

#endif /* CALENDAR_H */ 

set more off
set matsize 800

use nlswork.dta,clear

iis idcode
tis year

gen byte black = race==2
xtreg ln_wage grade age ttl_exp tenure black not_smsa south, re sa theta


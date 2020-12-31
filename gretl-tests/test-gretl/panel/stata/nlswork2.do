set mem 200m
set more off
set matsize 800

use nlswork.dta,clear

iis idcode
tis year

gen age2 = age^2
gen ttl_exp2 = ttl_exp^2
gen tenure2 = tenure^2
gen byte black = race==2

* GLS Random-effects model 
xtreg ln_w grade age* ttl_exp* tenure* black not_smsa south, re sa theta
estimates store re_est
quietly xtreg ln_w grade age* ttl_exp* tenure* black not_smsa south, fe
hausman . re_est


set mem 200m
set more off
set matsize 800

use nlswork.dta,clear

iis idcode
tis year

** describe the patterns of the data
* xtdes, patterns(30)

** estimate the model using 'xtreg'
* Dep variable = ln_wage
* Regressors = grade age ttl_exp tenure not_smsa south
*              And the square terms of age ttl_exp tenure are also included

gen age2 = age^2
gen ttl_exp2 = ttl_exp^2
gen tenure2 = tenure^2
gen byte black = race==2

* between-effects model
xtreg ln_wage grade age age2 ttl_exp ttl_exp2 tenure tenure2 black not_smsa south, be wls

* fixed-effects model
xtreg ln_wage grade age age2 ttl_exp ttl_exp2 tenure tenure2 black not_smsa south, fe 

* GLS Random-effects model 
xtreg ln_wage grade age age2 ttl_exp ttl_exp2 tenure tenure2 black not_smsa south, re sa theta
estimates store est1
xtreg ln_wage grade age ttl_exp tenure black not_smsa south, re sa theta
estimates store est2
hausman est1 est2


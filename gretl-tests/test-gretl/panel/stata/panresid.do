use nlswork.dta,clear

iis idcode
tis year

gen age2 = age^2
gen ttl_exp2 = ttl_exp^2
gen tenure2 = tenure^2
gen byte black = race==2

* fixed-effects model
xtreg ln_wage grade age age2 ttl_exp ttl_exp2 tenure tenure2 black not_smsa south, fe
predict e_fe, e
predict ue_fe, ue

* Random-effects model 
xtreg ln_wage grade age age2 ttl_exp ttl_exp2 tenure tenure2 black not_smsa south, re
predict e_re, e
predict ue_re, ue
predict u_re, u

summarize e_fe ue_fe e_re ue_re u_re
list e_fe ue_fe e_re ue_re u_re in 1/10



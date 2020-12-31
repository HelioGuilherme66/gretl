insheet using baltagi_Gasoline.csv
iis unit
tis year
xtreg gas y rp car, re sa theta
xtoverid
xtreg gas y rp car, re sa theta vce(cluster unit)
xtoverid



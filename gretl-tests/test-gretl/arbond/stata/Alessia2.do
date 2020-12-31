set mem 100m

use provagretl
egen pan=group(country)
gene panelid=isic+1000*pan
tsset panelid year
gene lprod = log(prod)
gene lemp = log(emp)
gene lva = log(va)

set matsize 800

xtabond lprod lemp lva, maxldep(2) twostep

*
* Labor Force Participation Trends
* pp 73-75
*
cal 1957 1 12
all 1995:3
open data lpart.prn
data(format=prn,org=columns)
*
* Graphs of women's and men's participation rates
*
graph(header='Figure 4.1 Labor Force Participation Rate for Females')
# lpartw
*
graph(header='Figure 4.3 Labor Force Participation Rate for Males')
# lpartm
*
* Regressions on linear trend, graphed with fitted values
*
set time = t
linreg lpartw
# constant time
prj fitw
linreg lpartm
# constant time
prj fitm
*
graph(header='Figure 4.4 Linear Trend: Labor Force Participation Rates for Females') 2
# lpartw
# fitw
graph(header='Figure 4.5 Linear Trend: Labor Force Participation Rates for Males') 2
# lpartm
# fitm


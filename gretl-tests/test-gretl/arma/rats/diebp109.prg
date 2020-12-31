*
* Housing starts example, pp 109-113
*
cal 1946 1 12
all 1994:11
open data hstarts.dat
data(format=prn,org=columns)
*
graph(header='Figure 5.4 Housing Starts, 1946.01-1994.11')
# hstarts
graph(header='Figure 5.5 Housing Starts, 1990.01-1994.11')
# hstarts 1990:1 *
*
* The seasonal instruction defines a dummy for the last period of
* the year (here the month of December). The first lead of this
* series (that is, the series created by taking the data for the
* next time period into the future) is November, the second lead
* is October, etc. The lag field seasonal{0 to -11} (in RATS, the
* -lags are leads) covers all 12 dummies. If you want all but one
* dummy, use {0 to -10}, which will leave out January.
*
seasonal seasons
*
* The extra parameters on the linreg (after the sample range)
* are for series for the residuals and for the coefficients. The
* coefficients are saved this way only in a very limited number
* of situations (such as this one) where they (in whole or in
* part) form an interesting sequence.
*
linreg hstarts * 1993:12 resids sfactors
# seasons{0 to -11}
*
graph(header='Figure 5.6 Residual Plot')
# resids
*
* This graphs the series of seasonal factors. The NODATES option is used
* to force the graph to be labeled with sequence numbers rather than dates.
* We've chosen to graph these as a bar rather than a line graph.
*
graph(header='Figure 5.7 Estimated Seasonal Factors: Housing Starts',nodates,style=bargraph)
# sfactors
*
* Generate forecasts with their standard errors, and compute 95% confidence bands.
*
prj(stderrs=stderrs) hfore 1994:1 1994:11
set upperf = hfore + 1.96*stderrs
set lowerf = hfore - 1.96*stderrs
*
* Create a series which is 1's in the period that we want to shade.
*
set forezone = %year(t)==1994
graph(header='Figure 5.8 Housing Starts: History and Forecast (1994:1-1994:11)',shading=forezone) 4
# hstarts 1990:1 1993:12
# hfore
# upperf
# lowerf
graph(header='Figure 5.9 Housing Starts: History. Forecast and Realization (1994:1-1994:11)',shading=forezone) 4
# hstarts 1990:1 1994:12
# hfore
# upperf
# lowerf


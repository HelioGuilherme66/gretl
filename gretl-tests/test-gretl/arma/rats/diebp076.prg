*
* Examples using NYSE volume data
* pp 76-80
*
cal 1947 1 12
all 1992:6
open data nyse_vol.dat
data(format=prn,org=columns)
*
* Simple graph
*
graph(header='Figure 4.6 Volume on the New York Stock Exchange')
# fsvol
*
* Fit to quadratic trend, graphed with fitted values
*
set time = t
set time2 = t**2
linreg fsvol
# constant time time2
prj qfit
graph(header='Figure 4.8 Quadratic Trend: Volume on the New York Stock Exchange',key=upleft) 2
# fsvol
# qfit
*
* Graph of logged series
*
set logvol = log(fsvol)
graph(header='Figure 4.9 Log Volume on the New York Stock Exchange')
# logvol
*
* Do linear trend in logs
*
linreg logvol
# constant time
prj efit
*
graph(header='Figure 4.11 Linear Trend: Log Volume on the New York Stock Exchange',key=upleft) 2
# logvol
# efit
*
* Exponential trend done by exponentiating the fitted values from the log linear trend
*
set efit = exp(efit)
*
graph(header='Figure 4.12 Exponential Trend: Volume on the New York Stock Exchange',key=upleft) 2
# fsvol
# efit


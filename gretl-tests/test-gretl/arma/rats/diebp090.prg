*
* Retail sales examples from pp 90-97
*
source(noecho) regcrits.src
cal 1954 1 12
all 1994:12
open data rsales.dat
data(format=prn,org=columns)
graph(header='Figure 4.14 Retail Sales')
# rsales
*
set time = t
set time2 = t**2
*
linreg rsales 1955:1 1993:12
# constant time
@regcrits
linreg rsales 1955:1 1993:12
# constant time time2
@regcrits
*
set logrsales = log(rsales)
linreg logrsales 1955:1 1993:12
# constant time
@regcrits
*
* Estimation of non-linear model using NLLS
*
nonlin b0 b1
frml etrend = b0*exp(b1*t)
nlls(frml=etrend) rsales 1955:1 1993:12
@regcrits
*
* Create series which is 1's where we want shading
*
set y1994 = %year(t)==1994
*
* Forecasts and 95% bands for quadratic trend forecasts
*
linreg rsales 1955:1 1993:12
# constant time time2
prj(stderrs=stderrs) qtrend 1994:1 1994:12
set upperq 1994:1 1994:12 = qtrend+1.96*stderrs
set lowerq 1994:1 1994:12 = qtrend-1.96*stderrs
*
graph(header='Figure 4.19 Retail Sales with Quadratic Trend Forecast',shading=y1994) 4
# rsales 1990:1 1993:12
# qtrend
# upperq
# lowerq
*
* Same thing, with actual sales through the forecast period
*
graph(header='Figure 4.20 Retail Sales with Quadratic Trend Forecast and Realization',shading=y1994) 4
# rsales 1990:1 1994:12
# qtrend
# upperq
# lowerq
*
* Forecasts and 95% bands for linear trend forecasts
*
linreg rsales 1955:1 1993:12
# constant time
prj(stderrs=stderrs) ltrend 1994:1 1994:12
set upperl 1994:1 1994:12 = ltrend+1.96*stderrs
set lowerl 1994:1 1994:12 = ltrend-1.96*stderrs
*
graph(header='Figure 4.21 Retail Sales with Linear Trend Forecast',shading=y1994) 4
# rsales 1990:1 1993:12
# ltrend
# upperl
# lowerl
*
* Same thing, with actual sales through the forecast period
*
graph(header='Figure 4.22 Retail Sales with Linear Trend Forecast and Realization',shading=y1994) 4
# rsales 1990:1 1994:12
# ltrend
# upperl
# lowerl


*
* Forecasting liquor sales, pp 210-227
*
source(noecho) regcorrs.src
source(noecho) uforecst.src
cal 1967 1 12
all 1994:12
open data liquor.dat
data(format=prn,org=columns)
graph(header='Figure 9.1 Liquor Sales, 1968.01-1993.12')
# liquor 1968:1 1993:12
set lsales = log(liquor)
graph(header='Figure 9.2 Log Liquor Sales, 1968.01-1993.12')
# lsales 1968:1 1993:12
*
set time = (t-1967:1)
set time2 = time**2
*
* The time trend variables are created in order to match up with the output
* given in the book. Ordinarily, you'd create these as
*
*  set time  = t   and
*  set time2 = t**2
*
* With those definitions, you'll get an identical fit, but the coefficients
* will be slightly different since the time2 used in the book is (t-1)**2 and
* time is t-1 so the regression function a*1+b*(t-1)+c*(t-1)**2 used in the
* text is equivalent to (a-b+c)+(b-2c)*t+c*t**2. Thus, the coefficient on
* time2 will be the same both ways, but the coefficients on "time" and the
* intercept need to adjust for the different definitions.
*
linreg lsales 1968:1 1993:12 resids
# constant time time2
*
@regcorrs(number=36,window='Quadratic Trend Regression') resids
*
* Generate seasonal dummy series. seasons{0 to -11} is used to
* include a complete set of dummies in the regression.
*
seasonal seasons
linreg lsales 1968:1 1993:12 resids
# time time2 seasons{0 to -11}
@regcorrs(number=36,window='Quadratic Trend + Seasonals') resids
*
* Now, quadratic trend + seasonals + AR(3) noise term. For a regression
* with ARMA errors, you use the BOXJENK instruction with the REGRESSORS
* option to include the variables other than CONSTANT and the ARMA
* parameters. You could also do this using a LINREG like the previous one
* with the addition of lsales{1 to 3} among the explanatory variables.
* This will give an identical fit, but with different coefficients. This
* works because:
*
*   (a) There are only AR parameters, not MA's
*   (b) The explanatory variables are simple functions of time.
*
* For the LINREG, the model being estimated is the straightforward regression:
*  lsales = g1 * lsales{1} + g2 * lsales{2} + g3*lsales{3} +
*       a0 + a1 * time + a2 * time**2 + b1 * D1 + .... + noise
*
* BOXJENK parameterizes it differently as
*   z = lsales - (A0 + A1 * time + A2 * time**2 + B1 * D1 + ....)
*   z = G1*z{1} + G2*z{2} + G3*z{3} + noise
*
* If you expand this out, you'll find that you get a lot of terms with lags
* of the time and time**2 variables and the seasonal dummies. But since lags of
* seasonal dummies are just seasonal dummies for other months, and lags of time
* and time**2, when expanded just generate other terms in time**2, time and constant,
* the end result is just a reshuffling of same set of regressors.
*
boxjenk(regressors,ar=3,noconstant) lsales 1968:1 1993:12 resids
# time time2 seasons{0 to -11}
@regcorrs(number=36,window='Quadratic Trend + Seasonals + AR(3)') resids
*
* In order to do forecasting, we need to define an equation for the model
*
boxjenk(regressors,ar=3,noconstant,define=ar3eq) lsales 1968:1 1993:12 resids
# time time2 seasons{0 to -11}
*
* Compute the forecasts, and their standard errors over the year 1994. Generate
* upper and lower 95% confidence bands
*
@uforecast(stderrs=stderrs,equation=ar3eq) fcst 1994:1 1994:12
set upper 1994:1 1994:12 = fcst+1.96*stderrs
set lower 1994:1 1994:12 = fcst-1.96*stderrs
*
* Create a dummy variable to be used in shading the forecast period. Since we're
* doing forecasts out to 1998:12 later, we create this over that entire period.
*
set forezone * 1998:12 = t>=1994:1
*
graph(header='Figure 9.10. Log Liquor Sales: History and 12-Month-Ahead Forecast',shading=forezone) 4
# lsales 1992:1 1993:12
# fcst 1994:1 1994:12
# upper 1994:1 1994:12
# lower 1994:1 1994:12
graph(header='Figure 9.11. Log Liquor Sales: History, 12-Month-Ahead Forecast and Realization',shading=forezone) 4
# lsales 1992:1 1994:12
# fcst 1994:1 1994:12
# upper 1994:1 1994:12
# lower 1994:1 1994:12
*
* To forecast out to 1998:12, we need to extend the time trend and dummy variables.
*
set time  1967:1 1998:12 = (t-1967:1)
set time2 1967:1 1998:12 = time**2
seasonal seasons 1967:1 1999:12
*
@uforecast(stderrs=stderrs,equation=ar3eq) fcst 1994:1 1998:12
set upper 1994:1 1998:12 = fcst+1.96*stderrs
set lower 1994:1 1998:12 = fcst-1.96*stderrs
graph(header='Figure 9.12. Log Liquor Sales: History and 60-Month-Ahead Forecast',shading=forezone) 4
# lsales 1990:1 1993:12
# fcst 1994:1 1998:12
# upper 1994:1 1998:12
# lower 1994:1 1998:12
graph(header='Figure 9.13. Log Liquor Sales: Long History and 60-Month-Ahead Forecast',shading=forezone) 2
# lsales 1968:1 1993:12
# fcst 1994:1 1998:12
*
* Transform log forecasts back to levels
*
set efcst 1994:1 1998:12 = exp(fcst)
graph(header='Figure 9.14. Liquor Sales: Long History and 60-Month-Ahead Forecast',shading=forezone) 2
# liquor 1968:1 1993:12
# efcst 1994:1 1998:12


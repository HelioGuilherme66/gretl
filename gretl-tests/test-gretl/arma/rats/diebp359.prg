*
* Exponential smoothing model, pp 359-360.
*
cal 1973 1 12
all 1996:7
open data exch.dat
data(format=prn,org=columns)
*
set logyen = log(yen)
*
* Holt-Winters smoothing is done using the ESMOOTH instruction in RATS.
* The first three parameters are series to be smoothed, followed by
* the estimation range. The fourth parameter (if included) is the
* series into which the forecasts are placed, while the fifth is the
* number of forecast periods.
*
* The specific form of smoothing is determined by the TREND and SEASONAL
* options. TREND can take the choices NONE, LINEAR and EXPONENTIAL.
* SEASONAL can be NONE, ADDITIVE and MULTIPLICATIVE. TREND=NONE,SEASONAL=NONE
* (which is the default) is the basic exponential smoothing process
* described on page 355. TREND=LINEAR,SEASONAL=NONE is on page 356, and
* TREND=LINEAR,SEASONAL=ADDITIVE is on page 357.
*
* The TREND and SEASONAL options can also take the choice SELECT, which
* has ESMOOTH try all the possibilities (all 3x3 if you use SELECT for both),
* and choose the one which minimizes the Schwarz Criterion.
*
* You can either have ESMOOTH estimate the parameters (option ESTIMATE) or
* set them directly. Since RATS comes up with different parameters than
* shown in the text (EViews uses a grid search procedure, RATS doesn't)
* when it estimates, we'll show both. Note, by the way, that RATS uses
* "gamma" for the trend smoothing parameter which is "beta" in the text.
*
set forezone * 2010:12 = t>=1995:1
esmooth(trend=linear,estimate) logyen * 1994:12 fsmooth 18
*
graph(header='Figure 12.18 Log Yen/Dollar Rate: History and Forecast',$
  subhead='Holt-Winters Smoothing',shading=forezone) 2
# logyen 1990:1 1994:12
# fsmooth 1995:1 1996:7
*
graph(header='Figure 12.20 Log Yen/Dollar Rate: History, Forecast and Realization',$
  subhead='Holt-Winters Smoothing',shading=forezone) 2
# logyen 1990:1 1996:7
# fsmooth 1995:1 1996:7
*
esmooth(trend=linear,estimate) logyen * 1994:12 fsmooth (2010:12-1994:12)
*
graph(header='Figure 12.19 Log Yen/Dollar Rate: History and Long-Horizon Forecast',$
  subhead='Holt-Winters Smoothing',shading=forezone) 2
# logyen 1990:1 1994:12
# fsmooth 1995:1 2010:12
*
*
*  Same long-horizon forecast, but with parameters set to the values used
*  in the text.
*
esmooth(trend=linear,alpha=1.0,gamma=.09) logyen * 1994:12 fsmooth (2010:12-1994:12)
*
graph(header='Figure 12.19 Log Yen/Dollar Rate: History and Long-Horizon Forecast',$
  subhead='Holt-Winters Smoothing',shading=forezone) 2
# logyen 1990:1 1994:12
# fsmooth 1995:1 2010:12


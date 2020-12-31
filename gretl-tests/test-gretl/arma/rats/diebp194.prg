*
* Forecasts from ARMA models for Canadian employment data
* pp 194-198
*
source(noecho) bjest.src
source(noecho) uforecst.src
cal 1961 1 4
all 1994:4
open data caemp.dat
data(format=prn,org=columns)
*
* Models with MA's will often produce somewhat different estimates from
* different software because there isn't a single "correct" way to
* deal with the unobservable pre-sample lagged residuals. Because the
* MA(4) isn't really a good model for this data set (in fact, this exercise
* is designed to demonstrate that), the differences are more pronounced
* than you typically see. RATS estimates a model with a somewhat lower
* intercept than the unconditional mean, so the long-run forecasts go to
* a value more like 96 than 100.2.
*
@bjest(ma=4,constant,define=ma4eq) caemp * 1993:4 ma4res
@uForecast(equation=ma4eq,stderr=ma4std) ma4fore 1994:1 1994:4
set upper = ma4fore+ma4std*1.96
set lower = ma4fore-ma4std*1.96
set forezone * 1996:4 = t>=1994:1
graph(header='Figure 8.1 Employment History and Forecast: MA(4) Model',shading=forezone) 4
# caemp 1990:1 1993:4
# ma4fore
# upper
# lower
*
@uForecast(equation=ma4eq,stderr=ma4std) ma4fore 1994:1 1996:4
set upper 1994:1 1996:12 = ma4fore+ma4std*1.96
set lower 1994:1 1996:12 = ma4fore-ma4std*1.96
graph(header='Figure 8.2 Employment History and Long-Horizon Forecast: MA(4) Model',shading=forezone) 4
# caemp 1990:1 1993:4
# ma4fore
# upper
# lower


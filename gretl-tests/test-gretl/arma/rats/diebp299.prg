*
* Shipping volume examples from pages 299-313
*
source(noecho) bjident.src
source(noecho) uforeerr.src
source(noecho) histogrm.src
*
cal(weekly) 1988 1 1
all 1997:7:18
open data data.dat
data(format=prn,org=columns)
*
graph(key=upleft,header='Figure 11.1 Shipping Volume: Quantitative Forecast and Realization') 2
# vol
# volq
graph(key=upleft,header='Figure 11.2 Shipping Volume: Judgemental Forecast and Realization') 2
# vol
# volj
set qerror = vol-volq
set jerror = vol-volj
*
* The min and max options are used to force a common scale to the
* graphs of the two errors. If each graph is allowed to choose its
* own scale, the judgemental errors will *look* bigger than they are
* when you compare graphs, because the graph instruction will try to
* spread the actual range (which is smaller for judgemental than for
* quantitative) across the full height of the graph.
*
graph(header='Figure 11.3 Quantitative Forecast Error',min=-6.0,max=6.0)
# qerror
graph(header='Figure 11.4 Judgemental Forecast Error',min=-6.0,max=6.0)
# jerror
*
*  The method=yule option is used to give the same results as shown
*  in the text. By default, @BJIDENT and the RATS instruction CORRELATE
*  use another algorithm, called "Burg," which gives very similar results
*  when the correlations are fairly small (as they are here), but differ
*  more substantially (with the Burg method considered to be more
*  accurate) when the correlations are much larger.
*
*  Note also that the standard error bands widen slightly with increasing
*  lags. This is because, unlike the Bartlett standard errors, which
*  assume white noise, the standard errors computed by RATS assume, when
*  computing the AC for lag K+1, that K+1 and up are zero, but 1,...K
*  aren't.
*
@bjident(number=6,print,method=yule) qerror
@bjident(number=6,print,method=yule) jerror
*
@histogram qerror
@histogram jerror
*
*
*  MA(1) regressions from page 307. Regressions with MA terms won't match
*  exactly because of differences in methods of handling pre-sample
*  residuals.
*
boxjenk(ma=1,constant) qerror
boxjenk(ma=1,constant) jerror
*
* Mincer-Zarnowitz Regressions
*
boxjenk(regressors,ma=1,constant) vol / resids
# volq
*
* Test joint restriction that the coefficient on the intercept is zero
* and that on the forecast is 1. For BOXJENK with regressors, the constant
* is the first coefficient, the ARMA terms come next, and regressors last,
* so our restriction is on coefficients 1 and 3.
*
test
# 1 3
# 0.0 1.0
*
* Same thing, for judgmental forecasts
*
boxjenk(regressors,ma=1,constant) vol
# volj
test
# 1 3
# 0.0 1.0
*
* Not in text. This handles the serial correlation by estimating a standard
* linear regression, then correcting the covariance matrix for the serially
* correlated residuals. This technique is commonly applied in such situations.
*
linreg(robusterrors,lags=1) vol
# constant volq
test
# 1 2
# 0.0 1.0
*
* Generate bias-adjusted forecasts from the regression
*
prj volqmz
*
linreg(robusterrors,lags=1) vol
# constant volj
test
# 1 2
# 0.0 1.0
prj voljmz
*
* Analysis of squared forecast errors
*
set qerrorsq = qerror**2
set jerrorsq = jerror**2
*
@histogram qerrorsq
@histogram jerrorsq
*
set dd = qerrorsq-jerrorsq
@bjident(number=6,print,method=yule) dd
boxjenk(constant,ma=1) dd
*
boxjenk(constant,regressors,ma=1) vol / resids
# volq volj
linreg vol
# constant volq volj
prj volc
*
@uForeErrors vol volq
@uForeErrors vol volj
@uForeErrors vol voljmz
@uForeErrors vol volqmz
@uForeErrors vol volc


source(noecho) bjident.src
source(noecho) regcorrs.src
source(noecho) regcrits.src
source(noecho) uforecst.src
source(noecho) uforeerr.src
source(noecho) dfunit.src
*
* Trend and difference models for Yen exchange rate, pp 341-352
*
cal 1973 1 12
all 1996:7
open data exch.dat
data(format=prn,org=columns)
*
set logyen = log(yen)
set dyen   = logyen-logyen{1}
graph(header='Figure 12.7 Log Yen Exchange Rate',window='Figure 12.7')
# logyen
graph(header='Figure 12.8 Change in Log Yen Exchange Rate',window='Figure 12.8')
# dyen
*
* The option DIFFS=1 causes BJIDENT to produce the autocorrelations
* for the series and its first difference. Examining these allows you
* to decide whether to difference the series or not. If you're interested
* in the possibility of seasonal differencing (item 3 in the Problems and
* Complements section on page 361), you can use the SDIFFS option as well.
* DIFFS=1,SDIFFS=1 will give you all combinations of 0 or 1 regular differences
* and 0 or 1 seasonal differences.
*
* METHOD=YULE uses the same method of calculation for the autocorrelations that
* is used in generating the graphs in the text. The default (METHOD=BURG)
* produces rather similar results for the first differenced series,
* but has quite different partial autocorrelation estimates for the
* undifferenced series. YULE works by estimating the autocorrelations
* directly, and the partial autocorrelations from them, while BURG
* estimates the partial autocorrelations directly, then "backs out"
* the autocorrelation estimates.
*
@bjident(number=12,diffs=1,method=yule) logyen
*
* The time trend is defined as t-1973:1 (which gives it a value of 0
* in 1973:1) rather than just t, so the regressions will match up. This
* has no effect on the fit.
*
set time  = t-1973:1
*
* This does some fancy table building to do the AIC and SBC (SIC) for
* all combinations of ar and ma up to 3. Because the subscripts on
* the arrays run from 1 to 4, you have to add one to the ar and ma
* numbers to get the right slot.
*
* Since we're probably not all that interested in all the estimation
* output for 16 models (most of which will be rejected out of hand
* once we see the information criteria), the BOXJENK instruction in
* the loop uses the option NOPRINT. The information criteria are computed
* using the procedure REGCRITS. This also uses NOPRINT. All that are
* used are the values %aic and %sbc.
*
* Note that the information criteria are slightly different from those
* in the text for any model with MA terms. When there are MA terms, the
* residuals have to be generated recursively, and different programs
* use different ways to do that, yielding somewhat different estimates
* and different residuals.
*
dec rect aic(4,4) sbc(4,4)
do mas=0,3
  do ars=0,3
     boxjenk(constant,ma=mas,ar=ars,regressors,noprint) logyen * 1994:12
     # time
     @regcrits(noprint)
     compute aic(ars+1,mas+1)=%aic
     compute sbc(ars+1,mas+1)=%sbc
  end do ars
end do mas
*
* The instruction MEDIT displays windows with the two arrays of information
* criteria. The hlabels and vlabels options are used to label the rows and columns
* the way we'd like. (Otherwise, they'd be 1 to 4, which would be a bit
* confusing). A less fancy way to handle this is to use the instructions
*
*   display 'AIC' aic
*   display 'SIC' sbc
*
* which will show unadorned arrays of numbers.
*
medit(hlabels=||'0','1','2','3'||,vlabels=||'0','1','2','3'||,noedit) aic sbc
*
* Pick the best model, estimate and forecast
*
disp 'Chosen model for estimation in levels'
boxjenk(constant,ar=2,regressors,define=levelseq) logyen * 1994:12 resids
# time
@uforecast(equation=levelseq,stderrs=sefore) fyen 1995:1 1996:7
*
* Compute upper and lower 2 standard error bands
*
set upper = fyen+2.0*sefore
set lower = fyen-2.0*sefore
*
* Create the dummy variable for the forecast shading zone. Since we're going to need this
* to stretch out to 2010:12 later on, we'll take care of the whole range now.
*
set forezone * 2010:12 = t>=1995:1
graph(header='Figure 12.11. Log Yen/Dollar Rate: History and Forecast',window='Figure 12.11',shading=forezone) 4
# logyen 1990:1 1994:12
# fyen 1995:1 1996:7
# upper 1995:1 1996:7 3
# lower 1995:1 1996:7 3
*
graph(header='Figure 12.13. Log Yen/Dollar Rate: History, Forecast and Realization',window='Figure 12.13',shading=forezone) 4
# logyen 1990:1 1996:7
# fyen 1995:1 1996:7
# upper 1995:1 1996:7 3
# lower 1995:1 1996:7 3
*
* uForeErrors computes forecast error statistics on the comparison between the
* actual data (logyen) and the forecasts (fyen).
*
@uForeErrors logyen fyen
*
* Longer term forecasts. We need to extend the time trend variable out
* to the end of the new forecast sample.
*
set time * 2010:12 = (t-1973:1)
@uforecast(equation=levelseq,stderrs=sefore) fyen 1995:1 2010:12
set upper 1995:1 2010:12 = fyen+2.0*sefore
set lower 1995:1 2010:12 = fyen-2.0*sefore
graph(header='Figure 12.12. Log Yen/Dollar Rate: History and Long-Horizon Forecast',window='Figure 12.12',shading=forezone) 4
# logyen 1990:1 1994:12
# fyen 1995:1 2010:12
# upper 1995:1 2010:12 3
# lower 1995:1 2010:12 3
*
* The Dickey-Fuller unit root test is done using the procedure DFUNIT.
* The main options for this are LAGS (number of added lags) and TREND,
* which puts a time trend into the regression.
*
@dfunit(lags=3,trend) logyen * 1994:12
*
* ARMA models for differenced series. Since we want forecasts of the
* undifferenced series, we use the DIFFS option on BOXJENK and use
* logyen as the dependent variable. The parameter estimates and any
* summary statistics based upon the residuals (such as the information
* criteria) will be identical to what you would get if you used the
* differenced series as the dependent variable and eliminated the diffs
* option. If you do the estimates both ways (and don't suppress the output),
* the one difference you note will be with any summary statistics which
* use information about the dependent variable. With DIFFS=1 and the level
* variable as the dependent variable, the R**2 will use the sum of squared
* deviations of the level variable in the denominator, while if you use
* the differenced variable directly, it will be the sum of squared deviations
* of the differenced series. Since the former is generally MUCH higher than
* the latter, while the sum of squared residuals are the same, the R**2 will
* generally be quite high when the differencing is done by BOXJENK and much
* lower when you do it yourself. This doesn't mean that the former is "better"
* than the latter - the two models are the same information stated in different
* forms. It's just an important reminder that the R**2 isn't very useful in
* comparing estimates with different dependent variable.
*
* Again, we use some fancy programming to generate the table of information
* criteria.
*
do mas=0,3
  do ars=0,3
     boxjenk(constant,diffs=1,ma=mas,ar=ars,iters=100,noprint) logyen 1973:5 1994:12
     @regcrits(noprint)
     compute aic(ars+1,mas+1)=%aic
     compute sbc(ars+1,mas+1)=%sbc
  end do ars
end do mas
medit(hlabels=||'0','1','2','3'||,vlabels=||'0','1','2','3'||,noedit) aic sbc
*
disp 'Chosen model for estimation in differences'
boxjenk(constant,diffs=1,ar=1,define=diffeq) logyen * 1994:12 resids
@uforecast(equation=diffeq,stderrs=sefore) fyen 1995:1 1996:7
*
set upper 1995:1 1996:7 = fyen+2.0*sefore
set lower 1995:1 1996:7 = fyen-2.0*sefore
graph(header='Figure 12.15. Log Yen/Dollar Rate: History and Forecast. Model in Differences',window='Figure 12.15',shading=forezone) 4
# logyen 1990:1 1994:12
# fyen 1995:1 1996:7
# upper 1995:1 1996:7 3
# lower 1995:1 1996:7 3
graph(header='Figure 12.17. Log Yen/Dollar Rate: History, Forecast and Realization',window='Figure 12.17',shading=forezone) 4
# logyen 1990:1 1996:12
# fyen 1995:1 1996:7
# upper 1995:1 1996:7 3
# lower 1995:1 1996:7 3
@uForeErrors logyen fyen
*
* Longer horizon forecasts
*
@uforecast(equation=diffeq,stderrs=sefore) fyen 1995:1 2010:12
set upper 1995:1 2010:12 = fyen+2.0*sefore
set lower 1995:1 2010:12 = fyen-2.0*sefore
graph(header='Figure 12.16. Log Yen/Dollar Rate: History and Long-Horizon Forecast',window='Figure 12.16',shading=forezone) 4
# logyen 1990:1 1994:12
# fyen 1995:1 2010:12
# upper 1995:1 2010:12 3
# lower 1995:1 2010:12 3


*
* VAR example from page 260-275
*
source(noecho) bjident.src
source(noecho) crosscor.src
*
open data house.dat
cal 1968 1 12
all 1996:6
data(format=prn,org=columns)
graph(header='Figure 10.2. U.S. Housing Starts and Completions, 1968:01-1996:06',key=upright) 2
# starts
# completions
@bjident(number=24) starts 1968:1 1991:12
@bjident(number=24) completions 1968:1 1991:12
@crosscorr(number=24) starts completions 1968:1 1991:12
*
* There's an error in the computations done in the book, so
* the AIC and BIC graph here will be somewhat different
*
do lags=1,36
   system(model=varmodel)
   variables starts completions
   lags 1 to lags
   det constant
   end(system)
   estimate(noprint) 1971:1 1991:12
*
*    These are the terms in the system log likelihood which
*
   compute baseaic=%nobs*%logdet
   compute aicpenalty=2*(2*lags+1)*2.0
   compute bicpenalty=2*(2*lags+1)*log(%nobs)
   compute fixedterms=%nobs*(2+2*log(2*%pi)+2*log(%nobs)-2*log(%nobs-2*lags-1))
   set varaic lags lags = (baseaic+aicpenalty+fixedterms)/%nobs
   set varbic lags lags = (baseaic+bicpenalty+fixedterms)/%nobs
end do lags
graph(header='Figure 10.6 VAR Order Selection with AIC and BIC',nodates,key=upleft,klabels=||'AIC','SIC'||) 2
# varaic
# varbic
*
* Work with the favored system
*
system(model=varmodel)
variables starts completions
lags 1 to 4
det constant
end(system)
*
estimate * 1991:12
errors(impulses,model=varmodel,results=vdecomp,window='Variance Decomposition') * 36 %sigma
spgraph(hfields=2,vfields=2,header='Figure 10.12. Variance Decomposition',$
  ylab=||'Starts','Completions'||,xlab=||'Starts','Completions'||)
do i=1,2
   do j=1,2
      graph(min=0.0,max=1.0,nodates)
      # vdecomp(i,j)
   end do j
end do i
spgraph(done)
spgraph(vfields=2)
do i=1,2
   list ivar = 1 to 2
   graph(header=%l(%eqndepvar(i)),style=stacked,nodates,max=1.0) 2
   cards vdecomp(i,ivar)
end do i
spgraph(done)
*
* Create the dummy series for the forecast period shading
*
set forezone = t>=1992:1
*
* The single forecast instruction does all the forecasts. The 54 is the number of period
* and the 1992:1 is the starting period of forecast. RATS will actually assume the 1992:1
* start because the estimation range ended at 1991:4. The forecasts go into the vector
* of series "forecasts" where forecasts(1) are the forecasts of the first of the variables
* in the system definition (that is, "starts") and forecasts(2) are those for "completions"
*
forecast(model=varmodel,results=forecasts) * 54 1992:1
graph(header='Figure 10.13 Housing Starts: History and Forecast',shading=forezone) 2
# starts * 1991:12
# forecasts(1) 1992:1 *
graph(header='Figure 10.14 Housing Starts: History, Forecast and Realization',shading=forezone) 2
# starts * 1996:12
# forecasts(1) 1992:1 *
graph(header='Figure 10.15 Housing Completions: History and Forecast',shading=forezone) 2
# completions * 1991:12
# forecasts(2) 1992:1 *
graph(header='Figure 10.16 Housing Starts: History, Forecast and Realization',shading=forezone) 2
# completions * 1996:12
# forecasts(2) 1992:1 *


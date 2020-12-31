*
* Examples of moving average processes
* pp 145-152
*
all 150
set u = %ran(1.0)
*
*  Moving average models. We use the same random shocks for each
*
set theta40 = u+.40*u{1}
set theta95 = u+.95*u{1}
graph(klabels=||'theta=.4','theta=.95'||,key=below,$
   header='Figure 7.1 Realizations of Two MA(1) Processes') 2
# theta40
# theta95
*
* The acf's are generated including 1.0 for the "0 lag", which is how they
* are produced using actual data using the RATS instruction CORRELATE. There
* are several reasons for this format:
*
*   1. It makes the autocorrelation and autocovariance functions have similar
*      structure
*   2. The 0 lag term is used in some calculations, such as turning
*      autocorrelations into partial autocorrelations.
*
* In order to make the correlation graphs look like those in the text, the
* series are graphed over entries 2 to 16, with the "number" option used
* to relabel them as starting with 1.
*
set ac4 1 16 = %if(t==1,1.0,%if(t==2,.40/(1+.40**2),0.0))
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.2 Population Autocorrelation Function: MA(1) Process, theta=.4')
# ac4 2 16
*
set ac95 1 16 = %if(t==1,1.0,%if(t==2,.95/(1+.95**2),0.0))
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.3 Population Autocorrelation Function: MA(1) Process, theta=.95')
# ac95 2 16
*
* The acf2pacf procedure generates a set of partial autocorrelations from an input
* autocorrelation function. With real data, you usually just use the CORRELATE
* instruction to compute both autocorrelations and partials, but here we have
* a theoetical ACF for which we want the PACF.
*
source(noecho) acf2pacf.src
@acf2pacf ac4 pac4
@acf2pacf ac95 pac95
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.4 Population Partial Autocorrelation Function: MA(1) Process, theta=.4')
# pac4 2 16
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.5 Population Partial Autocorrelation Function: MA(1) Process, theta=.95')
# pac95 2 16


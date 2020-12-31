*
* Examples of AR processes
* pp 152-160
*
all 150
set u = %ran(1.0)
*
* Generating AR processes is a bit trickier than MA's. An MA
* process only has a "memory" of the number of periods for which
* there are non-zero MA lags. For the MA(1) processes, we have
* the process exact for entries 2 and up by just taking the two
* term moving average of the u process. An AR, however, depends
* upon a lagged value of itself. Where do we get a time=0 value
* to compute time=1? There are two ways to handle this in practice:
*
* 1. Figure out the unconditional distribution for time=0 and
*    make a draw for that. That's what we do below: if u is N(0,s**2)
*    and y(t)=ay(t-1)+u(t), then y (unconditionally, that is, if we
*    have no information about earlier values of y) is N(0,s**2/(1-a**2))
* 2. Start with y(0)=0 (the unconditional) mean, and generate the process.
*    Discard or otherwise ignore the first few terms (50 is a good safe number),
*    which aren't really representative of the AR process. The initial draws
*    which you discard are known as the "burn-in."
*
set ar4 = %if(t==1,u/(1-.4**2),.4*ar4{1}+u)
set ar95 = %if(t==1,u/(1-.95**2),.95*ar95{1}+u)
graph(header='Figure 7.6 Realizations of Two AR(1) Processes',$
 klabels=||'phi=.4','phi=.95'||,key=below) 2
# ar4
# ar95
*
* Unlike the MA model, both the acf and pacf have closed form expressions
*
set ac4 1 16 = .4**(t-1)
set ac95 1 16 = .95**(t-1)
set pac4 1 16 = %if(t==1,1.0,%if(t==2,.4,0.0))
set pac95 1 16 = %if(t==1,1.0,%if(t==2,.95,0.0))
*
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.7 Population Autocorrelation Function: AR(1) Process, phi=.4')
# ac4 2 16
*
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.8 Population Autocorrelation Function: AR(1) Process, phi=.95')
# ac95 2 16
*
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.9 Population Partial Autocorrelation Function: AR(1) Process, phi=.4')
# pac4 2 16
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.10 Population Partial Autocorrelation Function: AR(1) Process, phi=.95')
# pac95 2 16




*
* Example of AR(2) process
* pp 159-160
*
all 200
*
* As described in diebp152.prg, this example uses the "burn-in"
* procedure for generating the AR process. t=1 and t=2 are set
* to zero, after which the AR process is put into effect. The first
* 50 terms generated are ignored in the graph of the realization.
*
set u = %ran(1.0)
set y = %if(t<=2,0.0,1.5*y{1}-.9*y{2}+u)
graph(header='Realization of y(t)=1.5y(t-1)-.9y(t-2)+u')
# y 51 200
*
* The acf has two "non-standard" terms followed by a formula. You
* could do this with a single SET and two %IF functions, but it's
* easier to read if you just handle the cases with separate SET
* instructions.
*
set acar2 1 1 = 1.0
set acar2 2 2 = 1.5/(1-(-.9))
set acar2 3 32 = 1.5*acar2{1}-.9*acar2{2}
graph(style=bargraph,max=1.0,min=-1.0,number=1,$
  header='Figure 7.11 Population Autocorrelation Function: AR(2) Process, Complex Roots')
# acar2 2 32




*
* ARMA models for Canadian employment data
* pp 163-173
*
source(noecho) bjest.src
cal 1961 1 4
all 1994:4
open data caemp.dat
data(format=prn,org=columns)
*
*  This loop runs MA(1) through MA(4) models. The BJEST procedure graphs
*  the residual autocorrelations and includes the Q stat, AIC and BIC (which
*  is called the SBC) on the graph.
*
do i=1,4
   @bjest(ma=i,constant,number=20) caemp 1962:1 1993:4
end do i
*
*  This loop runs AR(1) through AR(4)
*
do i=1,4
   @bjest(ar=i,constant,number=20) caemp 1962:1 1993:4
end do i
*
*  This estimates the favored model
*
@bjest(ar=3,ma=1,constant,number=20) caemp 1962:1 1993:4
*
*  The table on page 170 requires estimation of a lot of models. In
*  order to do this, we do a pair of nested loops over the ar and ma
*  parameters. These are put into a 5x5 array. Note that, because
*  matrix subscripts have to start at 1, we have to add 1 to the i
*  and j to get the proper slot.
*
dec rect aic(5,5) bic(5,5)
do i=0,4
   do j=0,4
      @bjest(ar=i,ma=j,constant,noprint,nograph) caemp 1962:1 1993:4
      compute aic(i+1,j+1)=%aic
      compute bic(i+1,j+1)=%sbc
   end do j
end do i
*
*  This does a fairly simple display of the matrix of information criteria
*
disp 'AIC Statistics' aic
disp 'BIC Statistics' bic
*
*  This displays them in a better formatted report window
*
medit(window='AIC',hlabels=||'0','1','2','3','4'||,vlabels=||'0','1','2','3','4'||,noedit) aic
medit(window='BIC',hlabels=||'0','1','2','3','4'||,vlabels=||'0','1','2','3','4'||,noedit) bic



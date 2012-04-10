#include "minpack.h"
#include <math.h>

/*
  subroutine fdjac2 

  this subroutine computes a forward-difference approximation 
  to the m by n jacobian matrix associated with a specified 
  problem of m functions in n variables. 

  the subroutine statement is 

  subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa) 

  where 

  fcn is the name of the user-supplied subroutine which 
  calculates the functions. fcn must be declared 
  in an external statement in the user calling 
  program, and should be written as follows. 

  subroutine fcn(m,n,x,fvec,iflag) 
  integer m,n,iflag 
  double precision x(n),fvec(m) 
  ---------- 
  calculate the functions at x and 
  return this vector in fvec. 
  ---------- 
  return 
  end 

  the value of iflag should not be changed by fcn unless 
  the user wants to terminate execution of fdjac2. 
  in this case set iflag to a negative integer. 

  m is a positive integer input variable set to the number 
  of functions. 

  n is a positive integer input variable set to the number 
  of variables. n must not exceed m. 

  x is an input array of length n. 

  fvec is an input array of length m which must contain the 
  functions evaluated at x. 

  fjac is an output m by n array which contains the 
  approximation to the jacobian matrix evaluated at x. 

  ldfjac is a positive integer input variable not less than m 
  which specifies the leading dimension of the array fjac. 

  iflag is an integer variable which can be used to terminate 
  the execution of fdjac2. see description of fcn. 

  epsfcn is an input variable used in determining a suitable 
  step length for the forward-difference approximation. this 
  approximation assumes that the relative errors in the 
  functions are of the order of epsfcn. if epsfcn is less 
  than the machine precision, it is assumed that the relative 
  errors in the functions are of the order of the machine 
  precision. 

  wa is a work array of length m. 

  p is a general-purpose pointer available in fcn. 

  argonne national laboratory. minpack project. march 1980. 
  burton s. garbow, kenneth e. hillstrom, jorge j. more 

*/

int fdjac2_(S_fp fcn, int m, int n, double *x, 
	    double *fvec, double *fjac, int ldfjac, int *iflag, 
	    double epsfcn, double *wa, void *p)
{
    /* System generated locals */
    int fjac_dim1, fjac_offset;

    double h;
    int i, j;
    double eps, temp, epsmch;

    /* Parameter adjustments */
    --wa;
    --fvec;
    --x;
    fjac_dim1 = ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;

    /* epsmch is the machine precision */
    epsmch = dpmpar_(1);
    eps = sqrt(max(epsfcn, epsmch));

    for (j = 1; j <= n; ++j) {
	temp = x[j];
	h = eps * fabs(temp);
	if (h == 0.0) {
	    h = eps;
	}
	x[j] = temp + h;
	(*fcn)(m, n, &x[1], &wa[1], iflag, p);
	if (*iflag < 0) {
	    return 0;
	}
	x[j] = temp;
	for (i = 1; i <= m; ++i) {
	    fjac[i + j * fjac_dim1] = (wa[i] - fvec[i]) / h;
	}
    }

    return 0;
}


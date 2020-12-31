*
* Autocorrelations of random walk and its differences, pp 334-5
*
source(noecho) bjident.src
all 300
set(first=1.0) y = y{1}+%ran(1.0)
*
* The BJIDENT procedure (which has been used before) also includes
* a DIFFS option to indicate the maximum number of differences to
* examine. With DIFFS=1, it will give (separate) graphs of
* correlations for both the levels (undifferenced) and first
* differenced series.
*
@bjident(diffs=1,number=12) y


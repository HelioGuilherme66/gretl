'''(Mostly time-series-related) functions needed and written by Sven Schreiber.

This is free but copyrighted software, distributed under the same license terms
(as of January 2007) as the 'gretl' program by Allin Cottrell and others, see
gretl.sf.net (in short: GPL v2, see www.gnu.org/copyleft/gpl.html).

(see end of this file for a changelog)
'''
from numpy import r_, c_, arange, diff, mean, sqrt, log, mat
from numpy import asarray, nan
from numpy.matlib import ones, zeros, rand, eye, empty
from numpy.linalg import eigh, cholesky, solve, lstsq
# (lstsq also as tool to determine rank)

# some constants/dictionaries first
quarter2month = {1: 1, 2: 4, 3: 7, 4: 10}
# in theory we only need the four months 1, 4, 7, 10, but well...
month2quarter = {1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 2, 7: 3, 8: 3, 9: 3, \
                 10: 4, 11: 4, 12: 4}
qNumber2qFloat = {1: 0.0, 2: 0.25, 3: 0.5, 4: 0.75}
mNumber2mFloat = {1: 0.0, 2: 0.0833, 3: 0.1666, 4: 0.2499, 5: 0.3332, \
                  6: 0.4165, 7: 0.4998, 8: 0.5831, 9: 0.6664, 10: 0.7497, \
                  11: 0.8330, 12: 0.9163}
qFracstring2qString = {0.0: 1, 0.25: 2, 0.5: 3, 0.75: 4}
mFloat2mNumber = {0.0: 1, 0.0833: 2, 0.1666: 3, 0.2499: 4, 0.3332: 5, \
                  0.4165: 6, 0.4998: 7, 0.5831: 8, 0.6664: 9, 0.7497: 10, \
                  0.8330: 11, 0.9163: 12}
# with onetwelfth == 0.0833 approx.


def vec(m):
    '''
    Returns all columns of the input as a stacked (column) vector.

    If m is a numpy-array, a 1d-array is returned. For a numpy-matrix m,
     the output has shape (n*m, 1).
    '''
    return m.T.ravel().T

from numpy import mat, asarray
def unvec(m, rows, cols):
    '''
    Turns (column) vector into matrix of shape rows, cols.

    Also accepts 1d-array input, but always returns numpy matrix.
    '''
    if type(m) == type(mat(m)):
        assert m.shape[1] == 1                            # col vector
        intype = 'matrix'
    else:
        assert len(m.shape) == 1                          # 1d array
        intype = 'array'
        m = mat(m).T
    assert cols * rows == m.shape[0]
    out = m.reshape(cols, rows).T
    if intype == 'array': return asarray(out)
    else: return out

from numpy import mat
def mat2gretlmatstring(m):
    '''
    Turns numpy matrix or array (or scalar!) m into gretl string representation.
    '''
    # mat(m) is necessary because if m is 1d-array, map() would fail
    out =  ';'.join( [  ','.join(map(str, row)) for row in mat(m).tolist() ] )
    return '{' + out + '}'

def startobs2obslist(startperiod, numofobs):
    '''
    Constructs list of observation labels following the input pattern.

    Example:
    startperiod = '1999q3', numofobs = 2 -> ['1999q3', '1999q4']
    Currently supports only annual (pure number), monthly, quarterly.
    Years must be in 4-digit format.
    '''
    if startperiod.isdigit():           # pure (integer) number
        startnumber = int(startperiod)
        return [ str(startnumber + ix) for ix in range(numofobs) ]
    elif startperiod[4] in 'qQ':        # quarterly dates
        wrap = 4
        period = int(startperiod[5])
    elif startperiod[4] in 'mM':
        wrap = 12
        period = int(startperiod[5:7])
    else: raise NotImplementedError

    year = int(startperiod[:4])
    out = [str(year) + startperiod[4] + str(period)]
    for ix in range(numofobs):
        if period == wrap:
            period = 1
            year += 1
        else: period += 1
        out.append(str(year) + startperiod[4] + str(period))

    return out

import csv
from numpy import mat
def writecsv(filename, data, orientation = 'cols', delim = ',', \
     varnames = [],  obslabels = [], comments = [], commentchar = '# '):
    '''
    Saves array or matrix <data> in csv format in file <filename> (path string).

    <comments> must be passed as a sequence of strings, one for each line,
     and will be written at the top of the file, each line starting with
     <commentchar>.
    <orientation> can be 'cols' or 'rows', determines whether the
     variable names will be used as column or row headers, and how to treat
     1d-input. (And observation labels will be written accordingly.)
    <varnames> and <obslabels> must be sequences of strings.
    '''
    data = mat(data)
    if orientation == 'rows':
        colheaders = obslabels
        rowheaders = varnames
        cell11 = 'var'
    else:                           # 'cols' orientation as fallback
        colheaders = varnames
        rowheaders = obslabels
        cell11 = 'obs'
        if data.shape[0] == 1: data = data.T    # make 1d-array a column vector
    if len(colheaders) > 0: assert len(colheaders) == data.shape[1]

    # start writing to the file
    target = csv.writer(open(filename, 'w'), delimiter = delim)
    target.writerows([ [commentchar + comment] for comment in comments])
    # (additional trivial list layer because otherwise the comment string itself
    #  would be split up with the delim character)
    if len(rowheaders) > 0:
        assert len(rowheaders) == data.shape[0]
        target.writerow(colheaders.insert(0, cell11))
    else: target.writerow(colheaders)
    temp = data.tolist()        # temp to have list conversion only once
    for ix in range(len(rowheaders)): temp[ix].insert(0, rowheaders[ix])
    target.writerows(temp)

    return 0            # success

import csv
from numpy import mat
def readcsv(filename, delim = ',', commentchar = '#', colheader = 'names', \
        rowheader = 'obs'):
    '''
    Read in a csv file (may contain comments starting with commentchar).

    The contents of the first non-comment row and column must be indicated in
     rowheader and colheader as one of 'names', 'obs' (labels), or None.
    The array (matrix) of the data is returned as is, i.e. w/o transpose, hence
     the caller must know whether variables are in rows or columns.
    If both colheader and rowheader are not None, the upper-left cell (header
     of the first row/col) is ignored (but must be non-empty).

    Returns a five-element tuple:
    0. numpy-matrix of the actual data as floats
    1. orientation of variables: 'cols', 'rows', or 'unknown'
    2. 1d-array of variable names (or None)
    3. 1d-array of observation labels (or None)
    4. the type/frequency of the data
        (currently one of 'a', 'q', 'm', guessed from the first date label)
        (if this deduction failed, 'unknown' is returned here)

    Easiest example with upper-left data cell in second row/second column:
    mydata = readcsv('myfile.csv')[0]
    '''
    read_from = csv.reader(open(filename, 'rb'), delimiter = delim, \
        skipinitialspace = True)
    tempnestedlist = [ line for line in read_from if not \
        line[0].strip().startswith(commentchar) ]
    data = mat(tempnestedlist, dtype = str)

    if colheader == 'names':
        orientation = 'cols'
        varnames, data = data[0, :].A1, data[1:, :]
        if rowheader == 'obs':
            obslabels, data = data[:, 0].A1, data[:, 1:]
            varnames = varnames[1:]
    elif rowheader == 'names':
        orientation = 'rows'
        varnames, data = data[:, 0].A1, data[:, 1:]
        if colheader == 'obs':
            obslabels, data = data[0, :].A1, data[1:, :]
            varnames = varnames[1:]
    elif colheader == 'obs':
        orientation = 'rows'
        obslabels, data = data[0, :].A1, data[1:, :]
        if rowheader == 'names':
            varnames, data = data[:, 0].A1, data[:, 1:]
            obslabels = obslabels[1:]
    elif rowheader == 'obs':
        orientation = 'cols'
        obslabels, data = data[:, 0].A1, data[:, 1:]
        if colheader == 'names':
            varnames, data = data[0, :].A1, data[1:, :]
            obslabels = obslabels[1:]
    else:
        assert colheader == None        # to catch typos, e.g. 'Names', 'OBS'
        assert rowheader == None
        orientation = 'unknown'
        varnames = None
        obslabels = None

    # detect the dataset type:
    # annual
    if len(obslabels[0]) == 4: freq = 'a'
    # quarterly
    elif len(obslabels[0]) == 6 and obslabels[0][4] in 'qQ': freq = 'q'
    # monthly
    elif len(obslabels[0]) == 7 and obslabels[0][4] in 'mM': freq = 'm'
    else: freq = 'unknown'

    return data.astype(float), orientation, varnames, obslabels, freq

from numpy import nan
def floatAndNanConverter(datapoint, nacode = 'na'):
    '''
    Converts nacode to numpy.nan value.

    Also returns other input as float (e.g. for matplotlib's load, asarray).
    '''
    if datapoint == nacode: return nan
    return float(datapoint)

def dateString2dateFloat(datestring):
    '''
    Converts '1999q2' -> 1999.25, '1999m2' -> 1999.0833, etc.

    So far only for quarterly and monthly.
    '''
    year, freq = float(datestring[:4]), datestring[4]
    assert freq in 'qQmM', 'sorry, only quarterly or monthly'
    if freq in 'qQ':    #quarterly
        result = year + qNumber2qFloat[int(datestring[5])]
    elif freq in 'mM':               #monthly
        result = year + mNumber2mFloat[int(datestring[5:7])]
    return result

from datetime import date, timedelta
def getQuarterlyDates(startyear, startquarter, t):
    '''
    Constructs a list of quarterly date labels for t obs.

    Algorithm to get a sequence of strings relating to quarterly dates:
     1. start with first day in the startquarter, e.g. 2006-04-01
     2. map the month to quarter and make string year + 'q' + quarter
     3. the longest quarters are 3rd and 4th (2*31 days + 30 days = 92 days),
        1st the shortest (90 or 91), so add a timedelta (in days,
        apparently default) of 100 days (anything between 92+1 and
        sum of shortest quarter plus one month = approx. 118)
     4. reset the day of that intermediate date to 1
     5. return to step 2
    '''
    try:
        y = int(startyear); q = int(startquarter); t = int(t)
    except: raise TypeError, 'need integers for year, quarter, t'
    if q not in range(1,5): raise ValueError, 'startquarter input out of range'
    # create list for date strings:
    datestrings = []
    # step 1.:
    d = date(y, quarter2month[startquarter], 1)
    for t in range(t):
        datestrings.append(str(d.year) + 'Q' + str(month2quarter[d.month]))
        d += timedelta(100)
        d = d.replace(day = 1)
    return datestrings

from numpy import mat, asarray
from numpy.matlib import empty, zeros, eye
from numpy.linalg import lstsq
def getOrthColumns(m):
    '''
    Constructs the orthogonally complementing columns of the input.

    Input of the form pxr is assumed to have r<=p,
    and have either full column rank r or rank 0 (scalar or matrix)
    Output is of the form px(p-r), except:
    a) if M square and full rank p, returns scalar 0
    b) if rank(M)=0 (zero matrix), returns I_p
    (Note you cannot pass scalar zero, because dimension info would be
    missing.)
    Return type is as input type.
    '''
    if type(m) == type(asarray(m)):
        m = mat(m)
        output = 'array'
    else: output = 'matrix'
    p, r = m.shape
    # first catch the stupid input case
    if p < r: raise ValueError, 'need at least as many rows as columns'
    # we use lstsq(M, ones) just to exploit its rank-finding algorithm,
    rk = lstsq(m, ones(p).T)[2]
    # first the square and full rank case:
    if rk == p: return 0
    # then the zero-matrix case (within machine precision):
    elif rk == 0: result = eye(p)
    # now the rank-deficient case:
    elif rk < r:
        raise ValueError, 'sorry, matrix does not have full column rank'
    # (what's left should be ok)
    else:
        # we have to watch out for zero rows in M,
        # if they are in the first p-r positions!
        # so the (probably inefficient) algorithm:
            # 1. check the rank of each row
            # 2. if zero, then also put a zero row in c
            # 3. if not, put the next unit vector in c-row
        idr = eye(r)
        idpr = eye(p-r)
        c = empty([0,r])    # starting point
        co = empty([0, p-r]) # will hold orth-compl.
        idrcount = 0
        for row in range(p):
            # (must be ones() instead of 1 because of 2d-requirement
            if lstsq( m[row,:], ones(1) )[2] == 0 or idrcount >= r:
                c = r_[ c, zeros(r) ]
                co = r_[ co, idpr[row-idrcount, :] ]
            else:     # row is non-zero, and we haven't used all unit vecs
                c = r_[ c, idr[idrcount, :] ]
                co = r_[ co, zeros(p-r) ]
                idrcount += 1
        # earlier non-general (=bug) line: c = mat(r_[eye(r), zeros((p-r, r))])
        # and:  co = mat( r_[zeros((r, p-r)), eye(p-r)] )
        # old:
        # result = ( eye(p) - c * (M.T * c).I * M.T ) * co
        result = co - c * solve(m.T * c, m.T * co)
    if output == 'array': return asarray(result)
    else: return result

from numpy import mat, asarray
def addLags(m, maxlag):
    '''
    Adds (contiguous) lags as additional columns to the TxN input.

    Early periods first. If maxlag is zero, original input is returned.
    maxlag rows are deleted (the matrix is shortened)
    '''
    if type(m) == type(asarray(m)):
        m = mat(m)
        output = 'array'
    else: output = 'matrix'
    T, N = m.shape
    if type(maxlag) != type(4):
        raise TypeError, 'addLags: need integer for lag order'
    if maxlag > m.shape[0]:
        raise ValueError, 'addLags: sample too short for this lag'
    temp = m[ maxlag: ,:]  # first maxlag periods must be dropped due to lags
    for lag in range(1, maxlag + 1) :
        temp = c_[ temp, m[(maxlag-lag):(T-lag) ,:] ]
    if output == 'array': return asarray(temp)
    else: return temp

from numpy.matlib import empty, ones, zeros
from numpy import mat, c_, r_
def getDeterministics(nobs, which = 'c', date = 0.5):
    '''
    Returns various useful deterministic terms for a given sample length T.

    Return object is a numpy-matrix-type of dimension Tx(len(which));
    (early periods first, where relevant).
    In the 'which' argument pass a string composed of the following letters,
    in arbitrary order:
    c - constant (=1) term
    t - trend (starting with 0)
    q - centered quarterly seasonal dummies (starting with 0.75, -0.25...)
    m - centered monthly seasonal dummies (starting with 11/12, -1/12, ...)
    l - level shift (date applies)
    s - slope shift (date applies)
    i - impulse dummy (date applies)

    If the date argument is a floating point number (between 0 and 1),
    it is treated as the fraction of the sample where the break occurs.
    If instead it is an integer between 0 and T, then that observation is
    treated as the shift date.
    '''
    # some input checks (as well as assignment of shiftperiod):
    if type(nobs) != type(4):  # is not an integer
        raise TypeError, 'need integer for sample length'
    if nobs <=0: raise ValueError, 'need positive sample length'
    if type(date) == type(0.5):     #is a float, treat as break fraction
        if date < 0 or date > 1:
            raise ValueError, 'need break fraction between 0 and 1'
        shiftperiod = int(date * nobs)
    elif type(date) == type(4):     # is integer, treat as period number
        if date not in range(1, nobs+1):
            raise ValueError, 'need period within sample range'
        shiftperiod = date
    else: raise TypeError, 'need float or integer input for date'
    if type(which) != type('a string'):
        raise TypeError, 'need string for case spec'
    # end input checks

    out = empty([nobs,0])   # create starting point
    if 'c' in which: out = c_[ out, ones(nobs).T ]
    if 't' in which: out = c_[ out, r_['c', :nobs] ]
    if 'l' in which:
        shift = r_[ zeros(shiftperiod).T, ones(nobs-shiftperiod).T ]
        out = c_[ out, shift ]
    if 's' in which:
        slopeshift = r_[ zeros(shiftperiod).T, r_['c', 1:(nobs - shiftperiod + 1)] ]
        out = c_[ out, slopeshift ]
    if 'i' in which:
        impulse = r_[ zeros(shiftperiod).T, ones(1), zeros(nobs-shiftperiod-1).T ]
        out = c_[ out, impulse ]
    if 'q' in which or 'Q' in which:
        # to end of next full year, thus need to slice at T below:
        q1 = [0.75, -0.25, -0.25, -0.25] * (1 + nobs/4)
        q2 = [-0.25, 0.75, -0.25, -0.25] * (1 + nobs/4)
        q3 = [-0.25, -0.25, 0.75, -0.25] * (1 + nobs/4)
        out = c_[ out, mat(q1[:nobs]).T, mat(q2[:nobs]).T, mat(q3[:nobs]).T ]
    if 'm' in which or 'M' in which:
        temp = [-1./12] * 11
        for month in range(11):
            temp.insert(month, 1-temp[0])
            # again, to end of next full year, thus need to slice at T below:
            monthly = temp * (1 + nobs/12)  # temp is still a list here!
            out = c_[ out, mat(monthly[:nobs]).T ]
    return out

from numpy.matlib import empty
def getImpulseDummies(sampledateslist, periodslist):
    '''
    Returns a (numpy-)matrix of impulse dummies for the specified periods.

    sampledateslist must consist of 1999.25 -style dates (quarterly or monthly).
    However, because periodslist is probably human-made, it expects strings
     such as '1999q3' or '1999M12'.
    Variables in columns.
    So far only for quarterly and monthly data.
    '''
    nobs = len(sampledateslist)
    result = empty([nobs,0])
    for periodstring in periodslist:
        period = dateString2dateFloat(periodstring)
        result = c_[result, getDeterministics(nobs, 'i', \
                            sampledateslist.index(period))]
    return result

from numpy import mat, asarray
from numpy.linalg import cholesky, eigh
def geneigsympos(A, B):
    ''' Solves symmetric-positive-def. generalized eigenvalue problem Az=lBz.

    Takes two real-valued symmetric matrices A and B (B must also be
    positive-definite) and returns the corresponding (also real-valued)
    eigenvalues and eigenvectors.

    Return format: as in scipy.linalg.eig, tuple (l, Z); l is taken from eigh
    output (a 1-dim array of length A.shape[0] ?) ordered ascending, and Z is
    an array or matrix (depending on type of input A) with the corresponding
    eigenvectors in columns (hopefully).

    Steps:
        1. get lower triang Choleski factor of B: L*L.T = B
         <=> A (LL^-1)' z = l LL' z
         <=> (L^-1 A L^-1') (L'z) = l (L'z)
        2. standard eig problem, with same eigvals l
        3. premultiply eigvecs L'z by L^-1' to get z
    '''
    output = 'matrix'
    if type(A) == type(asarray(A)):
        output = 'array'
        A, B = mat(A), mat(B)
    # step 1
    LI = cholesky(B).I
    # step 2
    evals, evecs = eigh(LI * A * LI.T)
    # sort
    evecs = evecs[:, evals.argsort()]
    evals.sort()        # in-place!
    # step 3
    evecs = LI.T * evecs
    if output == 'array': return evals, asarray(evecs)
    else:   return evals, evecs

from numpy.matlib import eye, c_
def vecm2varcoeffs(gammas, maxlag, alpha, beta):
    '''
    Converts Vecm coeffs to levels VAR representation.

    Gammas need to be coeffs in shape #endo x (maxlag-1)*#endo,
    such that contemp_diff = alpha*ect + Gammas * lagged_diffs
    is okay when contemp_diff is  #endo x 1.
    We expect matrix input!
    '''
    if alpha.shape != beta.shape:   # hope this computes for tuples
        raise ValueError, 'alpha and beta must have equal dim'
    N_y = alpha.shape[0]
    if beta.shape[0] != N_y:
        raise ValueError, "alpha or beta dim doesn't match"
    if gammas.shape[0] != N_y:
        raise ValueError, "alpha or gammas dim doesn't match"
    if gammas.shape[1] != (maxlag-1)*N_y:
        raise ValueError, "maxlag or gammas dim doesn't match"

    # starting point first lag:
    levelscoeffs = eye(N_y) + alpha * beta.T + gammas[ : , :N_y ]
    # intermediate lags:
    for lag in range(1, maxlag-1):
        levelscoeffs = c_[ levelscoeffs, gammas[:, N_y*lag : N_y*(lag+1)] - \
                          gammas[:,  N_y*(lag-1) : N_y*lag ] ]
    # last diff-lag, now this should be N_y x maxlags*N_y:
    return c_[ levelscoeffs, -gammas[:, -N_y: ] ]

def gammas2alternativegammas(gammas, alpha, beta):
    '''
    Converts Vecm-coeffs for ect at t-1 to the ones for ect at t-maxlag.

    The input gammas (shortrun coeffs) refer to a Vecm where the levels are
     lagged one period. In the alternative representation with the levels
     lagged maxlag periods the shortrun coeffs are different; the relation is:
         alt_gamma_i = alpha * beta' + gamma_i

    Actually with numpy's broadcasting the function is a one-liner so this here
     is mainly for documentation and reference purposes.
    In terms of the levels VAR coefficients A_i (i=1..maxlag) the gammas are
     defined as:
         gamma_i = - \sum_{j=i+1)^maxlag A_j for i=1..maxlag-1;
     and the alternative gammas (used py Proietti e.g.) are:
         alt_gamma_i = -I + \sum_{j=1}^i A_j for i=1..maxlag-1.
     (And \alpha \beta' = -I + \sum_{j=1}^maxlag A_j.)
    '''
    # use broadcasting to do the summation in one step:
    return  alpha * beta.T + gammas

################################
## now some more econometrically oriented helper functions
################################

from numpy import mat, asarray
from numpy.matlib import zeros
def autocovar(series, LagInput, Demeaned=False):
    '''
    Computes the autocovariance of a uni- or multivariate time series.

    Usage: autocovar(series, Lag [, Demeaned=False]) returns the NxN
    autocovariance matrix (even for N=1), where series is
    an TxN matrix holding the N-variable T-period data (early periods first),
    and Lag specifies the lag at which to compute the autocovariance.
    Specify Demeaned=True if passing ols-residuals to avoid double demeaning.
    Returns a numpy-matrix-type.
    '''
    if type(series) == type(asarray(series)):
        output = 'array'
        series = mat(series)
    else: output = 'matrix'
    t, n = series.shape
    try: Lag = int(LagInput)
    except: raise TypeError, 'autocovar: nonsense lag input type'
    if Demeaned == False:
        # axis=0 for columns (otherwise does overall-average):
        xbar = series.mean(axis=0)
    else: xbar = 0              # seems to broadcast to vector-0 ok (below)
    result = zeros([n,n])
    for tindex in range(Lag, t):
        xdev1 = series[tindex,:] - xbar
        xdev2 = series[tindex-Lag, :] - xbar
        result += xdev1.T * xdev2
    result /= t
    if output == 'array': return asarray(result)
    else: return result

from numpy import mat, asarray
from numpy.matlib import zeros
def longrunvar(series, Demeaned = False, LagTrunc = 4):
    '''
    Estimates the long-run variance (aka spectral density at frequency zero)
    of a uni- or multivariate time series.

    Usage: lrv = longrunvar(series [, Demeaned, LagTrunc]),
    where series is a TxN matrix holding
    the N-variable T-period data (early periods first).
    The Bartlett weighting function is used
    up to the specified lag truncation (default = 4).
    Specify Demeaned=True when passing Ols-residuals etc. (default False).
    Returns an NxN matrix (even for N=1).
    '''
    if type(series) == type(asarray(series)):
        output = 'array'
        series = mat(series)
    else: output = 'matrix'
    t, n = series.shape

    # set the lag window constant:
    try: Lag = int(LagTrunc)
    except: raise TypeError, 'longrunvar: nonsense lag input type'
    if Lag >= t-1:
        Lag = int(sqrt(t))
        print 'longrunvar warning: not enough data for chosen lag window'
        print '(was ', LagTrunc, ', reset to ', Lag, ')'

    result = zeros([n,n])
    for tau in range(1, Lag+1):
        Gamma = autocovar(series, tau, Demeaned)    # numpy-matrix here
        #the positive and negative range together:
        result += (1-tau/(Lag+1)) * (Gamma + Gamma.T)
    # add the tau=0 part:
    result +=  autocovar(series, 0, Demeaned)
    if output == 'array': return asarray(result)
    else: return result

from numpy import mat
from numpy.matlib import ones, zeros
from numpy.linalg import solve
def commontrendstest(series, LagTrunc=4, determ = 'c', breakpoint=0.5):
    '''
    The Nyblom&Harvey(2000)-type tests for K_0 against K>K_0
    common stochastic trends in time series.

    Usage: commontrendstest(series [, LagTrunc, Deterministics, breakpoint])
    returns a 1d N-array with the test statistics
    (partial sums of relevant eigenvalues),
    starting with the null hypothesis K_0=N-1 and ending with K_0=0.
    Pass a TxN array of data in series (early periods first).
    Optional:
    1)
    A number in LagTrunc to influence the lag window that is used
    for the implicit estimation of the long-run covariance matrix (default 4).
    2)
    Specify the type of deterministics; default is a constant mean,
    but you can set 't' to automatically de-trend the data (linearly),
    or use one of the following models with (one-time) deterministic shifts
    (see Busetti 2002):
    '1' (a string!) for a level shift w/o trend,
    '2' for a model with breaks in the mean and the trend slope,
    '2a' for a trend model where only the mean shifts.
    (Case 2b --broken trends with connected segments-- will not be implemented.)
    3) The relative breakpoint in the sample can be chosen for the shift-cases
    (otherwise it is ignored); it defaults to 0.5.
    '''
    series = mat(series)
    t, n = series.shape
    try: Lag = int(LagTrunc)
    except: raise TypeError, 'commontrendstest: nonsense lag input type'
    if Lag <= 0:
        print 'commontrendstest warning: lag trunc too small, set to default!'
        Lag = 4
    if type(breakpoint) != type(0.5):   # check for floating point input
        raise TypeError, 'commontrendstest: nonsense breakpoint input type'
    elif (breakpoint <= 0) or (breakpoint >= 1):
        raise ValueError, 'commontrendstest: breakpoint not in unit interval'
    # (end check input)
    if determ == 'c': D = ones(t).T
    elif determ == 't': D = getDeterministics(t, 'ct')
    elif determ == '1': D = getDeterministics(t, 'cl', breakpoint)
    elif determ == '2': D = getDeterministics(t, 'ctls', breakpoint)
    elif determ == '2a': D = getDeterministics(t, 'ctl', breakpoint)

    # okay, now remove the deterministics:
    # (by now, D should be Tx(1,2,3, or 4) )
    # this should do the projection:
    Resid = series - D * solve(D.T * D, D.T * series)
    Cmat = zeros([n,n])
    for i in range(t):
        temp = zeros((1,n))
        for tindex in range(i):
            temp += Resid[tindex,:]
        Cmat += temp.T * temp
    Cmat /= t**2
    Sm = longrunvar(Resid, True, Lag)
    # (True for data w/o deterministics, because everything removed)

    # generalized eigenvalues, corresponding to det(Cmat- lambda_j Sm)=0
    try: evals = geneigsympos(Cmat, Sm)[0]
    except:
        # most probably Sm wasn't pos-def, which can happen depending on lags,
        # then we try to switch over to scipy's more general eigenvalues
        print Sm    #to get some insight before everything dies
        from scipy import linalg as sl
        evals = sl.eigvals(Cmat, Sm)
        evals.sort()        # in-place!
    # default axis in cumsum works here:
    return evals.cumsum()

#############################################################

'''
Changelog:
15Jan2007:
    new unvec() function
11Jan2007:
    new writecsv() function,
    deleted writeGpl...,
    new startobs2obslist(),
    new vec() function
10Jan2007:
    fixed use of c_ / r_ due to change in numpy API,
    fix bug in readcsv
7Jan2007:
    rewrote input checks using assert,
    generalized readcsv (formerly known as readgretlcsv)
5Jan2007:
    explicit sorting of eigenvalues instead of relying on numpy implementation
3Jan2007:
    new and simpler readgretlcsv (after gretl cvs changes),
    converter for numpy matrix to gretl-type matrix string
17Aug2006:
    fixes for readGplFile,
    finished writeGplFile
16Aug2006:
    removed obsolete qString2qNumber(),
    rewrote readGplFile with csv module, scipy or matplotlib not required,
    started analogous writeGplFile
15Aug2006:
    minor cosmetics
12Aug2006:
    added readGplFile,
    added getImpulseDummies
11Aug2006:
    added helpers for use with matplotlib and or gpl-formatted csv files,
    renamed getDetermMatrix to getDeterministics
10Aug2006:
    commented out diagm, can use .diagonal() and diagflat() in numpy
21Jul2006:
    commented out zerosm, onesm, emptym, eyem, randm, which are obsoleted
     by the new numpy.matlib,
    what about diag?: still needs to be fixed in numpy,
    tried to avoid inefficient inverses (and use solve instead),
    replace asm/asmatrix by mat which now means the same in numpy,
    try to make makeNumpyMatrix redundant,
2Jun2006:
    added helpers zerosm, onesm, emptym, eyem, diagm to return numpy-matrices,
    added helper randm including workaround,
    switched to using ' instead of " where possible,
    don't add the replaced stuff like zeros etc. to the namespace anymore,
15May2006:
    kron is now in numpy, no need to substitute anymore;
    evals from geneigsympos are now always returned as array
1Mar2006:
    moved the Vecm class to file Vecmclass.py, and "stole" kron from scipy
28Feb2006:
    add Stock-Watson common-trends calculation
20Feb2006:
    work on deterministic adjustment of GG-decomp
14Feb2006:
    bugfixes and better treatment of S&L deterministics
12Feb2006:
    deterministics estimation a la S&L added to Vecm class,
    more use of numpy-lstsq-function 31Jan2006: all functions should
    return arrays or matrix-type according to the input, where that makes
    sense, i.e. whenever a data matrix is passed to a function(and where
    the purpose is not explicitly to produce matrices)
28Jan2006:
    bugfixing related to coeffs of restricted variables, wrote function
    for symmetric-def-gen.eigval problem to remove scipy-dependency
19Jan2006:
    work started on a vecm class 19Jan2006: switched over to
    raising exceptions instead of home-cooked string generation
19Jan2006:
    functions should all return numpy-matrix-type
6Jan2006:
    switched over to numpy/new-scipy
'''
'''
NumPy class for cointegration/VECM analysis written by Sven Schreiber.

This is free but copyrighted software, distributed under the same license terms
(as of January 2007) as the 'gretl' program by Allin Cottrell and others, see
gretl.sf.net (in short: GPL v2, see www.gnu.org/copyleft/gpl.html).

(see end of this file for a changelog)
'''
from numpy import ( r_, c_, fliplr, where,
 inf, log, diff, sqrt, diag, mat, kron, abs, trace )
from numpy.linalg import eigh, solve, lstsq, det
from numpy.matlib import ones, zeros, empty, rand, eye
import os

class Vecm:             # capitalized due to Python convention (?)
    '''
    Estimate a cointegrated VAR (aka Vector Error Correction Model).

    Usage:
    First specify your model,
    my = vecm(endo_vars, lags, [cirank, determ, sd, restr_exo, unrestr_exo])
    * endo_vars as a TxN array/matrix of the endogenous vector
    * lags refers to the lag order (in levels, only affects endo-vars)
    * cirank specifies rank(Pi)=rank(beta)=rank(alpha) = r
    * determ chooses "Eviews-cases" 1 through 5
    * sd requests centered seasonal dummies: pass "q" for quarterly and "m" for
       monthly, or pass True if info in the endo data file (else False)
    * restr-exo (optional) passes further variables, restricted to the ci-space
    * unrestr-exo (optional) for further variables "outside" the ci-space
        (these are included with lags, just pass the contemp level!)

    For all data matrices, note that early periods come first.
    Alternatively you can pass filenames as strings instead of numpy matrices,
     where the data are assumed to be in csv files (gretl export format).

    Then access the attributes you need:
     my.beta_star gets ( n_y + # of restricted variables) beta with
        coeffs for the restricted variables as well
     my.alpha
     my.beta_star_id (auto-identified by eye-matrix)
     my.alpha_id (corresponding)
     my.beta_o
     my.alpha_o
     etc.
    Returns numpy-matrices even if array is passed!

    The _id variants are immutable and come from the unrestricted estimation;
     any subsequent restrictions will potentially affect all other coefficients.
    So after creating an instance, all coeffs are for the unrestricted estimate.
    Then it depends...
    '''
    def __init__(self, endo, maxlag, cirank = 1, determcase = 3, sd = False, \
                 restr_exo = None, unrestr_exo = None):
        '''
        Prepares everything and estimates the unrestricted stuff.
        '''
        ## load the data
        if type(endo) == type('filename'):
            endo, orientation, varnames, obslabels, freq = readcsv(endo)
            if orientation == 'rows': endo = endo.T
        else: endo = mat(endo)
        if sd == True: sd = freq
        # (freq must be deductable from csv-file if True is chosen)
        if type(restr_exo) == type('filename'):
            restr_exo, orientation = readcsv(restr_exo)[0:2]
            if orientation == 'rows': restr_exo = restr_exo.T
        if type(unrestr_exo) == type('filename'):
            unrestr_exo, orientation = readcsv(unrestr_exo)[0:2]
            if orientation == 'rows': unrestr_exo = unrestr_exo.T

        ## some input checks
        if type(sd) == type('a'):
            try: assert freq == sd      # freq may not exist
            except: pass                # got error w/o this line
            assert sd in 'mqMQ', 'vecm: only "q" or "m" supported for seasonals'
        if restr_exo is not None:
            assert endo.shape[0] == restr_exo.shape[0], \
             'vecm: restricted exo var obs numbers unmatched'
        if unrestr_exo is not None:
            assert endo.shape[0] == unrestr_exo.shape[0], \
             'vecm: unrestricted exo var obs numbers unmatched'
        assert ( type(maxlag) == type(4) and type(cirank) == type(4) \
           and type(determcase) == type(4) ), \
             'vecm: bogus input where integer expected'
        assert maxlag in range(1, endo.shape[0]), 'vecm: bogus maxlag'
        assert cirank in range(1, endo.shape[1]), 'vecm: bogus cirank'
        assert determcase in range(1, 6), 'vecm: bogus determcase'

        ## calculate some useful numbers
        nobs, n_y = endo.shape      # nobs: full sample
        teff = nobs - maxlag        # teff: eff. sample
        if sd == False:     n_ud = 0    # num of unrestr. determ.
        elif sd in 'qQ':    n_ud = 3
        elif sd in 'mM':    n_ud = 11
        else:               n_ud = 0              # e.g. sd=='a', or 'unknown'
        if restr_exo == None: n_rx = 0       # to be able to use n_rx later
        else:
            restr_exo = mat(restr_exo)
            n_rx = restr_exo.shape[1]
        if unrestr_exo == None: n_ux = 0       # to be able to use n_ux later
        else:
            unrestr_exo = mat(unrestr_exo)
            n_ux = unrestr_exo.shape[1]
        if determcase == 2 or determcase == 4:  n_rd = 1
        else:                                   n_rd = 0
        if determcase == 3 or determcase == 4:  n_ud += 1
        elif determcase == 5:                   n_ud += 2
        n1 = n_y + n_rx + n_rd        # rows of beta_star, all ect-components
                                      # (==p1 in Boswijk/Doornik)

        ## construct the needed data matrices
        # starting points for unrestricted and restricted:
        dy = diff(endo, axis=0)
        dylags = addLags(dy, maxlag-1)      # addLags shortens the sample
        dy = dy[-teff: , :]                 # only preserve effective sample
        # discarding the contemp diff is starting point for unrestricted:
        unrestr = dylags[:, n_y:]
            # (if no lag-diffs specified, \
            #  then unrestr.shape == (teff, 0) should hold at this point)
        # lagged level is starting point for restricted, only effective sample:
        restr = endo[:-1, :][-teff:, :]

        # now the unrestricted exo- and d-stuff (only effective sample)
        if unrestr_exo is not None:
            unrestr = c_[unrestr, unrestr_exo[-teff:, :] ]
        if sd == 'q' or sd == 'Q' or sd == 'm' or sd == 'M':
            unrestr = c_[unrestr, getDeterministics(teff, sd)]
        if determcase >= 3:
            unrestr = c_[unrestr, getDeterministics(teff, 'c')]
        if determcase == 5:
            unrestr = c_[unrestr, getDeterministics(teff, 't')]
        # (If no data was added, unrestr.shape == (self.teff, 0) still !

        # now the restricted exo- and d-stuff
        #  (also adjust for lost starting obs due to endo-lags)
        if restr_exo is not None:
            # lag the restricted exog. vars by one period to match the y_{t-1}:
            restr = c_[restr, restr_exo[:-1, :][-teff:, :]]
        if determcase == 2: restr = c_[restr,  getDeterministics(teff, 'c')]
        # for determcase >= 3 the constant is already unrestricted, no need here
        if determcase == 4: restr = c_[restr,  getDeterministics(teff, 't')]
        # for determcase == 5 the trend is already unrestricted, no need here

        ## RRR in __init__, no use for vecm without it
        if unrestr.shape[1] > 0:    # actually some columns there
            R0 = dy - unrestr * lstsq(unrestr, dy)[0]
            R1 = restr - unrestr * lstsq(unrestr, restr)[0]
        else:
            R0 = dy
            R1 = restr
        S00 = R0.T * R0 / teff
        S01 = R0.T * R1 / teff
        S11 = R1.T * R1 / teff
        S = S01.T * solve(S00, S01)
        # generalized eigenvalues, corresponding to det(lambda_j S11 - S)=0
        # only the largest n_y eigenvalues are relevant,
        # the rest are zeroes due to the restricted variables;
        # need the evals ascending first for tracestats
        evals, evectors = geneigsympos(S, S11)
        evals = evals[-n_y:]
        tracestats = -teff * log(1-evals).cumsum()
        tracestats = tracestats[::-1]      # reversed
        evals = evals[::-1]                # accordingly descending
        logL = -teff/2 * ( log(det(S00)) + log(1-evals[:cirank]).sum() )

        ## first the most natural coeffs in the model with restricted vars:
        # since evals was ascending, we have to pick the _last_ vectors...!
        beta_star = evectors[:, -cirank:]
        # now offer the usual (trivial) identification for beta
        c_id = r_[eye(cirank), zeros((n1 - cirank, cirank))]
        c_id_o = r_[zeros((cirank, n1 - cirank)), eye(n1 - cirank)]
        beta_star_id = solve((c_id.T * beta_star).T, beta_star.T).T
        # old (delete later if no assertion errors):
        beta1 = beta_star[:cirank, :].T    # should now be rxr
        beta2 = beta_star[cirank:, :].T
        beta_star_id_old = r_[eye(cirank), solve(beta1, beta2).T]
        assert (beta_star_id == beta_star_id_old).all()

        ## estimate the unrestricted alpha:
        remainingVecmCoeffs = lstsq(c_[restr * beta_star_id, unrestr], dy)[0]
        # (Dimension is (cirank+N_unrestr x n_y), and the ordering:
        #  alpha.T, dylags coeffs --gammas--, unrestr_exo coeffs,
        #  seasonal coeffs, constant, trend)
        alpha_id = remainingVecmCoeffs[:cirank, :].T

        # (the cov matrix of beta_star_id and alpha_id is relegated to the end)

        ## make stuff available
        self.beta_star = beta_star
        self.beta_star_id = beta_star_id
        self.alpha_id = alpha_id
        self.R0 = R0            #e.g. for some restricted estimations
        self.R1 = R1
        self.S00 = S00
        self.S11 = S11
        self.S01 = S01
        self.cirank = cirank
        self.tracestats = tracestats
        self.evals = evals
        self.logL = logL
        self.maxlag = maxlag
        self.teff = teff
        self.endo = endo                # this has full sample!
        self.restr_exo = restr_exo      #  dito, may be None
        self.unrestr_exo = unrestr_exo  #  dito, may be None
        self.dy = dy                    #  only effective sample
        self.unrestr = unrestr          # effective sample
        self.restr = restr              #  dito
        self.sd = sd                    # seasonal spec. needed for S&L etc.
        self.determcase = determcase
        self.n1 = n1

        # outsourced all other coeffs
        #  (including omegamat, needed for alpha/beta standard errors)
        self.setOtherCoeffs(use_given_alpha = False)

        ## the cov-matrix of beta_star_id and alpha_id,
        #   based on eq 9 in Boswijk/Doornik
        # standard errors for alpha_id
        omI = self.omegamat.I; a = alpha_id; b = beta_star_id; co = c_id_o
        topright = kron(omI*a, b.T*S11*co)
        aBinfo = c_[kron(omI, b.T*S11*b), topright]
        aBinfo = teff * r_[aBinfo, c_[topright.T, kron(a.T*omI*a, co.T*S11*co)]]
        aBcov = aBinfo.I
        # first part refers to vec(a.T)
        self.alpha_id_se = sqrt(unvec(diag(aBcov)[:n_y*cirank], cirank, n_y).T)
        # and B is just the lower part of beta_star_id
        beta_star_id_se = zeros((cirank, cirank))
        self.beta_star_id_se = r_[beta_star_id_se, \
         sqrt(unvec(diag(aBcov)[n_y*cirank:], n1-cirank, cirank))]

    def setOtherCoeffs(self, use_given_alpha = True):
        '''
        Estimates the model given self.alpha_id and self.beta_star_id.

        This is useful because the user can specify restricted alphas and betas
         (on the already started vecm-instance), and then all methods of this
         class work with those estimates.

        If the new alpha is not restricted (use_given_alpha==False),
         for a given self.beta_star all other coefficients including
         alpha_id can be found by regression.
        '''
        n_y = self.endo.shape[1];   nobs = self.endo.shape[0]
        cirank = self.cirank;       beta_star = self.beta_star

        # then the beta-o estimate by matrix algebra:
        #  (only the rows relating to the endogenous variables)
        beta_o = getOrthColumns(beta_star[:n_y, :])

        if use_given_alpha:
            alpha = self.alpha
            shortruncoeffs = lstsq(self.unrestr, \
                self.dy - self.restr * beta_star * alpha.T)[0]
            # (Dimension is (N_unrestr x n_y), and the ordering:
            #  dylags coeffs --gammas--, unrestr_exo coeffs,
            #  seasonal coeffs, constant, trend)
            # determine alpha_orth by dumb mechanics
            alpha_o = getOrthColumns(alpha)

        else:
            # estimate the unrestricted alpha:
            remainingVecmCoeffs = lstsq(c_[self.restr * beta_star, \
             self.unrestr], self.dy)[0]
            alpha = remainingVecmCoeffs[:cirank, :].T
            shortruncoeffs = remainingVecmCoeffs[cirank:, :]
            # (Dimension is (cirank+N_unrestr x n_y), and the ordering:
            #  alpha.T, dylags coeffs --gammas--, unrestr_exo coeffs,
            #  seasonal coeffs, constant, trend)

            # estimate alpha_o (cf. Gonzalo-Granger eq 27)
            # (eigenvecs corresponding to the zero eigvals of the following,
            #  but this only allows a restricted beta, so can only use this
            #  formula for freely estimated alpha)
            evecs = geneigsympos(alpha * alpha.T, self.S00)[1]
            alpha_o = evecs[:, :n_y - cirank]

        # new: generalized cov matrix for the possibly restricted case
        #  (Boswijk/Doornik p. 455)
        omegamat = self.S00 - self.S01*beta_star*alpha.T \
         - alpha*beta_star.T*self.S01.T \
         + alpha*beta_star.T*self.S11*beta_star*alpha.T

        # gammas here are textbook-style for Nx1 convention, therefore .T
        gammas = shortruncoeffs[:(self.maxlag-1) * n_y, :].T
        # (skip the unrestr_exo coeffs and seasonal coeffs for now,
        #  don't need them?)
        if self.determcase == 5:             # unrestr. trend coeff is last
            mu_unr = shortruncoeffs[-2, :]
            tau_unr = shortruncoeffs[-1, :]
        elif self.determcase in range(3, 5):  # unrestr. constant is last
            mu_unr = shortruncoeffs[-1, :]
            tau_unr = zeros(n_y)
        else:
            mu_unr = zeros(n_y)
            tau_unr = zeros(n_y)
        # from the short-run gammas construct I-\sum G_i = Psi
        #  (don't forget our Nx1 convention gammas)
        psimat = eye(n_y)
        for lag in range(1, self.maxlag):
            psimat -= gammas[:, (lag-1) * n_y : lag * n_y]
        # and for convenience the C-matrix
        #  (long-run impact/common-trends, def. follows alpha, beta, etc)
        cmat = beta_o * solve(alpha_o.T * psimat * beta_o, alpha_o.T)
        # for the common trends the residuals are also useful:
        resids = self.dy - self.restr * beta_star * alpha.T \
            - self.unrestr * shortruncoeffs
        ## construct the full-sample ect:
        if self.restr_exo is None:
            temp = empty([nobs, 0])       # should do no harm
        else: temp = self.restr_exo
        if self.determcase == 2:          # need constant in cispace
            # as beta is a subset of beta_star, we pick the
            # relevant rows
            ect = c_[self.endo, ones(nobs).T, temp] * \
                  beta_star[:n_y + 1 + temp.shape[1], :]
        elif self.determcase == 4:      # need trend in cispace
            ect = c_[ self.endo, r_['c', 1:nobs+1], temp] * \
                  beta_star[:n_y + 1 + temp.shape[1], :]
        else:                           # just the exog. vars needed
            ect = c_[self.endo, temp] * \
                  beta_star[:n_y + temp.shape[1], :]

        ## Make the (possibly) changed stuff available:
        self.alpha = alpha
        self.alpha_o = alpha_o
        self.beta_o = beta_o
        self.omegamat = omegamat
        self.gammas = gammas
        self.psimat = psimat            # I - \sum G_i
        self.cmat = cmat                # bo (ao'Psi bo)^{-1} ao'
        self.mu_unr = mu_unr            # unrestr. const., may be zero vec
        self.tau_unr = tau_unr          # unrestr. trend, may be zero vec
        self.ect = ect                  # dito
        self.resids = resids            # dito

    def output2gretl(self, outfile, matnames = []):
        '''
        Writes a gretl script (genr and print statements) to transfer results.

        outfile should be a path string,
        matnames is a string list of wanted matrix names, e.g. ['beta'],
        '''
        out = open(outfile, 'w')
        # gretl needs double quotes for print!
        out.write('print "List of result objects (python side)):"' + os.linesep)
        for name in matnames:
            out.write('matrix ' + name + ' = ' + \
                mat2gretlmatstring(eval('self.' + name)))
            out.write(os.linesep + 'print "' + name + '"' + os.linesep)
        out.close()

    def getOrthogTrend(self, extraterms = None):
        '''
        Estimates linear trend in data under assumption b'mu_1 = 0.

        Reference is Saikkonen&Luetkepohl (2000, JoE). Only the first step
        is implemented, i.e. the coeff is bo(bo'bo)^{-1} tau_*, where
        tau_* is from regression bo' Delta(y_t) = tau_* + delta_* iota_t .
        Returns the (n_y x 1) coefficient (not the data component).
        This is the complete trend coeff of the levels under the assumption
        that b'tau = 0 (no trend in I(0) directions).
        '''
        # just an abbrev:
        bo = self.beta_o; y = self.endo; nobs, n_y = y.shape
        # prepare the data for aux regression (zero starting val, iota):
        dy = r_[zeros(n_y), y]
        dy = diff(dy, axis=0)      # should auto-adjust the sample
        rhs = ones(nobs).T    # the constant term
        iota = zeros(nobs).T
        iota[0, 0] = 1
        rhs = c_[rhs, iota]

        if extraterms is not None:
            ## input checks:
            extraterms = mat(extraterms)
            assert extraterms.shape[0] == nobs, 'vecm: dims unmatched'
            # difference those extra variables, but preserve their starting val:
            #  (first row values should actually be irrelevant due to iota)
            extraterms = r_[extraterms[0, :], diff(extraterms, axis=0)]
            rhs = c_[rhs, extraterms]
        # get the tau_* coeff (at first p-r in first tuple element)
        #  (and again watch out for TxN vs. NxT and pre-/postmulti)
        coeff = lstsq(rhs, dy * bo)[0]
            # (coeff should be (2 x n_y-cirank) )
        taustar = coeff[0, :].T
            # (and taustar now should be (n_y-cirank x 1) as in S&L 2000)
        # transform in the I(1) directions:
        return bo * solve(bo.T * bo, taustar)

    def getSW(self, withrestricted = True, pstyle = False):
        '''
        Returns the Stock-Watson permanent components (with drifts).

        Best for no trend in I(0)-directions.
        With respect to further restricted terms, we will (eventually) follow
         the approach as for Gonzalo-Granger adjustment below...
        The returned matrix holds data for the whole sample, but only the
         effective sample is filled non-trivially, the starting periods are
         set to y_0, i.e. first pre-eff-sample obs.

        If proiettistyle is chosen, then the permanent components are not
         calculated from the cumulated residuals (not robust w.r.t. outliers),
         but from distributed lags of the variables themselves:
             sw_proietti_t = (I-P)*(G(1)+a*b')^-1 * G(L) * X_t,
         where P = gabi * a * (b' * gabi * a)^-1 * b',
         and gabi = (G(1) + a*b')^-1,
         a = alpha, b = beta,
         and G(L) = I - G_1 L - G_2 L^2 -...-G_{maxlag-1) L^{maxlag-1},
         and the G_i are the alternative short-run gammas referring to the
         model with the ect lagged maxlag periods.
         (See Proietti 1997)
        Adjustment for unrestricted constant (drift) is automatically done,
         implicitly rendering the

        '''
        teff = self.teff; n_y = self.endo.shape[1];  maxlag = self.maxlag
        a = self.alpha; b = self.beta_star[:n_y,:]

        if pstyle:
            # transform our gammas into alternative gammas (with broadcasting)
            # (Both gammas and altgammas should have shape (n_y, n_y*(maxlag-1)),
            #  textbook style, used with gammas * dy)
            if maxlag == 1:         # no short run, because of kron-limitation
                altgamma1 = eye(n_y)
            else:
                altgammas = kron(ones(maxlag-1), a * b.T) + self.gammas
                # G(1), gabi, P
                altgamma1 = eye(n_y) - \
                 altgammas * kron(ones(maxlag-1).T, eye(n_y))
            gabi = (altgamma1 + a * b.T).I
            pmat = gabi * a * solve(b.T * gabi * a, b.T)

            # distributed lags of levels
            # (in endo: early periods first,
            #  but in (alt)gammas: low lags first!)
            # (implicitly set pre-sample values to zero)
            # (also remember (alt)gammas refer to gamma * y convention)
            glx = self.endo.copy()
            # (w/o copy the passed Y would be messed up!)
            for lag in range(1, maxlag):
                glx[lag:, :] -= self.endo[:-lag, :] * \
                    altgammas[:, (lag-1) * n_y : lag * n_y].T

            # the permanent component
            sw = glx * gabi.T * (eye(n_y) - pmat).T
            # (demeaning should come after dealing with restricted terms)

        else:                  # standard Stock-Watson residual-based components
            sw = zeros([maxlag + teff, n_y])
            # fill everything with the one initial value y_0
            sw += self.endo[maxlag - 1 , :] # hope for broadcasting here
            # add the drifts in I(1)-direction, but only for effective sample:
            # (watch out for different cmat Nx1 convention)
            sw[-teff:, :] += r_['c', 1:teff + 1] * self.mu_unr * self.cmat.T
            # add the random walks:
            sw[-teff:, :] += self.resids.cumsum(axis=0) * self.cmat.T

        # now the adjustment for restricted terms:
        # (this actually should also deal with restricted constant and trend)
        if withrestricted:
            # taken from GG, so this actually only deals with the PX_t part
            #  i.e. the difference of the restricted terms probably would
            #  show up in the psi2_t term (relating to a comb. of Delta X_t)
            sw += (self.endo * b - self.ect) * solve(a.T * b, a.T)

        if pstyle and (self.determcase in range(3,5)):
            # we shift the mean of the transitory comp. to the permanent one,
            #  but only calculated for the effective sample
            transeff = (self.endo - sw)[-teff:, :]
            sw[-teff:, :] += transeff.mean(axis=0)

        # use initial obs values for pre-eff sample starting periods
        sw[:maxlag, :] = self.endo[:maxlag, :]

        return sw

    def getGG(self, withrestricted = True, adjusted = True):
        '''
        Returns stuff coming from the Gonzalo-Granger decomposition.

        Returns 3 matrices (in tuple) w/
        1) permanent
        2) and transitory data components,
        3) and the common I(1) factors.

        All items are given for the full available sample.

        Treatment of adjustment for levels and other deterministics:
        (pass adjusted=False to leave all that out and return the
         plain gg-decomp, where the transitory comp is from beta'Y_t)
        1. For case 3 (unrestr. constant) we use the expectation of
         beta'Y_t:
         alpha-bar' (Psi C - I) mu, and so we de-mean the transitory component
         alpha (beta'alpha)^{-1} beta' Y_t by subtracting
         alpha (beta'alpha)^{-1} alpha-bar' (Psi C - I) mu.
         Then we add that mean to the permanent part to preserve Y = P + T.
        2. We don't adjust for other unrestricted terms (yet); i.e., we want
         to preserve the additivity of Y = P + T, so we don't want to adjust T
         without changing P accordingly. But then removing stuff from T would
         add noise to P, which we don't want. (But maybe we should do that for
         seasonals, because their effects belong in the P-part if the Y_t are
         not seasonally adjusted....)
        3. If there are restricted terms (constant --case 2--, trend --case 4--,
         or user-specified), we compare the ect's with and without them, i.e.
         beta'Y_t - beta_star'Y_star_t.
         The remaining deterministic terms should be removed from the transitory
         part, i.e., we de-mean it by alpha(beta'alpha)^{-1} (.-.), and as in 1.
         we add the same part to the permanent component.
         By passing 'withrestricted = False' you can leave that step out; the
         transitory component should then (more or less visibly) still contain
         the restricted variables.
        4. With restricted trend & unrestr. constant (case 4) we still do step
         1, with the proviso that the estimate is the expectation only after
         detrending.
        5. the ordering of steps 1. and 3. seems relevant, i.e. demeaning
         should come last (?)
        '''
        # just some abbrevs:
        n_y = self.endo.shape[1]
        bo = self.beta_o; ao = self.alpha_o
        b = self.beta_star[:n_y,:]; a = self.alpha

        postfix = solve(a.T * b, a.T)

        # (watch out for the TxN vs. NxT convention and required transpose)
        trans = self.endo * b * postfix
        perma = self.endo * solve(ao.T * bo, ao.T).T * bo.T
        factors = self.endo * ao
        # stop here if no adjustment wanted:
        if adjusted == False: return (perma, trans, factors)
        # do number 3 first:
        if withrestricted == True:
            transdiff = (self.endo * b - self.ect) * postfix
            trans -= transdiff
            perma += transdiff
        # now number 1:
        if self.determcase in range(3, 5):
            # (remember self.mu_unr comes from TxN convention)
            # (but transmean is textbook Nx1 style)
            transmean = postfix.T * solve(a.T * a, a.T) * \
                    (self.psimat * self.cmat - eye(n_y)) * self.mu_unr.T
            trans -= transmean.T
            perma += transmean.T
        return (perma, trans, factors)

    def restrictAlpha(self, amatinput = None):
        '''
        Estimates the model with alpha = amatinput * psi, with known amatinput.

        If no amatinput is given, a sequence of weak exogeneity tests will be
         carried out. (each with a row of zero in amatinput) In that case no
         coefficients of the vecm are changed.

        For given amatinput sets the resulting coefficients of the restricted
        model for further use, and returns the corresponding LR test.

        Returns a three-row matrix with typical column (Lr-stat, dof, pvalue)'.
         (pvalue will be -1 if importing scipy --which is necessary-- failed)

        See Johansen book p. 124-126, amat has dimension p x m2, m2 >= r.
         1: get amato = A_orthogonal
         2: define R^tilde_0t = R_0t - ..., where R_0t are data vectors in R0?
         3: same for R^tilde_1t, where restricted variables should be included
         4: define the analogues to S0, S01, S11
         5: solve the generalized eigenvalue problem:
             lambda S11_o - S10_Ao * Abar * (Abar.T*S00_Ao*Abar).I
                    * Abar.T * S01_Ao
         6: The "first" eigenvectors of this are the restricted beta_star
         7: estimate alpha by eq. (8.8), and other coefficients
        '''
        S00 = self.S00; S01 = self.S01; R0 = self.R0; R1 = self.R1
        n_y = self.endo.shape[1]

        if amatinput is not None:
            numoftests = 1
            assert amatinput.shape[0] == self.endo.shape[1], \
                'amat needs n_y rows'
            assert amatinput.shape[1] >= self.cirank, \
                'amat needs >= r (cirank) cols'
        else: numoftests = n_y

        output = empty([3, numoftests])
        for test in range(numoftests):
            if amatinput is None:
                amat = r_[eye(n_y-1)[:test, :], zeros(n_y-1), \
                 eye(n_y-1)[test:, :]]
            else: amat = amatinput
            #1:
            amato = getOrthColumns(amat)
            #2:
            temp =  amato * solve((amato.T * S00 * amato), amato.T * R0.T)
            Rt0 = R0 - temp.T * S00
            #3:
            Rt1 = R1 - temp.T * S01
            #4:
            S00ao = Rt0.T * Rt0 / self.teff
            S01ao = Rt0.T * Rt1 / self.teff
            S11ao = Rt1.T * Rt1 / self.teff
            #5:
            amatbar = solve((amat.T * amat), amat.T).T
            temp = solve((amatbar.T * S00ao * amatbar), amatbar.T * S01ao)
            evals, evecs = geneigsympos(S01ao.T * amatbar * temp, S11ao)
            # evals descending:
            evals = evals[::-1]
            # calculate the corresponding LR stat
            lr = log((1 - evals[:self.cirank]) / (1 - self.evals[:self.cirank]))
            output[0, test] = self.teff * lr.sum()          # is now a scalar
            output[1, test] = self.cirank * (n_y - amat.shape[1])
            try:
                from scipy import stats
                pval = float(stats.chi2.sf(lr, df))
                # sf (== survival function) is 1-cdf, so one-sided pvalue
            except: pval = -1.
            output[2, test] = pval

        # change vecm coefficients only if amat was provided
        if amatinput is not None:
            #6:
            self.beta_star = evecs[:, -self.cirank:]
            #7:
            temp = lstsq(c_[R1 * self.beta_star, R0 * amato], R0 * amatbar)[0]
            psi = temp[:self.cirank, :amat.shape[1]].T
            self.alpha = amat * psi
            self.setOtherCoeffs()

        return output

    def restrictAlphaBeta(self, G, H, h0 = None, maxiter = 1000, tol = 1e-7):
        '''
        Estimate the Vecm with general linear restrictions on alpha and beta.

        See Boswijk/Doornik 2004, section 4.4;
        Two input formats:
        1. either directly, corresponding to...
          vec(alpha.T) = gmat * psi
          vec(beta) = h0 + hmat * phi
        2. or specify h0 = None to trigger the following:
          gmat a pattern for alpha (not alpha.T!) consisting of 99s
           (free parameter) and 0s (literally zero);
          hmat a pattern of 99s (free parameter) and arbitrary other numbers;
           those other numbers will literally appear in beta!

        Returns a 2-tuple of:
        0. return code, integer:
         0 if everything went ok, results are definitely usable
         -1000: no convergence; None is returned as result
         -100: beta is not identified: test result is returned anyway,
          coeffs are updated, but standard errors are not set
         -10: information matrix is not pos. semi-def.: everything is set anyway
        1. actual results (matrix or None)

        If convergence was achieved, update the coefficients in the instance,
         and if identification holds:
        self.beta_star_se
        self.alpha_se
        '''
        p1 = self.beta_star.shape[0]
        n_y = self.endo.shape[1]
        r = self.cirank
        S00 = self.S00; S01 = self.S01; S10 = S01.T; S11 = self.S11
        vecPiT = vec(solve(S11, S01.T))

        if h0 == None:           # special input format, conversion needed
            assert G.shape == (n_y, r)
            assert H.shape == (p1, r)

            ## G first
            # how many (fuzzy) zeros are there, resulting in zero-rows in G
            gzerorows = where(G < 0.1, 1, 0).sum()
            # and how many 99s (free params), yielding unit vectors in G
            gfree = where(G > 98.9, 1, 0).sum()
            assert gzerorows + gfree == n_y * r
            # initialize new gmat
            gmat = empty((0, gfree))
            vecaTinput = vec(G.T)
            freecount = 0
            for row in range(n_y*r):
                if vecaTinput[row,0] < 0.1:         # restrict to zero
                    gmat = r_[gmat, zeros((1, gfree))]
                elif vecaTinput[row,0] > 98.9:      # free param
                    gmat = r_[gmat, eye(gfree)[freecount,:]]
                    freecount += 1

            ## now H
            # how many (fuzzy) free params
            bfree = where((H > 98.9) & (H < 99.1), 1, 0).sum()
            hmat = empty((0, bfree))
            h0 = empty((0, 1))
            vecbinput = vec(H)
            freecount = 0
            for row in range(p1*r):
                if 98.9 < vecbinput[row,0] < 99.1:  # free param
                    h0 = r_[h0, zeros(1)]
                    hmat = r_[hmat, eye(bfree)[freecount, :]]
                    freecount += 1
                else:                               # restrict number directly
                    h0 = r_[h0, vecbinput[row]]
                    hmat = r_[hmat, zeros((1,bfree))]
        else:               # G, H, h0 directly specified by user
            gmat = G
            hmat = H
        # input checks
        assert gmat.shape[0] == n_y * r, 'gmat must have n * r rows'
        assert (gmat.shape[1] <= gmat.shape[0] and \
         hmat.shape[1] <= hmat.shape[0]), 'gmat and hmat need rows >= cols'
        assert (h0.shape[0] == p1 * r and hmat.shape[0] == p1 * r), \
            'h0 and hmat must have (n_y+n_rexo)*r rows'
        assert h0.shape[1] == 1, 'h0 must be a column vector'

        # starting values (last used values in class):
        a = self.alpha; b = self.beta_star; omega = self.omegamat
        omI = omega.I; outstatus = 0
        iteration = 0; lik = -inf; ldiff = inf

        # switching algorithm
        while ldiff > tol and iteration < maxiter:
            # new phi
            phi = solve(hmat.T * kron(a.T * omI * a, S11) * hmat, hmat.T) \
                * kron(a.T * omI, S11) * (vecPiT - kron(a, eye(p1)) * h0)
            #print phi.shape
            # new beta
            vecb = h0 + hmat * phi
            #print vecb
            b = unvec(vecb, p1, r)
            #print b
            # new psi
            psi = solve(gmat.T * kron(omI, b.T * S11 * b) * gmat, gmat.T)\
                * kron(omI, b.T * S11) * vecPiT
            # new alpha
            vecaT = gmat * psi
            a = unvec(vecaT, r, n_y).T
            # new omega
            omega = S00 - S01*b*a.T - a*b.T*S10 + a*b.T*S11*b*a.T
            omI = omega.I
            # logL value eq (22), concentrated
            oldlik = +lik       # explicit copy not needed for python scalars
            lik = -self.teff/2 * log(det(omega))
            ldiff = lik - oldlik
            iteration += 1

        if iteration >= maxiter: return (-1000, None)   # no convergence

        # for generic identification: random stuff (eq 40)
        psirand = 2 * ( rand((gmat.shape[1], 1)) - 0.5 )
        phirand = 2 * ( rand((hmat.shape[1], 1)) - 0.5 )
        arand = unvec(gmat * psirand, r, n_y).T
        brand = unvec(h0 + hmat * phirand, p1, r)
        jacobi = c_[kron(eye(n_y), brand) * gmat, kron(arand, eye(p1)) * hmat]
        jrank = lstsq(jacobi, ones((jacobi.shape[0], 1)))[2]
        ## stuff for inference
        # LR stat eq (27)
        lrstat = 2 * (self.logL - lik)
        # dof eq (28)
        lrdof = (n_y + p1 - r) * r - jrank
        # (delegate the p-value calc e.g. to gretl...)
        if jrank < jacobi.shape[1]:
            outstatus = -100            # generic identification failed
            alpha_se = None
            beta_star_se = None
        else:
            # p. 455
            # the info matrix eq (41) (w.r.t. theta):
            bl11 = gmat.T * kron(omI, b.T * S11 * b) * gmat
            bl12 = gmat.T * kron(omI * a, b.T * S11) * hmat
            bl21 = hmat.T * kron(a.T * omI, S11 * b) * gmat
            bl22 = hmat.T * kron(a.T * omI * a, S11) * hmat
            infomattheta = self.teff * r_[c_[bl11, bl12], c_[bl21, bl22]]

            # check pos. semi-def'ness (possibly somewhat inefficient...)
            evals = eigh(infomattheta)[0]
            if (evals < 0).any():
                outstatus = -10     # inormation matrix is not pos. semi-def.
            # mapping from theta to alphaT/beta involves {G : 0 ; 0 : H}
            theta2abT = r_[c_[gmat, zeros((n_y * r, hmat.shape[1]))], \
                          c_[zeros((p1 * r, gmat.shape[1])), hmat]]
            covab = theta2abT * solve(infomattheta, theta2abT.T)
            # first n_y*r diag elements refer to alphaT
            alpha_se = sqrt(unvec(diag(covab)[:n_y * r], r, n_y).T)
            # and the rest to beta
            beta_star_se = sqrt(unvec(diag(covab)[n_y * r:], p1, r))

        self.alpha = a
        self.beta_star = b
        self.alpha_se = alpha_se
        self.beta_star_se = beta_star_se
        # omegamat will automatically be calculated in setOtherCoeffs()
        self.setOtherCoeffs()
        # and eventually return test result
        return (outstatus, c_[lrstat, lrdof])

    def testGGfactors(self, gmatinput = None, restrict_model = False):
        '''
        Tests linear hypotheses on the Gonzalo-Granger factors (alpha_o).

        If gmat is specified, it must be a n_y x m - matrix from the expression
         ao = gmat * theta, and m >= n_y - cirank must hold.
        If nothing is specified, all individual variables are tested separately
         if they can be excluded from the factors (sequence of zero rows in
         gmat).
        If restrict_model is True, the restriction of gmat is imposed
         model-wide for further estimation (this only makes sense for a
         specified gmat with m = n_y - cirank).
        Returned are LR test statistics, dof's, and p-values in a 3-row matrix.
         (p-values == -1. if scipy import fails)
        '''
        n_y = self.endo.shape[1]
        S01 = self.S01; S11 = self.S11; S00 = self.S00;
        R0 = self.R0; R1 = self.R1

        if gmatinput is None: numoftests = n_y
        else:                           # do sequence of individual tests
            numoftests = 1
            assert gmatinput.shape[0] == n_y, "gmat rows don't match model"
            assert gmatinput.shape[1] >= n_y - self.cirank, \
                'too few cols in gmat'
                # (previous line will raise exception for numpy-1d-array input)
        output = empty([3, numoftests])

        ## do a 1- or n_y-sequence of tests:
        for test in range(numoftests):
            if gmatinput is None:
                gmat = r_[eye(n_y-1)[:test, :], zeros(n_y-1), \
                    eye(n_y-1)[test:, :]]
            else: gmat = gmatinput
            sg = gmat.T * S01 * solve(S11, S01.T) * gmat
            # ascending:
            ggevals, ggevecs = geneigsympos(sg, gmat.T * S00 * gmat)
            # descending:
            ggevals = ggevals[::-1]
            lr = log((1 - ggevals[-(n_y - self.cirank):]) / \
                    (1 - self.evals[-(n_y - self.cirank):]))
            output[0, test] = -self.teff * lr.sum()  # is now a scalar
            output[1, test] = (n_y - self.cirank) * (n_y - gmat.shape[1])

            # calculate p-values
            try:
                from scipy import stats
                output[2, test] = float(stats.chi2.sf(output[0, test], \
                                    output[1, test]))
                # (sf (== survival function) is 1-cdf, so one-sided pvalue)
            except:
                output[2, test] = -1.

        if restrict_model:
            assert gmat.shape[1] == n_y - self.cirank, \
                'too many columns in gmat to restrict model'
            # convert restricted alpha_o into restricted alpha = amat*psi
            # (evecs refer to ascending evals):
            evals, amat = eigh(gmat*gmat.T)
            amat = amat[:, evals.argsort()]
            amat = amat[:, :n_y - gmat.shape[1]]
            # and estimate with that restriction
            self.restrictAlpha(amat)
            # then we estimate alpha_o as described in GG p. 30/31,
            #  overriding the alpha_o found by dumb algebra
            self.alpha_o = gmat * ggevecs[:, :n_y - self.cirank]

        return output

    def getSLdeterm(self, which = ''):
        '''
        Returns S&L-estimation of deterministics.

        W.r.t. the deterministics, the VECM spec is used;
         however, 1 and 5 are not implemented.
        Pass the setup in 'which' parameter:
        'u' in which: unrestricted exogenous variables will also be used.
        'r' in which: restricted exogenous variables will also be used.
         (Selecting subsets is not currently possible.)
        Furthermore, if 'r' (or 'ru') is specified and we are in determcase 3,
         then the estimation of the orthogonal trend also accounts for the
         restricted exogenous variables.
        So therefore the differences of the restricted exog. vars should
         normally be included in the unrestr. exog. vars by the user, just like
         a constant is always included unrestrictedly whenever a trend term is
         specified.
        Ordering is trend, const, seasonals, restr_exo terms, unrestr_exo terms.
        Returns a 3-tuple with data component (T x n_y),
         deterministic regressors (T x...), and the coefficients (n_y x ...).
        '''
        assert which in 'rur', 'bogus input for "which"'
        assert self.determcase in range(2,5), 'SL proc only for cases 2 to 4'

        # just some abbrevs:
        y = self.endo
        nobs, n_y = y.shape
        b = self.beta_star[:n_y,:]; bo = self.beta_o
        a = self.alpha; ao = self.alpha_o

        # case 3 = orthog. trend requires special treatment;
        # (but don't forget to add the trendcoeff at first position later!)
        if self.determcase == 3:
            # first adjust for that trend separately:
            # (and save the coeffs for output)
            if 'r' in which and self.restr_exo is not None:     #user choice!
                trendcoeff = self.getOrthogTrend(self.restr_exo)
            else: trendcoeff = self.getOrthogTrend()
            y = self.endo - r_['c', 1:nobs+1] * trendcoeff.T
            # and then implicitly apply the normal analysis w/o trend

        # build the appropriate auxiliary regressor matrix
        # starting point for determ-matrix
        determ = empty([nobs, 0])
        if self.determcase == 4:               # add trend
            # (try a subtle difference and have the trend start at 1...)
            determ = c_[determ, r_['c', 1:nobs+1]]
        # always add constant:
        determ = c_[determ, ones(nobs).T]
        # seasonals:
        if self.sd == 'q': determ = c_[determ, getDeterministics(nobs, 'q')]
        elif self.sd == 'm': determ = c_[determ, getDeterministics(nobs, 'm')]
        # other exogenous terms if requested:
        if 'r' in which and self.restr_exo is not None:
            determ = c_[determ, self.restr_exo]
        if 'u' in which and self.unrestr_exo is not None:
            determ = c_[determ, self.unrestr_exo]

        # first get the short-run dynamics:
        levelscoeffs = vecm2varcoeffs(self.gammas, self.maxlag, a, b)
        # the levels coeffs now should be of shape n_y x (maxlags*n_y)
        #  (to pre-multiply, as in diss)

        ## the GLS-transformation
        # calculate the Q-matrix:
        temp1 = (a.T * solve(self.omegamat, a)).I   # to be sqrt-ed
        temp2 = (ao.T * self.omegamat * ao).I       # dito
        evals1, evecs1 = eigh(temp1)
        evals2, evecs2 = eigh(temp2)
        Q = c_[solve(self.omegamat, a) * \
                evecs1 * sqrt(mat(diag(evals1))) * evecs1.T, \
                ao * (evecs2 * sqrt(mat(diag(evals2))) * evecs2.T)]
            # (Q should now have dim n_y x n_y)
        # transform A(L) with Q:
        QTAL = Q.T * levelscoeffs

        # fill the data matrices with zero starting values:
        y = r_[zeros([self.maxlag, n_y]), y]
        determ = r_[zeros([self.maxlag, determ.shape[1]]), determ]
            # (so now both y and determ have nobs + maxlag rows)

        ###### now we do the whole thing as we did earlier, to avoid mistakes;
        # no fancy stuff, where in contrast to the rest of this code we use
        # a stacked system (the endo-data is in a single column)
        # (it might be faster to use kronecker, but messy??)
        yAux = empty([0, 1])
            # (will have nobs*n_y rows)
        bigAuxMat = empty([0, n_y * determ.shape[1]])
            # (will also have nobs*n_y rows)
        for period in range(self.maxlag, nobs):
            tempy = y[period, :].T
            for lag in range(self.maxlag):
                currenttrafo = QTAL[:, lag * n_y : (lag+1) * n_y]
                tempy -= currenttrafo * (y[period-(lag+1), :]).T
            yAux = r_[yAux, tempy]
            temprowblock = empty([n_y, 0])  # will be n_y x num-of-exo*n_y
            for var in range(determ.shape[1]):
                tempexomat = eye(n_y) * determ[period, var]
                for lag in range(self.maxlag):
                    currenttrafo = QTAL[:, lag * n_y : (lag+1) * n_y]
                    tempexomat -= currenttrafo * determ[period-(lag+1), var]
                temprowblock = c_[temprowblock, tempexomat]
            bigAuxMat = r_[bigAuxMat, temprowblock]

        ### now we're ready for showtime
        determcoeff = lstsq(bigAuxMat, yAux)[0]
        # this should be n_y*num-of-exo x 1 !
        # so try elegant reshaping to conform to determ-shape (=(T,num-of-exo):
        # (careful here! my first direct attempt messed up columns and
        #   rows! reshape apparently fills row-wise first.
        #   Here determcoeff is for post-multi after transposing.)
        determcoeff = determcoeff.reshape(n_y, determ.shape[1])
        # discard the artificial zero starting values:
        determ = determ[-nobs:, :]
        # re-insert the trend and coeff if orthog. trend was specified:
        if self.determcase == 3:
            determ = c_[r_['c', 1:nobs+1], determ]
            determcoeff = c_[trendcoeff, determcoeff]
        # we return everything, to whom it may concern:
        return (determ * determcoeff.T, determ, determcoeff)

#### end of Vecm class ###############

'''
Changelog:
25Jan2007:
    bugfix for maxlag == 1 in getSW, pstyle True
18Jan2007:
    fixed and 'finished' restrictAlphaBeta, added matrix pattern input method
     (and made it the default)
    fixed the use of setOtherCoeffs
15Jan2007:
    switched sd API to one of False/True/'m'/'q',
    fixed c_/r_ according to new numpy API,
    started to extend restrictAlphaBeta()
11Jan2007:
    switched to cross-platform lineendings with os.linesep
10Jan2007:
    set a default =3 for determcase
9Jan2007:
    set a default =1 for cirank, for easy retrieval of tracestats/evals only,
    generalize restrictAlpha() to do automatic weak exog. tests
     (like testGGfactors)
7Jan2007:
    adapt to generalized readcsv,
    add output2gretl method
6Jan2007:
    allow for csv filenames instead of direct numpy matrices
5Jan2007:
    explicit sorting of eigenvals instead of relying on numpy's implementation
2Oct2006:
    small name-change fix,
7Sep2006:
    change estim_w_restr_alpha() to restrictAlpha(), add test,
    add p-values to alpha-o Test,
    added estimation under general linear restrictions on alpha and beta
26Aug2006:
    add Proietti-style Stock-Watson trends based on variables
15Aug2006:
    fix estimation with restricted alpha-o
11Aug2006:
    adapted to new name h.getDeterministics,
10Aug2006:
    estimation with restricted alpha-o (common trends) as
     part of .testGGfactors(),
25Jul2006:
    used new numpy.matlib,
    removed h.makeNumpymatrix calls,
    replaced most inverses by solve-commands
2Jun2006:
    removed extra imports, instead only depending on helpers.py,
    used new helpers zerosm, emptym etc.,
    removed a lot of unnecessary asm, discovered numpy.diff bug,
15May2006:
    adapted to evals from geneigsympos returned as array
8Mar2006:
    added setOtherCoeffs() so that user can specify own alphas and betas,
    added tests on alpha_o, ie., whether variables are not in GG-factors
1Mar2006:
    separated the VECM class into this file instead of ts.py
28Feb2006:
    add Stock-Watson common-trends calculation
20Feb2006:
    work on deterministic adjustment of GG-decomp
14Feb2006:
    bugfixes and better treatment of S&L deterministics
12Feb2006:
    deterministics estimation a la S&L added to Vecm class,
    more use of numpy-lstsq-function
31Jan2006:
    all functions should return arrays or matrix-type according to the input
     where that makes sense, i.e. whenever a data matrix is passed to a
     function (and where the purpose is not explicitly to produce matrices)
28Jan2006:
    bugfixing related to coeffs of restricted variables,
    function for symmetric-def-gen.eigval problem to remove scipy-dependency
19Jan2006:
    work started on a vecm class,
    switched over to exceptions instead of home-cooked string generation,
    functions should all return numpy-matrix-type
6Jan2006:
    switched over to numpy/new-scipy
'''
nv = Vecm('temp_p4g_endo.csv',  2, 2, 2, False, None, 'temp_p4g_uexo.csv')
filenames = ['temp_p4g_endo.csv', 'temp_p4g_uexo.csv']
import os
nv.output2gretl('temp_p4gvecmresults.inp', ['evals', 'tracestats'])
for f in filenames: os.remove(f)

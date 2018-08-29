
import logging
import numpy as np
import numpy.fft as fft

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def quadratic(p, x):
    """
    Short Summary
    -------------
    Calculate value of x at minimum or maximum value of y,
    (value of quadratic function at argument)

    Parameters
    ----------
    p: numpy array, 3 floats
    quadratic function: p[0]*x*x + p[1]*x + p[2]

    x: 1D float array
    arguments of p()

    Returns
    -------
    maxx: float
        value of x at minimum or maximum value of y

    maxy: float
        max y = -b^2/4a occurs at x = -b^2/2a

    fit_val: 1D float array
        values of quadratic function at arguments in x array
    """
    maxx = -p[1] / (2.0 * p[0])
    maxy = -p[1] * p[1] / (4.0 * p[0]) + p[2]
    fit_val = p[0] * x * x + p[1] * x + p[2]

    return maxx, maxy, fit_val


def makeA(nh):
    """
    Long Summary
    -------------
    Writes the 'NRM matrix' that gets pseudo-inverted to provide
    (arbitrarily constrained) zero-mean phases of the holes.
    Algorithm is taken verbatim from Anand's pseudoinverse.py

    Ax = b  where x are the nh hole phases, b the nh(nh-1)/2 fringe phases,
    and A the NRM matrix

    Solve for the hole phases:
        Apinv = np.linalg.pinv(A)
        Solution for unknown x's:
        x = np.dot(Apinv, b)

    Following Noah Gamper's convention of fringe phases,
    for holes 'a b c d e f g', rows of A are

        (-1 +1  0  0  ...)
        (0 -1 +1  0  ...)

    which is implemented in makeA() as:
        matrixA[row,h2] = -1
        matrixA[row,h1] = +1

    To change the convention just reverse the signs of the 'ones'.

    When tested against Alex'' NRM_Model.py 'piston_phase' text output
    of fringe phases, these signs appear to be correct -
    anand@stsci.edu 12 Nov 2014

    Parameters
    ----------
    nh: integer
        number of holes in NR mask

    Returns
    -------
    matrixA: 2D float array
         nh columns, nh(nh-1)/2 rows (eg 21 for nh=7)
    """
    log.debug('-------')
    log.debug(' makeA:')

    ncols = (nh * (nh - 1)) // 2
    nrows = nh
    matrixA = np.zeros((ncols, nrows))

    row = 0
    for h2 in range(nh):
        for h1 in range(h2 + 1, nh):
            if h1 >= nh:
                break
            else:
                log.debug(' row: %s, h1: %s, h2: %s', row, h1, h2)

                matrixA[row, h2] = -1
                matrixA[row, h1] = +1
                row += 1

    log.debug('matrixA:')
    log.debug(' %s', matrixA)

    return matrixA


def fringes2pistons(fringephases, nholes):
    """
    Short Summary
    -------------
    For NRM_Model.py to use to extract pistons out of fringes, given
    its hole bookkeeping, which apparently matches that of this module,
    and is the same as Noah Gamper's.

    Parameters
    ----------
    fringephases: 1D integer array
        fringe phases

    nholes: integer
        number of holes

    Returns
    -------
    np.dot(Apinv, fringephases): 1D integer array
        pistons in same units as fringe phases
    """
    Anrm = makeA(nholes)
    Apinv = np.linalg.pinv(Anrm)

    return -np.dot(Apinv, fringephases)


def rebin(a=None, rc=(2, 2)):
    """
    Short Summary
    -------------
    Perform simple-minded flux-conserving binning using specified binning
    kernel, clipping trailing size mismatch: eg a 10x3 array binned by
    3 results in a 3x1 array

    Parameters
    ----------
    a: 2D float array
        input array to bin

    rc: 2D float array
        binning kernel

    Returns
    -------
    binned_arr: float array
        binned array
    """
    binned_arr = krebin(a, (a.shape[0] // rc[0], a.shape[1] // rc[1]))

    return binned_arr


def krebin(a, shape):
    """
    Short Summary
    -------------
    Klaus P's fastrebin from web

    Parameters
    ----------
    a: 2D float array
        input array to rebin

    shape: tuple (integer, integer)
        dimensions of array 'a' binned down by dimensions of binning kernel

    Returns
    -------
    reshaped_a: 2D float array
        reshaped input array
    """
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    reshaped_a = a.reshape(sh).sum(-1).sum(1)

    return reshaped_a


def rcrosscorrelate(a=None, b=None):
    """
    Short Summary
    -------------
    Calculate cross correlation of two identically-shaped real arrays

    Parameters
    ----------
    a: 2D float array
        first input array

    b: 2D float array
        second input array

    Returns
    -------
    c.real.copy():
        real part of array that is the correlation of the two input arrays.
    """

    c = crosscorrelate(a=a, b=b)/(np.sqrt((a*a).sum())*np.sqrt((b*b).sum()))
    return c.real.copy()


def crosscorrelate(a=None, b=None):
    """
    Short Summary
    -------------
    Calculate cross correlation of two identically-shaped real or complex arrays

    Parameters
    ----------
    a: 2D complex float array
        first input array

    b: 2D complex float array
        second input array

    Returns
    -------
    fft.fftshift(c)
        complex array that is the correlation of the two input arrays.
    """
    if a.shape != b.shape:
        log.critical('crosscorrelate: need identical arrays')
        return None

    fac = np.sqrt(a.shape[0] * a.shape[1])

    A = fft.fft2(a) / fac
    B = fft.fft2(b) / fac
    c = fft.ifft2(A * B.conj()) * fac * fac

    log.debug('----------------')
    log.debug(' crosscorrelate:')
    log.debug(' a: %s:', a)
    log.debug(' A: %s:', A)
    log.debug(' b: %s:', b)
    log.debug(' B: %s:', B)
    log.debug(' c: %s:', c)
    log.debug(' a.sum(): %s:', a.sum())
    log.debug(' b.sum(): %s:', b.sum())
    log.debug(' c.sum(): %s:', c.sum())
    log.debug(' a.sum()*b.sum(): %s:', a.sum() * b.sum())
    log.debug(' c.sum().real: %s:', c.sum().real)
    log.debug(' a.sum()*b.sum()/c.sum().real: %s:', a.sum()*b.sum()/c.sum().real)

    return fft.fftshift(c)


def findmax(mag, vals, mid=1.0):
    """
    Short Summary
    -------------
    Fit a quadratic to the given input arrays mag and vals, and calculate the
    value of mag at the extreme value of vals.

    Parameters
    ----------
    mag: 1D float array
        array for abscissa

    vals: 1D float array
        array for ordinate

    mid: float
        midpoint of range

    Returns
    -------
    maxx: float
        value of mag at the extreme value of vals

    maxy: float
        value of vals corresponding to maxx
    """
    p = np.polyfit(mag, vals, 2)
    fitr = np.arange(0.95 * mid, 1.05 * mid, .01)
    maxx, maxy, fitc = quadratic(p, fitr)

    return maxx, maxy

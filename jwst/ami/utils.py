from jwst.datamodels import dqflags
from . import matrix_dft

import logging
import numpy as np
import numpy.fft as fft

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())


class Affine2d():
    """
    A class to help implement the Bracewell Fourier 2D affine transformation
    theorem to calculate appropriate coordinate grids in the Fourier domain (eg
    image space, momentum space), given an affine transformation of the
    original (eg pupil space, configuration space) space.  This class provides
    the required normalization for the Fourier transform, and provides for
    a single way to set pixel pitch (including independent x and y scales)
    in the image plane.

    The theorem states that if f(x,y) and F(u,v) are Fourier pairs, and
    g(x,y) = f(x',y'), where

        x' = mx * x  +  sx * y  +  xo
        y' = my * y  +  sy * x  +  yo,

    then G(u,v), the Fourier transform of g(x,y), is given by:

        G(u,v) = ( 1/|Delta| ) * exp { (2*Pi*i/Delta) *
                                          [ (my*xo - sx*yo) * u  +
                                            (mx*yo - sy*xo) * v  ] }  *
                                 F{ ( my*u - sy*v) / Delta,
                                    (-sx*u + mx*v) / Delta  }
    where:
                Delta = mx * my - sx * sy.

    The reverse transformation, from (x',y') to (x,y) is given by:

        x = (1/Delta) * ( my * x' - sx * y'  -  my * xo + sx * yo )
        y = (1/Delta) * ( mx * y' - sy * x'  -  mx * yo + sy * xo )

    For clarity we call (x,y) IDEAL coordinates, and (x',y') DISTORTED
    coordinates. We know the analytical form of F(u,v) from the literature,
    and need to calculate G(u,v) at a grid of points in the (u,v) space, with
    two lattice vectors a and b defining the grid.  These lattice vectors have
    components a=(au,av) and b=(bu,bv) along the u and v axes.

    Discussion with Randall Telfer (2018.05.18)  clarified that:

        These constants, properly applied to the analytical transform in a
        "pitch matrix" instead of a scalar "pitch" variable, provide the PSF
        sampled in radians on an imaginary detector that is perpendicular to
        the chief ray.  The actual detector might be tilted in a different
        manner, changing the x pitch and y pitch of the detector pixels
        independent of the effects of pupil distortion.

        We presume the main use of this object is to calculate intensity in the
        detector, so we include a DetectorTilt object in this class, although
        this object is constructed to have an 'identity' effect during the
        initial development and use of the Affine2d class in NRM data analysis.
        For most physical detector tilts we expect the DetectorTilt to have a
        small effect on an image simulated using the Fourier transform.  There
        are exceptions to this 'small effect' expectation (eg HST NICMOS 2 has
        a detector tilt of a few tens of degrees).  As long as the detector is
        small compared to the effective focal length (i.e. detector size <<<
        nominal f-ratio * primary diameter) of the system, detector tilts will
        change the pixel pitch (in radians) linearly across the field.

        There may be an ambiguity between the 'detector tilt effect' on pixel
        pitch and the diffractive effect (which results from pupil distortion
        between a pupil stop and the primary).  This might have to be broken
        using a pupil distortion from optical modelling such as ray tracing.
        Or it could be broken by requiring the detector tilt effect to be
        derived from optical models and known solid body models or metrology of
        the instrument/telescope, and the optical pupil distortion found from
        fitting on-sky data.

    Jean Baptiste Joseph Fourier 1768-1830
    Ron Bracewell 1921-2007
    Code by Anand Sivaramakrishnan 2018
    """

    def __init__(self, mx=None, my=None, sx=None, sy=None, xo=None, yo=None,
                 rotradccw=None, name="Affine"):
        """
        Short Summary
        -------------
        Initialize with transformation constants

        Parameters
        ----------
        mx: float
            dimensionless x-magnification

        my: float
            dimensionless y-magnification

        sx: float
            dimensionless x shear

        sy: float
            dimensionless y shear

        xo: float
            x-offset in pupil space

        yo: float
            y-offset in pupil space

        rotradccw: float
            a counter-clockwise rotation of *THE VECTOR FROM THE ORIGIN TO A
            POINT*, in a FIXED COORDINATE FRAME, by this angle (radians)
            (as viewed in ds9 or with fits NAXIS1 on X and NAXIS2 on Y);
            default is None

        name: string, optional

        Returns
        -------
        None
        """
        self.rotradccw = rotradccw
        if rotradccw is not None:
            mx = np.cos(rotradccw)
            my = np.cos(rotradccw)
            sx = -np.sin(rotradccw)
            sy = np.sin(rotradccw)
            xo = 0.0
            yo = 0.0

        self.mx = mx
        self.my = my
        self.sx = sx
        self.sy = sy
        self.xo = xo
        self.yo = yo
        self.determinant = mx * my - sx * sy
        self.absdeterminant = np.abs(self.determinant)
        self.name = name

        """
        numpy vector of length 2, (xprime,yprime) for use in manually writing
        the dot product needed for the exponent in the transform theorem.  Use
        this 2vec to dot with (x,y) in fromfunc to create the 'phase argument'
        Since this uses an offset xo yo in pixels of the affine transformation,
        these are *NOT* affected by the 'oversample' in image space.  The
        vector it is dotted with is in image space."""
        self.phase_2vector = np.array((my * xo - sx * yo,
                                       mx * yo - sy * xo)) / self.determinant

    def forward(self, point):
        """
        Short Summary
        -------------
        Create the forward affine transformation, in ideal-to-distorted coordinates.

        Parameters
        ----------
        point: float, float
            coordinates in ideal space, which to apply forward transform

        Returns
        -------
        trans_point: float, float
            coordinates in distorted space
        """
        trans_point = np.array((self.mx * point[0] + self.sx * point[1] + self.xo,
                                self.my * point[1] + self.sy * point[0] + self.yo))

        return trans_point

    def reverse(self, point):
        """
        Short Summary
        -------------
        Create the reverse affine transformation, in distorted-to-ideal coordinates.

        Parameters
        ----------
        point: float, float
            coordinates in distorted space, which to apply reverse transform

        Returns
        -------
        trans_point: float, float
            coordinates in ideal space
        """
        trans_point = np.array(
            (self.my * point[0] - self.sx * point[1] -
             self.my * self.xo + self.sx * self.yo,
             self.mx * point[1] - self.sy * point[0] -
             self.mx * self.yo + self.sy * self.xo)) * self.determinant

        return trans_point

    def distortFargs(self, u, v):
        """
        Short Summary
        -------------
        Implement the (u,v) to (u',v') change in arguments of F. See class
        documentation of Bracewell Fourier 2D affine transformation theorem.

        Parameters
        ----------
        u: float
            1st argument of F

        v: float
            2nd argument of F

        Returns
        -------
        uprime: float
            1st transformed argument of F
        vprime: float
            2nd transformed argument of F
        """
        uprime = (self.my * u - self.sy * v) / self.determinant
        vprime = (-self.sx * u + self.mx * v) / self.determinant
        return uprime, vprime

    def distortphase(self, u, v):
        """
        Short Summary
        -------------
        Calculate the phase term in the Bracewell Fourier 2D affine
        transformation theorem. The phase term is:

        1/|Delta| * exp{(2*Pi*i/Delta) * [(my*xo- x*yo) * u + (mx*yo-sy*xo)*v]}

        where u and v are in inverse length units.

        Parameters
        ----------
        u: float
            1st argument of F, in units of inverse length units

        v: float
            2nd argument of F, in units of inverse length units

        Returns
        -------
        phase: complex array
            phase term divided by the determinant.
        """
        phase = np.exp(2 * np.pi * 1j / self.determinant *
                       (self.phase_2vector[0] * u + self.phase_2vector[1] * v))

        return phase

    def get_rotd(self):
        """
        Short Summary
        -------------
        Calculate the rotation that was used to creat a pure rotation
        affine2d object.

        Parameters
        ----------
        None

        Returns
        -------
        rotd: float
            rotation used to creat a pure rotation affine2d
        """
        if self.rotradccw:
            rotd = 180.0 * self.rotradccw / np.pi
            return rotd
        else:
            return None


def affinepars2header(hdr, affine2d):
    """
    Short Summary
    -------------
    Write the affine2d parameters into fits header (will be modified or deleted
    in later build)

    Parameters
    ----------
    hdr: fits header
        fits header to write affine2d parameters into

    affine2d: Affine2d object

    Returns
    -------
    hdr: fits header
        fits header, updated with affine2d parameters
    """

    hdr['affine'] = (affine2d.name, 'Affine2d in pupil: name')
    hdr['aff_mx'] = (affine2d.mx, 'Affine2d in pupil: xmag')
    hdr['aff_my'] = (affine2d.my, 'Affine2d in pupil: ymag')
    hdr['aff_sx'] = (affine2d.sx, 'Affine2d in pupil: xshear')
    hdr['aff_sy'] = (affine2d.sx, 'Affine2d in pupil: yshear')
    hdr['aff_xo'] = (affine2d.xo, 'Affine2d in pupil: x offset')
    hdr['aff_yo'] = (affine2d.yo, 'Affine2d in pupil: y offset')
    hdr['aff_dev'] = ('analyticnrm2', 'dev_phasor')

    return hdr


def makedisk(N, R, ctr=(0, 0)):
    """
    Short Summary
    -------------
    Calculate a 'disk', an array whose values =1 in a circular region near
    the center of the array, and =0 elsewhere.

    Parameters
    ----------
    N: integer
        size of 1 dimension of the array to be returned

    R: integer
        radius of disk

    ctr: (integer, integer)
        center of disk

    Returns
    -------
    array: 2D integer array
        array whose values =1 in a circular region near the center of the
        array, and =0 elsewhere.
    """
    if N % 2 == 1:  # odd
        M = (N - 1) / 2
        xx = np.linspace(-M - ctr[0], M - ctr[0], N)
        yy = np.linspace(-M - ctr[1], M - ctr[1], N)
    if N % 2 == 0:  # even
        M = N / 2
        xx = np.linspace(-M - ctr[0], M - ctr[0] - 1, N)
        yy = np.linspace(-M - ctr[1], M - ctr[1] - 1, N)

    (x, y) = np.meshgrid(xx, yy.T)
    r = np.sqrt((x**2) + (y**2))
    array = np.zeros((N, N))
    array[r < R] = 1

    return array


def trim(m, s):
    """
    Short Summary
    -------------
    Remove the edge pixels from an index mask m.

    Parameters
    ----------
    m: (integer, integer) array
        2d index mask

    s: integer
        side of the parent array that was used to generate m.

    Returns
    -------
    m_masked: (integer, integer) array
        2d index mask with edge pixels trimmed
    """
    xl, yl = [], []  # trimmed lists
    for ii in range(len(m[0])):
        # Go through all indices in the mask:
        # the x & y lists test for any index being an edge index - if none are
        # on the edge, remember the indices in new list
        if (m[0][ii] == 0 or m[1][ii] == 0 or m[0][ii] == s - 1 or
                m[1][ii] == s - 1) is False:
            xl.append(m[0][ii])
            yl.append(m[1][ii])
    m_masked = (np.asarray(xl), np.asarray(yl))

    return m_masked


def avoidhexsingularity(rotation):
    """
    Short Summary
    -------------
    Avoid rotation of exact multiples of 15 degrees to avoid NaN's in
    hextransformee()

    Parameters
    ----------
    rotation: float
       rotation in degrees int or float

    Returns
    -------
    rotation_adjusted: float
        replacement value for rotation with epsilon = 1.0e-12 degrees added.
        Precondition before using rotationdegrees in Affine2d for hex geometries
    """
    diagnostic = rotation / 15.0 - int(rotation / 15.0)
    epsilon = 1.0e-12
    if abs(diagnostic) < epsilon / 2.0:
        rotation_adjusted = rotation + epsilon
    else:
        rotation_adjusted = rotation
    return rotation_adjusted


def center_imagepeak(img, r='default', cntrimg=True):
    """
    Short Summary
    -------------
    Calculate a cropped version of the input image centered on the peak pixel.

    Parameters
    ----------
    img: 2D float array
        input image array

    r: integer
        offset for center determination

    cntrimg: boolean
        If True, center on the peak pixel

    Returns
    -------
    cropped: 2D float array
        Cropped to place the brightest pixel at the center of the img array
    """
    peakx, peaky, h = min_distance_to_edge(img, cntrimg=cntrimg)
    log.debug(' peakx=%g, peaky=%g, distance to edge=%g', peakx, peaky, h)
    if r == 'default':
        r = h.copy()
    else:
        pass

    cropped = img[int(peakx - r):int(peakx + r + 1), int(peaky - r):int(peaky + r + 1)]

    return cropped


def centerpoint(s):
    """
    Short Summary
    -------------
    Calculate center of image, accounting for odd/even pixel size;
    used for jinc() and hex transform functions.

    Parameters
    ----------
    s: 2D integer or float tuple
        array shape

    Returns
    -------
    center: 2D integer or float tuple
        center of image
    """
    return (0.5 * s[0] - 0.5, 0.5 * s[1] - 0.5)


def min_distance_to_edge(img, cntrimg=True):
    """
    Short Summary
    -------------
    Calculate the coordinates of the brightest pixel, and the distance between
    the brightest pixel and the nearest edge of the input array.

    Parameters
    ----------
    img: 2D array
        input array

    cntrimg: boolean
        if True, only look for the peak pixel near the center of the image

    Returns
    -------
    peakx, peaky: integer, integer

    h: integer
        distance to the nearest image edge
    """
    if cntrimg is True:
        # Only look for the peak pixel at the center of the image
        ann = makedisk(img.shape[0], 31)  # search radius around array center
    else:
        # Peak of the image can be anywhere
        ann = np.ones((img.shape[0], img.shape[1]))

    peakmask = np.where(img == np.nanmax(np.ma.masked_invalid(img[ann == 1])))
    # following line takes care of peaks at two or more identical-value max
    #   pixel locations:
    peakx, peaky = peakmask[0][0], peakmask[1][0]

    dhigh = (img.shape[0] - peakx - 1, img.shape[1] - peaky - 1)
    dlow = (peakx, peaky)
    h0 = min((dhigh[0], dlow[0]))
    h1 = min((dhigh[1], dlow[1]))
    h = min(h0, h1)

    return peakx, peaky, h  # the 'half side' each way from the peak pixel


def find_centroid(a, thresh):
    """
    Short Summary
    -------------
    Calculate the centroid of the image

    Long Summary
    ------------
    Original domain a, Fourier domain CV
    sft square image a to CV array, no loss or oversampling - like an fft.
    Normalize peak of abs(CV) to unity
    Create 'live area' mask from abs(CV) with slight undersizing
        (allow for 1 pixel shifts to live data still)
        (splodges, or full image a la KP)
    Calculate phase slopes using CV.angle() phase array
    Calculate mean of phase slopes over mask
    Normalize phase slopes to reflect image centroid location in pixels

    XY conventions meshed to lg_model conventions:
    if you simulate a psf with pixel_offset = ( (0.2, 0.4), ) then blind
        application  centroid = utils.find_centroid()

    returns the image centroid (0.40036, 0.2000093) pixels in image space. To
    use this in lg_model, nrm_core,... you will want to calculate the new image
    center using:
    image_center = utils.centerpoint(s) + np.array((centroid[1], centroid[0])
    and everything holds together sensibly looking at DS9 images of a.

    Parameters
    ----------
    a: 2D square float
        input image array

    thresh: float
        Threshold for the absolute value of the FT(a). Normalize abs(CV = FT(a))
        for unity peak, and define the support of good CV when this is above
        threshold, then find the phase slope of the CV only over this support.

    Returns
    -------
    htilt, vtilt: float, float
        Centroid of a, as offset from array center, as calculated by the DFT's.
    """
    ft = matrix_dft.MatrixFourierTransform()

    cv = ft.perform(a, a.shape[0], a.shape[0])
    cvmod, cvpha = np.abs(cv), np.angle(cv)

    cvmod = cvmod / cvmod.max()  # normalize to unity peak

    cvmask = np.where(cvmod >= thresh)

    cvmask_edgetrim = trim(cvmask, a.shape[0])

    htilt, vtilt = findslope(cvpha, cvmask_edgetrim)

    return htilt, vtilt


def quadratic_extremum(p):
    """
    Short Summary
    -------------
    Calculate maximum of the quadratic

    Parameters
    ----------
    p: float, float, float
        quadratic coefficients

    Returns
    -------
    y_max: float
        maximum of the quadratic
    """
    y_max = -p[1] / (2.0 * p[0]), -p[1] * p[1] / (4.0 * p[0]) + p[2]

    return y_max


def findpeak_1d(yvec, xvec):
    """
    Short Summary
    -------------
    Calculate the fit function extreme for a given input vector

    Parameters
    ----------
    yvec: 1D float array
       function values for input vector

    xvec: 1D float array
       input vector

    Returns
    -------
    quad_ext: float
        fit function extreme for a given input vector
    """
    p = np.polyfit(np.array(xvec), np.array(yvec), 2)
    return quadratic_extremum(p)


def findslope(a, m):
    """
    Short Summary
    -------------
    Find slopes of an array

    Long Summary
    ------------
    Find slopes of an array, over pixels not bordering the edge of the array.
    There should be valid data on either side of every pixel selected by the mask
    m. a is in radians of phase (in Fourier domain) when used in NRM/KP
    applications. The std dev of the middle 9 pixels is used to further clean
    the mask 'm' of invalid slope data, where we're subtracting
    inside-mask-support from outside-mask-support. This mask is called newmask.
    Converting tilt in radians per Fourier Domain (eg pupil_ACF) pixel
    Original Domain (eg image intensity) pixels:

    If the tilt[0] is 2 pi radians per ODpixel you recover the same OD
    array you started with.  That means you shifted the ODarray one
    full lattice spacing, the input array size, so you moved it by
    OD.shape[0].

    2 pi/FDpixel of phase slope => ODarray.shape[0]
    1 rad/FDpixel of phase slope => ODarray.shape[0]/(2 pi) shift
    x rad/FDpixel of phase slope => x * ODarray.shape[0]/(2 pi) ODpixels shift

    Gain between rad/pix phase slope and original domin pixels is
         a.shape[0 or 1]/(2 pi)
    Multiply the measured phase slope by this gain for pixels of incoming array
         centroid shift away from array center.

    Parameters
    ----------
    a: 2D float array
        phase slope array

    m: 2D array, integer
        mask array

    Returns
    -------
    slopes: 2D float array
        slopes
    """
    a_up = np.zeros(a.shape)
    a_dn = np.zeros(a.shape)
    a_l = np.zeros(a.shape)
    a_r = np.zeros(a.shape)

    a_up[:, 1:] = a[:, :-1]
    a_dn[:, :-1] = a[:, 1:]

    a_r[1:, :] = a[:-1, :]
    a_l[:-1, :] = a[1:, :]

    offsetcube = np.zeros((4, a.shape[0], a.shape[1]))
    offsetcube[0, :, :] = a_up
    offsetcube[1, :, :] = a_dn
    offsetcube[2, :, :] = a_r
    offsetcube[3, :, :] = a_l

    tilt = np.zeros(a.shape), np.zeros(a.shape)
    tilt = (a_r - a_l) / 2.0, (a_up - a_dn) / 2.0  # raw estimate of phase slope
    c = centerpoint(a.shape)
    C = (int(c[0]), int(c[1]))
    sigh, sigv = tilt[0][C[0] - 1:C[0] + 1, C[1] - 1:C[1] + 1].std(), \
        tilt[1][C[0] - 1:C[0] + 1, C[1] - 1:C[1] + 1].std()
    avgh, avgv = tilt[0][C[0] - 1:C[0] + 1, C[1] - 1:C[1] + 1].mean(), \
        tilt[1][C[0] - 1:C[0] + 1, C[1] - 1:C[1] + 1].mean()

    # second stage mask cleaning: 5 sig rejection of mask
    newmaskh = np.where(np.abs(tilt[0] - avgh) < 5 * sigh)
    newmaskv = np.where(np.abs(tilt[1] - avgv) < 5 * sigv)

    th, tv = np.zeros(a.shape), np.zeros(a.shape)
    th[newmaskh] = tilt[0][newmaskh]
    tv[newmaskv] = tilt[1][newmaskv]

    # determine units of tilt -
    G = a.shape[0] / (2.0 * np.pi), a.shape[1] / (2.0 * np.pi)

    slopes = G[0] * tilt[0][newmaskh].mean(), G[1] * tilt[1][newmaskv].mean()
    return slopes


def quadratic(p, x):
    """
    Short Summary
    -------------
    Calculate value of x at minimum or maximum value of y,
    (value of quadratic function at argument)

    Parameters
    ----------
    p: float, float, float
        coefficients of quadratic function: p[0]*x*x + p[1]*x + p[2]

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

    When tested against Alex'' nrm_model.py 'piston_phase' text output
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
    For nrm_model.py to use to extract pistons out of fringes, given
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

    return np.dot(Apinv, fringephases)


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

    c = crosscorrelate(a=a, b=b) / (np.sqrt((a * a).sum()) * np.sqrt((b * b).sum()))
    return c.real.copy()


def lambdasteps(lam, frac_width, steps=4):
    """
    Short Summary
    -------------
    Create array of increments of lambda

    Parameters
    ---------

    lam: float
        lambda

    frac_width: float
        fractional bandwidth

    steps: integer
        With lam and frac, determines bin size of lambda array

    Returns
    -------
    lambda_array; 1D float array
        Array of increments of lambda
    """

    frac = frac_width / 2.0
    steps = steps / 2.0

    # add some very small number to the end to include the last number.
    lambda_array = np.arange(-1 * frac * lam + lam, frac * lam + lam + 10e-10,
                             frac * lam / steps)

    return lambda_array


def tophatfilter(lam_c, frac_width, npoints=10):
    """
    Short Summary
    -------------
    Create tophat filter list from array of lambda values

    Parameters
    ----------
    lam_c: float
        lambda

    frac_width: float
        fractional bandwidth

    npoints: integer
        number of bins in lambda array

    Returns
    -------
    filt: list
        tophat filter list
    """
    wllist = lambdasteps(lam_c, frac_width, steps=npoints)
    filt = []
    for ii in range(len(wllist)):
        filt.append(np.array([1.0, wllist[ii]]))
    return filt


def crosscorrelate(a=None, b=None):
    """
    Short Summary
    -------------
    Calculate cross correlation of two identically-shaped real or complex
    arrays

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
    log.debug(' a.sum: %s:', a.sum())
    log.debug(' b.sum: %s:', b.sum())
    log.debug(' c.sum: %s:', c.sum())
    log.debug(' a.sum*b.sum: %s:', a.sum() * b.sum())
    log.debug(' c.sum.real: %s:', c.sum().real)
    log.debug(' a.sum*b.sum/c.sum.real: %s:', a.sum() * b.sum() / c.sum().real)

    return fft.fftshift(c)


def rotate2dccw(vectors, thetarad):
    """
    Short Summary
    -------------
    Apply a CCW rotation to the given vectors. For a positive (CCW) rotation:
        x decreases under slight rotation
        y increases under slight rotation

    Parameters
    ----------
    vectors: list
       2D vectors

    thetarad: float
       rotation

    Returns
    -------
    rot_vectors: array of floats
        rotated vectors
    """
    c, s = (np.cos(thetarad), np.sin(thetarad))
    ctrs_rotated = []
    for vector in vectors:
        ctrs_rotated.append([c * vector[0] - s * vector[1],
                             s * vector[0] + c * vector[1]])
    rot_vectors = np.array(ctrs_rotated)

    return rot_vectors


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


def pix_median_fill_value(input_array, input_dq_array, bsize, xc, yc):
    """
    Short Summary
    -------------
    For the pixel specified by (xc, yc), calculate the median value of the
    good values within the box of size bsize neighboring pixels. If any of
    the box is outside the data, 0 will be returned.

    Parameters
    ----------
    input_array: ndarray
        2D input array to filter
    input_dq_array: ndarray
        2D input data quality array
    bsize: scalar
        square box size of the data to extract
    xc: scalar
        x position of the data extraction
    yc: scalar
        y position of the data extraction

    Returns
    -------
    median_value: float
        median value of good values within box of neighboring pixels

    """
    # set the half box size
    hbox = int(bsize / 2)

    # Extract the region of interest for the data
    try:
        data_array = input_array[yc - hbox:yc + hbox + 1, xc - hbox: xc + hbox + 1]
        dq_array = input_dq_array[yc - hbox:yc + hbox + 1, xc - hbox: xc + hbox + 1]
    except IndexError:
        # If the box is outside the data, return 0
        log.warning('Box for median filter is outside the data')
        return 0.

    # only keep pixels not flagged with DO_NOT_USE
    wh_good = np.where((np.bitwise_and(dq_array, dqflags.pixel['DO_NOT_USE'])
                        == 0))
    filtered_array = data_array[wh_good]

    # compute the median, excluding NaN's
    median_value = np.nanmedian(filtered_array)

    # check for bad result
    if np.isnan(median_value):
        log.warning('Median filter returned NaN; setting value to 0.')
        median_value = 0.

    return median_value


def mas2rad(mas):
    """
    Short Summary
    -------------
    Convert angle in milli arc-sec to radians

    Parameters
    ----------
    mas: float
        angle in milli arc-sec

    Returns
    -------
    rad: float
        angle in radians
    """
    rad = mas * (10**(-3)) / (3600 * 180 / np.pi)
    return rad


def img_median_replace(img_model, box_size):
    """
    Short Summary
    -------------
    Replace bad pixels (either due to a DQ value of DO_NOT_USE or having a
    value of NaN) with the median value of surrounding good pixels.

    Parameters
    ----------
    img_model: image model containing input array to filter.

    box_size: scalar
        box size for the median filter

    Returns
    -------
    img_model: input image model whose input array has its bad pixels replaced
        by the median of the surrounding good-value pixels.
    """
    input_data = img_model.data
    input_dq = img_model.dq

    num_nan = np.count_nonzero(np.isnan(input_data))
    num_dq_bad = np.count_nonzero(input_dq == dqflags.pixel['DO_NOT_USE'])

    # check to see if any of the pixels are bad
    if (num_nan + num_dq_bad > 0):

        log.info(f'Applying median filter for {num_nan} NaN and {num_dq_bad} DO_NOT_USE pixels')
        bad_locations = np.where(np.isnan(input_data) |
                                 np.equal(input_dq,
                                          dqflags.pixel['DO_NOT_USE']))

        # fill the bad pixel values with the median of the data in a box region
        for i_pos in range(len(bad_locations[0])):
            y_box_pos = bad_locations[0][i_pos]
            x_box_pos = bad_locations[1][i_pos]
            median_fill = pix_median_fill_value(input_data, input_dq,
                                                box_size, x_box_pos, y_box_pos)
            input_data[y_box_pos, x_box_pos] = median_fill

        img_model.data = input_data

    return img_model

from .matrix_dft import matrix_dft

import logging
import numpy as np
import numpy.fft as fft
from scipy.integrate import simpson
from astropy import units as u

import synphot

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())


class Affine2d:
    """
    Implement the Bracewell Fourier 2D affine transformation theorem.

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
    components a=(a_u,a_v) and b=(b_u,b_v) along the u and v axes.

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

    def __init__(
        self,
        mx=None,
        my=None,
        sx=None,
        sy=None,
        xo=None,
        yo=None,
        rotradccw=None,
        name="Affine",
    ):
        """
        Initialize with transformation constants.

        Parameters
        ----------
        mx : float
            Dimensionless x-magnification
        my : float
            Dimensionless y-magnification
        sx : float
            Dimensionless x shear
        sy : float
            Dimensionless y shear
        xo : float
            X-offset in pupil space
        yo : float
            Y-offset in pupil space
        rotradccw : float
            A counter-clockwise rotation of *THE VECTOR FROM THE ORIGIN TO A
            POINT*, in a FIXED COORDINATE FRAME, by this angle (radians)
            (as viewed in ds9 or with fits NAXIS1 on X and NAXIS2 on Y);
            default is None
        name : str, optional
            Name of the Affine2d object to store in name attribute
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

        # numpy vector of length 2, (xprime,yprime) for use in manually writing
        # the dot product needed for the exponent in the transform theorem.  Use
        # this 2vec to dot with (x,y) in fromfunc to create the 'phase argument'
        # Since this uses an offset xo yo in pixels of the affine transformation,
        # these are *NOT* affected by the 'oversample' in image space.  The
        # vector it is dotted with is in image space.
        self.phase_2vector = np.array((my * xo - sx * yo, mx * yo - sy * xo)) / self.determinant

    def forward(self, point):
        """
        Create the forward affine transformation, in ideal-to-distorted coordinates.

        Parameters
        ----------
        point : float, float
            Coordinates in ideal space, which to apply forward transform

        Returns
        -------
        trans_point : float, float
            Coordinates in distorted space
        """
        trans_point = np.array(
            (
                self.mx * point[0] + self.sx * point[1] + self.xo,
                self.my * point[1] + self.sy * point[0] + self.yo,
            )
        )

        return trans_point

    def reverse(self, point):
        """
        Create the reverse affine transformation, in distorted-to-ideal coordinates.

        Parameters
        ----------
        point : float, float
            Coordinates in distorted space, which to apply reverse transform

        Returns
        -------
        trans_point : float, float
            Coordinates in ideal space
        """
        trans_point = (
            np.array(
                (
                    self.my * point[0] - self.sx * point[1] - self.my * self.xo + self.sx * self.yo,
                    self.mx * point[1] - self.sy * point[0] - self.mx * self.yo + self.sy * self.xo,
                )
            )
            * self.determinant
        )

        return trans_point

    def distort_f_args(self, u, v):
        """
        Implement the (u,v) to (u',v') change in arguments of F.

        See class documentation of Bracewell Fourier 2D affine transformation theorem.

        Parameters
        ----------
        u : float
            1st argument of F
        v : float
            2nd argument of F

        Returns
        -------
        uprime : float
            1st transformed argument of F
        vprime : float
            2nd transformed argument of F
        """
        uprime = (self.my * u - self.sy * v) / self.determinant
        vprime = (-self.sx * u + self.mx * v) / self.determinant
        return uprime, vprime

    def distortphase(self, u, v):
        """
        Calculate the phase term in the Bracewell Fourier 2D affine transformation theorem.

        The phase term is:

        1/|Delta| * exp{(2*Pi*i/Delta) * [(my*xo- x*yo) * u + (mx*yo-sy*xo)*v]}

        where u and v are in inverse length units.

        Parameters
        ----------
        u : float
            1st argument of F, in units of inverse length units
        v : float
            2nd argument of F, in units of inverse length units

        Returns
        -------
        phase : complex array
            Phase term divided by the determinant.
        """
        phase = np.exp(
            2
            * np.pi
            * 1j
            / self.determinant
            * (self.phase_2vector[0] * u + self.phase_2vector[1] * v)
        )

        return phase

    def get_rotd(self):
        """
        Calculate the rotation that was used to creat a pure rotation affine2d object.

        Returns
        -------
        rotd : float
            Rotation used to creat a pure rotation affine2d
        """
        if self.rotradccw:
            rotd = 180.0 * self.rotradccw / np.pi
            return rotd
        else:
            return None


def makedisk(n, r):
    """
    Calculate a 'disk'.

    Disk is defined as an array whose values =1 in a circular region near
    the center of the array, and =0 elsewhere.

    Parameters
    ----------
    n : int
        Size of 1 dimension of the array to be returned
    r : int
        Radius of disk

    Returns
    -------
    array : 2D integer array
        Array whose values =1 in a circular region near the center of the
        array, and =0 elsewhere.
    """
    if n % 2 == 1:  # odd
        m = (n - 1) / 2
        xx = np.linspace(-m, m, n)
    if n % 2 == 0:  # even
        m = n / 2
        xx = np.linspace(-m + 0.5, m - 0.5, n)

    (x, y) = np.meshgrid(xx, xx.T)
    rad = np.sqrt((x**2) + (y**2))
    array = np.zeros((n, n))
    array[rad < r] = 1

    return array


def avoidhexsingularity(rotation):
    """
    Avoid rotation of exact multiples of 15 degrees to avoid NaNs in hextransformee().

    Parameters
    ----------
    rotation : float or int
       Rotation in degrees

    Returns
    -------
    rotation_adjusted : float
        Replacement value for rotation with epsilon = 1.0e-12 degrees added.
        Precondition before using rotationdegrees in Affine2d for hex geometries
    """
    diagnostic = rotation / 15.0 - int(rotation / 15.0)
    epsilon = 1.0e-12
    if abs(diagnostic) < epsilon / 2.0:
        rotation_adjusted = rotation + epsilon
    else:
        rotation_adjusted = rotation
    return rotation_adjusted


def centerpoint(s):
    """
    Calculate center of image, accounting for odd/even pixel size.

    Parameters
    ----------
    s : 2D int or float tuple
        Array shape

    Returns
    -------
    center : 2D int or float tuple
        Center of image
    """
    return (0.5 * s[0] - 0.5, 0.5 * s[1] - 0.5)


def min_distance_to_edge(img):
    """
    Calculate distance from the brightest pixel in img to the nearest edge of img.

    Parameters
    ----------
    img : 2D array
        Input array

    Returns
    -------
    peakx, peaky : integer, integer
        Coordinates of the peak pixel
    h : integer
        Distance to the nearest image edge
    """
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

    return peakx, peaky, h


def find_centroid(a):
    """
    Calculate the centroid of the image.

    Parameters
    ----------
    a : 2D square float
        Input image array

    Returns
    -------
    htilt, vtilt : float, float
        Centroid of a, as offset from array center, as calculated by the DFT's.

    Notes
    -----
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

    Returns the image centroid (0.40036, 0.2000093) pixels in image space. To
    use this in lg_model, nrm_core,... you will want to calculate the new image
    center using:
    image_center = utils.centerpoint(s) + np.array((centroid[1], centroid[0])
    and everything holds together sensibly looking at DS9 images of a.
    """
    cv = matrix_dft(a, a.shape[0], a.shape[0], centering="ADJUSTABLE")
    cvmod, cvpha = np.abs(cv), np.angle(cv)
    cvmod = cvmod / cvmod.max()  # normalize to unity peak
    htilt, vtilt = findslope(cvpha)
    return htilt, vtilt


def quadratic_extremum(p):
    """
    Calculate extremum of the quadratic.

    Parameters
    ----------
    p : float, float, float
        Quadratic coefficients p[0]*x*x + p[1]*x + p[2]

    Returns
    -------
    float
        Extremum of the quadratic
    """
    return -p[1] / (2.0 * p[0]), -p[1] * p[1] / (4.0 * p[0]) + p[2]


def findpeak_1d(xvec, yvec):
    """
    Calculate the fit function extreme for a given input vector.

    Parameters
    ----------
    xvec : 1D float array
       Input vector
    yvec : 1D float array
       Function values for input vector

    Returns
    -------
    quad_ext : float
        Fit function extreme for a given input vector
    """
    p = np.polyfit(np.array(xvec), np.array(yvec), 2)
    return quadratic_extremum(p)


def findslope(a):
    """
    Find slopes of an array.

    Parameters
    ----------
    a : 2D float array
        Phase slope array

    Returns
    -------
    slopes : 2D float array
        Slopes

    Notes
    -----
    Find slopes of an array, over pixels not bordering the edge of the array.
    There should be valid data on either side of every pixel selected by the mask
    m. a is in radians of phase (in Fourier domain) when used in NRM/KP
    applications. The std dev of the middle 9 pixels is used to clean
    the mask of invalid slope data, where we're subtracting
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

    Gain between rad/pix phase slope and original domain pixels is
         a.shape[0 or 1]/(2 pi)
    Multiply the measured phase slope by this gain for pixels of incoming array
         centroid shift away from array center.
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
    c = (int(c[0]), int(c[1]))
    sigh, sigv = (
        tilt[0][c[0] - 1 : c[0] + 1, c[1] - 1 : c[1] + 1].std(),
        tilt[1][c[0] - 1 : c[0] + 1, c[1] - 1 : c[1] + 1].std(),
    )
    avgh, avgv = (
        tilt[0][c[0] - 1 : c[0] + 1, c[1] - 1 : c[1] + 1].mean(),
        tilt[1][c[0] - 1 : c[0] + 1, c[1] - 1 : c[1] + 1].mean(),
    )

    # second stage mask cleaning: 5 sig rejection of mask
    newmaskh = np.where(np.abs(tilt[0] - avgh) < 5 * sigh)
    newmaskv = np.where(np.abs(tilt[1] - avgv) < 5 * sigv)

    th, tv = np.zeros(a.shape), np.zeros(a.shape)
    th[newmaskh] = tilt[0][newmaskh]
    tv[newmaskv] = tilt[1][newmaskv]

    # determine units of tilt -
    G = a.shape[0] / (2.0 * np.pi), a.shape[1] / (2.0 * np.pi)  # noqa: N806

    slopes = G[0] * tilt[0][newmaskh].mean(), G[1] * tilt[1][newmaskv].mean()
    return slopes


def make_a(nh):
    """
    Write the 'NRM matrix'.

    The NRM matrix later (?) gets pseudo-inverted to provide (arbitrarily constrained)
    zero-mean phases of the holes.
    Algorithm is taken verbatim from Anand's pseudoinverse.py

    Parameters
    ----------
    nh : int
        Number of holes in NR mask

    Returns
    -------
    matrix_a: 2D float array
        Shape nh columns, nh(nh-1)/2 rows (eg 21 for nh=7)

    Notes
    -----
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

    which is implemented in make_a() as:
        matrix_a[row,h2] = -1
        matrix_a[row,h1] = +1

    To change the convention just reverse the signs of the 'ones'.

    When tested against Alex'' nrm_model.py 'piston_phase' text output
    of fringe phases, these signs appear to be correct -
    anand@stsci.edu 12 Nov 2014
    """
    log.debug("-------")
    log.debug(" make_a:")

    ncols = (nh * (nh - 1)) // 2
    nrows = nh
    matrix_a = np.zeros((ncols, nrows))

    row = 0
    for h2 in range(nh):
        for h1 in range(h2 + 1, nh):
            if h1 >= nh:
                break
            else:
                log.debug(" row: %s, h1: %s, h2: %s", row, h1, h2)

                matrix_a[row, h2] = -1
                matrix_a[row, h1] = +1
                row += 1

    log.debug("matrix_a:")
    log.debug(" %s", matrix_a)

    return matrix_a


def fringes2pistons(fringephases, nholes):
    """
    Extract pistons out of fringes.

    For nrm_model.py to use to extract pistons out of fringes, given
    its hole bookkeeping, which apparently matches that of this module,
    and is the same as Noah Gamper's.

    Parameters
    ----------
    fringephases : 1D int array
        Fringe phases
    nholes : int
        Number of holes

    Returns
    -------
    np.dot(Apinv, fringephases) : 1D int array
        Pistons in same units as fringe phases
    """
    a_nrm = make_a(nholes)
    a_p_inv = np.linalg.pinv(a_nrm)

    return np.dot(a_p_inv, fringephases)


def rebin(a=None, rc=(2, 2)):
    """
    Perform simple flux-conserving binning using specified binning kernel.

    Trailing size mismatch is clipped: eg a 10x3 array binned by
    3 results in a 3x1 array

    Parameters
    ----------
    a : 2D float array
        Input array to bin
    rc : 2D float array
        Binning kernel

    Returns
    -------
    binned_arr : float array
        Binned array
    """
    binned_arr = krebin(a, (a.shape[0] // rc[0], a.shape[1] // rc[1]))

    return binned_arr


def krebin(a, shape):
    """
    Rebin, applying Klaus P's fastrebin from web.

    Parameters
    ----------
    a : 2D float array
        Input array to rebin
    shape : tuple (int, int)
        Dimensions of array 'a' binned down by dimensions of binning kernel

    Returns
    -------
    reshaped_a: 2D float array
        Reshaped input array
    """
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    reshaped_a = a.reshape(sh).sum(-1).sum(1)

    return reshaped_a


def rcrosscorrelate(a=None, b=None):
    """
    Calculate cross correlation of two identically-shaped real arrays.

    Parameters
    ----------
    a : 2D float array
        First input array
    b : 2D float array
        Second input array

    Returns
    -------
    c.real.copy() :
        Real part of array that is the correlation of the two input arrays.
    """
    c = crosscorrelate(a=a, b=b) / (np.sqrt((a * a).sum()) * np.sqrt((b * b).sum()))
    return c.real.copy()


def crosscorrelate(a=None, b=None):
    """
    Calculate cross correlation of two identically-shaped real or complex arrays.

    Parameters
    ----------
    a : 2D complex float array
        First input array
    b : 2D complex float array
        Second input array

    Returns
    -------
    fft.fftshift(c)
        Complex array that is the correlation of the two input arrays.
    """
    if a.shape != b.shape:
        log.critical("crosscorrelate: need identical arrays")
        return None

    fac = np.sqrt(a.shape[0] * a.shape[1])

    A = fft.fft2(a) / fac  # noqa: N806
    B = fft.fft2(b) / fac  # noqa: N806
    c = fft.ifft2(A * B.conj()) * fac * fac

    log.debug("----------------")
    log.debug(" crosscorrelate:")
    log.debug(" a: %s:", a)
    log.debug(" A: %s:", A)
    log.debug(" b: %s:", b)
    log.debug(" B: %s:", B)
    log.debug(" c: %s:", c)
    log.debug(" a.sum: %s:", a.sum())
    log.debug(" b.sum: %s:", b.sum())
    log.debug(" c.sum: %s:", c.sum())
    log.debug(" a.sum*b.sum: %s:", a.sum() * b.sum())
    log.debug(" c.sum.real: %s:", c.sum().real)
    log.debug(" a.sum*b.sum/c.sum.real: %s:", a.sum() * b.sum() / c.sum().real)

    return fft.fftshift(c)


def rotate2dccw(vectors, thetarad):
    """
    Apply a CCW rotation to the given vectors.

    For a positive (CCW) rotation:
        x decreases under slight rotation
        y increases under slight rotation

    Parameters
    ----------
    vectors : np.ndarray[float]
       List of 2-D vectors, so the shape is Nx2
    thetarad : float
       Rotation to apply in radians

    Returns
    -------
    np.ndarray[float]
        Rotated vectors
    """
    c, s = (np.cos(thetarad), np.sin(thetarad))
    ctrs_rotated = []
    for vector in vectors:
        ctrs_rotated.append([c * vector[0] - s * vector[1], s * vector[0] + c * vector[1]])
    rot_vectors = np.array(ctrs_rotated)
    return rot_vectors


def get_filt_spec(throughput_model):
    """
    Load filter throughput data into synphot spectrum object.

    Parameters
    ----------
    throughput_model : ThroughputModel
        Datamodel containing normalized fractional throughput
        data for one of the four AMI filters

    Returns
    -------
    band : synphot Spectrum object
        Filter bandpass
    """
    thruput = throughput_model.filter_table
    wl_list = np.asarray([tup[0] for tup in thruput])  # angstroms
    tr_list = np.asarray([tup[1] for tup in thruput])
    band = synphot.spectrum.SpectralElement(
        synphot.models.Empirical1D, points=wl_list, lookup_table=tr_list, keep_neg=False
    )
    return band


def get_flat_spec():
    """
    Produce a synphot spectrum object with constant (unity) flux.

    Returns
    -------
    synphot Spectrum object
        Spectrum with constant flux
    """
    return synphot.SourceSpectrum(synphot.models.ConstFlux1D, amplitude=1)


def combine_src_filt(bandpass, srcspec, trim=0.01, nlambda=19):
    """
    Get the observed spectrum through a filter.

    Largely copied from Poppy instrument.py
    Define nlambda bins of wavelengths, calculate effstim for each, normalize by effstim total.
    nlambda should be calculated so there are ~10 wavelengths per resolution element
    (19 should work)

    Parameters
    ----------
    bandpass : synphot Spectrum
        Filter bandpass (from get_filt_spec)
    srcspec : synphot Spectrum
         Source spectrum (from get_src_spec)
    trim : float, None
        If not None, trim bandpass to where throughput greater than trim
    nlambda : int
        Number of wavelengths across filter to return

    Returns
    -------
    finalsrc : numpy array
        Array of shape (nlambda,2) containing wavelengths, final throughputs
    """
    wl_filt, th_filt = bandpass._get_arrays(bandpass.waveset)  # noqa: SLF001

    if trim:
        log.debug(f"Trimming bandpass to above {trim:.1e} throughput")
        goodthru = np.where(np.asarray(th_filt) > trim)
        low_idx, high_idx = goodthru[0][0], goodthru[0][-1]
        wl_filt, th_filt = wl_filt[low_idx:high_idx], th_filt[low_idx:high_idx]
    ptsin = len(wl_filt)
    if nlambda is None:
        nlambda = ptsin  # Don't bin throughput
    # get effstim for bins of wavelengths
    minwave, maxwave = wl_filt.min(), wl_filt.max()  # trimmed or not
    wave_bin_edges = np.linspace(minwave, maxwave, nlambda + 1)
    wavesteps = (wave_bin_edges[:-1] + wave_bin_edges[1:]) / 2
    deltawave = wave_bin_edges[1] - wave_bin_edges[0]
    area = 1 * (u.m * u.m)
    effstims = []

    binfac = ptsin // nlambda
    log.debug(f"Binning spectrum by {binfac:d} from {ptsin:d} points to {nlambda:d} points")
    for wave in wavesteps:
        log.debug(
            f"\t Integrating across band centered at {wave.to(u.micron):.2f} "
            f"with width {deltawave.to(u.micron):.2f}"
        )
        box = (
            synphot.spectrum.SpectralElement(
                synphot.models.Box1D, amplitude=1, x_0=wave, width=deltawave
            )
            * bandpass
        )

        binset = np.linspace(wave - deltawave, wave + deltawave, 30)
        binset = binset[binset >= 0]  # remove any negative values
        result = synphot.observation.Observation(srcspec, box, binset=binset).effstim(
            "count", area=area
        )
        effstims.append(result)

    effstims = u.Quantity(effstims)
    effstims /= effstims.sum()  # Normalized count rate (total=1) is unitless
    wave_m = wavesteps.to_value(u.m)  # convert to meters
    effstims = effstims.to_value()  # strip units

    finalsrc = np.array((effstims, wave_m)).T  # this is the order expected by InstrumentData

    return finalsrc


def get_cw_beta(bandpass):
    """
    Convert input bandpass array into format expected by code.

    Format is:
    Weighted mean wavelength of each wavelength bin, fractional bandpass

    Parameters
    ----------
    bandpass : array
        Array of weights, wavelengths

    Returns
    -------
    bandpass : array
        Weighted mean wavelength in meters, fractional bandpass
    """
    wt = bandpass[:, 0]
    wl = bandpass[:, 1]
    cw = (wl * wt).sum() / wt.sum()  # Weighted mean wavelength in meters "central wavelength"
    area = simpson(wt, x=wl)
    ew = area / wt.max()  # equivalent width
    beta = ew / cw  # fractional bandpass
    return cw, beta


def handle_bandpass(bandpass, throughput_model):
    """
    Determine what to do with the input bandpass.

    If user-provided, return in format expected by code. If none,
    fetch filter throughput and combine with flat spectrum to
    produce appropriate array.

    Parameters
    ----------
    bandpass : Synphot spectrum or array, or None
        User-defined bandpass to override filter/source.
        If bandpass is an array, wavelengths must have units of meters.
    throughput_model : ThroughputModel
        Datamodel containing filter throughput info.
        Wavelengths are in Angstroms.
        Will not be used if bandpass is not None.

    Returns
    -------
    bandpass : array
        Array of weights, wavelengths used to generate model.
        Wavelengths are in meters.
    """
    # user-defined bandpass can be synphot object or appropriate array
    if bandpass is not None:
        if isinstance(bandpass, synphot.spectrum.SpectralElement):
            log.info("User-defined synphot spectrum provided")
            wl, wt = bandpass._get_arrays(bandpass.waveset.to(u.m))  # noqa: SLF001
            return np.array((wt, wl)).T
        return np.array(bandpass)

    # Default behavior: get the filter and source spectrum
    log.info(f"Reading throughput model data for {throughput_model.meta.instrument.filter}.")
    filt_spec = get_filt_spec(throughput_model)
    log.info("Using flat spectrum model.")
    flat_spec = get_flat_spec()
    nspecbin = 19  # how many wavelngth bins used across bandpass -- affects runtime
    bandpass = combine_src_filt(
        filt_spec,
        flat_spec,
        trim=0.01,
        nlambda=nspecbin,
    )
    return bandpass


def _cdmatrix_to_sky(vec, cd11, cd12, cd21, cd22):
    """
    Convert the CD matrix into RA, Dec pixel scale.

    Parameters
    ----------
    vec : array
        2d, units of pixels
    cd11 : float
        Linear transform matrix element (axis 1 w.r.t x)
    cd12 : float
        Linear transform matrix element (axis 1 w.r.t y)
    cd21 : float
        Linear transform matrix element (axis 2 w.r.t x)
    cd22 : float
        Linear transform matrix element (axis 2 w.r.t y)

    Returns
    -------
    np.ndarray
        Array containing x pixel scale vector, y pixel scale vector

    Notes
    -----
    Use the global header values explicitly, for clarity.
    CD inputs are 4 scalars, conceptually 2x2 array in units degrees/pixel
    """
    return np.array((cd11 * vec[0] + cd12 * vec[1], cd21 * vec[0] + cd22 * vec[1]))


def degrees_per_pixel(datamodel):
    """
    Get pixel scale info from data model.

    If it fails to find the right keywords, use 0.0656 as/pixel

    Parameters
    ----------
    datamodel : datamodel object
        NIRISS data model containing WCS information

    Returns
    -------
    float
        Pixel scale in degrees/pixel
    """
    wcsinfo = datamodel.meta.wcsinfo._instance  # noqa: SLF001
    if "cd1_1" in wcsinfo and "cd1_2" in wcsinfo and "cd2_1" in wcsinfo and "cd2_2" in wcsinfo:
        cd11 = datamodel.meta.wcsinfo.cd1_1
        cd12 = datamodel.meta.wcsinfo.cd1_2
        cd21 = datamodel.meta.wcsinfo.cd2_1
        cd22 = datamodel.meta.wcsinfo.cd2_2
        # Create unit vectors in detector pixel X and Y directions, units: detector pixels
        dxpix = np.array((1.0, 0.0))  # axis 1 step
        dypix = np.array((0.0, 1.0))  # axis 2 step
        # transform pixel x and y steps to RA-tan, Dec-tan degrees
        dxsky = _cdmatrix_to_sky(dxpix, cd11, cd12, cd21, cd22)
        dysky = _cdmatrix_to_sky(dypix, cd11, cd12, cd21, cd22)
        log.debug("Used CD matrix for pixel scales")
        return np.linalg.norm(dxsky, ord=2), np.linalg.norm(dysky, ord=2)
    elif "cdelt1" in wcsinfo and "cdelt2" in wcsinfo:
        log.debug("Used CDELT[12] for pixel scales")
        return datamodel.meta.wcsinfo.cdelt1, datamodel.meta.wcsinfo.cdelt2
    else:
        log.warning("WARNING: NIRISS pixel scales not in header.  Using 65.6 mas in deg/pix")
        return 65.6 / (60.0 * 60.0 * 1000), 65.6 / (60.0 * 60.0 * 1000)

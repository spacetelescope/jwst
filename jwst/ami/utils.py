from stdatamodels.jwst.datamodels import dqflags

from . import matrix_dft

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

        """
        numpy vector of length 2, (xprime,yprime) for use in manually writing
        the dot product needed for the exponent in the transform theorem.  Use
        this 2vec to dot with (x,y) in fromfunc to create the 'phase argument'
        Since this uses an offset xo yo in pixels of the affine transformation,
        these are *NOT* affected by the 'oversample' in image space.  The
        vector it is dotted with is in image space."""
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


def affinepars2header(hdr, affine2d):
    """
    Write the affine2d parameters into fits header (will be modified or deleted in later build).

    Parameters
    ----------
    hdr : fits header
        FITS header to write affine2d parameters into

    affine2d : Affine2d object
        The affine2d object to write into

    Returns
    -------
    hdr : fits header
        FITS header, updated with affine2d parameters
    """
    hdr["affine"] = (affine2d.name, "Affine2d in pupil: name")
    hdr["aff_mx"] = (affine2d.mx, "Affine2d in pupil: xmag")
    hdr["aff_my"] = (affine2d.my, "Affine2d in pupil: ymag")
    hdr["aff_sx"] = (affine2d.sx, "Affine2d in pupil: xshear")
    hdr["aff_sy"] = (affine2d.sx, "Affine2d in pupil: yshear")
    hdr["aff_xo"] = (affine2d.xo, "Affine2d in pupil: x offset")
    hdr["aff_yo"] = (affine2d.yo, "Affine2d in pupil: y offset")
    hdr["aff_dev"] = ("analyticnrm2", "dev_phasor")

    return hdr


def makedisk(n, r, ctr=(0, 0)):
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

    ctr : (int, int)
        Center of disk

    Returns
    -------
    array : 2D integer array
        Array whose values =1 in a circular region near the center of the
        array, and =0 elsewhere.
    """
    if n % 2 == 1:  # odd
        m = (n - 1) / 2
        xx = np.linspace(-m - ctr[0], m - ctr[0], n)
        yy = np.linspace(-m - ctr[1], m - ctr[1], n)
    if n % 2 == 0:  # even
        m = n / 2
        xx = np.linspace(-m - ctr[0], m - ctr[0] - 1, n)
        yy = np.linspace(-m - ctr[1], m - ctr[1] - 1, n)

    (x, y) = np.meshgrid(xx, yy.T)
    rad = np.sqrt((x**2) + (y**2))
    array = np.zeros((n, n))
    array[rad < r] = 1

    return array


def trim(m, s):
    """
    Remove the edge pixels from an index mask m.

    Parameters
    ----------
    m : (int, int) array
        2d index mask

    s : int
        Side of the parent array that was used to generate m.

    Returns
    -------
    m_masked : (integer, integer) array
        2d index mask with edge pixels trimmed
    """
    xl, yl = [], []  # trimmed lists
    for ii in range(len(m[0])):
        # Go through all indices in the mask:
        # the x & y lists test for any index being an edge index - if none are
        # on the edge, remember the indices in new list
        if (m[0][ii] == 0 or m[1][ii] == 0 or m[0][ii] == s - 1 or m[1][ii] == s - 1) is False:
            xl.append(m[0][ii])
            yl.append(m[1][ii])
    m_masked = (np.asarray(xl), np.asarray(yl))

    return m_masked


def avoidhexsingularity(rotation):
    """
    Avoid rotation of exact multiples of 15 degrees to avoid NaNs in hextransformee().

    Parameters
    ----------
    rotation : float
       Rotation in degrees int or float

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


def center_imagepeak(img, r="default", cntrimg=True):
    """
    Calculate a cropped version of the input image centered on the peak pixel.

    Parameters
    ----------
    img : 2D float array
        Input image array

    r : int
        Offset for center determination

    cntrimg : bool
        If True, center on the peak pixel

    Returns
    -------
    cropped : 2D float array
        Cropped to place the brightest pixel at the center of the img array
    """
    peakx, peaky, h = min_distance_to_edge(img, cntrimg=cntrimg)
    log.debug(" peakx=%g, peaky=%g, distance to edge=%g", peakx, peaky, h)
    if r == "default":
        r = h.copy()
    else:
        pass

    cropped = img[int(peakx - r) : int(peakx + r + 1), int(peaky - r) : int(peaky + r + 1)]

    return cropped


def centerpoint(s):
    """
    Calculate center of image, accounting for odd/even pixel size.

    Used for jinc() and hex transform functions.

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


def min_distance_to_edge(img, cntrimg=False):
    """
    Calculate distance from the brightest pixel in img to the nearest edge of img.

    Parameters
    ----------
    img : 2D array
        Input array

    cntrimg : bool
        If True, only look for the peak pixel near the center of the image

    Returns
    -------
    peakx, peaky : integer, integer
        Coordinates of the peak pixel

    h : integer
        Distance to the nearest image edge
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

    returns the image centroid (0.40036, 0.2000093) pixels in image space. To
    use this in lg_model, nrm_core,... you will want to calculate the new image
    center using:
    image_center = utils.centerpoint(s) + np.array((centroid[1], centroid[0])
    and everything holds together sensibly looking at DS9 images of a.
    """
    ft = matrix_dft.MatrixFourierTransform()

    cv = ft.perform(a, a.shape[0], a.shape[0])
    cvmod, cvpha = np.abs(cv), np.angle(cv)

    cvmod = cvmod / cvmod.max()  # normalize to unity peak

    htilt, vtilt = findslope(cvpha)

    return htilt, vtilt


def quadratic_extremum(p):
    """
    Calculate maximum of the quadratic.

    Parameters
    ----------
    p : float, float, float
        Quadratic coefficients

    Returns
    -------
    y_max : float
        Maximum of the quadratic
    """
    y_max = -p[1] / (2.0 * p[0]), -p[1] * p[1] / (4.0 * p[0]) + p[2]

    return y_max


def findpeak_1d(yvec, xvec):
    """
    Calculate the fit function extreme for a given input vector.

    Parameters
    ----------
    yvec : 1D float array
       Function values for input vector

    xvec : 1D float array
       Input vector

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


def quadratic(p, x):
    """
    Calculate value of x at min or max value of y given coefficients of quadratic function.

    Parameters
    ----------
    p : float, float, float
        Coefficients of quadratic function: p[0]*x*x + p[1]*x + p[2]

    x : 1D float array
        Arguments of p()

    Returns
    -------
    maxx : float
        Value of x at minimum or maximum value of y

    maxy : float
        Max y = -b^2/4a occurs at x = -b^2/2a

    fit_val : 1D float array
        Values of quadratic function at arguments in x array
    """
    maxx = -p[1] / (2.0 * p[0])
    maxy = -p[1] * p[1] / (4.0 * p[0]) + p[2]
    fit_val = p[0] * x * x + p[1] * x + p[2]

    return maxx, maxy, fit_val


def make_a(nh):
    """
    Write the 'NRM matrix'.

    The NRM matrix later (?) gets pseudo-inverted to provide (arbitrarily constrained)
    zero-mean phases of the holes. TODO: check with Rachel the inversion happens outside
    the function!

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


def lambdasteps(lam, frac_width, steps=4):
    """
    Create array of increments of lambda.

    Parameters
    ----------
    lam : float
        Lambda

    frac_width : float
        Fractional bandwidth

    steps : int
        With lam and frac, determines bin size of lambda array

    Returns
    -------
    lambda_array : 1D float array
        Array of increments of lambda
    """
    frac = frac_width / 2.0
    steps = steps / 2.0

    # add some very small number to the end to include the last number.
    lambda_array = np.arange(-1 * frac * lam + lam, frac * lam + lam + 10e-10, frac * lam / steps)

    return lambda_array


def tophatfilter(lam_c, frac_width, npoints=10):
    """
    Create tophat filter list from array of lambda values.

    Parameters
    ----------
    lam_c : float
        Lambda

    frac_width : float
        Fractional bandwidth

    npoints : int
        Number of bins in lambda array

    Returns
    -------
    filt : list
        Tophat filter list
    """
    wllist = lambdasteps(lam_c, frac_width, steps=npoints)
    filt = []
    for ii in range(len(wllist)):
        filt.append(np.array([1.0, wllist[ii]]))
    return filt


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
    vectors : list
       2D vectors

    thetarad : float
       Rotation

    Returns
    -------
    rot_vectors : array of floats
        Rotated vectors
    """
    c, s = (np.cos(thetarad), np.sin(thetarad))
    ctrs_rotated = []
    for vector in vectors:
        ctrs_rotated.append([c * vector[0] - s * vector[1], s * vector[0] + c * vector[1]])
    rot_vectors = np.array(ctrs_rotated)

    return rot_vectors


def findmax(mag, vals, mid=1.0):
    """
    Find values at extrema of input vals using quadratic fit to mags, vals.

    Parameters
    ----------
    mag : 1D float array
        Array for abscissa

    vals : 1D float array
        Array for ordinate

    mid : float
        Midpoint of range

    Returns
    -------
    maxx : float
        Value of mag at the extreme value of vals

    maxy : float
        Value of vals corresponding to maxx
    """
    p = np.polyfit(mag, vals, 2)
    fitr = np.arange(0.95 * mid, 1.05 * mid, 0.01)
    maxx, maxy, fitc = quadratic(p, fitr)

    return maxx, maxy


def pix_median_fill_value(input_array, input_dq_array, bsize, xc, yc):
    """
    Calculate the median value of good values within the box of neighboring pixels.

    For the pixel specified by (xc, yc), calculate the median value of the
    good values within the box of size bsize neighboring pixels. If any of
    the box is outside the data, 0 will be returned.

    Parameters
    ----------
    input_array : ndarray
        2D input array to filter
    input_dq_array : ndarray
        2D input data quality array
    bsize : scalar
        Square box size of the data to extract
    xc : scalar
        X position of the data extraction
    yc : scalar
        Y position of the data extraction

    Returns
    -------
    median_value : float
        Median value of good values within box of neighboring pixels
    """
    # set the half box size
    hbox = int(bsize / 2)

    # Extract the region of interest for the data
    try:
        data_array = input_array[yc - hbox : yc + hbox + 1, xc - hbox : xc + hbox + 1]
        dq_array = input_dq_array[yc - hbox : yc + hbox + 1, xc - hbox : xc + hbox + 1]
    except IndexError:
        # If the box is outside the data, return 0
        log.warning("Box for median filter is outside the data")
        return 0.0

    # only keep pixels not flagged with DO_NOT_USE
    wh_good = np.where(np.bitwise_and(dq_array, dqflags.pixel["DO_NOT_USE"]) == 0)
    filtered_array = data_array[wh_good]

    # compute the median, excluding NaN's
    median_value = np.nanmedian(filtered_array)

    # check for bad result
    if np.isnan(median_value):
        log.warning("Median filter returned NaN; setting value to 0.")
        median_value = 0.0

    return median_value


def mas2rad(mas):
    """
    Convert angle in milli arc-sec to radians.

    Parameters
    ----------
    mas : float
        Angle in milli arc-sec

    Returns
    -------
    rad : float
        Angle in radians
    """
    rad = mas * (10 ** (-3)) / (3600 * 180 / np.pi)
    return rad


def img_median_replace(img_model, box_size):
    """
    Replace bad pixels with the median value of surrounding good pixels.

    Bad pixels may arise here either due to a DQ value of DO_NOT_USE or having a
    value of NaN.

    Parameters
    ----------
    img_model : image model
        Image model containing input array to filter.

    box_size : scalar
        Box size for the median filter

    Returns
    -------
    img_model : datamodel
        Input image model whose input array has its bad pixels replaced
        by the median of the surrounding good-value pixels.
    """
    input_data = img_model.data
    input_dq = img_model.dq

    num_nan = np.count_nonzero(np.isnan(input_data))
    num_dq_bad = np.count_nonzero(input_dq == dqflags.pixel["DO_NOT_USE"])

    # check to see if any of the pixels are bad
    if num_nan + num_dq_bad > 0:
        log.info(f"Applying median filter for {num_nan} NaN and {num_dq_bad} DO_NOT_USE pixels")
        bad_locations = np.where(
            np.isnan(input_data) | np.equal(input_dq, dqflags.pixel["DO_NOT_USE"])
        )

        # fill the bad pixel values with the median of the data in a box region
        for i_pos in range(len(bad_locations[0])):
            y_box_pos = bad_locations[0][i_pos]
            x_box_pos = bad_locations[1][i_pos]
            median_fill = pix_median_fill_value(
                input_data, input_dq, box_size, x_box_pos, y_box_pos
            )
            input_data[y_box_pos, x_box_pos] = median_fill

        img_model.data = input_data

    return img_model


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
    flatspec : synphot Spectrum object
        Spectrum with constant flux
    """
    flatspec = synphot.SourceSpectrum(synphot.models.ConstFlux1D, amplitude=1)

    return flatspec


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
        User-defined bandpass to override filter/source
    throughput_model : ThroughputModel
        Datamodel containing filter throughput info.
        Will not be used if bandpass is not None.

    Returns
    -------
    bandpass : array
        Array of weights, wavelengths used to generate model
    """
    if bandpass is not None:
        # bandpass can be user-defined synphot object or appropriate array
        if isinstance(bandpass, synphot.spectrum.SpectralElement):
            log.info("User-defined synphot spectrum provided")
            wl, wt = bandpass._get_arrays(bandpass.waveset)  # noqa: SLF001
            bandpass = np.array((wt, wl)).T
        else:
            bandpass = np.array(bandpass)

    else:
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

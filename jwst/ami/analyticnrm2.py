# Heritage mathematica nb from Alex & Laurent
# Python by Alex Greenbaum & Anand Sivaramakrishnan Jan 2013
# updated May 2013 to include hexagonal envelope

import logging
import warnings
import numpy as np
import scipy.special
from . import leastsqnrm
from . import utils
from . import hextransformee

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def jinc(x, y, d, lam, pitch, offx=0.0, offy=0.0):
    """
    Compute 2d Jinc for given coordinates.

    Parameters
    ----------
    x, y : float, float
        Input coordinates
    d : float
        Hole diameter
    lam : float
        Wavelength
    pitch : float
        Sampling pitch in radians in image plane
    offx, offy : float, float
        Offsets from image center in detector pixels

    Returns
    -------
    jinc_2d : float array
        2d jinc at the given coordinates, with NaNs replaced by pi/4.
    """
    r = (d / lam) * pitch * np.sqrt((x - offx) ** 2 + (y - offy) ** 2)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=RuntimeWarning, message="invalid value encountered in divide"
        )
        jinc_2d = leastsqnrm.replacenan(scipy.special.jv(1, np.pi * r) / (2.0 * r))
    return jinc_2d


def ffc(kx, ky, ko, baseline, lam, pitch, affine2d):
    """
    Calculate cosine terms of analytic model.

    Parameters
    ----------
    kx, ky : float, float
        X-component and y-component of image plane (spatial frequency) vector
    ko : array-like[float, float]
        Center of PSF.
    baseline : 2D float array
        Hole centers
    lam : float
        Wavelength
    pitch : float
        Sampling pitch in radians in image plane
    affine2d : Affine2d object
        Distortion object.

    Returns
    -------
    cos_array : 2D float array
        Cosine terms of analytic model
    """
    kxprime, kyprime = affine2d.distort_f_args(kx - ko[0], ky - ko[1])

    cos_array = 2 * np.cos(
        2 * np.pi * pitch * (kxprime * baseline[0] + kyprime * baseline[1]) / lam
    )
    return cos_array


def ffs(kx, ky, ko, baseline, lam, pitch, affine2d):
    """
    Calculate sine terms of analytic model.

    Parameters
    ----------
    kx, ky : float, float
        X-component and y-component of image plane (spatial frequency) vector
    ko : array-like[float, float]
        Center of PSF.
    baseline : 2D float array
        Hole centers
    lam : float
        Wavelength
    pitch : float
        Sampling pitch in radians in image plane
    affine2d : Affine2d object
        Distortion object.

    Returns
    -------
    sin_array : 2D float array
        Sine terms of analytic model
    """
    kxprime, kyprime = affine2d.distort_f_args(kx - ko[0], ky - ko[1])

    sin_array = 2 * np.sin(
        2 * np.pi * pitch * (kxprime * baseline[0] + kyprime * baseline[1]) / lam
    )
    return sin_array


def harmonicfringes(fov, pitch, baseline, lam, oversample, affine2d, psf_offset=(0, 0)):
    """
    Calculate the sine and cosine fringes.

    This is in image space and, for later
    versions, this works in the oversampled space that is each slice of the model.

    Parameters
    ----------
    fov : int, default=None
        Number of detector pixels on a side
    pitch : float
        Sampling pitch in radians in image plane
    baseline : 2D float array
        Hole centers vector, units of meters.
    lam : float
        Wavelength in meters.
    oversample : int
        Number of samples per detector pixel pitch
    affine2d : Affine2d object
        The affine2d object
    psf_offset : 2D float array, optional
        Offset from image center in detector pixels, default is (0,0).

    Returns
    -------
    (cosine_fringes, sine_fringes) : tuple
        Sine and cosine fringes: float arrays
    """
    cpitch = pitch / oversample
    im_ctr = image_center(fov, oversample, psf_offset)

    cosine_fringes = np.fromfunction(
        ffc,
        (fov * oversample, fov * oversample),
        ko=im_ctr,
        baseline=baseline,
        lam=lam,
        pitch=cpitch,
        affine2d=affine2d,
    )
    sine_fringes = np.fromfunction(
        ffs,
        (fov * oversample, fov * oversample),
        ko=im_ctr,
        baseline=baseline,
        lam=lam,
        pitch=cpitch,
        affine2d=affine2d,
    )
    return cosine_fringes, sine_fringes


def phasor(kx, ky, hx, hy, lam, phi_m, pitch, affine2d):
    """
    Calculate the wavefront for a single hole.

    This routine returns the complex
    amplitude array of fringes phi to units of meters, which is more physical for
    broadband simulations.

    Parameters
    ----------
    kx, ky : float
        Image plane coords in units of sampling pitch (oversampled, or not)
    hx, hy : float
        Hole centers in meters
    lam : float
        Wavelength
    phi_m : float
        Distance of fringe from hole center in units of meters
    pitch : float
        Sampling pitch in radians in image plane
    affine2d : Affine2d object
        The affine2d object

    Returns
    -------
    phasor : complex
        Calculate wavefront for a single hole
    """
    kxprime, kyprime = affine2d.distort_f_args(kx, ky)
    return np.exp(
        -2 * np.pi * 1j * ((pitch * hx * kxprime + pitch * hy * kyprime) / lam + phi_m / lam)
    ) * affine2d.distortphase(kx, ky)


def image_center(fov, oversample, psf_offset):
    """
    Calculate the Image center location in oversampled pixels.

    Parameters
    ----------
    fov : int
        Number of detector pixels of field of view
    oversample : int
        Number of samples per detector pixel pitch
    psf_offset : 2D int array
        Offset from image center in detector pixels

    Returns
    -------
    offsets_from_center : 2D int array
        Offset of the psf center from the array center.
    """
    offsets_from_center = (
        np.array(utils.centerpoint((oversample * fov, oversample * fov)))
        + np.array((psf_offset[1], psf_offset[0])) * oversample
    )

    return offsets_from_center


def interf(kx, ky, c, ctrs, phi, lam, pitch, affine2d):
    """
    Calculate the complex amplitudes for all holes.

    Parameters
    ----------
    kx, ky : float, float radians
        X-component and y-component of image plane (spatial frequency) vector
    c : 2D float array
        Center of PSF, in simulation pixels (i.e. oversampled)
    ctrs : 2D float array
        Centers of holes
    phi : float
        Distance of fringe from hole center in units of waves
    lam : float
        Wavelength
    pitch : float
        Sampling pitch in radians in image plane
    affine2d : Affine2d object
        The affine2d object

    Returns
    -------
    fringe_complexamp : 2D complex array
        Interference for all holes
    """
    fringe_complexamp = 0j
    for hole, ctr in enumerate(ctrs):
        fringe_complexamp += phasor(
            (kx - c[0]),
            (ky - c[1]),
            ctr[0],
            ctr[1],
            lam,
            phi[hole],
            pitch,
            affine2d,
        )

    return fringe_complexamp


def model_array(
    ctrs,
    lam,
    oversample,
    pitch,
    fov,
    d,
    psf_offset=(0, 0),
    phi=None,
    shape="circ",
    affine2d=None,
):
    """
    Create a model using the specified wavelength.

    Parameters
    ----------
    ctrs : 2D float array
        Centers of holes
    lam : float
        Wavelength in the bandpass for this particular model
    oversample : int
        Oversampling factor
    pitch : float
        Sampling pitch in radians in image plane
    fov : int
        Number of detector pixels on a side.
    d : float
        Hole diameter for 'circ'; flat to flat distance for 'hex
    psf_offset : 2D int array
        Offset from image center in detector pixels
    phi : float
        Distance of fringe from hole center in units of waves
    shape : str
        Shape of hole; possible values are 'circ', 'hex', and 'fringe'
    affine2d : Affine2d object
        The affine2d object

    Returns
    -------
    primary_beam : float 2D array
        Array of primary beam,
    ffmodel : list of fringe arrays
        List of fringe arrays
    """
    nholes = ctrs.shape[0]
    if phi is None:
        np.zeros((nholes,))  # no phase errors in the model slices...
    modelshape = (
        fov * oversample,
        fov * oversample,
    )  # spatial extent of image model - the oversampled array

    # calculate primary beam envelope (non-negative real)
    if shape == "circ":
        asf_pb = asf(pitch, fov, oversample, d, lam, affine2d)
    elif shape == "hex":
        asf_pb = asf_hex(pitch, fov, oversample, d, lam, psf_offset, affine2d)
    else:
        raise KeyError(
            "Must provide a valid hole shape. Current supported shapes are 'circ' and 'hex'."
        )

    primary_beam = (asf_pb * asf_pb.conj()).real

    alist = []
    for i in range(nholes - 1):
        for j in range(nholes - 1):
            if j + i + 1 < nholes:
                alist = np.append(alist, i)
                alist = np.append(alist, j + i + 1)
    alist = alist.reshape((len(alist) // 2, 2))

    ffmodel = []
    ffmodel.append(nholes * np.ones(modelshape))
    for basepair in alist:
        baseline = ctrs[int(basepair[0])] - ctrs[int(basepair[1])]
        cosfringe, sinfringe = harmonicfringes(
            fov=fov,
            pitch=pitch,
            psf_offset=psf_offset,
            baseline=baseline,
            oversample=oversample,
            lam=lam,
            affine2d=affine2d,
        )
        ffmodel.append(cosfringe)
        ffmodel.append(sinfringe)

    return primary_beam, ffmodel


def asf(detpixel, fov, oversample, d, lam, psf_offset):
    """
    Calculate the Amplitude Spread Function for a circular aperture.

    Amplitude Spread Function is also know as image plane complex amplitude.

    Parameters
    ----------
    detpixel : float
        Pixel scale
    fov : int
        Number of detector pixels on a side
    oversample : int
        Oversampling factor
    d : float
        Hole diameter
    lam : float
        Wavelength
    psf_offset : 2D float array
        Offset from image center in detector pixels

    Returns
    -------
    asf : 2D real array
        Amplitude Spread Function (a.k.a. image plane complex amplitude) for
        a circular aperture
    """
    pitch = detpixel / float(oversample)
    im_ctr = image_center(fov, oversample, psf_offset)

    return np.fromfunction(
        jinc,
        (oversample * fov, oversample * fov),
        d=d,
        lam=lam,
        pitch=pitch,
        offx=im_ctr[0],
        offy=im_ctr[1],
    )


def asffringe(detpixel, fov, oversample, ctrs, lam, phi, psf_offset, affine2d):
    """
    Amplitude Spread Function (a.k.a. image plane complex amplitude) for a fringe.

    Parameters
    ----------
    detpixel : float
        Pixel scale
    fov : int
        Number of detector pixels on a side
    oversample : int
        Oversampling factor
    ctrs : 2D float array
        Centers of holes
    lam : float
        Wavelength
    phi : float
        Distance of fringe from hole center in units of waves
    psf_offset : 2D float array
        Offset from image center in detector pixels
    affine2d : Affine2d object
        The affine2d object

    Returns
    -------
    fringing : 2D complex array
        Amplitude Spread Function (a.k.a. image plane complex amplitude) for
        a fringe
    """
    pitch = detpixel / float(oversample)
    im_ctr = image_center(fov, oversample, psf_offset)

    return np.fromfunction(
        interf,
        (oversample * fov, oversample * fov),
        c=im_ctr,
        ctrs=ctrs,
        phi=phi,
        lam=lam,
        pitch=pitch,
        affine2d=affine2d,
    )


def asf_hex(detpixel, fov, oversample, d, lam, psf_offset, affine2d):
    """
    Amplitude Spread Function (a.k.a. image plane complex amplitude) for a hexagonal aperture.

    Parameters
    ----------
    detpixel : float
        Pixel scale
    fov : int
        Number of detector pixels on a side
    oversample : int
        Oversampling factor
    d : float
        Flat-to-flat distance across hexagon
    lam : float
        Wavelength
    psf_offset : 2D float array
        Offset from image center in detector pixels
    affine2d : Affine2d object
        The affine2d object

    Returns
    -------
    asf : 2D complex array
        Amplitude Spread Function (a.k.a. image plane complex amplitude) for
        a hexagonal aperture
    """
    pitch = detpixel / float(oversample)

    im_ctr = image_center(fov, oversample, psf_offset)

    return hextransformee.hextransform(
        s=(oversample * fov, oversample * fov),
        c=im_ctr,
        d=d,
        lam=lam,
        pitch=pitch,
        affine2d=affine2d,
    )


def psf(detpixel, fov, oversample, ctrs, d, lam, phi, psf_offset, affine2d, shape="circ"):
    """
    Calculate the PSF for the requested aperture shape.

    Parameters
    ----------
    detpixel : float
        Pixel scale
    fov : int
        Number of detector pixels on a side
    oversample : int
        Oversampling factor
    ctrs : 2D float array
        Centers of holes
    d : float
        Hole diameter for 'circ'; flat-to-flat distance across for 'hex'
    lam : float
        Wavelength
    phi : float
        Distance of fringe from hole center in units of waves
    psf_offset : 2D float array
        Offset from image center in detector pixels
    affine2d : Affine2d object
        The affine2d object
    shape : str
        Shape of hole; possible values are 'circ', 'circonly', 'hex', 'hexonly', 'fringeonly'

    Returns
    -------
    PSF : 2D float array
        The point-spread function
    """
    if shape == "circ":
        asf_fringe = asffringe(detpixel, fov, oversample, ctrs, lam, phi, psf_offset, affine2d)
        asf_2d = asf(detpixel, fov, oversample, d, lam, psf_offset) * asf_fringe

    elif shape == "circonly":
        asf_2d = asf(detpixel, fov, oversample, d, lam, psf_offset)

    elif shape == "hex":
        asf_fringe = asffringe(detpixel, fov, oversample, ctrs, lam, phi, psf_offset, affine2d)
        asf_2d = asf_hex(detpixel, fov, oversample, d, lam, psf_offset, affine2d) * asf_fringe

    elif shape == "hexonly":
        asf_2d = asf_hex(detpixel, fov, oversample, d, lam, psf_offset, affine2d)

    elif shape == "fringeonly":
        asf_2d = asffringe(detpixel, fov, oversample, ctrs, lam, phi, psf_offset, affine2d)

    else:
        raise ValueError(
            f"pupil shape {shape} not supported - choices: "
            "'circ', 'circonly', 'hex', 'hexonly', and 'fringeonly'"
        )

    return (asf_2d * asf_2d.conj()).real

"""
Matrix-based discrete Fourier transforms for computing PSFs.

MatrixDFT: Matrix-based discrete Fourier transforms for computing PSFs.
Internally this will call one of several subfunctions depending on the
specified centering type. These have to do with where the (0, 0) element of
the Fourier transform is located, i.e. where the PSF center ends up.
- 'FFTSTYLE' centered on one pixel
- 'SYMMETRIC' centered on crosshairs between middle pixel
- 'ADJUSTABLE', always centered in output array depending on
whether it is even or odd
'ADJUSTABLE' is the default.
This module was originally called "Slow Fourier Transform", and this
terminology still appears in some places in the code.  Note that this is
'slow' only in the sense that if you perform the exact same calculation as
an FFT, the FFT algorithm is much faster. However this algorithm gives you
much more flexibility in choosing array sizes and sampling, and often lets
you replace "fast calculations on very large arrays" with "relatively slow
calculations on much smaller ones".

Code originally by A. Sivaramakrishnan
2010-11-05 Revised normalizations for flux conservation consistent
with Soummer et al. 2007. Updated documentation.  -- M. Perrin
2011-2012: Various enhancements, detailed history not kept, sorry.
2012-05-18: module renamed SFT.py -> matrixDFT.py
2012-09-26: minor big fixes
2015-01-21: Eliminate redundant code paths, correct parity flip,
PEP8 formatting pass (except variable names)-- J. Long

References
----------
Soummer et al. 2007, Opt. Express  15, 15935-15951 (2007)
https://doi.org/10.1364/OE.15.015935

Examples
--------
result = matrix_dft.matrix_dft(pupilArray, focalplane_size, focalplane_npix)
"""

__all__ = ["matrix_dft", "matrix_idft"]

import numpy as np

FFTSTYLE = "FFTSTYLE"
FFTRECT = "FFTRECT"
SYMMETRIC = "SYMMETRIC"
ADJUSTABLE = "ADJUSTABLE"
CENTERING_CHOICES = (FFTSTYLE, SYMMETRIC, ADJUSTABLE, FFTRECT)


def matrix_dft(plane, nlam_d, npix, offset=None, inverse=False, centering=FFTSTYLE):
    """
    Perform a matrix discrete Fourier transform with selectable output sampling and centering.

    Where parameters can be supplied as either
    scalars or 2-tuples, the first element of the 2-tuple is used for the
    Y dimension and the second for the X dimension. This ordering matches
    that of numpy.ndarray.shape attributes and that of Python indexing.
    To achieve exact correspondence to the FFT set nlam_d and npix to the size
    of the input array in pixels and use 'FFTSTYLE' centering. (n.b. When
    using `numpy.fft.fft2` you must `numpy.fft.fftshift` the input pupil both
    before and after applying fft2 or else it will introduce a checkerboard
    pattern in the signs of alternating pixels!)

    Parameters
    ----------
    plane : 2D ndarray
        2D array (either real or complex) representing the input image plane or
        pupil plane to transform.
    nlam_d : float or 2-tuple of floats (nlam_dy, nlam_dx)
        Size of desired output region in lambda / D units, assuming that the
        pupil fills the input array (corresponds to 'm' in
        Soummer et al. 2007 4.2). This is in units of the spatial frequency
        that is just Nyquist sampled by the input array.) If given as a tuple,
        interpreted as (nlam_dy, nlam_dx).
    npix : int or 2-tuple of ints (npix_y, npix_x)
        Number of pixels per side side of destination plane array (corresponds
        to 'N_B' in Soummer et al. 2007 4.2). This will be the # of pixels in
        the image plane for a forward transformation, in the pupil plane for an
        inverse. If given as a tuple, interpreted as (npix_y, npix_x).
    offset : 2-tuple of floats (offset_y, offset_x)
        For ADJUSTABLE-style transforms, an offset in pixels by which the PSF
        will be displaced from the central pixel (or cross). Given as
        (offset_y, offset_x).
    inverse : bool, optional
        Is this a forward or inverse transformation? (Default is False,
        implying a forward transformation.)
    centering : {'FFTSTYLE', 'SYMMETRIC', 'ADJUSTABLE'}, optional
        What type of centering convention should be used for this FFT?
        * ADJUSTABLE (the default) For an output array with ODD size n,
          the PSF center will be at the center of pixel (n-1)/2. For an output
          array with EVEN size n, the PSF center will be in the corner between
          pixel (n/2-1, n/2-1) and (n/2, n/2)
        * FFTSTYLE puts the zero-order term in a single pixel.
        * SYMMETRIC spreads the zero-order term evenly between the center
          four pixels

    Returns
    -------
    norm_coeff * t2; float, ndarray
        Normalized FT coeffs
    """
    npup_y, npup_x = plane.shape

    if np.isscalar(npix):
        npix_y, npix_x = npix, npix
    else:
        try:
            npix_y, npix_x = npix
        except ValueError:
            raise ValueError(
                "'npix' must be supplied as a scalar (for square arrays) or as"
                "a 2-tuple of ints (npix_y, npix_x)"
            ) from None

    if np.isscalar(nlam_d):
        nlam_dy, nlam_dx = nlam_d, nlam_d
    else:
        try:
            nlam_dy, nlam_dx = nlam_d
        except ValueError:
            raise ValueError(
                "'nlam_d' must be supplied as a scalar (for square arrays) or"
                " as a 2-tuple of floats (nlam_dy, nlam_dx)"
            ) from None

    centering = centering.upper()

    # In the following: X and Y are coordinates in the input plane
    #                   U and V are coordinates in the output plane
    if inverse:
        dx = nlam_dx / float(npup_x)
        dy = nlam_dy / float(npup_y)
        du = 1.0 / float(npix_x)
        dv = 1.0 / float(npix_y)
    else:
        du = nlam_dx / float(npix_x)
        dv = nlam_dy / float(npix_y)
        dx = 1.0 / float(npup_x)
        dy = 1.0 / float(npup_y)

    if centering == FFTSTYLE:
        xs = (np.arange(npup_x) - (npup_x / 2)) * dx
        ys = (np.arange(npup_y) - (npup_y / 2)) * dy

        us = (np.arange(npix_x) - npix_x / 2) * du
        vs = (np.arange(npix_y) - npix_y / 2) * dv
    elif centering == ADJUSTABLE:
        if offset is None:
            offset_y, offset_x = 0.0, 0.0
        else:
            try:
                offset_y, offset_x = offset
            except (ValueError, TypeError) as e:
                raise ValueError(
                    "'offset' must be supplied as a 2-tuple with "
                    "(y_offset, x_offset) as floating point values"
                ) from e
        xs = (np.arange(npup_x) - float(npup_x) / 2.0 - offset_x + 0.5) * dx
        ys = (np.arange(npup_y) - float(npup_y) / 2.0 - offset_y + 0.5) * dy

        us = (np.arange(npix_x) - float(npix_x) / 2.0 - offset_x + 0.5) * du
        vs = (np.arange(npix_y) - float(npix_y) / 2.0 - offset_y + 0.5) * dv
    elif centering == SYMMETRIC:
        xs = (np.arange(npup_x) - float(npup_x) / 2.0 + 0.5) * dx
        ys = (np.arange(npup_y) - float(npup_y) / 2.0 + 0.5) * dy

        us = (np.arange(npix_x) - float(npix_x) / 2.0 + 0.5) * du
        vs = (np.arange(npix_y) - float(npix_y) / 2.0 + 0.5) * dv
    else:
        raise ValueError("Invalid centering style")

    xu = np.outer(xs, us)
    yv = np.outer(ys, vs)

    if inverse:
        exp_yv = np.exp(-2.0 * np.pi * -1j * yv).T
        exp_xu = np.exp(-2.0 * np.pi * -1j * xu)
        t1 = np.dot(exp_yv, plane)
        t2 = np.dot(t1, exp_xu)
    else:
        exp_xu = np.exp(-2.0 * np.pi * 1j * xu)
        exp_yv = np.exp(-2.0 * np.pi * 1j * yv).T
        t1 = np.dot(exp_yv, plane)
        t2 = np.dot(t1, exp_xu)

    norm_coeff = np.sqrt((nlam_dy * nlam_dx) / (npup_y * npup_x * npix_y * npix_x))

    return norm_coeff * t2


def matrix_idft(*args, **kwargs):  # noqa: D103
    kwargs["inverse"] = True
    return matrix_dft(*args, **kwargs)


matrix_idft.__doc__ = matrix_dft.__doc__.replace(  # type: ignore[union-attr]
    "Perform a matrix discrete Fourier transform",
    "Perform an inverse matrix discrete Fourier transform",
)

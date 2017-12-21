"""
A module that provides utilities for working with WCS.

:Authors: Mihai Cara (contact: help@stsci.edu)

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

# STDLIB
import logging

# THIRD PARTY
import numpy as np

from . import linearfit

# LOCAL
from . import __version__
from . import __vdate__


__all__ = ['apply_affine_to_wcs']


__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


def linearize(wcsim, wcsima, wcsref, imcrpix, f, shift, hx=1.0, hy=1.0):
    # NOTE: we increase the accuracy of intermediate results in order
    #       to avoid possible loss of accuracy when computing
    #       small differences.
    # linearization using 5-point formula for first order derivative
    x0 = imcrpix[0]
    y0 = imcrpix[1]
    p = np.asarray([[x0, y0],
                    [x0 - hx, y0],
                    [x0 - hx * 0.5, y0],
                    [x0 + hx * 0.5, y0],
                    [x0 + hx, y0],
                    [x0, y0 - hy],
                    [x0, y0 - hy * 0.5],
                    [x0, y0 + hy * 0.5],
                    [x0, y0 + hy]],
                   dtype=np.float)

    # convert image coordinates to reference image coordinates:
    p = wcsref.wcs_world2pix(*wcsim.wcs_pix2world(p[:, 0], p[:, 1]))
    p = np.asanyarray(p, dtype=ndfloat128).T
    # apply linear fit transformation:
    p = (np.dot(f, (p - shift).T).T).astype(np.float)
    # convert back to image coordinate system:
    p = wcsima.wcs_world2pix(*wcsref.wcs_pix2world(p[:, 0], p[:, 1]))
    p = np.asanyarray(p, dtype=ndfloat128).T

    # derivative with regard to x:
    u1 = ((p[1] - p[4]) + 8 * (p[3] - p[2])) / (6 * hx)
    # derivative with regard to y:
    u2 = ((p[5] - p[8]) + 8 * (p[7] - p[6])) / (6 * hy)

    return (np.array([u1, u2], dtype=np.float).T,
            np.asanyarray(p[0], dtype=np.float))



def apply_affine_to_wcs(imwcs, rot=0.0, scale=1.0, xsh=0.0, ysh=0.0,
                        matrix=None):
    """
    Apply affine transformation given in the reference image to WCS of an
    image/catalog by adjusting this WCS's parameters.

    Parameters
    ----------
    imwcs : gwcs.WCS
        WCS of the image to be adjusted. On return, this WCS will contain
        modified parameters.

    refwcs : gwcs.WCS
        WCS of the reference coordinate system in which affine transformations
        were defined.

    imshape : tuple, optional
        A tuple indicating the shape of the image's data array. Must follow
        the same convention as the shape of the `~numpy.ndarray`.

    rot : tuple, float, optional
        Rotation angle in degrees. Two values (one for each axis) can be
        provided as a tuple. `rot` is ignored if `matrix` is provided.

    scale : tuple, float, optional
        Scale of the liniar transformation. Two values (one for each axis)
        can be provided as a tuple. `scale` is ignored if `matrix` is provided.

    xsh : float, optional
        Shift (pixels) along the X-axis. Shifts are applied before shifts and
        rotations.

    ysh : float, optional
        Shift (pixels) along the Y-axis. Shifts are applied before shifts and
        rotations.

    matrix : numpy.ndarray, None, optional
        A 2x2 matrix providing the coefficients of the liniar transformation.

    """
    # compute the matrix for the scale and rotation correction
    if matrix is None:
        matrix = linearfit.build_fit_matrix(rot, scale)

    matrix = matrix.T
    shift = -np.dot(np.linalg.inv(matrix), [xsh, ysh])
    imwcs.set_correction(matrix, shift)

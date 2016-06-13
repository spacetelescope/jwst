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
from .simplewcs import SimpleWCS


__all__ = ['create_ref_wcs', 'apply_affine_to_wcs']

__version__ = '0.1.0'
__vdate__ = '17-April-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


def create_ref_wcs(images):
    """
    Create a reference WCS from a list of input images.

    Parameters
    ----------
    images : list of WCSImageCatalog or WCSGroupCatalog
        A list of (groups of) image catalogs to be used to derive an
        undistorted reference WCS with similar pixel scales and orientation
        as the input images.

    Returns
    -------
    wcs : gwcs.WCS
        A `gwcs.WCS` object.

    shape : tuple
        A tuple showing the size of the "reference image" would have had if it
        were large enough to include all the input images.

    """
    # TODO: Later we can implement a more sophisticated and comprehensive
    # algorithm, hopefully one used across multiple products (such as drizzle).
    # For now simply rectify the WCS of the first image in the first group.

    if hasattr(images[0], '__iter__'):
        img = images[0][0]
    else:
        img = images[0]

    swcs = SimpleWCS(img.wcs, copy=True)

    # find bounding rectangle for all input images.
    # TODO: in the future implement an algorithm to find
    #       minimum-area enclosing rectangle and its orientation.
    pscale = None
    bb_ra = []
    bb_dec = []
    for img in images:
        if hasattr(img, '__iter__'):
            for i in img:
                ra, dec = i.bb_radec
                bb_ra += ra.tolist()
                bb_dec += dec.tolist()
                if pscale is None or pscale > i.pscale:
                    pscale = i.pscale
        else:
            ra, dec = img.bb_radec
            bb_ra += ra.tolist()
            bb_dec += dec.tolist()
            if pscale is None or pscale > img.pscale:
                pscale = img.pscale

    # CD matrix and parity:
    cd = np.dot(np.diag([swcs.cdelt1, swcs.cdelt2]), swcs.pc)
    det = np.linalg.det(cd)
    par = -1 if det < 0.0 else 1

    # find Y-axis orientation:
    rot = np.arctan2(cd[0, 1], cd[1, 1])  # angle of the Y-axis

    # create a perfectly square, orthogonal WCS
    sn = np.sin(rot)
    cs = np.cos(rot)
    orthogonal_pc = np.array([[par * cs, sn], [-par * sn, cs]])

    swcs.begin_update()
    swcs.cdelt = [pscale / 3600.0, pscale / 3600.0]
    swcs.pc = orthogonal_pc
    swcs.end_update()

    # find coordinates of the lower-bottom corner of the bounding box:
    bbx, bby = swcs.wcs_world2pix(bb_ra, bb_dec)
    shx = np.amin(bbx)
    shy = np.amin(bby)
    mx = np.amax(bbx)
    my = np.amax(bby)
    shape = (int(np.ceil(my)) - int(np.floor(shy)),
             int(np.ceil(mx)) - int(np.floor(shx)))

    # shift CRPIX so that blc is at (0, 0):
    swcs.begin_update()
    swcs.crpix1 -= shx
    swcs.crpix2 -= shy
    swcs.end_update()

    wcs = swcs.get_std_wcs()
    wcs.name = 'ReferenceWCS'

    return wcs, shape


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


def _invmat(mat):
    if mat.shape != (2, 2):
        return np.linalg.inv(mat)

    # use own code for 2x2 matrices
    inv = mat.astype(ndfloat128)

    det = inv[0, 0] * inv[1, 1] - inv[0, 1] * inv[1, 0]
    if np.abs(det) < np.finfo(np.float64).tiny:
        raise ArithmeticError('Singular matrix.')

    d = inv[1, 1]
    inv[1, 0] *= -1.0
    inv[0, 1] *= -1.0
    inv[0, 0] = d
    inv[1, 1] = d
    inv /= det
    inv = inv.astype(np.float64)
    if not np.all(np.isfinite(inv)):
        raise ArithmeticError('Singular matrix.')

    return inv


def _cd_to_pc(cd):
    cdelt = np.sqrt((np.asanyarray(cd)**2).sum(axis=0, dtype=np.float))
    pc = np.dot(np.linalg.inv(np.diag(cdelt)), cd)
    return pc, cdelt


def apply_affine_to_wcs(imwcs, refwcs, imshape=None, rot=0.0, scale=1.0,
                        xsh=0.0, ysh=0.0, matrix=None):
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

    refwcs = SimpleWCS(refwcs)
    imwcs = SimpleWCS(imwcs, copy=False)
    cwcs = imwcs.copy()

    shift = np.asanyarray([xsh, ysh]) - np.dot(refwcs.crpix, matrix) + \
        refwcs.crpix
    matrix = _invmat(matrix).T

    # estimate step for numerical differentiation. We need a step
    # large enough to avoid rounding errors and small enough to get a
    # better precision for numerical differentiation.
    # TODO: The logic below should be revised at a later time so that it
    # better takes into account the two competing requirements.
    if imshape is None:
        hx = 2.0
        hy = 2.0

    else:
        hx = max(1.0, min(2.0, np.fabs(imwcs.crpix1 - 1.0) / 100.0,
                          np.fabs(imshape[1] - imwcs.crpix1) / 100.0))
        hy = max(1.0, min(2.0, np.fabs(imwcs.crpix2 - 1.0) / 100.0,
                          np.fabs(imshape[0] - imwcs.crpix2) / 100.0))

    # compute new CRVAL for the image WCS:
    crpixinref = np.asanyarray(refwcs.wcs_world2pix(
        *imwcs.wcs_pix2world(imwcs.crpix1, imwcs.crpix2)))
    crpixinref = np.dot(matrix, (crpixinref - shift).T).T
    imwcs.crval = refwcs.wcs_pix2world(*crpixinref)

    # compute new CD matrix of the image WCS:
    (U, u) = linearize(cwcs, imwcs, refwcs, imwcs.crpix,
                       matrix, shift, hx=hx, hy=hy)

    cd = np.dot(imwcs.cd.astype(ndfloat128), U).astype(np.float64)
    pc, cdelt = _cd_to_pc(cd)
    imwcs.begin_update()
    imwcs.pc = pc
    imwcs.cdelt = cdelt
    imwcs.end_update()

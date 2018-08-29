"""
A module that provides algorithms for performing linear fits between
sets of 2D points.

:Authors: Mihai Cara, Warren Hack (contact: help@stsci.edu)


"""
import logging
import numpy as np


__all__ = ['iter_linear_fit', 'build_fit_matrix']

__author__ = 'Mihai Cara, Warren Hack'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


def iter_linear_fit(xy, uv, xyindx=None, uvindx=None, xyorig=None, uvorig=None,
                    fitgeom='general', nclip=3, sigma=3.0, center=None):
    """
    Compute iteratively using sigma-clipping linear transformation parameters
    that fit `xy` sources to `uv` sources.

    """

    minobj_per_fitgeom = {'shift': 1, 'rscale': 2, 'general': 3}
    minobj = minobj_per_fitgeom[fitgeom]

    xy = np.asanyarray(xy, dtype=np.float64)
    uv = np.asanyarray(uv, dtype=np.float64)

    if xy.shape[0] < nclip:
        log.warning("The number of sources for the fit < number of clipping "
                    "iterations.")
        log.warning("Resetting number of clipping iterations to 0.")
        nclip = 0

    if center is None:
        xcen = uv[:, 0].mean(dtype=np.float64)
        ycen = uv[:, 1].mean(dtype=np.float64)
        center = [xcen, ycen]

    xy -= center
    uv -= center

    fit = linear_fit(xy, uv, fitgeom=fitgeom, verbose=True)

    npts = xy.shape[0]
    npts0 = 0

    if nclip is None:
        nclip = 0

    # define index to initially include all points
    for n in range(nclip):
        resids = fit['resids']

        # redefine what pixels will be included in next iteration
        whtfrac = npts / (npts - npts0 - 1.0)
        cutx = sigma * (fit['rms'][0] * whtfrac)
        cuty = sigma * (fit['rms'][1] * whtfrac)

        goodpix = (np.abs(resids[:, 0]) < cutx) & (np.abs(resids[:, 1]) < cuty)
        ngoodpix = np.count_nonzero(goodpix)

        if ngoodpix < minobj:
            break

        npts0 = npts - goodpix.shape[0]
        xy = xy[goodpix]
        uv = uv[goodpix]

        if xyindx is not None:
            xyindx = xyindx[goodpix]
        if uvindx is not None:
            uvindx = uvindx[goodpix]

        if xyorig is not None:
            xyorig = xyorig[goodpix]
        if uvorig is not None:
            uvorig = uvorig[goodpix]

        fit = linear_fit(xy, uv, fitgeom=fitgeom, verbose=False)

    fit['img_coords'] = xy
    fit['ref_coords'] = uv
    fit['img_indx'] = xyindx
    fit['ref_indx'] = uvindx
    fit['img_orig_xy'] = xyorig
    fit['ref_orig_xy'] = uvorig

    return fit


def linear_fit(xy, uv, fitgeom='rscale', verbose=False):
    """ Performs an 'rscale' fit between matched lists of pixel positions
        xy and uv
    """
    fitgeom = fitgeom.lower()

    xy = np.asanyarray(xy)
    uv = np.asanyarray(uv)

    if verbose:
        log.info("Performing '{:s}' fit".format(fitgeom))

    if fitgeom == 'general':
        result = fit_general(xy, uv)
    elif fitgeom == 'rscale':
        result = fit_rscale(xy, uv)
    elif fitgeom == 'shift':
        result = fit_shifts(xy, uv)
    else:
        raise ValueError("Unsupported 'fitgeom' value: '{}'".format(fitgeom))

    return result


def fit_shifts(xy, uv):
    """ Performs a simple fit for the shift only between
        matched lists of positions 'xy' and 'uv'.

        =================================
        DEVELOPMENT NOTE:
            Checks need to be put in place to verify that
            enough objects are available for a fit.
        =================================

        Output:
           (Xo,Yo),Rot,(Scale,Sx,Sy)
           where
                Xo,Yo:  offset,
                Rot:    rotation,
                Scale:  average scale change, and
                Sx,Sy:  scale changes in X and Y separately.

        Algorithm and nomenclature provided by: Colin Cox (11 Nov 2004)

    """
    if len(xy) < 1:
        raise ValueError("At least one point is required to find shifts.")

    diff_pts = xy - uv
    meanx = (diff_pts[:, 0].mean(dtype=np.float64)).astype(np.float64)
    meany = (diff_pts[:, 1].mean(dtype=np.float64)).astype(np.float64)
    Pcoeffs = np.array([1.0, 0.0, meanx])
    Qcoeffs = np.array([0.0, 1.0, meany])

    fit = build_fit(Pcoeffs, Qcoeffs, 'shift')
    resids = diff_pts - fit['offset']
    rms = [resids[:, 0].std(dtype=np.float64),
           resids[:, 1].std(dtype=np.float64)]
    fit['resids'] = resids
    fit['rms'] = rms

    return fit


# Implementation of geomap 'rscale' fitting based on 'lib/geofit.x'
# by Warren Hack. Support for axis flips added by Mihai Cara.
def fit_rscale(xyin, xyref):
    """
    Set up the products used for computing the fit derived using the code from
    lib/geofit.x for the function 'geo_fmagnify()'. Comparisons with results
    from geomap (no additional clipping) were made and produced the same
    results out to 5 decimal places.

    Output
    ------
    fit: dict
        Dictionary containing full solution for fit.
    """
    if len(xyin) < 2:
        raise ValueError("At least two points are required to find "
                         "shifts, rotation, and scale.")

    dx = xyref[:, 0].astype(ndfloat128)
    dy = xyref[:, 1].astype(ndfloat128)
    du = xyin[:, 0].astype(ndfloat128)
    dv = xyin[:, 1].astype(ndfloat128)

    n = xyref.shape[0]
    Sx = dx.sum()
    Sy = dy.sum()
    Su = du.sum()
    Sv = dv.sum()
    xr0 = Sx / n
    yr0 = Sy / n
    xi0 = Su / n
    yi0 = Sv / n
    Sxrxr = np.power((dx - xr0), 2).sum()
    Syryr = np.power((dy - yr0), 2).sum()
    Syrxi = ((dy - yr0) * (du - xi0)).sum()
    Sxryi = ((dx - xr0) * (dv - yi0)).sum()
    Sxrxi = ((dx - xr0) * (du - xi0)).sum()
    Syryi = ((dy - yr0) * (dv - yi0)).sum()

    rot_num = Sxrxi * Syryi
    rot_denom = Syrxi * Sxryi

    if rot_num == rot_denom:
        det = 0.0
    else:
        det = rot_num - rot_denom

    if (det < 0):
        rot_num = Syrxi + Sxryi
        rot_denom = Sxrxi - Syryi
    else:
        rot_num = Syrxi - Sxryi
        rot_denom = Sxrxi + Syryi

    if rot_num == rot_denom:
        theta = 0.0
    else:
        theta = np.rad2deg(np.arctan2(rot_num, rot_denom))
        if theta < 0:
            theta += 360.0

    ctheta = np.cos(np.deg2rad(theta))
    stheta = np.sin(np.deg2rad(theta))
    s_num = rot_denom * ctheta + rot_num * stheta
    s_denom = Sxrxr + Syryr

    if s_denom < 0.0:
        mag = 1.0
    elif s_denom > 0.0:
        mag = s_num / s_denom
    else:
        raise ArithmeticError("Singular matrix.")

    if det < 0:
        # "flip" y-axis (reflection about x-axis *after* rotation)
        # NOTE: keep in mind that 'fit_matrix'
        #       is the transposed rotation matrix.
        sthetax = -mag * stheta
        cthetay = -mag * ctheta
    else:
        sthetax = mag * stheta
        cthetay = mag * ctheta

    cthetax = mag * ctheta
    sthetay = mag * stheta

    sdet = np.sign(det)
    xshift = (xi0 - (xr0 * cthetax + sdet * yr0 * sthetax)).astype(np.float64)
    yshift = (yi0 - (-sdet * xr0 * sthetay + yr0 * cthetay)).astype(np.float64)

    P = np.array([cthetax, sthetay, xshift], dtype=np.float64)
    Q = np.array([-sthetax, cthetay, yshift], dtype=np.float64)

    # Return the shift, rotation, and scale changes
    result = build_fit(P, Q, fitgeom='rscale')
    resids = xyin - np.dot((xyref), result['fit_matrix']) - result['offset']
    rms = [resids[:, 0].std(dtype=np.float64),
           resids[:, 1].std(dtype=np.float64)]
    result['rms'] = rms
    result['resids'] = resids

    return result


def _inv3x3(x):
    """ Return inverse of a 3-by-3 matrix.
    """
    x0 = x[:, 0]
    x1 = x[:, 1]
    x2 = x[:, 2]
    m = np.array([np.cross(x1, x2), np.cross(x2, x0), np.cross(x0, x1)])
    d = np.dot(x0, np.cross(x1, x2))
    if np.abs(d) < np.finfo(np.float64).tiny:
        raise ArithmeticError("Singular matrix.")
    return (m / d)


def fit_general(xy, uv):
    """ Performs a simple fit for the shift only between
        matched lists of positions 'xy' and 'uv'.

        =================================
        DEVELOPMENT NOTE:
            Checks need to be put in place to verify that
            enough objects are available for a fit.
        =================================

        Output:
           (Xo,Yo),Rot,(Scale,Sx,Sy)
           where
                Xo,Yo:  offset,
                Rot:    rotation,
                Scale:  average scale change, and
                Sx,Sy:  scale changes in X and Y separately.

        Algorithm and nomenclature provided by: Colin Cox (11 Nov 2004)

    """
    if len(xy) < 3:
        raise ValueError("At least three points are required to find "
                         "6-parameter linear affine transformations.")

    # Set up products used for computing the fit
    gxy = xy.astype(ndfloat128)
    guv = uv.astype(ndfloat128)

    Sx = gxy[:, 0].sum()
    Sy = gxy[:, 1].sum()
    Su = guv[:, 0].sum()
    Sv = guv[:, 1].sum()

    Sxu = np.dot(gxy[:, 0], guv[:, 0])
    Syu = np.dot(gxy[:, 1], guv[:, 0])
    Sxv = np.dot(gxy[:, 0], guv[:, 1])
    Syv = np.dot(gxy[:, 1], guv[:, 1])
    Suu = np.dot(guv[:, 0], guv[:, 0])
    Svv = np.dot(guv[:, 1], guv[:, 1])
    Suv = np.dot(guv[:, 0], guv[:, 1])

    n = len(xy[:, 0])
    M = np.array([[Su, Sv, n], [Suu, Suv, Su], [Suv, Svv, Sv]])
    U = np.array([Sx, Sxu, Sxv])
    V = np.array([Sy, Syu, Syv])

    # The fit solutioN...
    # where
    #   u = P0 + P1*x + P2*y
    #   v = Q0 + Q1*x + Q2*y
    #
    invM = _inv3x3(M)
    P = np.dot(invM, U).astype(np.float64)
    Q = np.dot(invM, V).astype(np.float64)
    if not (np.all(np.isfinite(P)) and np.all(np.isfinite(Q))):
        raise ArithmeticError("Singular matrix.")

    # Return the shift, rotation, and scale changes
    result = build_fit(P, Q, 'general')
    resids = xy - np.dot(uv, result['fit_matrix']) - result['offset']
    rms = [resids[:, 0].std(dtype=np.float64),
           resids[:, 1].std(dtype=np.float64)]
    result['rms'] = rms
    result['resids'] = resids

    return result


def build_fit(P, Q, fitgeom):
    # Build fit matrix:
    fit_matrix = np.dstack((P[:2], Q[:2]))[0]

    # determinant of the transformation
    det = P[0] * Q[1] - P[1] * Q[0]
    sdet = np.sign(det)
    proper = (sdet >= 0)

    # Create a working copy (no reflections) for computing transformation
    # parameters (scale, rotation angle, skew):
    wfit = fit_matrix.copy()

    # Default skew:
    skew = 0.0

    if fitgeom == 'shift':
        return {'offset': (P[2], Q[2]),
                'fit_matrix': fit_matrix,
                'rot': 0.0,
                'rotxy': (0.0, 0.0, 0.0, skew),
                'scale': (1.0, 1.0, 1.0),
                'coeffs': (P, Q),
                'skew': skew,
                'proper': proper,
                'fitgeom': fitgeom}

    # Compute average scale:
    s = np.sqrt(np.abs(det))
    # Compute scales for each axis:
    if fitgeom == 'general':
        sx = np.sqrt(P[0]**2 + Q[0]**2)
        sy = np.sqrt(P[1]**2 + Q[1]**2)
    else:
        sx = s
        sy = s

    # Remove scale from the transformation matrix:
    wfit[0, :] /= sx
    wfit[1, :] /= sy

    # Compute rotation angle as if we have a proper rotation.
    # This will also act as *some sort* of "average rotation" even for
    # transformations with different rot_x and rot_y:
    prop_rot = np.rad2deg(np.arctan2(wfit[1, 0] - sdet * wfit[0, 1],
                                     wfit[0, 0] + sdet * wfit[1, 1])) % 360.0

    if proper and fitgeom == 'rscale':
        rotx = prop_rot
        roty = prop_rot
        rot = prop_rot
        skew = 0.0
    else:
        rotx = np.rad2deg(np.arctan2(-wfit[0, 1], wfit[0, 0])) % 360.0
        roty = np.rad2deg(np.arctan2(wfit[1, 0], wfit[1, 1])) % 360.0
        rot = 0.5 * (rotx + roty)
        skew = roty - rotx

    return {'offset': (P[2], Q[2]),
            'fit_matrix': fit_matrix,
            'rot': prop_rot,
            'rotxy': (rotx, roty, rot, skew),
            'scale': (s, sx, sy),
            'coeffs': (P, Q),
            'skew': skew,
            'proper': proper,
            'fitgeom': fitgeom}


def build_fit_matrix(rot, scale=1):
    """
    Create an affine transformation matrix (2x2) from the provided rotation
    and scale transformations.

    Parameters
    ----------
    rot : tuple, float, optional
        Rotation angle in degrees. Two values (one for each axis) can be
        provided as a tuple.

    scale : tuple, float, optional
        Scale of the liniar transformation. Two values (one for each axis)
        can be provided as a tuple.

    Returns
    -------
    matrix : numpy.ndarray
       A 2x2 `numpy.ndarray` containing coefficients of a liniear
       transformation.

    """
    if hasattr(rot, '__iter__'):
        rx = rot[0]
        ry = rot[1]
    else:
        rx = float(rot)
        ry = rx

    if hasattr(scale, '__iter__'):
        sx = scale[0]
        sy = scale[1]
    else:
        sx = float(scale)
        sy = sx

    matrix = np.array(
        [
            [sx * np.cos(np.deg2rad(rx)), -sx * np.sin(np.deg2rad(rx))],
            [sy * np.sin(np.deg2rad(ry)), sy * np.cos(np.deg2rad(ry))]
        ]
    )

    return matrix

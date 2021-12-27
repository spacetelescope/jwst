""" A set of routines to assist in the WCS transforms used in the
cube_build step
"""
import numpy as np
import math
# _______________________________________________________________________


def radec2std(crval1, crval2, ra, dec, rot_angle=None):
    """ Compute the tangent projection coordinates (xi,eta) from ra,dec
    using crval1 and crval2 (the tangent point).

    Parameters
    ____________
    crval1 : float
      RA value of tangent point
    crval2 : float
      DEC value of tangent point
    rot_angle: float
      rotation angle given in degrees
    ra : numpy.ndarray or float
      A list (or single value) of ra points to convert
    dec : numpy.ndarray  or float
      A list (or single value) of ra points to convert

    Return Values
    _____________
    xi, eta - rectangular coordinates of tangent plane projected ra,dec

    """

    if np.isscalar(ra):
        ra = np.asarray([ra])
        dec = np.asarray([dec])

    rad2arcsec = (180.0 * 3600.0) / math.pi
    deg2rad = math.pi / 180.0

    ra0 = crval1 * deg2rad
    dec0 = crval2 * deg2rad
    radiff = ra * deg2rad - ra0
    decr = dec * deg2rad

    h = np.sin(decr) * math.sin(dec0) + np.cos(decr) * math.cos(dec0) * np.cos(radiff)

    xi = np.cos(decr) * np.sin(radiff) / h
    eta = (np.sin(decr) * math.cos(dec0) -
           np.cos(decr) * math.sin(dec0) * np.cos(radiff)) / h

    xi = -xi
    xi = xi * rad2arcsec
    eta = eta * rad2arcsec

    # xi is made negative so it increases in the opposite direction
    # of ra to match the images the Parity of the ifu_cube.

    if rot_angle is not None:
        temp1 = xi * np.cos(-rot_angle * deg2rad) - eta * np.sin(-rot_angle * deg2rad)
        temp2 = xi * np.sin(-rot_angle * deg2rad) + eta * np.cos(-rot_angle * deg2rad)
        xi = temp1
        eta = temp2

    return xi, eta
# ________________________________________________________________________________


def std2radec(crval1, crval2, xi, eta):
    """ Compute ra,dec from the tangent plane rectangular coordinates

    Compute the ra,dec values of  tangent plane rectangular coordinates using
    crval1, crval2(the tangent point). This routine takes the rectangular
    plane and projects it to the spherical plane using crval1, crval2 as
    the tangent plane.

    Parameters
    ____________
    crval1 : float
      RA value of tangent point
    crval2 : float
      DEC value of tangent point
    xi : float
      xi rectangular coordinates of tangent plane projected ra,dec
    eta  : float
      eta rectangular coordinates of tangent plane projected ra,dec

    Return Values
    _____________
    ra : float
      list (or single value) of ra computed values
    dec : float
      list (or single value) of dec computed values
    """

    if np.isscalar(xi):
        eta = np.asarray([eta])
        xi = np.asarray([xi])

    rad2arcsec = (180.0 * 3600.0) / math.pi
    deg2rad = math.pi / 180.0

    ra0 = crval1 * deg2rad
    dec0 = crval2 * deg2rad

    # tangent projection
    xi = xi / rad2arcsec
    eta = eta / rad2arcsec
    xi = -xi  # xi is made negative so it increases in the opposite direction
    # of ra to match the images the Parity of the ifu_cube is
    # for ra is PC1_1 = -1

    ra0 = crval1 * deg2rad
    dec0 = crval2 * deg2rad

    beta = math.cos(dec0) - eta * math.sin(dec0)
    angle = xi / beta

    ra = np.arctan(angle) + ra0
    gamma = np.sqrt(xi * xi + beta * beta)
    angle = eta * math.cos(dec0) + math.sin(dec0)
    angle = angle / gamma
    dec = np.arctan(angle)
    ra = ra / deg2rad
    dec = dec / deg2rad

    mask = ra < 0
    ra[mask] += 360.0
    mask = ra > 360.0
    ra[mask] -= 360.0
    return ra, dec

# _______________________________________________________________________


def v2v32radec_estimate(ra_ref, dec_ref, roll_ref, v2_ref, v3_ref, v2, v3):
    """ Estimation of ra and dec from the v2, v3 coordinates

    This routine is used for debugging purposes. It is not actually used
    in the cube_build step for routine IFU cube building.
    The conversion from V2,V3 to ra,dec is handled more accurately by
    the transforms provided by assign_wcs.

    Parameters
    ----------
    ra_ref : float
       ra of reference point given in arc seconds
    dec_ref : float
       dec of reference point given in arc seconds
    roll_ref : float
       roll angle given in degrees
    v2_ref : float
       v2 coordinate of reference point given in arc seconds
    v3_ref : float
       v3 coordinate of reference point given in arc seconds
    v2 : float
       v2 coordinate given in arc seconds
    v3 :  float
       v3 coordinate given in arc seconds

    Returns
    -------
    ra : float
    dec : float

    Notes
    ----
    it is assumed that the v2,v3 coordinates have the effects of dithering included
    """
    d2r = math.pi / 180.0

    v2deg = v2.copy() / 3600.0   # convert to degrees
    v3deg = v3.copy() / 3600.0   # convert to degrees

    v2_ref = v2_ref / 3600.0  # covert to degrees
    v3_ref = v3_ref / 3600.0  # convert to degrees
    v3_ref_rad = v3_ref * d2r
    roll_ref_rad = roll_ref * d2r

    delta_v2 = (v2deg - v2_ref) * math.cos(v3_ref_rad)
    delta_v3 = (v3deg - v3_ref)
    delta_ra = delta_v2 * math.cos(roll_ref_rad) + delta_v3 * math.sin(roll_ref_rad)
    delta_dec = -delta_v2 * math.sin(roll_ref_rad) + delta_v3 * math.cos(roll_ref_rad)

    ra = ra_ref + delta_ra / math.cos(dec_ref * d2r)
    dec = delta_dec + dec_ref

    return ra, dec
# _______________________________________________________________________


def radec2v2v3_estimate(ra_ref, dec_ref, roll_ref, v2_ref, v3_ref, ra, dec):
    """ Convert ra,dec to v2, v3

    This routine is used for debugging purposes. It is not actually used
    in the cube_build step for routine IFU cube building.
    The conversion from Ra,Dec to V2,V3 is handled more accurately by
    the transforms provided by assign_wcs.

    Parameters
    ----------
    ra_ref : float
       ra of reference point given in degrees
    dec_ref : float
       dec of reference point given in degrees
    roll_ref : float
       roll angle given in degrees
    v2_ref : float
       v2 coordinate of reference point given in arc seconds
    v3_ref : float
       v3 coordinate of reference point given in arc seconds
    ra : float
       ra coordinate given in degrees
    dec :  float
       dec coordinate given in degrees

    Returns
    -------
    v2 : float
       v2 coordinate in arc seconds
    v3 : float
       v2 coordinate in arc seconds
     """

    d2r = math.pi / 180.0
    r2d = 180.0 / math.pi

    v3_ref_rad = (v3_ref) / 3600.0 * d2r  # convert from arc seconds to radians
    v2_ref_rad = (v2_ref) / 3600.0 * d2r  # convert from arc seconds to radians

    roll_ref_rad = roll_ref * d2r

    ra_ref_rad = ra_ref * d2r
    dec_ref_rad = dec_ref * d2r
    this_ra = ra * d2r
    this_dec = dec * d2r

    delta_ra = (this_ra - ra_ref_rad) * math.cos(dec_ref_rad)
    delta_dec = (this_dec - dec_ref_rad)

    dv2 = delta_ra * math.cos(roll_ref_rad) - delta_dec * math.sin(roll_ref_rad)
    dv3 = delta_ra * math.sin(roll_ref_rad) + delta_dec * math.cos(roll_ref_rad)

    v2 = v2_ref_rad + (dv2 / math.cos(v3_ref_rad))
    v3 = v3_ref_rad + dv3

    v2 = v2 * r2d * 3600.0  # convert from radians to arc seconds
    v3 = v3 * r2d * 3600.0  # convert from radians to arc seconds

    return v2, v3

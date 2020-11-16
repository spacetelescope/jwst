'''
This script adds velocity aberration correction information to the FITS
files provided to it on the command line (one or more).

It assumes the following keywords are present in the file header:

JWST_DX (km/sec)
JWST_DY (km/sec)
JWST_DZ (km/sec)
RA_REF (deg)
DEC_REF (deg)

The keywords added are:

VA_SCALE (dimensionless scale factor)

It does not currently place the new keywords in any particular location
in the header other than what is required by the standard.
'''

import astropy.io.fits as fits
import logging
import math


# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

SPEED_OF_LIGHT = 299792.458     # km / s
d_to_r = math.pi / 180.


def aberration_scale(velocity_x, velocity_y, velocity_z,
                     targ_ra, targ_dec):
    """Compute the scale factor due to velocity aberration.

    Parameters
    ----------
    velocity_x, velocity_y, velocity_z: float
        The components of the velocity of JWST, in km / s with respect to
        the Sun.  These are celestial coordinates, with x toward the
        vernal equinox, y toward right ascension 90 degrees and declination
        0, z toward the north celestial pole.

    targ_ra, targ_dec: float
        The right ascension and declination of the target (or some other
        point, such as the center of a detector).  The equator and equinox
        should be the same as the coordinate system for the velocity.

    Returns
    -------
    scale_factor: float
        Multiply the nominal image scale (e.g. in degrees per pixel) by
        this value to obtain the image scale corrected for the "aberration
        of starlight" due to the velocity of JWST with respect to the Sun.
    """

    speed = math.sqrt(velocity_x**2 + velocity_y**2 + velocity_z**2)
    if speed == 0.0:
        logger.warning('Speed is zero. Forcing scale to 1.0')
        return 1.0

    beta = speed / SPEED_OF_LIGHT
    gamma = 1. / math.sqrt(1. - beta**2)

    # [targ_x, targ_y, targ_z] is a unit vector.
    r_xy = math.cos(targ_dec * d_to_r)          # radial distance in xy-plane
    targ_x = r_xy * math.cos(targ_ra * d_to_r)
    targ_y = r_xy * math.sin(targ_ra * d_to_r)
    targ_z = math.sin(targ_dec * d_to_r)

    dot_prod = (velocity_x * targ_x +
                velocity_y * targ_y +
                velocity_z * targ_z)
    cos_theta = dot_prod / speed
    # This sin_theta is only valid over the range [0, pi], but so is the
    # angle between the velocity vector and the direction toward the target.
    sin_theta = math.sqrt(1. - cos_theta**2)

    tan_theta_p = sin_theta / (gamma * (cos_theta + beta))
    theta_p = math.atan(tan_theta_p)

    scale_factor = (gamma * (cos_theta + beta)**2 /
                    (math.cos(theta_p)**2 * (1. + beta * cos_theta)))

    return scale_factor


def aberration_offset(velocity_x, velocity_y, velocity_z,
                      targ_ra, targ_dec):
    """Compute the RA/Dec offsets due to velocity aberration.

    Parameters
    ----------
    velocity_x, velocity_y, velocity_z: float
        The components of the velocity of JWST, in km / s with respect to
        the Sun.  These are celestial coordinates, with x toward the
        vernal equinox, y toward right ascension 90 degrees and declination
        0, z toward the north celestial pole.

    targ_ra, targ_dec: float
        The right ascension and declination of the target (or some other
        point, such as the center of a detector).  The equator and equinox
        should be the same as the coordinate system for the velocity.

    Returns
    -------
    delta_ra, delta_dec: float
        The offset to be added to the input RA/Dec, in units of radians.
    """

    xdot = velocity_x / SPEED_OF_LIGHT
    ydot = velocity_y / SPEED_OF_LIGHT
    zdot = velocity_z / SPEED_OF_LIGHT

    sin_alpha = math.sin(targ_ra * d_to_r)
    cos_alpha = math.cos(targ_ra * d_to_r)
    sin_delta = math.sin(targ_dec * d_to_r)
    cos_delta = math.cos(targ_dec * d_to_r)

    delta_ra = (-xdot * sin_alpha + ydot * cos_alpha) / cos_delta
    delta_dec = (-xdot * cos_alpha * sin_delta -
                 ydot * sin_alpha * sin_delta +
                 zdot * cos_delta)

    return delta_ra, delta_dec


def add_dva(filename):
    '''
    Given the name of a valid partially populated level 1b JWST file,
    determine the velocity aberration scale factor.

    It presumes all the accessed keywords are present (see first block).
    '''
    hdulist = fits.open(filename, 'update')
    pheader = hdulist[0].header
    sheader = hdulist['SCI'].header
    jwst_dx = float(pheader['JWST_DX'])
    jwst_dy = float(pheader['JWST_DY'])
    jwst_dz = float(pheader['JWST_DZ'])
    ra_ref = float(sheader['RA_REF'])
    dec_ref = float(sheader['DEC_REF'])

    # compute the velocity aberration information
    scale_factor = aberration_scale(jwst_dx, jwst_dy, jwst_dz,
                                    ra_ref, dec_ref)
    ra_off, dec_off = aberration_offset(jwst_dx, jwst_dy, jwst_dz,
                                        ra_ref, dec_ref)

    # update header
    pheader['DVA_RA'] = ra_off
    pheader['DVA_DEC'] = dec_off
    sheader['VA_SCALE'] = scale_factor
    hdulist.flush()
    hdulist.close()

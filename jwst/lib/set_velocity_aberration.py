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

import logging
import numpy as np
import astropy.io.fits as fits
from gwcs.geometry import SphericalToCartesian
from scipy.constants import speed_of_light
import math

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

SPEED_OF_LIGHT = speed_of_light / 1000     # km / s
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
    beta = np.array([velocity_x, velocity_y, velocity_z]) / SPEED_OF_LIGHT
    beta2 = np.dot(beta, beta)
    u = SphericalToCartesian()(targ_ra, targ_dec)
    beta_u =  np.dot(beta, u)
    igamma = np.sqrt(1.0 - beta2)  # inverse of usual gamma
    scale_factor = (1.0 + beta_u) /  igamma
    return scale_factor


def aberration_offset(velocity_x, velocity_y, velocity_z, targ_ra, targ_dec):
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
    # Algorithm from Colin Cox notebook
    beta = np.array([velocity_x, velocity_y, velocity_z]) / SPEED_OF_LIGHT
    u = np.array(SphericalToCartesian()(targ_ra, targ_dec))
    beta2 = np.dot(beta, beta)
    if beta2 == 0.0:
        return 0.0, 0.0
    igamma = np.sqrt(1.0 - beta2)  # inverse of usual gamma
    beta_u = np.dot(beta, u)
    u_corr = (igamma * u + beta * (1.0 + (1.0 - igamma) * beta_u / beta2)) / (1.0 + beta_u)
    ra_corr = np.arctan2(u_corr[1], u_corr[0])
    dec_corr = np.arcsin(u_corr[2])
    return ra_corr - np.deg2rad(targ_ra), dec_corr -  np.deg2rad(targ_dec)


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

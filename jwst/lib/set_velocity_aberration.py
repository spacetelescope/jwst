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
from gwcs.geometry import SphericalToCartesian, CartesianToSpherical
from scipy.constants import speed_of_light

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

SPEED_OF_LIGHT = speed_of_light / 1000  # km / s


def compute_va_effects(velocity_x, velocity_y, velocity_z, targ_ra, targ_dec):
    """ Computes constant scale factor due to velocity aberration as well as
    corrected ``RA`` and ``DEC`` values.

    Parameters
    ----------
    velocity_x, velocity_y, velocity_z: float
        The components of the velocity of JWST, in km / s with respect to
        the Sun.  These are celestial coordinates, with x toward the
        vernal equinox, y toward right ascension 90 degrees and declination
        0, z toward the north celestial pole.

    targ_ra, targ_dec: float
        The right ascension and declination of the target (or some other
        point, such as the center of a detector) in the barycentric coordinate
        system.  The equator and equinox should be the same as the coordinate
        system for the velocity.

    Returns
    -------
    scale_factor: float
        Multiply the nominal image scale (e.g. in degrees per pixel) by
        this value to obtain the image scale corrected for the "aberration
        of starlight" due to the velocity of JWST with respect to the Sun.

    apparent_ra: float
        Aparent star position in the moving telescope frame.

    apparent_dec: float
        Aparent star position in the moving telescope frame.

    """
    beta = np.array([velocity_x, velocity_y, velocity_z]) / SPEED_OF_LIGHT
    beta2 = np.dot(beta, beta)
    if beta2 == 0.0:
        return 1.0, targ_ra, targ_dec

    u = np.asanyarray(SphericalToCartesian()(targ_ra, targ_dec))
    beta_u =  np.dot(beta, u)
    igamma = np.sqrt(1.0 - beta2)  # inverse of usual gamma
    scale_factor = (1.0 + beta_u) / igamma

    # Algorithm below is from Colin Cox notebook.
    # Also see: Instrument Science Report OSG-CAL-97-06 by Colin Cox (1997).
    u_corr = (igamma * u + beta * (1.0 + (1.0 - igamma) * beta_u / beta2)) / (1.0 + beta_u)

    apparent_ra, apparent_dec = CartesianToSpherical()(*u_corr)
    return scale_factor, apparent_ra, apparent_dec


def add_dva(filename):
    """
    Given the name of a valid partially populated level 1b JWST file,
    determine the velocity aberration scale factor and apparent target position
    in the moving (telescope) frame.

    It presumes all the accessed keywords are present (see first block).
    """
    hdulist = fits.open(filename, 'update')
    pheader = hdulist[0].header
    sheader = hdulist['SCI'].header
    jwst_dx = float(pheader['JWST_DX'])
    jwst_dy = float(pheader['JWST_DY'])
    jwst_dz = float(pheader['JWST_DZ'])
    ra_ref = float(sheader['RA_REF'])
    dec_ref = float(sheader['DEC_REF'])

    # compute the velocity aberration information
    scale_factor, apparent_ra, apparent_dec = compute_va_effects(
        jwst_dx, jwst_dy, jwst_dz, ra_ref, dec_ref
    )

    # update header
    pheader['VA_RA'] = apparent_ra
    pheader['VA_DEC'] = apparent_dec
    sheader['VA_SCALE'] = scale_factor
    hdulist.flush()
    hdulist.close()

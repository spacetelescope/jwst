"""
Utilities for velocity aberration correction.

Script to add velocity aberration correction information to the FITS
files provided to it on the command line (one or more).

It assumes the following keywords are present in the file header:

* JWST_DX (km/sec)
* JWST_DY (km/sec)
* JWST_DZ (km/sec)
* RA_REF (deg)
* DEC_REF (deg)

The keywords added are:

* VA_SCALE (dimensionless scale factor)

It does not currently place the new keywords in any particular location
in the header other than what is required by the standard.
"""

import logging

import stcal.velocity_aberration as va

import jwst.datamodels as dm
from jwst.datamodels import Level1bModel  # type: ignore[attr-defined]

# Configure logging
logger = logging.getLogger(__name__)

__all__ = ["add_dva"]


def add_dva(filename, force_level1bmodel=True):
    """
    Determine velocity aberration.

    Given the name of a valid partially populated level 1b JWST file,
    determine the velocity aberration scale factor and apparent target position
    in the moving (telescope) frame.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    filename : str
        The name of the file to be updated.
    force_level1bmodel : bool, optional
        If True, the input file will be force-opened as a Level1bModel.  If False,
        the file will be opened using the generic DataModel.  The default is True.
    """
    if force_level1bmodel:
        model = Level1bModel(filename)
    else:
        model = dm.open(filename)
    scale_factor, apparent_ra, apparent_dec = va.compute_va_effects(
        velocity_x=model.meta.ephemeris.velocity_x_bary,
        velocity_y=model.meta.ephemeris.velocity_y_bary,
        velocity_z=model.meta.ephemeris.velocity_z_bary,
        ra=model.meta.wcsinfo.ra_ref,
        dec=model.meta.wcsinfo.dec_ref,
    )

    # update header
    model.meta.velocity_aberration.scale_factor = scale_factor
    model.meta.velocity_aberration.va_ra_ref = apparent_ra
    model.meta.velocity_aberration.va_dec_ref = apparent_dec
    model.save(filename)

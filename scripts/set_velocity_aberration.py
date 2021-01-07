#!/usr/bin/env python

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
import sys
import logging

from jwst.lib.set_velocity_aberration import add_dva

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        raise ValueError('missing filename argument(s)')
    for filename in sys.argv[1:]:
        add_dva(filename)

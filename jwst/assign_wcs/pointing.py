from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import numpy as np
from ..datamodels import fits_support
from gwcs import utils as gwutils


##TODO: populate the core model with basiic WCS keywords?


def create_fitswcs_transform(input_model):
    ff = fits_support.to_fits(input_model._instance, input_model._schema)
    hdu = fits_support.get_hdu(ff._hdulist, "PRIMARY")
    header = hdu.header
    transform = gwutils.make_fitswcs_transform(header)
    return transform


def fitswcs_transform_from_model(input_model):
    """
    Create a fits wcs from a datamodel.meta.wcsinfo.
    """
    wcsinfo = {}
    wcsaxes = input_model._instance['meta']['wcsinfo']['wcsaxes']
    wcsinfo['WCSAXES'] = wcsaxes
    for key in ['CTYPE', 'CRPIX', 'CRVAL', 'CDELT', 'CUNIT']:
        val = []
        for ax in range(1, wcsaxes + 1):
            try:
                val.append(input_model._instance['meta']['wcsinfo'][(key + "{0}".format(ax)).lower()])
            except KeyError:
                pass
        wcsinfo[key.upper()] = val

    pc = np.zeros((wcsaxes, 2))

    for i in range(1, wcsaxes + 1):
        for j in range(1, 3):
            pc[i-1, j-1] = input_model._instance['meta']['wcsinfo']['pc{0}_{1}'.format(i, j)]
    wcsinfo['PC'] = pc
    # getting "coordinates" from the schema is currently broken.
    wcsinfo['RADESYS'] = 'ICRS' #input_model._instance['meta']['coordinates']['reference_frame']
    wcsinfo['has_cd'] = False
    transform = gwutils.make_fitswcs_transform(wcsinfo)
    return transform

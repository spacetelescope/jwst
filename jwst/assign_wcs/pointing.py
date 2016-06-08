from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from jwst.datamodels import fits_support
from gwcs import utils as gwutils


##TODO: populate the core model with basiic WCS keywords?


def create_fitswcs_transform(input_model):
    ff = fits_support.to_fits(input_model._instance, input_model._schema)
    hdu = fits_support.get_hdu(ff._hdulist, "PRIMARY")
    header = hdu.header
    transform = gwutils.make_fitswcs_transform(header)
    return transform



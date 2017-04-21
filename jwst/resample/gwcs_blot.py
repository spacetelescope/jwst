from __future__ import (division, print_function, unicode_literals,
    absolute_import)

import os
import os.path

import numpy as np
from astropy import wcs
from astropy.io import fits

from astropy import wcs as fitswcs
from gwcs import wcs
from gwcs import wcstools

from drizzle import util
from drizzle import doblot
from drizzle import cdrizzle
from . import resample_utils

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class GWCSBlot(object):
    """
    Combine images using the drizzle algorithm
    """
    def __init__(self, product):
        """
        Create new blotted output objects and set the blot parameters.

        All parameters are optional, but either infile or outwcs must be supplied.
        If infile initializes the object from a file written after a
        previous run of drizzle. Results from the previous run will be combined
        with new results. The value passed in outwcs will be ignored. If infile is
        not set, outwcs will be used to initilize a new run of drizzle.

        Parameters
        ----------

        product : datamodel
            A data model containing results from a previous run. The three
            extensions SCI, WHT, and CTX contain the combined image, total counts
            and image id bitmap, repectively. The WCS of the combined image is
            also read from the SCI extension.

        """

        # Initialize the object fields
        self.source_model = product
        self.source_wcs = product.meta.wcs
        self.source = product.data

    def extract_image(self, blot_img, interp='poly5', sinscl=1.0):
        """
        Resample the output/resampled image to recreate an input image based on
        the input image's world coordinate system

        Parameters
        ----------

        blot_img : datamodel
            Datamodel containing header and WCS to define the 'blotted' image

        interp : str, optional
            The type of interpolation used in the resampling. The
            possible values are "nearest" (nearest neighbor interpolation),
            "linear" (bilinear interpolation), "poly3" (cubic polynomial
            interpolation), "poly5" (quintic polynomial interpolation),
            "sinc" (sinc interpolation), "lan3" (3rd order Lanczos
            interpolation), and "lan5" (5th order Lanczos interpolation).

        sincscl : float, optional
            The scaling factor for sinc interpolation.
        """
        blot_wcs = blot_img.meta.wcs
        blot_shape = blot_img.shape
        outsci = np.zeros((blot_shape[0], blot_shape[1]), dtype=np.float32)

        # Compute the mapping between the input and output pixel coordinates
        #log.info("Creating PIXMAP for blotted image...")
        pixmap = resample_utils.calc_gwcs_pixmap(self.source_wcs, blot_wcs)

        source_pscale = self.source_model.meta.wcsinfo.cdelt1
        blot_pscale = blot_img.meta.wcsinfo.cdelt1

        pix_ratio = source_pscale / blot_pscale
        log.info("Starting blot...")
        cdrizzle.tblot(self.source, pixmap, outsci, scale=pix_ratio, kscale=1.0,
                       interp=interp, exptime=1.0, misval=0.0, sinscl=sinscl)

        return outsci

#! /usr/bin/env python
"""
JWST pipeline step for image alignment.

:Authors: Mihai Cara

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

import os
import logging
import numpy as np
from astropy.table import Table

from ..stpipe import Step, cmdline
from .. import datamodels

from . imalign import align
from . wcsimage import *

__all__ = ['TweakRegStep']


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class TweakRegStep(Step):
    """
    TweakRegStep: Image alignment based on catalogs of sources detected in
    input images.
    """

    spec = """
        # Optimize alignment order:
        enforce_user_order = boolean(default=True) # Align images in user specified order?

        # Reference Catalog parameters:
        expand_refcat = boolean(default=False) # Expand reference catalog with new sources?

        # Object matching parameters:
        minobj = integer(default=15) # Minimum number of objects acceptable for matching
        searchrad = float(default=1.0) # The search radius for a match
        searchunits = option('arcseconds', 'pixels', default='arcseconds') # Units for search radius
        use2dhist = boolean(default=True) # Use 2d histogram to find initial offset?
        separation = float(default=0.5) # Minimum object separation in pixels
        tolerance = float(default=1.0) # Matching tolerance for xyxymatch in pixels
        xoffset = float(default=0.0), # Initial guess for X offset in pixels
        yoffset = float(default=0.0) # Initial guess for Y offset in pixels

        # Catalog fitting parameters:
        fitgeometry = option('shift', 'rscale', 'general', default='general') # Fitting geometry
        nclip = integer(min=0, default=3) # Number of clipping iterations in fit
        sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units
    """

    reference_file_types = []

    def process(self, input):

        img = datamodels.ModelContainer(input, persist=True)

        if len(img) == 0:
            raise ValueError("Input must contain at least one image model.")

        # group images by their "group id":
        grp_img = img.models_grouped

        if len(grp_img) == 1:
            # we need at least two images/groups to perform image alignment
            log.info("At least two images or groups are required for image "
                     "alinment.")
            log.info("Nothing to do. Exiting 'TweakRegStep'...")
            return input

        # find maximum length of the filename:
        max_name_len = None
        for im in img:
            fnlen = _get_filename(im)
            if max_name_len is None or max_name_len < fnlen:
                max_name_len = fnlen

        # create a list of WCS-Catalog-Images Info and/or their Groups:
        images = []

        for g in grp_img:
            if len(g) == 0:
                raise AssertionError("Logical error in the pipeline code.")

            wcsimlist = list(map(self._imodel2wcsim, g))
            wgroup = WCSGroupCatalog(wcsimlist, name=wcsimlist[0].name)
            images.append(wgroup)

        # align images:
        align(
            imcat=images,
            refcat=None,
            enforce_user_order=self.enforce_user_order,
            expand_refcat=self.expand_refcat,
            minobj=self.minobj,
            searchrad=self.searchrad,
            searchunits=self.searchunits,
            use2dhist=self.use2dhist,
            separation=self.separation,
            tolerance=self.tolerance,
            xoffset=self.xoffset,
            yoffset=self.yoffset,
            fitgeom=self.fitgeometry,
            nclip=self.nclip,
            sigma=self.sigma
        )

        return img

    def _imodel2wcsim(self, image_model):
        # make sure that we have a catalog:
        if hasattr(image_model, 'catalog'):
            catalog = image_model.catalog
        else:
            catalog = image_model.meta.tweakreg_catalog.filename

        if not isinstance(catalog, Table):
            catalog = Table.read(catalog, format='ascii.ecsv')

        if 'xcentroid' in catalog.colnames:
            catalog.rename_column('xcentroid', 'x')
            catalog.rename_column('ycentroid', 'y')

        # create WCSImageCatalog object:
        im = WCSImageCatalog(
            shape=image_model.data.shape,
            wcs=image_model.meta.wcs,
            catalog=catalog,
            name=_get_filename(image_model),
            meta={'image_model': image_model}
        )

        return im


def _get_filename(image_model):
    if image_model.meta.filename is None:
        return None
    else:
        return os.path.splitext(os.path.basename(image_model.meta.filename))[0]


if __name__ == '__main__':
    cmdline.step_script(TweakRegStep)

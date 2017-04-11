#! /usr/bin/env python
"""
JWST pipeline step for sky matching.

:Authors: Mihai Cara

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

import collections
import numpy as np
import logging

from ..stpipe import Step, cmdline
from .. import datamodels

from stsci.tools import bitmask

#LOCAL:
from . skymatch import match
from . skyimage import SkyImage, SkyGroup
from . skystatistics import SkyStats

__all__ = ['SkyMatchStep']


class SkyMatchStep(Step):
    """
    SkyMatchStep: Subtraction or equalization of sky background in science images.
    """

    spec = """
        # General sky matching parameters:
        skymethod = option('local', 'global', 'match', 'global+match', default='global+match') # sky computation method
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?

        # Image's bounding polygon parameters:
        stepsize = integer(default=None) # Max vertex separation

        # Sky statistics parameters:
        skystat = option('median', 'midpt', 'mean', 'mode', default='mode') # sky statistics
        dqbits = string(default=None) # "good" DQ bits
        lower = float(default=None) # Lower limit of "good" pixel values
        upper = float(default=None) # Upper limit of "good" pixel values
        nclip = integer(min=0, default=5) # number of sky clipping iterations
        lsigma = float(min=0.0, default=4.0) # Lower clipping limit, in sigma
        usigma = float(min=0.0, default=4.0) # Upper clipping limit, in sigma
        binwidth = float(min=0.0, default=0.1) # Bin width for 'mode' and 'midpt' `skystat`, in sigma
    """

    reference_file_types = []

    def process(self, input):
        self.log.setLevel(logging.DEBUG)
        img = datamodels.ModelContainer(input)

        self._dqbits = bitmask.interpret_bit_flags(self.dqbits)

        # set sky stattistics:
        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth
        )

        # group images by their "group id":
        grp_img = img.models_grouped

        # create a list of "Sky" Images and/or Groups:
        images = []
        grp_id = 1

        for g in grp_img:
            if len(g) > 1:
                images.append(
                    SkyGroup(
                        list(map(self._imodel2skyim, g)),
                        id=grp_id
                    )
                )
                grp_id += 1
            elif len(g) == 1:
                images.append(self._imodel2skyim(g[0]))
            else:
                raise AssertionError("Logical error in the pipeline code.")

        # match/compute sky values:
        match(images, skymethod=self.skymethod, match_down=self.match_down,
              subtract=self.subtract)

        # set sky background value in each image's meta:
        for im in images:
            if isinstance(im, SkyImage):
                self._set_sky_background(im.meta['imagemodel'], im.sky)
            else:
                for gim in im:
                    self._set_sky_background(gim.meta['imagemodel'], gim.sky)

        return img

    #DEBUG: remove _group_images_by_id() once we are past
    # debugging stage.
    def _group_images_by_id(self, image_list):

        # group images by their "group id":
        imdict = collections.OrderedDict()

        for im in image_list:
            groupid = im.meta.observation.observation_label

            if groupid in imdict:
                imdict[groupid].append(im)
            else:
                imdict[groupid] = [im]

        return imdict.values()

    def _imodel2skyim(self, image_model):

        # create
        if self._dqbits is None:
            dqmask = None
        else:
            dqmask = bitmask.bitmask2mask(
                bitmask=image_model.dq,
                ignore_bits=self._dqbits,
                good_mask_value=1,
                dtype=np.uint8
            )

        sky_im = SkyImage(
            image=image_model.data,
            wcs_fwd=image_model.meta.wcs.__call__,
            wcs_inv=image_model.meta.wcs.invert,
            pix_area=1.0, #TODO: pixel area
            convf=1.0,    #TODO: conv. factor to brightness
            mask=dqmask,
            id=image_model.meta.filename, # file name?
            skystat=self._skystat,
            stepsize=self.stepsize,
            meta={'imagemodel': image_model}
        )

        if self.subtract:
            if hasattr(image_model.meta, 'bkglevel'):
                sky_im.sky = image_model.meta.bkglevel

        return sky_im

    def _set_sky_background(self, image, sky):
        #skybg_schema = {
            #"meta.bkglevel":
            #{
                #"title": "Computed sky background value",
                #"type": "float",
                #"fits_keyword": "BKGLEVEL"
            #}
        #}
        #image.extend_schema(skybg_schema)
        image.meta.bkglevel = sky
        image.meta.bkgsub = self.subtract

if __name__ == '__main__':
    cmdline.step_script(SkyMatchStep)

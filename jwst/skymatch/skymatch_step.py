#! /usr/bin/env python
"""
JWST pipeline step for sky matching.

:Authors: Mihai Cara


"""

from copy import deepcopy
import logging

import numpy as np

from ..stpipe import Step
from .. import datamodels
from ..datamodels.dqflags import pixel
from ..datamodels.util import open as datamodel_open


from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)

# LOCAL:
from .skymatch import match
from .skyimage import SkyImage, SkyGroup
from .skystatistics import SkyStats


__all__ = ['SkyMatchStep']


class SkyMatchStep(Step):
    """
    SkyMatchStep: Subtraction or equalization of sky background in science images.
    """

    class_alias = "skymatch"

    spec = """
        # General sky matching parameters:
        skymethod = option('local', 'global', 'match', 'global+match', default='match') # sky computation method
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?

        # Image's bounding polygon parameters:
        stepsize = integer(default=None) # Max vertex separation

        # Sky statistics parameters:
        skystat = option('median', 'midpt', 'mean', 'mode', default='mode') # sky statistics
        dqbits = string(default='~DO_NOT_USE+NON_SCIENCE') # "good" DQ bits
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
        # for now turn off memory optimization until we have better machinery
        # to handle outputs in a consistent way.
        self._is_asn = False
        # self._is_asn = datamodels.util.is_association(input) or isinstance(input, str)

        img = datamodels.ModelContainer(
            input,
            save_open=not self._is_asn,
            return_open=not self._is_asn
        )

        self._dqbits = interpret_bit_flags(self.dqbits, flag_name_map=pixel)

        # set sky statistics:
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
                self._set_sky_background(
                    im,
                    "COMPLETE" if im.is_sky_valid else "SKIPPED"
                )
            else:
                for gim in im:
                    self._set_sky_background(
                        gim,
                        "COMPLETE" if gim.is_sky_valid else "SKIPPED"
                    )

        return input if self._is_asn else img

    def _imodel2skyim(self, image_model):
        input_image_model = image_model
        if self._is_asn:
            image_model = datamodel_open(image_model)

        if self._dqbits is None:
            dqmask = np.isfinite(image_model.data).astype(dtype=np.uint8)
        else:
            dqmask = bitfield_to_boolean_mask(
                image_model.dq,
                self._dqbits,
                good_mask_value=1,
                dtype=np.uint8
            ) * np.isfinite(image_model.data)

        # see if 'skymatch' was previously run and raise an exception
        # if 'subtract' mode has changed compared to the previous pass:
        if image_model.meta.background.subtracted is None:
            if image_model.meta.background.level is not None:
                if self._is_asn:
                    image_model.close()

                # report inconsistency:
                raise ValueError("Background level was set but the "
                                 "'subtracted' property is undefined (None).")
            level = 0.0

        else:
            level = image_model.meta.background.level
            if level is None:
                # NOTE: In principle we could assume that level is 0 and
                # possibly add a log entry documenting this, however,
                # at this moment I think it is saver to quit and...
                #
                # report inconsistency:
                if self._is_asn:
                    image_model.close()

                raise ValueError("Background level was subtracted but the "
                                 "'level' property is undefined (None).")

            if image_model.meta.background.subtracted != self.subtract:
                # cannot run 'skymatch' step on already "skymatched" images
                # when 'subtract' spec is inconsistent with
                # meta.background.subtracted:
                if self._is_asn:
                    image_model.close()

                raise ValueError("'subtract' step's specification is "
                                 "inconsistent with background info already "
                                 "present in image '{:s}' meta."
                                 .format(image_model.meta.filename))

        wcs = deepcopy(image_model.meta.wcs)

        sky_im = SkyImage(
            image=image_model.data,
            wcs_fwd=wcs.__call__,
            wcs_inv=wcs.invert,
            pix_area=1.0,  # TODO: pixel area
            convf=1.0,  # TODO: conv. factor to brightness
            mask=dqmask,
            id=image_model.meta.filename,  # file name?
            skystat=self._skystat,
            stepsize=self.stepsize,
            reduce_memory_usage=self._is_asn,
            meta={'image_model': input_image_model}
        )

        if self._is_asn:
            image_model.close()

        if self.subtract:
            sky_im.sky = level

        return sky_im

    def _set_sky_background(self, sky_image, step_status):
        image = sky_image.meta['image_model']
        sky = sky_image.sky

        if self._is_asn:
            dm = datamodel_open(image)
        else:
            dm = image

        if step_status == "COMPLETE":
            dm.meta.background.method = str(self.skymethod)
            dm.meta.background.level = sky
            dm.meta.background.subtracted = self.subtract
            if self.subtract:
                dm.data[...] = sky_image.image[...]

        dm.meta.cal_step.skymatch = step_status

        if self._is_asn:
            dm.save(image)
            dm.close()

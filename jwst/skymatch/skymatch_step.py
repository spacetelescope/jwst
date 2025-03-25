#! /usr/bin/env python
"""
JWST pipeline step for sky matching.

:Authors: Mihai Cara
"""

from copy import deepcopy
import logging

import numpy as np

from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)

from stcal.skymatch import skymatch, SkyImage, SkyGroup, SkyStats

from stdatamodels.jwst.datamodels.dqflags import pixel

from jwst.datamodels import ModelLibrary
from jwst.lib.suffix import remove_suffix
from pathlib import Path

from jwst.stpipe import Step

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["SkyMatchStep"]


class SkyMatchStep(Step):
    """SkyMatchStep: Subtraction or equalization of sky background in science images."""

    class_alias = "skymatch"

    spec = """
        # General sky matching parameters:
        skymethod = option('local', 'global', 'match', 'global+match', 'user', default='match') # sky computation method
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?
        skylist = string(default=None) # Filename pointing to list of (imagename skyval) pairs

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

        # Memory management:
        in_memory = boolean(default=True) # If False, preserve memory using temporary files
    """  # noqa: E501

    reference_file_types: list = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def process(self, input_models):
        """
        Run the step.

        Parameters
        ----------
        input_models : Any data type readable into a ModelLibrary, e.g. an asn file
            An association of datamodels to input.

        Returns
        -------
        ModelLibrary
            A library of datamodels with the skymatch step applied.
        """
        self.log.setLevel(logging.DEBUG)

        if isinstance(input_models, ModelLibrary):
            library = input_models
        else:
            library = ModelLibrary(input_models, on_disk=not self.in_memory)

        # Method: "user". Use user-provided sky values, and bypass skymatch() altogether.
        if self.skymethod == "user":
            return self._user_sky(library)

        self._dqbits = interpret_bit_flags(self.dqbits, flag_name_map=pixel)

        # set sky statistics:
        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth,
        )

        images = []
        with library:
            for group_index, (_group_id, group_inds) in enumerate(library.group_indices.items()):
                sky_images = []
                for index in group_inds:
                    model = library.borrow(index)
                    try:
                        sky_images.append(self._imodel2skyim(model, index))
                    finally:
                        library.shelve(model, index, modify=False)
                if len(sky_images) == 1:
                    images.extend(sky_images)
                else:
                    images.append(SkyGroup(sky_images, sky_id=group_index))

        # match/compute sky values:
        skymatch(
            images, skymethod=self.skymethod, match_down=self.match_down, subtract=self.subtract
        )

        # set sky background value in each image's meta:
        with library:
            for im in images:
                if isinstance(im, SkyImage):
                    self._set_sky_background(
                        im, library, "COMPLETE" if im.is_sky_valid else "SKIPPED"
                    )
                else:
                    for gim in im:
                        self._set_sky_background(
                            gim, library, "COMPLETE" if gim.is_sky_valid else "SKIPPED"
                        )

        return library

    def _imodel2skyim(self, image_model, index):
        if self._dqbits is None:
            dqmask = np.isfinite(image_model.data).astype(dtype=np.uint8)
        else:
            dqmask = bitfield_to_boolean_mask(
                image_model.dq, self._dqbits, good_mask_value=1, dtype=np.uint8
            ) * np.isfinite(image_model.data)

        # see if 'skymatch' was previously run and raise an exception
        # if 'subtract' mode has changed compared to the previous pass:
        if image_model.meta.background.subtracted is None:
            if image_model.meta.background.level is not None:
                # report inconsistency:
                raise ValueError(
                    "Background level was set but the 'subtracted' property is undefined (None)."
                )
            level = 0.0

        else:
            level = image_model.meta.background.level
            if level is None:
                # NOTE: In principle we could assume that level is 0 and
                # possibly add a log entry documenting this, however,
                # at this moment I think it is saver to quit and...
                #
                # report inconsistency:
                raise ValueError(
                    "Background level was subtracted but the 'level' property is undefined (None)."
                )

            if image_model.meta.background.subtracted != self.subtract:
                # cannot run 'skymatch' step on already "skymatched" images
                # when 'subtract' spec is inconsistent with
                # meta.background.subtracted:
                raise ValueError(
                    "'subtract' step's specification is "
                    "inconsistent with background info already "
                    f"present in image '{image_model.meta.filename:s}' meta."
                )

        wcs = deepcopy(image_model.meta.wcs)

        sky_im = SkyImage(
            image=image_model.data,
            wcs_fwd=wcs.__call__,
            wcs_inv=wcs.invert,
            pix_area=1.0,  # TODO: pixel area
            convf=1.0,  # TODO: conv. factor to brightness
            mask=dqmask,
            sky_id=image_model.meta.filename,
            skystat=self._skystat,
            stepsize=self.stepsize,
            reduce_memory_usage=False,  # this overwrote input files
            meta={"index": index},
        )

        if self.subtract:
            sky_im.sky = level

        return sky_im

    def _set_sky_background(self, sky_image, library, step_status):
        """
        Set sky background values in the image's metadata.

        Parameters
        ----------
        sky_image : SkyImage
            SkyImage object containing sky image data and metadata.

        library : ModelLibrary
            Library of input data models, must be open

        step_status : str
            Status of the sky subtraction step. Must be one of the following:
            'COMPLETE', 'SKIPPED'.
        """
        index = sky_image.meta["index"]
        dm = library.borrow(index)
        sky = sky_image.sky

        if step_status == "COMPLETE":
            dm.meta.background.method = str(self.skymethod)
            dm.meta.background.level = sky
            dm.meta.background.subtracted = self.subtract
            if self.subtract:
                dm.data[...] = sky_image.image[...]

        dm.meta.cal_step.skymatch = step_status
        library.shelve(dm, index)

    def _user_sky(self, library):
        """
        Handle user-provided sky values for each image.

        Parameters
        ----------
        library : ModelLibrary
            Library of input data models.

        Returns
        -------
        ModelLibrary
            Library of input data models with sky background values set to user-provided values.
        """
        if self.skylist is None:
            raise ValueError('skymethod set to "user", but no sky value file provided.')

        log.info(" ")
        log.info(
            "Setting sky background of input images to user-provided values "
            f"from `skylist` ({self.skylist})."
        )

        # read the comma separated file and get just the stem of the filename
        skylist = np.genfromtxt(
            self.skylist,
            dtype=[("fname", "<S128"), ("sky", "f")],
        )
        skyfnames, skyvals = skylist["fname"], skylist["sky"]
        skyfnames = skyfnames.astype(str)
        skyfnames = [remove_suffix(Path(fname).stem)[0] for fname in skyfnames]
        skyfnames = np.array(skyfnames)

        if len(skyvals) != len(library):
            raise ValueError(
                f"Number of entries in skylist ({len(self.skylist)}) does not match "
                f"number of input images ({len(library)})."
            )

        with library:
            for model in library:
                fname, _ = remove_suffix(Path(model.meta.filename).stem)
                sky = skyvals[np.where(skyfnames == fname)]
                if len(sky) == 0:
                    raise ValueError(f"Image with stem '{fname}' not found in the skylist.")
                if len(sky) > 1:
                    raise ValueError(
                        f"Image with stem '{fname}' found multiple times in the skylist."
                    )

                sky = float(sky[0])
                log.debug(f"Setting sky background of image '{model.meta.filename}' to {sky}.")

                model.meta.background.level = sky
                model.meta.background.subtracted = self.subtract
                model.meta.background.method = self.skymethod
                if self.subtract:
                    model.data -= sky
                model.meta.cal_step.skymatch = "COMPLETE"
                library.shelve(model)

        return library

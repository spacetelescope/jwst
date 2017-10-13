#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import background_sub


class BackgroundStep(Step):
    """
    BackgroundStep:  Subtract background exposures from target exposures.
    """

    spec = """
    """

    # These reference files are only used for WFSS/GRISM data.
    # xxx can't be used yet:
    # xxx reference_file_types = ["wfssbkg", "wavelengthrange"]

    def process(self, input, bkg_list):

        """
        Subtract the background signal from target exposures by subtracting
        designated background images from them.

        Parameters
        ----------
        input: JWST data model
            input target data model to which background subtraction is applied

        bkg_list: filename list
            list of background exposure file names

        Returns
        -------
        result: JWST data model
            the background-subtracted target data model
        """

        # Load the input data model
        with datamodels.open(input) as input_model:

            if input_model.meta.exposure.type in ["NIS_WFSS", "NRC_GRISM"]:
                # Delete these four lines when reference files are available.
                self.log.info("No 'wfssbkg' reference file available; "
                              "skipping background subtraction.")
                result = input_model.copy()
                result.meta.cal_step.back_sub = 'SKIPPED'

                # Uncomment this block when reference files are available.
                # bkg_filename = self.get_reference_file(
                #                 input_model, "wfssbkg")
                # wl_range_name = self.get_reference_file(
                #                 input_model, "wavelengthrange")
                # Do the background subtraction for WFSS/GRISM data
                # result = background_sub.subtract_wfss_bkg(
                #                 input_model, bkg_filename, wl_range_name)
                # result.meta.cal_step.back_sub = 'COMPLETE'
            else:
                # Do the background subtraction
                result = background_sub.background_sub(input_model, bkg_list)
                result.meta.cal_step.back_sub = 'COMPLETE'

        return result

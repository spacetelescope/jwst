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

    # This reference file is only used for WFSS/GRISM data.
    # xxx can't be used yet xxx reference_file_types = ["wfssbkg"]

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
                # xxx self.bkg_filename = self.get_reference_file(
                # xxx                 input_model, "wfssbkg")
                # xxx bkg_list = [self.bkg_filename]
                bkg_list = []

            # Do the background subtraction
            result = background_sub.background_sub(input_model, bkg_list)
            if result.meta.cal_step.back_sub != 'SKIPPED':
                result.meta.cal_step.back_sub = 'COMPLETE'

        return result

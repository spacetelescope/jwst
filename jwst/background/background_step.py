#! /usr/bin/env python

from jwst.stpipe import Step
from jwst import datamodels
from . import background_sub

class BackgroundStep(Step):
    """
    BackgroundStep:  Subtract background exposures from target exposures.
    """

    spec = """
    """

    def process(self, input, bkg_list):

        """
        Short Summary
        -------------
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
        with models.open(input) as input_model:

            # Do the background subtraction
            result = background_sub.background_sub(input_model, bkg_list)
            result.meta.cal_step.back_sub = 'COMPLETE'

        return result


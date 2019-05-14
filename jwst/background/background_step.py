#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import background_sub
import numpy as np
__all__ = ["BackgroundStep"]


class BackgroundStep(Step):
    """
    BackgroundStep:  Subtract background exposures from target exposures.
    """

    spec = """
        sigma = float(default=3.0)  # Clipping threshold
        maxiters = integer(default=None)  # Number of clipping iterations
    """

    # These reference files are only used for WFSS/GRISM data.
    reference_file_types = ["wfssbkg", "wavelengthrange"]

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

            if input_model.meta.exposure.type in ["NIS_WFSS", "NRC_WFSS"]:

                # Get the reference file names
                bkg_name = self.get_reference_file(input_model, "wfssbkg")
                wlrange_name = self.get_reference_file(input_model,
                                                       "wavelengthrange")
                self.log.info('Using WFSSBKG reference file %s', bkg_name)
                self.log.info('Using WavelengthRange reference file %s',
                              wlrange_name)

                # Do the background subtraction for WFSS/GRISM data
                result = background_sub.subtract_wfss_bkg(
                                input_model, bkg_name, wlrange_name)
                result.meta.cal_step.back_sub = 'COMPLETE'
            else:
                # check if input data is NRS_IFU
                tolerance = 1.0e-8
                do_sub = True
                if input_model.meta.instrument.name in ["NIRSPEC"]:
                    # check if GWA_XTIL & GWA_YTIL values of source
                    # background are the same. If not skip step
                    input_xtilt = input_model.meta.instrument.gwa_xtilt
                    input_ytilt = input_model.meta.instrument.gwa_ytilt
                    for bkg_file in bkg_list:
                        with datamodels.open(bkg_file) as bkg_model:
                            bkg_xtilt = bkg_model.meta.instrument.gwa_xtilt
                            bkg_ytilt = bkg_model.meta.instrument.gwa_ytilt
                            if np.allclose((input_xtilt, input_ytilt),
                                           (bkg_xtilt, bkg_ytilt), atol=tolerance, rtol=0):
                                pass
                            else:
                                do_sub  = False
                                break
                # Do the background subtraction
                if do_sub:
                    result = background_sub.background_sub(input_model,
                                                           bkg_list,
                                                           self.sigma,
                                                           self.maxiters)
                    result.meta.cal_step.back_sub = 'COMPLETE'
                else:
                    result = input_model.copy()
                    result.meta.cal_step.back_sub = 'SKIPPED'
                    self.log.warning('Skipping background subtraction')
                    self.log.warning('GWA_XTIL and GWA_YTIL source values '
                                     'are not the same as bkg values')

        return result

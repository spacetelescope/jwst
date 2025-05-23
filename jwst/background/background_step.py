#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from .background_sub import background_sub
from .background_sub_wfss import subtract_wfss_bkg
from .asn_intake import asn_get_data

import numpy as np


__all__ = ["BackgroundStep"]


WFSS_TYPES = ["NIS_WFSS", "NRC_GRISM", "NRC_WFSS"]


class BackgroundStep(Step):
    """Subtract background exposures from target exposures."""

    class_alias = "background"

    spec = """
        bkg_list = force_list(default=None)  # List of background files
        save_combined_background = boolean(default=False)  # Save combined background image
        sigma = float(default=3.0)  # Clipping threshold
        maxiters = integer(default=None)  # Number of clipping iterations
        wfss_mmag_extract = float(default=None)  # WFSS minimum abmag to extract
        wfss_maxiter = integer(default=5)  # WFSS iterative outlier rejection max iterations
        wfss_rms_stop = float(default=0)  # WFSS iterative outlier rejection RMS improvement threshold (percent)
        wfss_outlier_percent = float(default=1)  # WFSS outlier percentile to reject per iteration
    """  # noqa: E501

    # These reference files are only used for WFSS/GRISM data.
    reference_file_types = ["wfssbkg", "wavelengthrange"]

    # Define a suffix for optional saved output of the combined background
    bkg_suffix = "combinedbackground"

    def process(self, step_input, input_bkg_list=None):
        """
        Subtract designated background images from target exposures.

        Parameters
        ----------
        step_input : str, ImageModel or IFUImageModel
            Input target data model to which background subtraction is applied or asn file

        input_bkg_list : list, optional
            File name list of background exposures. If an association file is provided as
            input, this parameter can be used to provide the name of the file that has
            been processed through the assign_wcs step.

        Returns
        -------
        result : ImageModel or IFUImageModel
            The background-subtracted target data model
        """
        # If the input is an asn get all the info, otherwise assume input is
        # a fits file or datamodel
        if isinstance(step_input, str) and ".fits" not in step_input:
            asn = self.load_as_level2_asn(step_input)
            # The input_bkg_list variable can be used to provide the name of the
            # file processed through the assign_wcs step
            if input_bkg_list is not None:
                model_with_wcs = input_bkg_list
            else:
                model_with_wcs = None
            input_model, members_by_type = asn_get_data(asn, model_with_wcs=model_with_wcs)
            bkg_list = members_by_type["background"]
        else:
            input_model = datamodels.open(step_input)
            if input_bkg_list is not None:
                bkg_list = input_bkg_list
            else:
                bkg_list = self.bkg_list

        if input_model.meta.exposure.type in ["NIS_WFSS", "NRC_WFSS"]:
            # Get the reference file names
            bkg_name = self.get_reference_file(input_model, "wfssbkg")
            wlrange_name = self.get_reference_file(input_model, "wavelengthrange")
            self.log.info("Using WFSSBKG reference file %s", bkg_name)
            self.log.info("Using WavelengthRange reference file %s", wlrange_name)

            # Do the background subtraction for WFSS/GRISM data
            rescaler_kwargs = {
                "p": self.wfss_outlier_percent,
                "maxiter": self.wfss_maxiter,
                "delta_rms_thresh": self.wfss_rms_stop / 100,
            }
            result = subtract_wfss_bkg(
                input_model,
                bkg_name,
                wlrange_name,
                self.wfss_mmag_extract,
                rescaler_kwargs=rescaler_kwargs,
            )
            if result is None:
                result = input_model
                result.meta.cal_step.back_sub = "SKIPPED"
            else:
                result.meta.cal_step.back_sub = "COMPLETE"

        else:
            # Make sure that the background list is not empty for this case,
            # or report and skip the step
            if bkg_list is None:
                self.log.warning("* No background list provided * Skipping step.")
                input_model.meta.cal_step.back_sub = "SKIPPED"
                return input_model

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
                        if np.allclose(
                            (input_xtilt, input_ytilt),
                            (bkg_xtilt, bkg_ytilt),
                            atol=tolerance,
                            rtol=0,
                        ):
                            pass
                        else:
                            do_sub = False
                            break
            # Do the background subtraction
            if do_sub:
                bkg_model, result = background_sub(input_model, bkg_list, self.sigma, self.maxiters)
                result.meta.cal_step.back_sub = "COMPLETE"
                if self.save_combined_background:
                    comb_bkg_path = self.save_model(bkg_model, suffix=self.bkg_suffix, force=True)
                    self.log.info(f"Combined background written to {comb_bkg_path}.")

            else:
                result = input_model.copy()
                result.meta.cal_step.back_sub = "SKIPPED"
                self.log.warning("Skipping background subtraction")
                self.log.warning(
                    "GWA_XTIL and GWA_YTIL source values are not the same as bkg values"
                )

        # Cleanup
        del input_model

        return result

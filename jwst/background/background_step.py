#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from .background_sub import background_sub
from .background_sub_wfss import subtract_wfss_bkg
from .asn_intake import asn_get_data
from jwst.background.background_sub_soss import subtract_soss_bkg

import numpy as np


__all__ = ["BackgroundStep"]


class BackgroundStep(Step):
    """Subtract background exposures from target exposures."""

    class_alias = "bkg_subtract"

    spec = """
        bkg_list = force_list(default=None)  # List of background files. Ignored for WFSS or if asn is provided
        save_combined_background = boolean(default=False)  # Save combined background image
        sigma = float(default=3.0)  # Clipping threshold
        maxiters = integer(default=None)  # Number of clipping iterations
        soss_source_percentile = float(default=35.0) # Threshold flux percentile to mask out source pixels
        soss_bkg_percentile = float_list(min=2, max=2, default=None) # Background percentiles to use; default is [25.0, 50.0]
        wfss_mmag_extract = float(default=None)  # WFSS minimum abmag to extract
        wfss_maxiter = integer(default=5)  # WFSS iterative outlier rejection max iterations
        wfss_rms_stop = float(default=0)  # WFSS iterative outlier rejection RMS improvement threshold (percent)
        wfss_outlier_percent = float(default=1)  # WFSS outlier percentile to reject per iteration
    """  # noqa: E501

    # These reference files are only used for WFSS/GRISM or SOSS data.
    reference_file_types = ["bkg", "wavelengthrange"]

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
            File name list of background exposures.

        Returns
        -------
        result : ImageModel or IFUImageModel
            The background-subtracted target data model
        """
        asn = self.load_as_level2_asn(step_input)
        input_model, members_by_type = asn_get_data(asn)
        result = input_model.copy()

        if result.meta.exposure.type in ["NIS_WFSS", "NRC_WFSS"]:
            # Get the reference file names
            bkg_name = self.get_reference_file(result, "bkg")
            wlrange_name = self.get_reference_file(result, "wavelengthrange")

            if bkg_name == "N/A":
                self.log.warning("No BKG reference file found. Skipping background subtraction.")
                result.meta.cal_step.bkg_subtract = "SKIPPED"
                input_model.close()
                return result

            self.log.info("Using BKG reference file %s", bkg_name)
            self.log.info("Using WavelengthRange reference file %s", wlrange_name)

            # Do the background subtraction for WFSS/GRISM data
            rescaler_kwargs = {
                "p": self.wfss_outlier_percent,
                "maxiter": self.wfss_maxiter,
                "delta_rms_thresh": self.wfss_rms_stop / 100,
            }
            result = subtract_wfss_bkg(
                result,
                bkg_name,
                wlrange_name,
                self.wfss_mmag_extract,
                rescaler_kwargs=rescaler_kwargs,
            )
            if result is None:
                result = input_model.copy()
                result.meta.cal_step.bkg_subtract = "SKIPPED"
            else:
                result.meta.cal_step.bkg_subtract = "COMPLETE"

        elif input_model.meta.exposure.type == "NIS_SOSS":
            # Fetch the background reference filename
            bkg_name = self.get_reference_file(input_model, "bkg")
            self.log.info("Using BKG reference file %s", bkg_name)

            if self.soss_bkg_percentile is None:
                soss_bkg_percentile = [25.0, 50.0]
            else:
                soss_bkg_percentile = self.soss_bkg_percentile

            result = subtract_soss_bkg(
                input_model, bkg_name, self.soss_source_percentile, soss_bkg_percentile
            )
            if result is None:
                result = input_model
                result.meta.cal_step.back_sub = "SKIPPED"
            else:
                result.meta.cal_step.back_sub = "COMPLETE"

        else:
            # Get the background files to be subtracted
            if len(members_by_type["background"]) >= 1:
                bkg_list = members_by_type["background"]
            elif input_bkg_list is not None:
                if isinstance(input_bkg_list, str):
                    if "," in input_bkg_list:
                        bkg_list = input_bkg_list.split(sep=",")
                    else:
                        bkg_list = [input_bkg_list]
                else:
                    bkg_list = input_bkg_list
            else:
                bkg_list = self.bkg_list

            # Make sure to catch a trailing comma and ignore it
            if bkg_list is not None:
                bkg_list = [bg for bg in bkg_list if bg]

            # Make sure that the background list is not empty for this case,
            # or report and skip the step
            if bkg_list is None or len(bkg_list) == 0:
                self.log.warning("* No background list provided * Skipping step.")
                result.meta.cal_step.bkg_subtract = "SKIPPED"
                return result

            # check if input data is NRS_IFU
            tolerance = 1.0e-8
            do_sub = True
            if result.meta.instrument.name in ["NIRSPEC"]:
                # check if GWA_XTIL & GWA_YTIL values of source
                # background are the same. If not skip step
                input_xtilt = result.meta.instrument.gwa_xtilt
                input_ytilt = result.meta.instrument.gwa_ytilt
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
                bkg_model, result = background_sub(result, bkg_list, self.sigma, self.maxiters)
                result.meta.cal_step.bkg_subtract = "COMPLETE"
                if self.save_combined_background:
                    comb_bkg_path = self.save_model(bkg_model, suffix=self.bkg_suffix, force=True)
                    self.log.info(f"Combined background written to {comb_bkg_path}.")

            else:
                result.meta.cal_step.bkg_subtract = "SKIPPED"
                self.log.warning("Skipping background subtraction")
                self.log.warning(
                    "GWA_XTIL and GWA_YTIL source values are not the same as bkg values"
                )

        input_model.close()
        return result

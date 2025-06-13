from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from jwst.lib import pipe_utils
from . import reference_pixels
from . import irs2_subtract_reference


__all__ = ["RefPixStep"]


class RefPixStep(Step):
    """Use reference pixels to correct bias drifts."""

    class_alias = "refpix"

    spec = """
        odd_even_columns = boolean(default=True) # Compute reference signal separately for even/odd columns
        use_side_ref_pixels = boolean(default=True) # Use side reference pixels for reference signal for each row
        side_smoothing_length = integer(default=11) # Median window smoothing height for side reference signal
        side_gain = float(default=1.0) # Multiplicative factor for side reference signal before subtracting from rows
        odd_even_rows = boolean(default=True) # Compute reference signal separately for even- and odd-numbered rows
        ovr_corr_mitigation_ftr = float(default=3.0) # Factor to avoid overcorrection of bad reference pixels for IRS2
        preserve_irs2_refpix = boolean(default=False) # Preserve reference pixels in output
        irs2_mean_subtraction = boolean(default=False) # Apply a mean offset subtraction before IRS2 correction
        refpix_algorithm = option("median", "sirs", default="median") # NIR full-frame side pixel algorithm
        sigreject = float(default=4.0) # Number of sigmas to reject as outliers
        gaussmooth = float(default=1.0) # Width of Gaussian smoothing kernel to use as a low-pass filter
        halfwidth = integer(default=30) # Half-width of convolution kernel to build
    """  # noqa: E501

    reference_file_types = ["refpix", "sirskernel"]

    def process(self, step_input):
        """
        Execute the reference pixel correction step.

        Parameters
        ----------
        step_input : DataModel
            Input datamodel to the step

        Returns
        -------
        result : DataModel
            Result of applying the reference pixel correction step
        """
        conv_kernel_params = {
            "refpix_algorithm": self.refpix_algorithm,
            "sirs_kernel_model": None,
            "sigreject": self.sigreject,
            "gaussmooth": self.gaussmooth,
            "halfwidth": self.halfwidth,
        }

        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # Work on a copy
            result = input_model.copy()

            if pipe_utils.is_irs2(result):
                # Flag bad reference pixels first
                irs2_subtract_reference.flag_bad_refpix(
                    result, n_sigma=self.ovr_corr_mitigation_ftr, flag_only=True
                )

                # If desired, do the normal refpix correction before IRS2, without
                # side pixel handling
                if self.irs2_mean_subtraction:
                    if self.use_side_ref_pixels:
                        self.log.info("Turning off side pixel correction for IRS2")
                        self.use_side_ref_pixels = False
                    reference_pixels.correct_model(
                        result,
                        self.odd_even_columns,
                        self.use_side_ref_pixels,
                        self.side_smoothing_length,
                        self.side_gain,
                        self.odd_even_rows,
                        conv_kernel_params,
                    )

                # Now that values are updated, replace bad reference pixels
                irs2_subtract_reference.flag_bad_refpix(result, replace_only=True)

                # Get the necessary refpix reference file for IRS2 correction
                self.irs2_name = self.get_reference_file(result, "refpix")
                self.log.info(f"Using refpix reference file: {self.irs2_name}")

                # Check for a valid reference file
                if self.irs2_name == "N/A":
                    self.log.warning("No refpix reference file found")
                    self.log.warning("RefPix step will be skipped")
                    result.meta.cal_step.refpix = "SKIPPED"
                    return result

                # Load the reference file into a datamodel
                irs2_model = datamodels.IRS2Model(self.irs2_name)

                # Apply the IRS2 correction scheme
                result = irs2_subtract_reference.correct_model(
                    result, irs2_model, preserve_refpix=self.preserve_irs2_refpix
                )

                if result.meta.cal_step.refpix != "SKIPPED":
                    result.meta.cal_step.refpix = "COMPLETE"
                del irs2_model
                return result

            else:
                # Not an NRS IRS2 exposure. Do the normal refpix correction.

                # Get the reference file from CRDS or use user-supplied one
                # only for NIR full-frame data
                if self.refpix_algorithm == "sirs":
                    if (
                        input_model.meta.instrument.name != "MIRI"
                        and "FULL" in input_model.meta.subarray.name
                    ):
                        sirs_ref_filename = self.get_reference_file(result, "sirskernel")
                        if sirs_ref_filename == "N/A":
                            self.log.warning(
                                "No reference file found for the optimized convolution kernel."
                            )
                            self.log.warning(
                                "REFPIX step will use the running median algorithm for side pixels."
                            )
                        else:
                            self.log.info(f"Using SIRS reference file: {sirs_ref_filename}")
                            sirs_kernel_model = datamodels.SIRSKernelModel(sirs_ref_filename)
                            conv_kernel_params["sirs_kernel_model"] = sirs_kernel_model
                    elif input_model.meta.instrument.name == "MIRI":
                        self.log.info(
                            "Simple Improved Reference Subtraction (SIRS) "
                            "not applied for MIRI data."
                        )
                    elif "FULL" not in input_model.meta.subarray.name:
                        self.log.info(
                            "Simple Improved Reference Subtraction (SIRS) "
                            "not applied for subarray data."
                        )

                status = reference_pixels.correct_model(
                    result,
                    self.odd_even_columns,
                    self.use_side_ref_pixels,
                    self.side_smoothing_length,
                    self.side_gain,
                    self.odd_even_rows,
                    conv_kernel_params,
                )

                if status == reference_pixels.REFPIX_OK:
                    result.meta.cal_step.refpix = "COMPLETE"
                elif status == reference_pixels.SUBARRAY_DOESNTFIT:
                    self.log.warning("Subarray doesn't fit in full-sized array")
                    result.meta.cal_step.refpix = "SKIPPED"
                elif status == reference_pixels.BAD_REFERENCE_PIXELS:
                    self.log.warning("No valid reference pixels, refpix step skipped")
                    result.meta.cal_step.refpix = "SKIPPED"
                elif status == reference_pixels.SUBARRAY_SKIPPED:
                    result.meta.cal_step.refpix = "SKIPPED"

                return result

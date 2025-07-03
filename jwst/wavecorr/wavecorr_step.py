from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import wavecorr

__all__ = ["WavecorrStep"]


class WavecorrStep(Step):
    """
    Apply wavelength offsets to off-center NIRSpec sources.

    Wavelength corrections are applied only to point sources in NIRSpec MOS
    and FS data.

    The algorithm uses a reference file which is a look-up table of
    wavelength_correction as a function of slit_x_position and wavelength.
    The x direction is the one parallel to dispersion/wavelength for
    both MOS and FS slits.

    The slit_x_position is read from the `source_xpos` attribute in the input
    slit metadata.  For MOS data, the x position is read from the msa_metadata_file
    in the assign_wcs step.  For FS data, the x position is calculated from
    the dither `x_offset` value in the extract_2d step.

    The wavelength value used to look up the wavelength correction at each dispersion
    element is an average of the wavelength values in the cross-dispersion direction
    at that element.
    """

    class_alias = "wavecorr"

    spec = """
    """  # noqa: E501

    reference_file_types = ["wavecorr"]

    def process(self, step_input):
        """
        Apply the wavelength correction to the input data.

        Parameters
        ----------
        step_input : DataModel or str
            Input data to correct.

        Returns
        -------
        DataModel
            The corrected data.
        """
        wavecorr_supported_modes = ["NRS_FIXEDSLIT", "NRS_MSASPEC", "NRS_BRIGHTOBJ", "NRS_AUTOFLAT"]

        # Open the input
        with datamodels.open(step_input) as input_model:
            # Check for valid exposure type
            exp_type = input_model.meta.exposure.type.upper()
            if exp_type not in wavecorr_supported_modes:
                self.log.info(f"Skipping wavecorr correction for EXP_TYPE {exp_type}")
                input_model.meta.cal_step.wavecorr = "SKIPPED"
                return input_model

            # Check for prerequisites
            if (
                hasattr(input_model.meta.cal_step, "assign_wcs")
                and input_model.meta.cal_step.assign_wcs == "SKIPPED"
            ):
                self.log.warning("assign_wcs was skipped")
                self.log.warning("Wavecorr step will be skipped")
                input_model.meta.cal_step.wavecorr = "SKIPPED"
                return input_model

            # Check for existence of WCS
            if isinstance(input_model, datamodels.SlitModel):
                _check_slit_metadata_attributes(input_model)
            elif isinstance(input_model, datamodels.MultiSlitModel):
                _check_slit_metadata_attributes(input_model.slits[0])
            else:
                raise TypeError(
                    f"Input model must be a SlitModel or MultiSlitModel, not {type(input_model)}"
                )

            # Get the reference file
            reffile = self.get_reference_file(input_model, "wavecorr")
            self.log.info(f"Using WAVECORR reference file {reffile}")
            if reffile == "N/A":
                self.log.warning("No WAVECORR reference file found")
                self.log.warning("Wavecorr step will be skipped")
                input_model.meta.cal_step.wavecorr = "SKIPPED"
                return input_model

            # Apply the correction
            output_model = wavecorr.do_correction(input_model, reffile)

        return output_model


def _check_slit_metadata_attributes(slit):
    if not (hasattr(slit.meta, "wcs") and slit.meta.wcs is not None):
        raise AttributeError(
            "Input model does not have a WCS object; assign_wcs should be run before wavecorr."
        )

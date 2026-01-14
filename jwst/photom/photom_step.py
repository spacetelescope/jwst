import logging

from stdatamodels.jwst import datamodels

from jwst.photom import photom
from jwst.stpipe import Step

__all__ = ["PhotomStep"]

log = logging.getLogger(__name__)


class PhotomStep(Step):
    """Apply photometric calibration."""

    class_alias = "photom"

    spec = """
        inverse = boolean(default=False)    # Invert the operation
        source_type = string(default=None)  # Process as specified source type
        apply_time_correction = boolean(default=True) # Apply time dependent corrections if available
    """  # noqa: E501

    reference_file_types = ["photom", "area"]

    def process(self, input_data):
        """
        Execute the photom calibration step.

        This step loads photometric conversion information from
        reference files and attaches or applies them to the input science
        data model.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel`
            Input science data file to be calibrated.

        Returns
        -------
        result : DataModel
            DataModel with the photom correction applied
        """
        output_model = self.prepare_output(input_data)

        # Report the detected type of input model
        log.debug(f"Input is {str(output_model)}")
        if not isinstance(
            output_model,
            datamodels.CubeModel
            | datamodels.SlitModel
            | datamodels.ImageModel
            | datamodels.IFUImageModel
            | datamodels.MultiSlitModel
            | datamodels.TSOMultiSpecModel,
        ):
            log.warning(
                "Input is not one of the supported model types: "
                "CubeModel, ImageModel, IFUImageModel, "
                "SlitModel, MultiSlitModel, or TSOMultiSpecModel."
            )
            log.warning("Photom step will be skipped")
            output_model.meta.cal_step.photom = "SKIPPED"
            return output_model

        # Setup reference files and whether previous correction information
        # should be used.
        if self.use_correction_pars and self.correction_pars:
            log.info("Using previously specified correction parameters.")
            correction_pars = self.correction_pars
            phot_filename = correction_pars["refs"]["photom"]
            area_filename = correction_pars["refs"]["area"]
        else:
            correction_pars = None
            phot_filename = self.get_reference_file(output_model, "photom")
            area_filename = self.get_reference_file(output_model, "area")

        log.info("Using photom reference file: %s", phot_filename)
        log.info("Using area reference file: %s", area_filename)

        # Check for a valid photom reference file
        if phot_filename == "N/A":
            log.warning("No PHOTOM reference file found")
            log.warning("Photom step will be skipped")
            output_model.meta.cal_step.photom = "SKIPPED"
            return output_model

        try:
            # Do the correction
            phot = photom.DataSet(
                output_model,
                self.inverse,
                self.source_type,
                self.apply_time_correction,
                correction_pars,
            )
            result = phot.apply_photom(phot_filename, area_filename)
            result.meta.cal_step.photom = "COMPLETE"
            self.correction_pars = phot.correction_pars
            self.correction_pars["refs"] = {"photom": phot_filename, "area": area_filename}

        except photom.DataModelTypeError:
            # should trip e.g. for NIRISS SOSS data in FULL subarray
            log.error(
                f"Unexpected data model type {str(output_model)} for "
                f"{output_model.meta.instrument.name.upper()}. "
                "Photom will be skipped."
            )
            output_model.meta.cal_step.photom = "SKIPPED"
            return output_model

        return result

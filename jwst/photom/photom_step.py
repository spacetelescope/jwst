#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import photom

__all__ = ["PhotomStep"]


class PhotomStep(Step):
    """
    Class for applying photom calibration step.

    This step loads photometric conversion information from
    reference files and attaches or applies them to the input science
    data model
    """

    class_alias = "photom"

    spec = """
        inverse = boolean(default=False)    # Invert the operation
        source_type = string(default=None)  # Process as specified source type.
        mrs_time_correction = boolean(default=True) # Apply the MIRI MRS time dependent correction
    """  # noqa: E501

    reference_file_types = ["photom", "area"]

    def process(self, input_data):
        """
        Execute the photom calibration step.

        Parameters
        ----------
        input_data : str
            Input science data file to be calibrated

        Returns
        -------
        result : DataModel
            DataModel with the photom correction applied
        """
        try:
            input_model = datamodels.open(input_data)
        except OSError:
            self.log.error("Input can not be opened as a Model.")

        # Report the detected type of input model
        self.log.debug(f"Input is {str(input_model)}")
        if not isinstance(
            input_model,
            datamodels.CubeModel
            | datamodels.SlitModel
            | datamodels.ImageModel
            | datamodels.IFUImageModel
            | datamodels.MultiSlitModel
            | datamodels.TSOMultiSpecModel,
        ):
            self.log.warning(
                "Input is not one of the supported model types: "
                "CubeModel, ImageModel, IFUImageModel, "
                "SlitModel, MultiSlitModel, or TSOMultiSpecModel."
            )

        # Setup reference files and whether previous correction information
        # should be used.
        if self.use_correction_pars and self.correction_pars:
            self.log.info("Using previously specified correction parameters.")
            correction_pars = self.correction_pars
            phot_filename = correction_pars["refs"]["photom"]
            area_filename = correction_pars["refs"]["area"]
        else:
            correction_pars = None
            phot_filename = self.get_reference_file(input_model, "photom")
            area_filename = self.get_reference_file(input_model, "area")

        self.log.info("Using photom reference file: %s", phot_filename)
        self.log.info("Using area reference file: %s", area_filename)

        # Check for a valid photom reference file
        if phot_filename == "N/A":
            self.log.warning("No PHOTOM reference file found")
            self.log.warning("Photom step will be skipped")
            result = input_model.copy()
            result.meta.cal_step.photom = "SKIPPED"
            return result

        try:
            # Do the correction
            phot = photom.DataSet(
                input_model,
                self.inverse,
                self.source_type,
                self.mrs_time_correction,
                correction_pars,
            )
            result = phot.apply_photom(phot_filename, area_filename)
            result.meta.cal_step.photom = "COMPLETE"
            self.correction_pars = phot.correction_pars
            self.correction_pars["refs"] = {"photom": phot_filename, "area": area_filename}

        except photom.DataModelTypeError:
            # should trip e.g. for NIRISS SOSS data in FULL subarray
            self.log.error(
                f"Unexpected data model type {str(input_model)} for "
                f"{input_model.meta.instrument.name.upper()}. "
                "Photom will be skipped."
            )
            result = input_model.copy()
            result.meta.cal_step.photom = "SKIPPED"
            return result
        return result

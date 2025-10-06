import logging

from stdatamodels.jwst import datamodels

from jwst.linearity import linearity
from jwst.stpipe import Step

__all__ = ["LinearityStep"]

log = logging.getLogger(__name__)


class LinearityStep(Step):
    """Perform a correction for non-linear detector response, using the polynomial method."""

    class_alias = "linearity"

    spec = """
    """  # noqa: E501

    reference_file_types = ["linearity"]

    def process(self, step_input):
        """
        Read in linearity correction and apply it to science data.

        Parameters
        ----------
        step_input : `~stdatamodels.jwst.datamodels.RampModel`
            The input ramp model.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            The output ramp model with linearity correction applied.
        """
        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        # Get the name of the linearity reference file to use
        lin_name = self.get_reference_file(result, "linearity")
        log.info("Using Linearity reference file %s", lin_name)

        # Check for a valid reference file
        if lin_name == "N/A":
            log.warning("No Linearity reference file found")
            log.warning("Linearity step will be skipped")
            result.meta.cal_step.linearity = "SKIPPED"
            return result

        # Open the linearity reference file data model
        lin_model = datamodels.LinearityModel(lin_name)

        # Do the linearity correction
        result = linearity.do_correction(result, lin_model)
        result.meta.cal_step.linearity = "COMPLETE"

        # Cleanup
        del lin_model

        return result

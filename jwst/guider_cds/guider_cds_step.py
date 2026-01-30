import logging

from stdatamodels.jwst import datamodels

from jwst.guider_cds import guider_cds
from jwst.stpipe import Step

__all__ = ["GuiderCdsStep"]

log = logging.getLogger(__name__)


class GuiderCdsStep(Step):
    """Calculate the count rate for each pixel for FGS modes."""

    class_alias = "guider_cds"

    reference_file_types = ["gain", "readnoise"]

    def process(self, input_data):
        """
        Execute the step.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.GuiderRawModel`
            The input file name or datamodel.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.GuiderRawModel` or \
                       `~stdatamodels.jwst.datamodels.GuiderCalModel`
            This will be a GuiderRawModel if the step was skipped; otherwise,
            it will be a GuiderCalModel containing calibrated count rates.
        """
        output_model = self.prepare_output(input_data, open_as_type=datamodels.GuiderRawModel)

        # Get the gain reference file
        gain_filename = self.get_reference_file(output_model, "gain")
        if gain_filename == "N/A":
            log.warning("No GAIN reference file found!")
            log.warning("guider_cds step will be skipped.")
            output_model.meta.cal_step.guider_cds = "SKIPPED"
            return output_model

        log.info("Using GAIN reference file: %s", gain_filename)
        gain_model = datamodels.GainModel(gain_filename)

        # Get the readnoise reference file
        readnoise_filename = self.get_reference_file(output_model, "readnoise")
        if readnoise_filename == "N/A":
            log.warning("No READNOISE reference file found!")
            log.warning("guider_cds step will be skipped.")
            output_model.meta.cal_step.guider_cds = "SKIPPED"
            return output_model

        log.info("Using READNOISE reference file: %s", readnoise_filename)
        readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)

        result = guider_cds.guider_cds(output_model, gain_model, readnoise_model)
        result.meta.cal_step.guider_cds = "COMPLETE"

        # Output is a new model, so close the input if it was opened here
        if output_model is not input_data:
            output_model.close()

        return result

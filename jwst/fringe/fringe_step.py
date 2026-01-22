import logging

from stdatamodels.jwst import datamodels

from jwst.fringe import fringe
from jwst.stpipe import Step

__all__ = ["FringeStep"]

log = logging.getLogger(__name__)


class FringeStep(Step):
    """Apply fringe correction to a science image using a fringe reference image."""

    class_alias = "fringe"

    reference_file_types = ["fringe"]

    def process(self, input_data):
        """
        Apply fringe correction to a science image using a fringe reference image.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.IFUImageModel`
            Input MIRI MRS science file name or datamodel.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
            Fringe corrected MIRI MRS science data.
        """
        output_model = self.prepare_output(input_data)

        # Open the reference file
        fringe_filename = self.get_reference_file(output_model, "fringe")
        log.info("Using FRINGE reference file: %s", fringe_filename)

        # Check for a valid reference file
        if fringe_filename == "N/A":
            log.warning("No FRINGE reference file found")
            log.warning("Fringe step will be skipped")
            output_model.meta.cal_step.fringe = "SKIPPED"
            return output_model

        # Load the fringe reference file
        with datamodels.FringeModel(fringe_filename) as fringe_model:
            # Do the correction
            output_model = fringe.apply_fringe(output_model, fringe_model)

        output_model.meta.cal_step.fringe = "COMPLETE"
        return output_model

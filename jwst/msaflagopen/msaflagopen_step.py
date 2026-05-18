"""Flag pixels affected by open MSA shutters in NIRSpec exposures."""

import logging

from jwst.assign_wcs.util import get_wcs_reference_files
from jwst.msaflagopen import msaflag_open
from jwst.stpipe import Step

__all__ = ["MSAFlagOpenStep"]

log = logging.getLogger(__name__)


class MSAFlagOpenStep(Step):
    """Flag pixels affected by MSA failed open shutters."""

    class_alias = "msa_flagging"

    spec = """
    """  # noqa: E501

    reference_file_types = ["msaoper"]

    def process(self, input_data):
        """
        Flag data affected by open shutters.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.ImageModel`
            Science data to be corrected.

        Returns
        -------
        `~stdatamodels.jwst.datamodels.ImageModel`
            Science data with DQ array modified.
        """
        # Open the input data model
        output_model = self.prepare_output(input_data)

        self.reference_name = self.get_reference_file(output_model, "msaoper")
        log.info("Using reference file %s", self.reference_name)

        # Check for a valid reference file
        if self.reference_name == "N/A":
            log.warning("No reference file found")
            log.warning("Step will be skipped")
            output_model.meta.cal_step.msa_flagging = "SKIPPED"
            return output_model

        # Get the reference file names for constructing the WCS pipeline
        wcs_reffile_names = get_wcs_reference_files(output_model)

        # Do the DQ flagging
        output_model = msaflag_open.do_correction(
            output_model, self.reference_name, wcs_reffile_names
        )

        # set the step status to complete
        output_model.meta.cal_step.msa_flagging = "COMPLETE"

        return output_model

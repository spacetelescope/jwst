"""Initialize the data with data quality information."""

import logging

from astropy.io import fits
from stdatamodels.jwst import datamodels

from jwst.dq_init import dq_initialization
from jwst.stpipe import Step

__all__ = ["DQInitStep"]

log = logging.getLogger(__name__)


class DQInitStep(Step):
    """
    Initialize the Data Quality (DQ) extension from the mask reference file.

    The dq_init step initializes the pixeldq attribute of the
    input datamodel using the MASK reference file.  For some
    FGS exp_types, initialize the dq attribute of the input model
    instead.  The dq attribute of the MASK model is bitwise OR'd
    with the pixeldq (or dq) attribute of the input model.
    """

    class_alias = "dq_init"

    spec = """
        user_supplied_dq = string(default=None)  # FITS user-supplied DQ file
    """  # noqa: E501
    reference_file_types = ["mask"]

    def process(self, step_input):
        """
        Perform the ``dq_init`` calibration step.

        Parameters
        ----------
        step_input : str or `~stdatamodels.jwst.datamodels.RampModel`
            Input JWST filename or datamodel.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.RampModel` or \
                       `~stdatamodels.jwst.datamodels.SuperstripeRampModel` or \
                       `~stdatamodels.jwst.datamodels.GuiderRawModel`
            Result JWST datamodel.
        """
        # Open the input as a ramp
        result = self.prepare_output(step_input, open_as_ramp=True)

        # Retrieve the mask reference file name
        mask_filename = self.get_reference_file(result, "mask")
        log.info("Using MASK reference file %s", mask_filename)

        # Check for a valid reference file
        if mask_filename == "N/A":
            log.warning("No MASK reference file found")
            log.warning("DQ initialization step will be skipped")
            result.meta.cal_step.dq_init = "SKIPPED"
            return result

        # Load the reference file
        mask_model = datamodels.MaskModel(mask_filename)

        # Load the user-supplied DQ file, if applicable
        if self.user_supplied_dq is not None:
            log.info(
                f"User-supplied DQ {self.user_supplied_dq} given. Adding it to DQ initialization."
            )
            user_dq = fits.getdata(self.user_supplied_dq)
        else:
            user_dq = None

        # Apply the step
        result = dq_initialization.do_dqinit(result, mask_model, user_dq=user_dq)

        # Cleanup
        del mask_model
        if user_dq is not None:
            del user_dq

        return result

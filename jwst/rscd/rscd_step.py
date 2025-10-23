import logging

from stdatamodels.jwst import datamodels

from jwst.rscd import rscd_sub
from jwst.stpipe import Step

__all__ = ["RscdStep"]

log = logging.getLogger(__name__)


class RscdStep(Step):
    """
    Flag the first N groups of MIRI data to 'DO_NOT_USE' in the 2nd and later integrations.

    The number of groups, N, for which to
    set the GROUPDQ flag to 'DO_NOT_USE' is read in from the RSCD reference file. This number
    depends on the readout model and subarray size. The step checks that the total number of groups
    in an integration is greater than N+3 before flagging the GROUPDQ array. If the number of groups
    is less than N+3 then no flagging is performed, because doing so would leave too few groups
    to work with in later steps.
    """

    class_alias = "rscd"

    spec = """
    """  # noqa: E501

    reference_file_types = ["rscd"]

    def process(self, step_input):
        """
        Flag the initial groups to 'DO_NOT_USE' in the 2nd and later integrations.

        The number of initial groups to flag is read in from the RSCD reference file. This number
        varies based on readout mode and subarray size.

        Parameters
        ----------
        step_input : `~stdatamodels.jwst.datamodels.RampModel`
            Ramp datamodel to be corrected, or the path to the ramp file.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            Ramp datamodel with initial groups in an integration flagged as DO_NOT_USE.
            Flags are only set of integration 2 and higher.
        """
        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        # check the data is MIRI data
        detector = result.meta.instrument.detector
        if not detector.startswith("MIR"):
            log.warning("RSCD correction is only for MIRI data")
            log.warning("RSCD step will be skipped")
            result.meta.cal_step.rscd = "SKIPPED"
            return result

        # Get the name of the rscd reference file to use
        rscd_name = self.get_reference_file(result, "rscd")
        log.info("Using RSCD reference file %s", rscd_name)

        # Check for a valid reference file
        if rscd_name == "N/A":
            log.warning("No RSCD reference file found")
            log.warning("RSCD step will be skipped")
            result.meta.cal_step.rscd = "SKIPPED"
            return result

        # Load the rscd ref file data model
        rscd_model = datamodels.RSCDModel(rscd_name)

        # Do the rscd correction
        result = rscd_sub.do_correction(result, rscd_model)

        # Cleanup
        del rscd_model

        return result

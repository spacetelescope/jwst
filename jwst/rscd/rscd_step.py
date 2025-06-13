from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import rscd_sub

__all__ = ["RscdStep"]


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
        step_input : RampModel
            Ramp datamodel to be corrected, or the path to the ramp file.

        Returns
        -------
        result : RampModel
            Ramp datamodel with initial groups in an integration flagged as DO_NOT_USE.
            Flags are only set of integration 2 and higher.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if not detector.startswith("MIR"):
                self.log.warning("RSCD correction is only for MIRI data")
                self.log.warning("RSCD step will be skipped")
                input_model.meta.cal_step.rscd = "SKIPPED"
                return input_model

            # Get the name of the rscd reference file to use
            self.rscd_name = self.get_reference_file(input_model, "rscd")
            self.log.info("Using RSCD reference file %s", self.rscd_name)

            # Check for a valid reference file
            if self.rscd_name == "N/A":
                self.log.warning("No RSCD reference file found")
                self.log.warning("RSCD step will be skipped")
                input_model.meta.cal_step.rscd = "SKIPPED"
                return input_model

            # Load the rscd ref file data model
            rscd_model = datamodels.RSCDModel(self.rscd_name)

            # Work on a copy
            result = input_model.copy()

            # Do the rscd correction
            result = rscd_sub.do_correction(result, rscd_model)

            # Cleanup
            del rscd_model

        return result

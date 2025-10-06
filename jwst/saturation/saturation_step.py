import logging

from stdatamodels.jwst import datamodels

from jwst.lib import pipe_utils, reffile_utils
from jwst.saturation import saturation
from jwst.stpipe import Step

__all__ = ["SaturationStep"]

log = logging.getLogger(__name__)


class SaturationStep(Step):
    """Set saturation flags."""

    class_alias = "saturation"

    spec = """
        n_pix_grow_sat = integer(default=1) # number of layers adjacent pixels to flag
        use_readpatt = boolean(default=True) # Use grouped read pattern information to assist with flagging
    """  # noqa: E501

    reference_file_types = ["saturation", "superbias"]

    def process(self, step_input):
        """
        Set the saturation flags.

        Parameters
        ----------
        step_input : `~stdatamodels.jwst.datamodels.RampModel` or str
            Input datamodel or string name of the fits file.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            Output datamodel with saturation flags set for saturated pixels.
        """
        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        # Get the name of the saturation reference file
        ref_name = self.get_reference_file(result, "saturation")
        bias_name = self.get_reference_file(result, "superbias")
        log.info("Using SATURATION reference file %s", ref_name)
        log.info("Using SUPERBIAS reference file %s", bias_name)

        # Check for a valid reference file
        if ref_name == "N/A":
            log.warning("No SATURATION reference file found")
            log.warning("Saturation step will be skipped")
            result.meta.cal_step.saturation = "SKIPPED"
            return result

        # Open the reference file data model
        ref_model = datamodels.SaturationModel(ref_name)

        # Open the superbias if one is available
        bias_model = None
        if bias_name != "N/A":
            bias_model = datamodels.SuperBiasModel(bias_name)
            # Check for subarray mode and extract subarray from the
            # bias reference data if necessary
            if not reffile_utils.ref_matches_sci(result, bias_model):
                bias_model = reffile_utils.get_subarray_model(result, bias_model)

        # Do the saturation check
        if pipe_utils.is_irs2(result):
            result = saturation.irs2_flag_saturation(
                result, ref_model, self.n_pix_grow_sat, self.use_readpatt, bias_model=bias_model
            )
        else:
            result = saturation.flag_saturation(
                result, ref_model, self.n_pix_grow_sat, self.use_readpatt, bias_model=bias_model
            )
        result.meta.cal_step.saturation = "COMPLETE"

        # Cleanup
        del ref_model

        return result

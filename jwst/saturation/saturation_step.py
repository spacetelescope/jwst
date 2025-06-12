#! /usr/bin/env python
from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from jwst.lib import pipe_utils, reffile_utils
from . import saturation


__all__ = ["SaturationStep"]


class SaturationStep(Step):
    """Sets saturation flags."""

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
        step_input : DataModel or str
            Input datamodel or string name of the fits file.

        Returns
        -------
        result : DataModel
            Output datamodel with saturation flags set for saturated pixels.
        """
        # Open the input data model
        with datamodels.open(step_input) as input_model:
            # Get the name of the saturation reference file
            self.ref_name = self.get_reference_file(input_model, "saturation")
            self.bias_name = self.get_reference_file(input_model, "superbias")
            self.log.info("Using SATURATION reference file %s", self.ref_name)
            self.log.info("Using SUPERBIAS reference file %s", self.bias_name)

            # Check for a valid reference file
            if self.ref_name == "N/A":
                self.log.warning("No SATURATION reference file found")
                self.log.warning("Saturation step will be skipped")
                input_model.meta.cal_step.saturation = "SKIPPED"
                return input_model

            # Open the reference file data model
            ref_model = datamodels.SaturationModel(self.ref_name)

            # Open the superbias if one is available
            bias_model = None
            if self.bias_name != "N/A":
                bias_model = datamodels.SuperBiasModel(self.bias_name)
                # Check for subarray mode and extract subarray from the
                # bias reference data if necessary
                if not reffile_utils.ref_matches_sci(input_model, bias_model):
                    bias_model = reffile_utils.get_subarray_model(input_model, bias_model)

            # Work on a copy
            result = input_model.copy()

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

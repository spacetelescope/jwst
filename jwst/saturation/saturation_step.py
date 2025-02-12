#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from ..lib import pipe_utils
from . import saturation


__all__ = ["SaturationStep"]


class SaturationStep(Step):
    """
    This Step sets saturation flags.
    """

    class_alias = "saturation"

    spec = """
        n_pix_grow_sat = integer(default=1) # number of layers adjacent pixels to flag
        use_readpatt = boolean(default=True) # Use grouped read pattern information to assist with flagging
    """ # noqa: E501

    reference_file_types = ['saturation']

    def process(self, step_input):

        # Open the input data model
        with datamodels.open(step_input) as input_model:

            # Get the name of the saturation reference file
            self.ref_name = self.get_reference_file(input_model, 'saturation')
            self.log.info('Using SATURATION reference file %s', self.ref_name)

            # Check for a valid reference file
            if self.ref_name == 'N/A':
                self.log.warning('No SATURATION reference file found')
                self.log.warning('Saturation step will be skipped')
                input_model.meta.cal_step.saturation = 'SKIPPED'
                return input_model

            # Open the reference file data model
            ref_model = datamodels.SaturationModel(self.ref_name)

            # Work on a copy
            result = input_model.copy()

            # Do the saturation check
            if pipe_utils.is_irs2(result):
                result = saturation.irs2_flag_saturation(result, ref_model, self.n_pix_grow_sat, self.use_readpatt)
            else:
                result = saturation.flag_saturation(result, ref_model, self.n_pix_grow_sat, self.use_readpatt)
            result.meta.cal_step.saturation = 'COMPLETE'

            # Cleanup
            del ref_model

        return result

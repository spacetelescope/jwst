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
    """

    reference_file_types = ['saturation']

    def process(self, input):

        # Open the input data model
        with datamodels.RampModel(input) as input_model:

            # Get the name of the saturation reference file
            self.ref_name = self.get_reference_file(input_model, 'saturation')
            self.log.info('Using SATURATION reference file %s', self.ref_name)

            # Check for a valid reference file
            if self.ref_name == 'N/A':
                self.log.warning('No SATURATION reference file found')
                self.log.warning('Saturation step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.saturation = 'SKIPPED'
                return result

            # Open the reference file data model
            ref_model = datamodels.SaturationModel(self.ref_name)

            # Do the saturation check
            if pipe_utils.is_irs2(input_model):
                sat = saturation.irs2_flag_saturation(input_model, ref_model, self.n_pix_grow_sat)
            else:
                sat = saturation.flag_saturation(input_model, ref_model, self.n_pix_grow_sat)

            # Close the reference file and update the step status
            ref_model.close()
            sat.meta.cal_step.saturation = 'COMPLETE'

        return sat

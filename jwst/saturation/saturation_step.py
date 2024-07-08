#! /usr/bin/env python
import gc
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from ..lib import pipe_utils
from . import saturation
from jwst.lib.basic_utils import use_datamodel, copy_datamodel


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

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model, model_class=datamodels.RampModel)

        result, input_model = copy_datamodel(input_model, self.parent)

        # Get the name of the saturation reference file
        self.ref_name = self.get_reference_file(result, 'saturation')
        self.log.info('Using SATURATION reference file %s', self.ref_name)

        # Check for a valid reference file
        if self.ref_name == 'N/A':
            self.log.warning('No SATURATION reference file found')
            self.log.warning('Saturation step will be skipped')
            result.meta.cal_step.saturation = 'SKIPPED'
            gc.collect()
            return result

        # Open the reference file data model
        ref_model = datamodels.SaturationModel(self.ref_name)

        # Do the saturation check
        if pipe_utils.is_irs2(result):
            sat = saturation.irs2_flag_saturation(result, ref_model, self.n_pix_grow_sat)
        else:
            sat = saturation.flag_saturation(result, ref_model, self.n_pix_grow_sat)

        # Close the reference file and update the step status
        del ref_model
        sat.meta.cal_step.saturation = 'COMPLETE'

        del result
        gc.collect()
        return sat

#! /usr/bin/env python
from ..stpipe import Step
from .. import datamodels
from . import fringe


__all__ = ["FringeStep"]


class FringeStep(Step):
    """
    FringeStep: Apply fringe correction to a science image using a fringe
    reference image.
    """

    reference_file_types = ['fringe']

    def process(self, input):
        with datamodels.open(input) as input_model:

            # Open the reference file
            self.fringe_filename = self.get_reference_file(input_model,
                                                           'fringe')
            self.log.info('Using FRINGE reference file: %s',
                          self.fringe_filename)

            # Check for a valid reference file
            if self.fringe_filename == 'N/A':
                self.log.warning('No FRINGE reference file found')
                self.log.warning('Fringe step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.fringe = 'SKIPPED'
                return result

            # Load the fringe reference file
            fringe_model = datamodels.FringeModel(self.fringe_filename)

            # Do the correction
            output_model = fringe.do_correction(input_model, fringe_model)

            fringe_model.close()

        return output_model

from __future__ import (division, print_function, unicode_literals,
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import resample

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResampleStep(Step):
    """
    ResampleStep: Uses the drizzle process to resample (geometric correction)
    a single 2D image or resample and combine a set of 2D images specified as
    an association.

    Parameters
    -----------
    input : str or model
        Single filename for either a single image or an association table.
    """

    spec = """
        single = boolean(default=False)
        wht_type = option('exptime','error',None,default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        good_bits = integer(default=4)
        blendheaders = boolean(default=True)
    """
    reference_file_types = ['drizpars']

    def process(self, input):

        input_models = datamodels.open(input)
        # If single input, insert into a ModelContainer
        if input_models.__class__ is not datamodels.ModelContainer:
            s = datamodels.ModelContainer()
            s.append(input_models)
            s.assign_group_ids()
            s.meta.resample.output = s.group_names[0].replace('_resamp', '_drz')

            self.input_models = s
        else:
            self.input_models = input_models

        # identify what reference file has been associated with these inputs
        try:
            self.ref_filename = self.get_reference_file(self.input_models[0],
                self.reference_file_types[0])
            log.info('{} reffile is {}'.format(self.reference_file_types[0],
                self.ref_filename))
        except:
            log.error('{} reffile is not found.'.format(
                self.reference_file_types[0]))

        # Call the resampling routine
        resamp = resample.ResampleData(self.input_models,
            ref_filename=self.ref_filename,
            single=self.single, wht_type=self.wht_type, pixfrac=self.pixfrac,
            kernel=self.kernel, fillval=self.fillval, good_bits=self.good_bits,
            blendheaders=self.blendheaders)
        resamp.do_drizzle()

        if len(resamp.output_models) == 1:
            output_model = resamp.output_models[0]
        else:
            output_model = resamp.output_models


        output_model.meta.cal_step.resample = "COMPLETE"

        return output_model


if __name__ == '__main__':
    cmdline.step_script(ResampleStep)

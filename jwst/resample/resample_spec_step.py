from __future__ import (division, print_function, unicode_literals, 
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import resample_spec
from ..exp_to_source import multislit_to_container

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResampleSpecStep(Step):
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
        wht_type = option('exptime', 'error', None, default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        good_bits = integer(default=4)
    """
    reference_file_types = ['drizpars']

    def process(self, input):
        
        with datamodels.open(input) as self.input_models:
            
            # Single input model, single resample output
            if not isinstance(self.input_models, datamodels.ModelContainer):
                s = datamodels.ModelContainer()
                s.append(self.input_models)
                self.input_models = s

            self.driz_filename = self.get_reference_file(self.input_models[0], 'drizpars')
            
            # Multislits get converted to a ModelContainer per slit
            if all([isinstance(i, datamodels.MultiSlitModel) for i in self.input_models]):
                log.info('Converting MultiSlit to ModelContainer')
                container_dict = multislit_to_container(self.input_models)
                output_product = datamodels.MultiProductModel()
                output_product.update(self.input_models[0])
                for k, v in container_dict.items():
                    self.input_models = v

                    # Set up the resampling object as part of this step
                    self.step = resample_spec.ResampleSpecData(self.input_models,
                        self.driz_filename, single=self.single,
                        wht_type=self.wht_type, pixfrac=self.pixfrac,
                        kernel=self.kernel, fillval=self.fillval,
                        good_bits=self.good_bits)
                    # Do the resampling
                    self.step.do_drizzle()
                    if len(self.step.output_models) == 1:
                        out_slit = self.step.output_models[0]
                        output_product.products.append(out_slit)
                    else:
                        out_slit = self.step.output_models
                result = output_product
            else:
                # Set up the resampling object as part of this step
                self.step = resample_spec.ResampleSpecData(self.input_models,
                    self.driz_filename, single=self.single, wht_type=self.wht_type,
                    pixfrac=self.pixfrac, kernel=self.kernel,
                    fillval=self.fillval, good_bits=self.good_bits)
                # Do the resampling
                self.step.do_drizzle()

                # Return either the single resampled datamodel, or the container
                # of datamodels.
                if len(self.step.output_models) == 1:
                    result = self.step.output_models[0]
                else:
                    result = self.step.output_models
        
        result.meta.cal_step.resample = 'COMPLETE'
        return result


if __name__ == '__main__':
    cmdline.step_script(ResampleSpecStep)

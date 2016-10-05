#! /usr/bin/env python

from __future__ import division

from ..stpipe import Step, cmdline
from .. import datamodels
from . import ramp_fit

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class RampFitStep (Step):

    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    spec = """
        int_name = string(default='')
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        algorithm = option('OLS', 'GLS', default='OLS') # 'OLS' or 'GLS'
        weighting = option('unweighted', 'optimal', default='unweighted') \
        # 'unweighted' or 'optimal'
    """

    reference_file_types = ['readnoise', 'gain']


    def process(self, input):

        with datamodels.open(input) as input_model:

            readnoise_filename = self.get_reference_file(input_model,
                                                          'readnoise')
            gain_filename = self.get_reference_file(input_model,
                                                     'gain')

            log.info('Using READNOISE reference file: %s', readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)
            log.info('Using GAIN reference file: %s', gain_filename)
            gain_model = datamodels.GainModel(gain_filename)

            log.info('Using algorithm = %s' % self.algorithm)
            log.info('Using weighting = %s' % self.weighting)

            buffsize = ramp_fit.BUFSIZE
            if self.algorithm == "GLS":
                buffsize //= 10
            out_model, int_model, opt_model, gls_opt_model = \
                        ramp_fit.ramp_fit(input_model,
                                           buffsize, self.save_opt,
                                           readnoise_model, gain_model,
                                           self.algorithm, self.weighting)

            readnoise_model.close()

        if int_model is not None:
            if self.int_name != '':
                int_model.save(self.int_name)
            else:
                self.save_model(int_model, 'rateints')
        if opt_model is not None:
            if self.opt_name != '':
                opt_model.save(self.opt_name)
            else:
                self.save_model(opt_model, 'fitopt')
        if gls_opt_model is not None:
            if self.opt_name != '':
                gls_opt_model.save(self.opt_name)
            else:
                self.save_model(gls_opt_model, 'fitoptgls')

        out_model.meta.cal_step.ramp_fit = 'COMPLETE'

        return out_model

if __name__ == '__main__':
    cmdline.step_script(ramp_fit_step)

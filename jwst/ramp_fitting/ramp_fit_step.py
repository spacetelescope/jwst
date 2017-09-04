#! /usr/bin/env python

from __future__ import division

from ..stpipe import Step, cmdline
from .. import datamodels
from ..gain_scale import gain_scale
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
    """

    # Prior to 04/26/17, the following were also in the spec above:
    #      algorithm = option('OLS', 'GLS', default='OLS') # 'OLS' or 'GLS'
    #      weighting = option('unweighted', 'optimal', default='unweighted') \
    #      # 'unweighted' or 'optimal'
    # As of 04/26/17, the only allowed algorithm is 'ols', and the 
    #      only allowed weighting is 'optimal'.
    algorithm = 'ols'      # Only algorithm allowed for Build 7.1
    weighting = 'optimal'  # Only weighting allowed for Build 7.1

    reference_file_types = ['readnoise', 'gain']

    def process(self, input):

        with datamodels.open(input) as input_model:

            readnoise_filename = self.get_reference_file(input_model,
                                                          'readnoise')
            gain_filename = self.get_reference_file(input_model, 'gain')

            log.info('Using READNOISE reference file: %s', readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)
            log.info('Using GAIN reference file: %s', gain_filename)
            gain_model = datamodels.GainModel(gain_filename)

            # Try to retrieve the gain factor from the gain reference file
            if gain_model.meta.gain_factor is not None:
                gain_factor = gain_model.meta.gain_factor
                input_model.meta.exposure.gain_factor = gain_factor
            else:
                self.log.warning('GAINFACT not found in gain reference file')
                input_model.meta.exposure.gain_factor = None
                gain_factor = None

            log.info('Using algorithm = %s' % self.algorithm)
            log.info('Using weighting = %s' % self.weighting)

            buffsize = ramp_fit.BUFSIZE
            if self.algorithm == "GLS":
                buffsize //= 10

            out_model, int_model, opt_model, gls_opt_model, \
               var_slope_r, var_slope_p, var_r_s_4d, var_p_s_4d = \
                            ramp_fit.ramp_fit(input_model,
                                               buffsize, self.save_opt,
                                               readnoise_model, gain_model,
                                               self.algorithm, self.weighting) 

            out_model.instance['extra_fits']['PoissonNoise'] = {'data': var_p_s_4d} 
            out_model.instance['extra_fits']['ReadNoise'] = {'data': var_r_s_4d}  

            readnoise_model.close()

        # Save the multi-integration product, if it exists
        if int_model is not None:

            # Apply the gain scale to the multi-integration product
            if gain_factor is not None:
                int_model = gain_scale.do_correction(int_model, gain_factor)
                int_model.meta.cal_step.gain_scale = 'COMPLETE'
            else:
                int_model.meta.cal_step.gain_scale = 'SKIPPED'

            # Save the model
            if self.int_name != '':
                int_model.save(self.int_name)
            else:
                self.save_model(int_model, 'rateints')

        # Save the OLS optional fit product, if it exists
        if self.save_opt:
            if self.opt_name != '':
                opt_model.save(self.opt_name)
            else:
                self.save_model(opt_model, 'fitopt')

        # Save the GLS optional fit product, if it exists
        if gls_opt_model is not None:
            if self.opt_name != '':
                gls_opt_model.save(self.opt_name)
            else:
                self.save_model(gls_opt_model, 'fitoptgls')

        out_model.meta.cal_step.ramp_fit = 'COMPLETE'

        return out_model

if __name__ == '__main__':
    cmdline.step_script(ramp_fit_step)

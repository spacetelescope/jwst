#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import ramp_fit

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["RampFitStep"]


class RampFitStep (Step):

    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    spec = """
        int_name = string(default='')
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
    """

    # Prior to 04/26/17, the following were also in the spec above:
    #      algorithm = option('OLS', 'GLS', default='OLS') # 'OLS' or 'GLS'
    #      weighting = option('unweighted', 'optimal', default='unweighted') \
    #      # 'unweighted' or 'optimal'
    # As of 04/26/17, the only allowed algorithm is 'ols', and the
    #      only allowed weighting is 'optimal'.

    algorithm = 'ols'      # Only algorithm allowed for Build 7.1
#    algorithm = 'gls'       # 032520

    weighting = 'optimal'  # Only weighting allowed for Build 7.1

    reference_file_types = ['readnoise', 'gain']

    def process(self, input):

        with datamodels.RampModel(input) as input_model:
            max_cores = self.maximum_cores
            readnoise_filename = self.get_reference_file(input_model, 'readnoise')
            gain_filename = self.get_reference_file(input_model, 'gain')

            log.info('Using READNOISE reference file: %s', readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)
            log.info('Using GAIN reference file: %s', gain_filename)
            gain_model = datamodels.GainModel(gain_filename)

            # Try to retrieve the gain factor from the gain reference file.
            # If found, store it in the science model meta data, so that it's
            # available later in the gain_scale step, which avoids having to
            # load the gain ref file again in that step.
            if gain_model.meta.exposure.gain_factor is not None:
                input_model.meta.exposure.gain_factor = \
                    gain_model.meta.exposure.gain_factor

            log.info('Using algorithm = %s' % self.algorithm)
            log.info('Using weighting = %s' % self.weighting)

            buffsize = ramp_fit.BUFSIZE
            if self.algorithm == "GLS":
                buffsize //= 10

            out_model, int_model, opt_model, gls_opt_model = ramp_fit.ramp_fit(
                input_model, buffsize,
                self.save_opt, readnoise_model, gain_model, self.algorithm,
                self.weighting, max_cores
            )

            readnoise_model.close()
            gain_model.close()

        # Save the OLS optional fit product, if it exists
        if opt_model is not None:
            self.save_model(opt_model, 'fitopt', output_file=self.opt_name)

        # Save the GLS optional fit product, if it exists
        if gls_opt_model is not None:
            self.save_model(
                gls_opt_model, 'fitoptgls', output_file=self.opt_name
            )

        if out_model is not None:
            out_model.meta.cal_step.ramp_fit = 'COMPLETE'
            if (input_model.meta.exposure.type in ['NRS_IFU', 'MIR_MRS']) or (
                input_model.meta.exposure.type in ['NRS_AUTOWAVE', 'NRS_LAMP'] and
                    input_model.meta.instrument.lamp_mode == 'IFU'):

                out_model = datamodels.IFUImageModel(out_model)

        if int_model is not None:
            int_model.meta.cal_step.ramp_fit = 'COMPLETE'

        return out_model, int_model

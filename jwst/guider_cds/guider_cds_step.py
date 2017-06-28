#! /usr/bin/env python

from __future__ import division

from ..stpipe import Step, cmdline
from .. import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class GuiderCdsStep (Step):

    """

    """

    spec = """

    """
    def process(self, input):

        with datamodels.open(input) as input_model:

            buffsize = ramp_fit.BUFSIZE
          #  out_model, int_model, opt_model, gls_opt_model = \
          #              ramp_fit.ramp_fit(input_model,
          #                                 buffsize, self.save_opt,
          #                                 readnoise_model, gain_model,
          #                                 self.algorithm, self.weighting)

        out_model.meta.cal_step.guider_cds = 'COMPLETE'

        return out_model

if __name__ == '__main__':
    cmdline.step_script(guider_cds_step)
    

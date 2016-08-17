#! /usr/bin/env python

import os

from ..stpipe import Step, cmdline
from .. import datamodels

from . import klip

class KlipStep(Step):

    """
    KlipStep: Performs KLIP processing on a science target coronagraphic
    exposure. The input science exposure is assumed to be a fully calibrated
    level-2b image. The processing is performed using a set of reference PSF
    images observed in the same coronagraphic mode.
    """

    spec = """
        truncate = integer(default=50,min=0) # The number of KL transform rows to keep
    """

    def process(self, target, psfrefs):

        with datamodels.CubeModel(target) as target_model:

            # Retrieve the parameter values
            truncate = self.truncate
            self.log.info('KL transform truncation = %d', truncate)

            # Get the PSF reference images
            refs_model = datamodels.QuadModel(psfrefs)

            # Call the KLIP routine
            (target, psf) = klip.klip(target_model, refs_model, truncate)

            refs_model.close()

        result.meta.cal_step.klip = 'COMPLETE'

        return (target, psf)


if __name__ == '__main__':
    cmdline.step_script(KlipStep)

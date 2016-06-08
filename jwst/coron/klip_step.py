#! /usr/bin/env python

import os

from jwst.stpipe import Step, cmdline
from jwst import datamodels

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

        with models.ImageModel(target) as target_model:

            # Retrieve the parameter values
            truncate = self.truncate
            self.log.info('KL transform truncation = %d', truncate)

            # Get the PSF reference images
            refs_model = models.CubeModel(psfrefs)

            # Call the KLIP routine
            (target, psf) = klip.klip(target_model, refs_model, truncate)

            refs_model.close()

        # Save the output target and psf images
        root = os.path.abspath(os.path.splitext(target.meta.filename)[0])
        target_file = root + "_klip.fits"
        psf_file = root + "_psf.fits"
        target.save(target_file)
        target.close()
        psf.save(psf_file)
        psf.close()

        #result.meta.cal_step.klip = 'COMPLETE'

        return


if __name__ == '__main__':
    cmdline.step_script(KlipStep)

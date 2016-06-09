#! /usr/bin/env python

import os

from ..stpipe import Step, cmdline
from .. import datamodels

from . import hlsp

class HlspStep(Step):

    """
    HlspStep: Make High-Level Science Products (HLSP's) from the results of
    coronagraphic exposure that's had KLIP processing applied to it.
    """

    spec = """
        annuli_width = integer(default=2, min=1) # Width of contrast annuli
    """

    def process(self, target):

        width = self.annuli_width

        # Open the input target image model
        with datamodels.ImageModel(target) as target_model:

            # Create a signal-to-noise ratio image
            self.log.info ('Creating SNR image')
            snr = hlsp.snr_image(target_model)

            # Create a contrast curve
            self.log.info ('Creating contrast curve')
            contrast = hlsp.contrast_curve(target_model, width)

        # Save the SNR output file
        root = os.path.abspath(os.path.splitext(snr.meta.filename)[0])
        snr.meta.cal_step.hlsp = 'COMPLETE'
        snr_file = root + "_snr.fits"
        snr.save(snr_file)
        snr.close()

        # Save the Contrast curve file
        contrast.meta.cal_step.hlsp = 'COMPLETE'
        contrast_file = root + "_contrast.fits"
        contrast.save(contrast_file)
        contrast.close()

        return


if __name__ == '__main__':
    cmdline.step_script(HlspStep)


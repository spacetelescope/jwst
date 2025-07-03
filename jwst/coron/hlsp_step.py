#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from . import hlsp
from jwst.stpipe import Step

__all__ = ["HlspStep"]


class HlspStep(Step):
    """Make High-Level Science Products (HLSP's) after KLIP processing has been applied."""

    class_alias = "hlsp"

    spec = """
        annuli_width = integer(default=2, min=1) # Width of contrast annuli
        save_results = boolean(default=true) # Save results
    """  # noqa: E501

    def process(self, target):
        """Execute the HLSP calibration step."""
        width = self.annuli_width

        # Open the input target image model
        with datamodels.ImageModel(target) as target_model:
            # Create a signal-to-noise ratio image
            self.log.info("Creating SNR image")
            snr = hlsp.snr_image(target_model)

            # Create a contrast curve
            self.log.info("Creating contrast curve")
            contrast = hlsp.contrast_curve(target_model, width)

        # Save the SNR output file
        if self.output_file is None:
            self.output_file = target_model.meta.filename
        snr.meta.cal_step.hlsp = "COMPLETE"
        self.save_model(snr, suffix="snr")
        snr.close()

        # Save the Contrast curve file
        contrast.meta.cal_step.hlsp = "COMPLETE"
        self.save_model(contrast, "contrast")
        contrast.close()

        return

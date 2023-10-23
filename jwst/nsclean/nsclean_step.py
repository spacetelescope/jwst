from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import nsclean

__all__ = ["NSCleanStep"]


class NSCleanStep(Step):
    """
    NSCleanStep: This step performs 1/f noise correction ("cleaning")
    of NIRSpec images, using the "NSClean" method.
    """

    class_alias = "nsclean"

    spec = """
        n_sigma = float(default=5.0)
        save_mask = boolean(default=False)
    """

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # For now, only allow operations on IFU and MOS mode.
            # Skip for fixed-slit and BOTS.
            if input_model.meta.exposure.type.lower() not in ['nrs_ifu', 'nrs_msaspec']:
                self.log.warning("Step can only be applied to IFU and MOS images.")
                self.log.warning("Step will be skipped.")
                output_model = input_model.copy()
                output_model.meta.cal_step.nsclean = 'SKIPPED'
                return output_model

            # Do the NSClean correction
            result = nsclean.do_correction(input_model, self.n_sigma, self.save_mask)
            output_model, mask_model = result

            # Save the mask, if requested
            if self.save_mask:
                mask_path = self.make_output_path(basepath=input_model.meta.filename, suffix='mask')
                self.log.info(f"Saving mask file {mask_path}")
                mask_model.save(mask_path)
                mask_model.close()

        return output_model

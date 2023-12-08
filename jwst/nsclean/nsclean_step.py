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
        mask_spectral_regions = boolean(default=True)  # Mask WCS-defined regions
        n_sigma = float(default=5.0)  # Clipping level for outliers
        save_mask = boolean(default=False)  # Save the created mask
        user_mask = string(default=None)  # Path to user-supplied mask
        skip = boolean(default=True)  # By default, skip the step
    """

    def process(self, input):
        """
        Fit and subtract 1/f background noise from a NIRSpec image

        Parameters
        ----------
        input : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`
            Input datamodel to be corrected

        n_sigma : float, optional
            Sigma clipping threshold to be used in detecting outliers in the image

        save_mask : bool, optional
            Save the computed mask image

        user_mask : None, str, or `~jwst.datamodels.ImageModel`
            Optional user-supplied mask image; path to file or opened datamodel

        mask_spectral_regions : bool, optional
            Mask regions of the image defined by WCS bounding boxes for slits/slices

        Returns
        -------
        output_model : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`
            The 1/f corrected datamodel
        """

        # Open the input data model
        with datamodels.open(input) as input_model:

            # Do the NSClean correction
            result = nsclean.do_correction(input_model, self.mask_spectral_regions, self.n_sigma,
                                           self.save_mask, self.user_mask)
            output_model, mask_model = result

            # Save the mask, if requested
            if self.save_mask and mask_model is not None:
                mask_path = self.make_output_path(basepath=input_model.meta.filename, suffix='mask')
                self.log.info(f"Saving mask file {mask_path}")
                mask_model.save(mask_path)
                mask_model.close()

        return output_model

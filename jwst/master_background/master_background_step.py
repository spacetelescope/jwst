from ..stpipe import Step
from .. import datamodels


__all__ = ["MasterBackgroundStep"]


class MasterBackgroundStep(Step):
    """
    MasterBackgroundStep:  Compute and subtract master background from spectra
    """

    spec = """
        user_background = str(default=None) # Path to user-supplied master background
        save_background = boolean(default=False) # Save computed master background
    """


    def process(self, input):

        """
        Compute a master background from N > 1 dedicated background exposures
        and subtract it from the target science exposure(s)

        Parameters
        ----------
        input: JWST data model, ModelContainer or association
            input target data model(s) to which background subtraction is applied

        user_background: string or SpecModel
            Path to user-supplied master background 1D spectrum, or SpecModel

        Returns
        -------
        result: JWST data model
            the background-subtracted target data model(s)
        """

        # Load the input data model
        with datamodels.open(input) as input_data:

            # Handle individual NIRSpec FS, NIRSpec MOS
            if isinstance(input_data, datamodels.MultiSlitModel):
                pass

            # Handle associations, or input ModelContainers
            elif isinstance(input_data, datamodels.ModelContainer):
                pass

            # Handle MIRI LRS
            elif isinstance(input_data, datamodels.ImageModel):
                pass

            # Handle MIRI MRS and NIRSpec IFU
            elif isinstance(input_data, datamodels.IFUImageModel):
                pass

            else:
                result = input_data
                self.log.warning(
                    "Input %s cannot be handled.  Doing nothing.", input
                    )

            # Check if user has supplied a master background spectrum.
            if self.user_background is None:
                pass
            else:
                pass

            # Return input as dummy result for now
            result = input_data
            result.meta.cal_step.master_back_sub = 'COMPLETE'

        return result

from ..stpipe import Step
from .. import datamodels


__all__ = ["MasterBackgroundStep"]


class MasterBackgroundStep(Step):
    """
    MasterBackgroundStep:  Compute and subtract master background from spectra
    """

    spec = """
        user_background = string(default=None) # Path to user-supplied master background
        save_background = boolean(default=False) # Save computed master background
    """


    def process(self, input):
        """
        Compute and subtract a master background spectrum

        Parameters
        ----------
        input : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.ModelContainer`, association
            Input target data model(s) to which master background subtraction is
            to be applied

        user_background : None, string, or `~jwst.datamodels.MultiSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel

        save_background : bool, optional
            Save master background.

        Returns
        -------
        result: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.ModelContainer`
            The background-subtracted target data model(s)
        """

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
                    "Input %s of type %s cannot be handled.  Step skipped.",
                    input, type(input)
                    )
                try:
                    result.meta.cal_step.master_back_sub = 'SKIPPED'
                except AttributeError:
                    for model in result:
                        model.meta.cal_step.master_back_sub = 'SKIPPED'

                return result

            # Check if user has supplied a master background spectrum.
            if self.user_background is None:
                pass
            else:
                pass

            # Save the computed background if requested by user
            if self.save_background and self.user_background is None:
                pass
                # TODO, figure out step suffix, use stpipe tools to generate
                # the filename from the input filename(s) and/or association.
                # Notice this is a custom `output_file` type in the spec above

            # Return input as dummy result for now
            result = input_data
            
            try:
                result.meta.cal_step.master_back_sub = 'COMPLETE'
            except AttributeError:
                for model in result:
                    model.meta.cal_step.master_back_sub = 'COMPLETE'

        return result

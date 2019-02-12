from os.path import basename

from ..stpipe import Step
from .. import datamodels

from .expand_to_2d import expand_to_2d


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

        # Get association info if available
        # asn_id = ???

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
                result = input_data.copy()
                self.log.warning(
                    "Input %s of type %s cannot be handled.  Step skipped.",
                    input, type(input)
                    )
                self.record_step_status(result, 'master_back_sub', success=False)

                return result

            # Check if user has supplied a master background spectrum.
            if self.user_background is None:
                # TODO: 1. compute master background from asn, 2. subtract it
                # Return input as dummy result for now
                result = input_data.copy()
            else:
                background_2d = expand_to_2d(input_data, self.user_background)
                result = subtract_2d_background(input_data, background_2d)

                # Record name of user-supplied master background spectrum
                if isinstance(result, datamodels.ModelContainer):
                    for model in result:
                        model.meta.master_background = basename(self.user_background)
                else:
                    result.meta.master_background = basename(self.user_background)

            # Save the computed background if requested by user
            if self.save_background and self.user_background is None:
                # self.save_model(background, suffix='masterbg', asn_id=asn_id)
                pass
            
            self.record_step_status(result, 'master_back_sub')

        return result


def subtract_2d_background(source, background):
    """Subtract a 2D background

    Parameters
    ----------
    source : `~jwst.datamodels.DataModel` or `~jwst.datamodels.ModelContainer`
        The input science data.

    background : `~jwst.datamodels.DataModel`
        The input background data.  Must be the same datamodel type as `source`.
        For a `~jwst.datamodels.ModelContainer`, the source and background
        models in the input containers must match one-to-one.

    Returns
    -------
    `~jwst.datamodels.DataModel`
        Background subtracted from source.
    """
    def _subtract_2d_background(model, background):
        result = model.copy()
        # Handle individual NIRSpec FS, NIRSpec MOS
        if isinstance(model, datamodels.MultiSlitModel):
            for slit, slitbg in zip(result.slits, background.slits):
                slit.data -= slitbg.data

        # Handle MIRI LRS, MIRI MRS and NIRSpec IFU
        elif isinstance(model, (datamodels.ImageModel, datamodels.IFUImageModel)):
            result.data -= background.data

        else:
            # Shouldn't get here.
            raise RuntimeError("Input type {} is not supported."
                               .format(type(model)))
        return result

    # Handle containers of many datamodels
    if isinstance(source, datamodels.ModelContainer):
        result = datamodels.ModelContainer()
        result.update(source)
        for model, bg in zip(source, background):
            result.append(_subtract_2d_background(model, bg))

    # Handle single datamodels
    elif isinstance(source, (datamodels.ImageModel, datamodels.IFUImageModel, datamodels.MultiSlitModel)):
        result = _subtract_2d_background(source, background)

    else:
        # Shouldn't get here.
        raise RuntimeError("Input type {} is not supported."
                           .format(type(source)))

    return result

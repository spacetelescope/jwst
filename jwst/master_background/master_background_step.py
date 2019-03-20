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
        force_subtract = boolean(default=False) # Force subtracting master background
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

        force_subtract : bool, optional
            Optional user-supplied flag which subtracts the master background overriding the checks
            on whether the background in calspec2 has already been subtracted.
            Default is set to false. The step logic determines if the master background should be
            subtracted. If set to true then the step logic is bypassed and the master background is
            subtracted.

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
                self.record_step_status(result, 'master_background', success=False)

                return result
            # Check if background was subtracted in calspec 2 and force_subtract is False
            # first check if input is a model container

            do_sub = True  # flag if set to False then no master background subtraction is done
            if not self.force_subtract:  # default mode

                # check if the input data is a model container. If it is then loop over
                # container and see if the background was subtracted in calspec2.
                # If all data was  background subtracted. Print log.info and skip master bgk subtrction
                # If there is a mixture of some being background subtracted print, warning
                # and message to user to use force_subtract to force the master background
                # to be subtracted.
                if isinstance(input_data, datamodels.ModelContainer):
                    isub = 0
                    for indata in input_data:
                        if indata.meta.cal_step.back_sub == 'COMPLETE':
                            do_sub = False
                            isub += 1

                    if not do_sub and isub == len(input_data):
                        self.log.info(
                            "Not subtracting master background, background was subtracted in calspec2")
                        self.log.info("To force the master background to be subtracted from this data, "
                            "run again and set force_subtract = True.")

                    if not do_sub and isub != len(input_data):
                        self.log.warning("Not subtracting master background.")
                        self.log.warning("Input data contains a mixture of data with and without "
                            "background subtraction done in calspec2.")
                        self.log.warning("To force the master background to be subtracted from this data, "
                            "run again and set force_subtract = True.")
                # input data is a single file
                else:
                    if input_data.meta.cal_step.back_sub == 'COMPLETE':
                        do_sub = False
                        self.log.info(
                            "Not subtracting master background, background was subtracted in calspec2")
                        self.log.info("To force the master background to be subtracted from this data, "
                            "run again and set force_subtract = True.")

            # checked all the input data - do we subtract (continue) or not (return)
            if not do_sub:
                result = input_data.copy()
                self.record_step_status(result, 'master_background', success=False)
                return result

            # various tests have passed and now we want to subtract the master background
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
                        model.meta.background.master_background_file = basename(self.user_background)
                else:
                    result.meta.background.master_background_file = basename(self.user_background)

            # Save the computed background if requested by user
            if self.save_background and self.user_background is None:
                # self.save_model(background, suffix='masterbg', asn_id=asn_id)
                pass

            self.record_step_status(result, 'master_background', success=True)

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

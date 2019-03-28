from os.path import basename
from ..stpipe import Step
from .. import datamodels
from jwst.combine_1d.combine1d import combine_1d_spectra

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
        output_use_model = boolean(default=True)
    """

    def process(self, input):
        """
        Compute and subtract a master background spectrum

        Parameters
        ----------
        input : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.ModelContainer`, association
            Input target datamodel(s) to which master background subtraction is
            to be applied

        user_background : None, string, or `~jwst.datamodels.MultiSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel

        save_background : bool, optional
            Save computed master background.

        force_subtract : bool, optional
            Optional user-supplied flag that overrides step logic to force subtraction of the
            master background.
            Default is False, in which case the step logic determines if the calspec2 background step
            has already been applied and, if so, the master background step is skipped.
            If set to True, the step logic is bypassed and the master background is subtracted.

        Returns
        -------
        result: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.ModelContainer`
            The background-subtracted science datamodel(s)
        """

        # Get association info if available
        # asn_id = ???

        with datamodels.open(input) as input_data:
            background_data = None

            # Handle individual NIRSpec FS, NIRSpec MOS
            if isinstance(input_data, datamodels.MultiSlitModel):
                pass

            # Handle associations, or input ModelContainers
            elif isinstance(input_data, datamodels.ModelContainer):
                input_data, background_data = split_container(input_data)
                asn_id = input_data.meta.asn_table.asn_id

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

            # Check if background was subtracted in calspec2 and force_subtract is False
            do_sub = True  # If False, then no master background subtraction is done
            if not self.force_subtract:  # default mode

                # Check if the input data is a model container. If it is, then loop over the
                # container and see if the background was subtracted in calspec2.
                # If all data was background subtracted, print log info and skip master bkg subtraction.
                # If there is a mixture of some being background subtracted, print warning
                # message to user to use force_subtract to force the master background
                # to be subtracted.
                if isinstance(input_data, datamodels.ModelContainer):
                    isub = 0
                    for indata in input_data:
                        if indata.meta.cal_step.back_sub == 'COMPLETE':
                            do_sub = False
                            isub += 1

                    if not do_sub and isub == len(input_data):
                        self.log.info(
                            "Not subtracting master background, because background was subtracted in calspec2")
                        self.log.info("To force master background subtraction for these data, "
                            "run again and set force_subtract = True.")

                    if not do_sub and isub != len(input_data):
                        self.log.warning("Not subtracting master background, because")
                        self.log.warning("input data contains a mixture of data with and without "
                            "background subtraction done in calspec2.")
                        self.log.info("To force master background subtraction for these data, "
                            "run again and set force_subtract = True.")

                # Input is a single model
                else:
                    if input_data.meta.cal_step.back_sub == 'COMPLETE':
                        do_sub = False
                        self.log.info(
                            "Not subtracting master background, because background was subtracted in calspec2")
                        self.log.info("To force master background subtraction for these data, "
                            "run again and set force_subtract = True.")

            # Checked all the input data - do we subtract (continue) or not (return)?
            if not do_sub:
                result = input_data.copy()
                self.record_step_status(result, 'master_background', success=False)
                return result

            # Compute master background and subtract it
            if self.user_background is None and background_data is not None:
                master_background = combine_1d_spectra(
                    background_data,
                    exptime_key='exposure_time',
                    background=True,
                )
                if isinstance(input_data, datamodels.ModelContainer):
                    input_data, background = split_container(input_data)
                    asn_id = input_data.meta.asn_table.asn_id

                    result = datamodels.ModelContainer()
                    for model in input_data:
                        background_2d = expand_to_2d(model, master_background)
                        result.append(subtract_2d_background(model, background_2d))
                else:
                    # Is it possible that we will not have a container and not
                    # have a list of background 1D spectra?  Probably not.
                    pass

            # Use user-supplied master background and subtract it
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
                self.save_model(master_background, suffix='masterbg', asn_id=asn_id)
            
            self.record_step_status(result, 'master_background', success=True)

        return result


def split_container(container):
    """Divide a ModelContainer with science and background into one of each
    """
    asn = container.meta.asn_table.instance
    background = datamodels.ModelContainer()
    science = datamodels.ModelContainer()
    for product in asn['products']:
        for member in product['members']:
            if member['exptype'] == 'science':
                science.append(datamodels.open(member['expname']))
            if member['exptype'] == 'background':
                background.append(datamodels.open(member['expname']))

    science.meta.asn_table = {}
    datamodels.model_base.properties.merge_tree(
        science.meta.asn_table._instance, asn
    )
    for p in science.meta.asn_table.instance['products']:
        p['members'] = [m for m in p['members'] if m['exptype'] != 'background']
    return science, background


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

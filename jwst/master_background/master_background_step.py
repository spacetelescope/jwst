from pathlib import Path

import numpy as np
from stdatamodels.properties import merge_tree
from stdatamodels.jwst import datamodels
from scipy.signal import medfilt

from jwst.datamodels import ModelContainer
from jwst.stpipe import record_step_status

from jwst.stpipe import Step
from jwst.combine_1d.combine1d import combine_1d_spectra
from .expand_to_2d import expand_to_2d

__all__ = ["MasterBackgroundStep"]


class MasterBackgroundStep(Step):
    """
    Compute and subtract master background from spectra.

    Attributes
    ----------
    median_kernel : int
        Optional user-supplied kernel with which to moving-median boxcar
        filter the master background spectrum.  Must be an odd integer; even integers
        will be rounded down to the nearest odd integer.
    user_background : None, str, or `~jwst.datamodels.MultiSpecModel`
        Optional user-supplied master background 1D spectrum, path to file
        or opened datamodel
    save_background : bool, optional
        Save computed master background.
    force_subtract : bool, optional
        Optional user-supplied flag that overrides step logic to force subtraction of the
        master background. Default is False, in which case the step logic determines
        if the calspec2 background step has already been applied and, if so, the master
        background step is skipped. If set to True, the step logic is bypassed and the
        master background is subtracted.
    """

    class_alias = "master_background"

    spec = """
        median_kernel = integer(default=1) # Moving-median boxcar size to filter the background
        user_background = string(default=None) # Path to user-supplied master background
        save_background = boolean(default=False) # Save computed master background
        force_subtract = boolean(default=False) # Force subtracting master background
        output_use_model = boolean(default=True)
    """  # noqa: E501

    def process(self, input_data):
        """
        Compute and subtract a master background spectrum.

        Parameters
        ----------
        input_data : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
                     `~jwst.datamodels.ModelContainer`, str
            Input target datamodel(s) or association file to which master background
            subtraction is to be applied.

        Returns
        -------
        result : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
                 `~jwst.datamodels.ModelContainer`
            The background-subtracted science datamodel(s)
        """
        with datamodels.open(input_data) as input_model:
            # Make the input data available to self
            self.input_model = input_model

            # First check if we should even do the subtraction.  If not, bail.
            if not self._do_sub:
                result = input_model.copy()
                record_step_status(result, "master_background", success=False)
                return result

            # Check that data is a supported datamodel. If not, bail.
            if not isinstance(
                input_model,
                (
                    ModelContainer,
                    datamodels.MultiSlitModel,
                    datamodels.ImageModel,
                    datamodels.IFUImageModel,
                ),
            ):
                result = input_model.copy()
                self.log.warning(
                    f"Input {input_data} of type {type(input_data)} cannot be handled.  "
                    f"Step skipped."
                )
                record_step_status(result, "master_background", success=False)
                return result

            # If user-supplied master background, subtract it
            if self.user_background:
                if isinstance(input_model, ModelContainer):
                    input_model, _ = split_container(input_model)
                    asn_id = input_model.asn_table["asn_id"]
                    del _
                    result = ModelContainer()
                    background_2d_collection = ModelContainer()
                    for model in input_model:
                        background_2d = expand_to_2d(model, self.user_background)
                        result.append(subtract_2d_background(model, background_2d))
                        background_2d_collection.append(background_2d)
                    # Record name of user-supplied master background spectrum
                    for model in result:
                        model.meta.background.master_background_file = Path(
                            self.user_background
                        ).name
                # Use user-supplied master background and subtract it
                else:
                    asn_id = None
                    background_2d = expand_to_2d(input_model, self.user_background)
                    background_2d_collection = ModelContainer([background_2d])
                    result = subtract_2d_background(input_model, background_2d)
                    # Record name of user-supplied master background spectrum
                    result.meta.background.master_background_file = Path(self.user_background).name

                # Save the computed 2d background if requested by user. The user has supplied
                # the master background so just save the expanded 2d background
                if self.save_background:
                    self.save_container(
                        background_2d_collection, suffix="masterbg2d", force=True, asn_id=asn_id
                    )

            # Compute master background and subtract it
            else:
                if isinstance(input_model, ModelContainer):
                    input_model, background_data = split_container(input_model)
                    if len(background_data) == 0:
                        msg = (
                            "No background data found in input container, "
                            "and no user-supplied background provided.  Skipping step."
                        )
                        self.log.warning(msg)
                        result = input_model.copy()
                        record_step_status(result, "master_background", success=False)
                        return result
                    asn_id = input_model.asn_table["asn_id"]

                    for model in background_data:
                        # Check if the background members are nodded x1d extractions
                        # or background from dedicated background exposures.
                        # Use "bkgdtarg == False" so we don't also get None cases
                        # for simulated data that didn't bother populating this
                        # keyword.
                        this_is_ifu_extended = False
                        if (
                            model.meta.exposure.type == "NRS_IFU"
                            and model.spec[0].source_type == "EXTENDED"
                        ):
                            this_is_ifu_extended = True
                        if model.meta.exposure.type == "MIR_MRS":
                            # always treat as extended for MIRI MRS
                            this_is_ifu_extended = True

                        if model.meta.observation.bkgdtarg is False or this_is_ifu_extended:
                            self.log.debug("Copying BACKGROUND column to SURF_BRIGHT")
                            copy_background_to_surf_bright(model)

                    master_background = combine_1d_spectra(
                        background_data,
                        exptime_key="exposure_time",
                    )

                    # If requested, apply a moving-median boxcar filter to the
                    # master background spectrum.
                    # Round down even kernel sizes because only odd kernel sizes are supported.
                    if self.median_kernel % 2 == 0:
                        self.median_kernel -= 1
                        self.log.info(
                            "Even median filter kernels are not supported."
                            f" Rounding the median kernel size down to {self.median_kernel}."
                        )

                    if self.median_kernel > 1:
                        self.log.info(
                            f"Applying moving-median boxcar of width {self.median_kernel}."
                        )
                        master_background.spec[0].spec_table["surf_bright"] = medfilt(
                            master_background.spec[0].spec_table["surf_bright"],
                            kernel_size=[self.median_kernel],
                        )
                        master_background.spec[0].spec_table["flux"] = medfilt(
                            master_background.spec[0].spec_table["flux"],
                            kernel_size=[self.median_kernel],
                        )

                    background_data.close()

                    result = ModelContainer()
                    background_2d_collection = ModelContainer()
                    for model in input_model:
                        background_2d = expand_to_2d(model, master_background)
                        result.append(subtract_2d_background(model, background_2d))
                        background_2d_collection.append(background_2d)

                    input_model.close()

                else:
                    result = input_model.copy()
                    input_model.close()
                    self.log.warning(
                        f"Input {input_data} of type {type(input_data)} cannot be "
                        "handled without user-supplied background.  Step skipped."
                    )
                    record_step_status(result, "master_background", success=False)
                    return result

                # Save the computed background if requested by user
                if self.save_background:
                    self.save_model(
                        master_background, suffix="masterbg1d", force=True, asn_id=asn_id
                    )
                    self.save_container(
                        background_2d_collection, suffix="masterbg2d", force=True, asn_id=asn_id
                    )

            record_step_status(result, "master_background", success=True)

        return result

    @property
    def _do_sub(self):
        """
        Decide if subtraction is to be done.

        Encapsulates logic that checks if background step has already been run
        on the data, or if the user has selected to force_subtract regardless.

        Returns
        -------
        do_sub : bool
            If ``True``, do the subtraction
        """
        do_sub = True
        if not self.force_subtract:
            input_model = self.input_model
            # check if the input data is a model container. If it is then loop over
            # container and see if the background was subtracted in calspec2.
            # If all data was background subtracted, skip master bgk subtraction.
            # If there is a mixture of some being background subtracted, don't
            # subtract and print warning message
            if isinstance(input_model, ModelContainer):
                isub = 0
                for indata in input_model:
                    if (
                        indata.meta.cal_step.bkg_subtract == "COMPLETE"
                        or indata.meta.cal_step.master_background == "COMPLETE"
                    ):
                        do_sub = False
                        isub += 1

                if not do_sub and isub == len(input_model):
                    self.log.info(
                        "Not subtracting master background, background was subtracted in calspec2"
                    )
                    self.log.info(
                        "To force the master background to be subtracted from this data, "
                        "run again and set force_subtract = True."
                    )

                if not do_sub and isub != len(input_model):
                    self.log.warning("Not subtracting master background.")
                    self.log.warning(
                        "Input data contains a mixture of data with and without "
                        "background subtraction done in calspec2."
                    )
                    self.log.warning(
                        "To force the master background to be subtracted from this data, "
                        "run again and set force_subtract = True."
                    )
            # input data is a single file
            else:
                if (
                    input_model.meta.cal_step.bkg_subtract == "COMPLETE"
                    or input_model.meta.cal_step.master_background == "COMPLETE"
                ):
                    do_sub = False
                    self.log.info(
                        "Not subtracting master background, background was subtracted in calspec2"
                    )
                    self.log.info(
                        "To force the master background to be subtracted from this data, "
                        "run again and set force_subtract = True."
                    )

        return do_sub

    def save_container(self, container, suffix="", asn_id="", force=True):
        """Save all models in container for intermediate background subtraction."""
        for i, model in enumerate(container):
            self.save_model(model, suffix=suffix, force=force, asn_id=asn_id, idx=i)


def copy_background_to_surf_bright(spectrum):
    """Copy the background column to the surf_bright column in a MultiSpecModel in-place."""
    for spec in spectrum.spec:
        spec.spec_table["SURF_BRIGHT"][:] = spec.spec_table["BACKGROUND"].copy()
        spec.spec_table["SB_ERROR"][:] = spec.spec_table["BKGD_ERROR"].copy()
        # Zero out the background column for safety
        spec.spec_table["BACKGROUND"][:] = 0
        # Set BERROR to dummy val of 0.0, as in extract_1d currently
        spec.spec_table["BKGD_ERROR"][:] = 0.0


def split_container(container):
    """
    Divide a ModelContainer with science and background into one of each.

    Parameters
    ----------
    container : ModelContainer
        Input model container

    Returns
    -------
    science : ModelContainer
        Container for science data.
    background : ModelContainer
        Container for background data.
    """
    background = ModelContainer()
    science = ModelContainer()

    for ind_science in container.ind_asn_type("science"):
        science.append(container[ind_science])

    for ind_bkgd in container.ind_asn_type("background"):
        background.append(container[ind_bkgd])

    # Pass along the association table to the output science container
    science.asn_pool_name = container.asn_pool_name
    science.asn_table_name = container.asn_table_name
    merge_tree(science.asn_table, container.asn_table)
    # Prune the background members from the table
    for p in science.asn_table["products"]:
        p["members"] = [m for m in p["members"] if m["exptype"].lower() != "background"]
    return science, background


def subtract_2d_background(source, background):
    """
    Subtract a 2D background.

    Parameters
    ----------
    source : `~jwst.datamodels.JwstDataModel` or `~jwst.datamodels.ModelContainer`
        The input science data.
    background : `~jwst.datamodels.JwstDataModel`
        The input background data.  Must be the same datamodel type as `source`.
        For a `~jwst.datamodels.ModelContainer`, the source and background
        models in the input containers must match one-to-one.

    Returns
    -------
    `~jwst.datamodels.JwstDataModel`
        Background subtracted from source.
    """

    def _subtract_2d_background(model, background):
        result = model.copy()
        # Handle individual NIRSpec FS, NIRSpec MOS
        if isinstance(model, datamodels.MultiSlitModel):
            for slit, slitbg in zip(result.slits, background.slits, strict=False):
                slit.data -= slitbg.data
                slit.dq = np.bitwise_or(slit.dq, slitbg.dq)

        # Handle MIRI LRS, MIRI MRS and NIRSpec IFU
        elif isinstance(model, (datamodels.ImageModel, datamodels.IFUImageModel)):
            result.data -= background.data
            result.dq = np.bitwise_or(result.dq, background.dq)
        else:
            # Shouldn't get here.
            raise TypeError(f"Input type {type(model)} is not supported.")
        return result

    # Handle containers of many datamodels
    if isinstance(source, ModelContainer):
        result = ModelContainer()
        result.update(source)
        for model, bg in zip(source, background, strict=False):
            result.append(_subtract_2d_background(model, bg))

    # Handle single datamodels
    elif isinstance(
        source, (datamodels.ImageModel, datamodels.IFUImageModel, datamodels.MultiSlitModel)
    ):
        result = _subtract_2d_background(source, background)

    else:
        # Shouldn't get here.
        raise TypeError(f"Input type {type(source)} is not supported.")

    return result

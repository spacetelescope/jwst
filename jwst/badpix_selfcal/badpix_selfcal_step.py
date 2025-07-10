import warnings
from jwst.stpipe import Step
from . import badpix_selfcal
import numpy as np
from jwst import datamodels as dm

__all__ = ["BadpixSelfcalStep"]

# Define the regions between MRS channels with which to estimate
# the pedestal dark.  Format is (Y1, Y2, X1, X2)
PEDESTAL_INDX_MIRIFULONG = (50, 974, 474, 507)
PEDESTAL_INDX_MIRIFUSHORT = (50, 974, 503, 516)


class BadpixSelfcalStep(Step):
    """
    Flag residual artifacts as bad pixels using a median filter and percentile cutoffs.

    All input exposures in the association file (or manually-provided bkg_list) are combined
    into a single background model using a MIN operation. The bad pixels are then identified
    using a median filter and percentile cutoffs, and applied to the science data by setting
    the flagged pixels, errors, and variances to NaN,
    and the DQ flag to DO_NOT_USE + OTHER_BAD_PIXEL.
    """

    class_alias = "badpix_selfcal"

    spec = """
    flagfrac_lower = float(default=0.001, min=0.0, max=0.5)  # fraction of pixels to flag on the low-flux end
    flagfrac_upper = float(default=0.001, min=0.0, max=0.5)  # fraction of pixels to flag on the high-flux end
    kernel_size = integer(default=15, min=1)  # size of kernel for median filter
    force_single = boolean(default=False)  # force single input exposure
    save_flagged_bkg = boolean(default=False)  # save flagged background exposures to file
    skip = boolean(default=True)
    """  # noqa: E501

    def save_model(self, model, *args, **kwargs):
        """
        Override save_model to suppress index 0 when save_model is True.

        Parameters
        ----------
        model : JWST data model
            Data model to save
        *args : tuple
            Additional arguments to pass to Step.save_model
        **kwargs : dict
            Additional keyword arguments to pass to Step.save_model

        Returns
        -------
        list[str]
            List of output paths for the saved models
        """
        kwargs["idx"] = None
        return Step.save_model(self, model, *args, **kwargs)

    def save_bkg(self, bkg_list, suffix="badpix_selfcal_bkg"):
        """
        Save the background exposures to file with correct indexing.

        Parameters
        ----------
        bkg_list : list of ImageModels
            Background exposures to save
        suffix : str
            Suffix to append to the filename
        """
        for i, bkg_model in enumerate(bkg_list):
            self.save_model(bkg_model, suffix=suffix + f"_{str(i)}")

    def process(self, input_data, selfcal_list=None, bkg_list=None):
        """
        Flag residual artifacts as bad pixels in the DQ array of a JWST exposure.

        Parameters
        ----------
        input_data : JWST data model or association
            Input science data to be corrected, or tuple of (sci, bkg, selfcal)
        selfcal_list : list of ImageModels or filenames, default None
            Exposures to include as part of median background model used to find bad pixels,
            but that are not flagged and returned as background exposures.
        bkg_list : list of ImageModels or filenames, default None
            Exposures to include as part of median background model used to find bad pixels,
            and that are flagged and returned as background exposures.

        Returns
        -------
        output : JWST data model or association
            Data model with CRs flagged

        Notes
        -----
        If an association file is read in, all exposures in the
        association file, including science, background, and selfcal exposures,
        are included in the MIN frame from which outliers are detected.
        If selfcal_list and/or bkg_list are specified manually, they are appended to any
        selfcal or background exposures found in the input association file.
        If selfcal_list and bkg_list are both set to None and input is
        a single science exposure, the step will be skipped with a warning unless
        the force_single parameter is set True.
        In that case, the input exposure will be used as the sole background exposure,
        i.e., true self-calibration.
        """
        input_sci, selfcal_list, bkg_list = _parse_inputs(input_data, selfcal_list, bkg_list)

        # ensure that there are background exposures to use, otherwise skip the step
        # unless forced
        if (len(selfcal_list + bkg_list) == 0) and (not self.force_single):
            self.log.warning(
                "No selfcal or background exposures provided for self-calibration. Skipping step."
            )
            self.log.info(
                "If you want to force self-calibration with the science "
                "exposure alone (generally not recommended), set force_single=True."
            )
            input_sci.meta.cal_step.badpix_selfcal = "SKIPPED"
            return input_sci, bkg_list

        # get the dispersion axis
        try:
            dispaxis = input_sci.meta.wcsinfo.dispersion_direction
        except AttributeError:
            self.log.warning(
                "Dispersion axis not found in input science image metadata. "
                "Kernel for self-calibration will be two-dimensional."
            )
            dispaxis = None

        # collapse all selfcal exposures into a single background model
        # note that selfcal_list includes the science exposure. This is expected.
        # all exposures are combined into a single background model using a MIN operation.
        selfcal_list = [input_sci] + selfcal_list
        selfcal_3d = []
        for selfcal_model in selfcal_list:
            # If working with MIRI MRS data, subtract any pedestal residual dark
            # using the median of count rates from between the channels
            if input_sci.meta.instrument.detector.upper() == "MIRIFUSHORT":
                selfcal_3d.append(
                    selfcal_model.data
                    - np.nanmedian(
                        selfcal_model.data[
                            PEDESTAL_INDX_MIRIFUSHORT[0] : PEDESTAL_INDX_MIRIFUSHORT[1],
                            PEDESTAL_INDX_MIRIFUSHORT[2] : PEDESTAL_INDX_MIRIFUSHORT[3],
                        ]
                    )
                )
            elif input_sci.meta.instrument.detector.upper() == "MIRIFULONG":
                selfcal_3d.append(
                    selfcal_model.data
                    - np.nanmedian(
                        selfcal_model.data[
                            PEDESTAL_INDX_MIRIFULONG[0] : PEDESTAL_INDX_MIRIFULONG[1],
                            PEDESTAL_INDX_MIRIFULONG[2] : PEDESTAL_INDX_MIRIFULONG[3],
                        ]
                    )
                )
            else:
                selfcal_3d.append(selfcal_model.data)

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")
            minimg = np.nanmin(np.asarray(selfcal_3d), axis=0)
        bad_indices = badpix_selfcal.badpix_selfcal(
            minimg, self.flagfrac_lower, self.flagfrac_upper, self.kernel_size, dispaxis
        )

        # apply the flags to the science data
        input_sci = badpix_selfcal.apply_flags(input_sci, bad_indices)

        self.log.info(f"Number of new bad pixels flagged: {len(bad_indices[0])}")
        # apply the flags to the background data to be passed to background sub step
        if len(bkg_list) > 0:
            for i, background_model in enumerate(bkg_list):
                bkg_list[i] = badpix_selfcal.apply_flags(dm.open(background_model), bad_indices)

        if self.save_flagged_bkg:
            self.save_bkg(bkg_list)

        input_sci.meta.cal_step.badpix_selfcal = "COMPLETE"
        return input_sci, bkg_list


def _parse_inputs(input_data, selfcal_list, bkg_list):
    """
    Parse the input to the step.

    Parameters
    ----------
    input_data : JWSTDataModel, filename, or association
        Input exposures to be split into science, background, and selfcal lists.
    selfcal_list : list of ImageModels or filenames, default None
        Exposures to include as part of median background model used to find bad pixels,
        but that are not flagged and returned as background exposures
    bkg_list : list of ImageModels or filenames, default None
        Exposures to include as part of median background model used to find bad pixels,
        and that are flagged and returned as background exposures

    Returns
    -------
    input_sci : JWSTDataModel
        Input science data to be corrected. Will be a single datamodel regardless of input type.
    selfcal_list : list[JWSTDataModel]
        Images to use for self-calibration.
    bkg_list : list[JWSTDataModel]
        Images to use as background exposures.
    """
    if selfcal_list is None:
        selfcal_list = []
    selfcal_list = [dm.open(selfcal) for selfcal in selfcal_list]
    if bkg_list is None:
        bkg_list = []
    bkg_list = [dm.open(bkg) for bkg in bkg_list]
    selfcal_list = selfcal_list + bkg_list

    with dm.open(input_data) as input_dm:
        # find science and background exposures in association file
        if isinstance(input_dm, dm.ModelContainer):
            sci_models, bkg_list_asn, selfcal_list_asn = split_container_by_asn_exptype(
                input_dm, exptypes=["science", "background", "selfcal"]
            )

            selfcal_list = selfcal_list + list(bkg_list_asn) + list(selfcal_list_asn)
            bkg_list += list(bkg_list_asn)

            if len(sci_models) > 1:
                raise ValueError(
                    "Input data contains multiple science exposures. "
                    "This is not supported in calwebb_spec2 steps."
                )
            input_sci = sci_models[0]

        elif isinstance(input_dm, dm.IFUImageModel) or isinstance(input_dm, dm.ImageModel):
            input_sci = input_dm

        else:
            raise TypeError(
                "Input data is not a ModelContainer, ImageModel, or IFUImageModel. Cannot continue."
            )

    return input_sci, selfcal_list, bkg_list


def split_container_by_asn_exptype(container: dm.ModelContainer, exptypes: list) -> list:
    """
    Split a ModelContainer into a list of ImageModels based on exposure type.

    Parameters
    ----------
    container : ModelContainer
        The input association.
    exptypes : list[str]
        List of exposure types to split on.

    Returns
    -------
    split_list : list of lists
        Lists of ImageModels, where the outer list is indexed by the input exptypes and the
        inner list contains the ImageModels of that type.
    """
    split_list = []
    for exptype in exptypes:
        good_models = [
            container[j] for j in range(len(container)) if container[j].meta.asn.exptype == exptype
        ]
        split_list.append(good_models)

    return split_list

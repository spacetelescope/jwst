import warnings
from ..stpipe import Step
from . import badpix_selfcal
import numpy as np
from jwst import datamodels as dm

__all__ = ["BadpixSelfcalStep"]


class BadpixSelfcalStep(Step):
    """
    BadpixSelfcalStep: Flags residual artifacts as bad pixels in the DQ array
    of a JWST exposure using a median filter and percentile cutoffs.

    All input exposures in the association file (or manually-provided bkg_list) are combined
    into a single background model using a MIN operation. The bad pixels are then identified
    using a median filter and percentile cutoffs, and applied to the science data by setting
    the flagged pixels, errors, and variances to NaN
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
    """

    def save_model(self, model, *args, **kwargs):
        """Override save_model to suppress index 0 when save_model is True
        """
        kwargs["idx"] = None
        return Step.save_model(self, model, *args, **kwargs)

    def save_bkg(self, bkg_list, suffix="badpix_selfcal_bkg"):
        """Save the background exposures to file with correct indexing
        """
        for i, bkg_model in enumerate(bkg_list):
            self.save_model(bkg_model, suffix=suffix + f"_{str(i)}")

    def process(self, input, selfcal_list=None, bkg_list=None):
        """
        Flag residual artifacts as bad pixels in the DQ array of a JWST exposure

        Parameters
        ----------
        input : JWST data model or association
            input science data to be corrected, or tuple of (sci, bkg, selfcal)

        selfcal_list : list of ImageModels or filenames, default None
            exposures to include as part of median background model used to find bad pixels,
            but that are not flagged and returned as background exposures

        bkg_list : list of ImageModels or filenames, default None
            exposures to include as part of median background model used to find bad pixels,
            and that are flagged and returned as background exposures

        Returns
        -------
        output : JWST data model or association
            data model with CRs flagged

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
        input_sci, selfcal_list, bkg_list = _parse_inputs(input, selfcal_list, bkg_list)

        # ensure that there are background exposures to use, otherwise skip the step
        # unless forced
        if (len(selfcal_list + bkg_list) == 0) and (not self.force_single):
            self.log.warning("No selfcal or background exposures provided for self-calibration. "
                        "Skipping step.")
            self.log.info("If you want to force self-calibration with the science "
                        "exposure alone (generally not recommended), set force_single=True.")
            input_sci.meta.cal_step.badpix_selfcal = 'SKIPPED'
            return input_sci, bkg_list

        # get the dispersion axis
        try:
            dispaxis = input_sci.meta.wcsinfo.dispersion_direction
        except AttributeError:
            self.log.warning("Dispersion axis not found in input science image metadata. "
                        "Kernel for self-calibration will be two-dimensional.")
            dispaxis = None

        # collapse all selfcal exposures into a single background model
        # note that selfcal_list includes the science exposure. This is expected.
        # all exposures are combined into a single background model using a MIN operation.
        selfcal_list = [input_sci] + selfcal_list
        selfcal_3d = []
        for i, selfcal_model in enumerate(selfcal_list):
            selfcal_3d.append(selfcal_model.data)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")
            minimg = np.nanmin(np.asarray(selfcal_3d), axis=0)
        bad_indices = badpix_selfcal.badpix_selfcal(minimg, self.flagfrac_lower, self.flagfrac_upper, self.kernel_size, dispaxis)

        # apply the flags to the science data
        input_sci = badpix_selfcal.apply_flags(input_sci, bad_indices)

        # apply the flags to the background data to be passed to background sub step
        if len(bkg_list) > 0:
            for i, background_model in enumerate(bkg_list):
                bkg_list[i] = badpix_selfcal.apply_flags(dm.open(background_model), bad_indices)

        if self.save_flagged_bkg:
            self.save_bkg(bkg_list)

        input_sci.meta.cal_step.badpix_selfcal = 'COMPLETE'
        return input_sci, bkg_list


def _parse_inputs(input, selfcal_list, bkg_list):
    """
    Parse the input to the step. This is a helper function that is used in the
    command line interface to the step.

    Parameters
    ----------
    input : str
        input science exposure or association

    selfcal_list : list of ImageModels or filenames, default None
        exposures to include as part of median background model used to find bad pixels,
        but that are not flagged and returned as background exposures

    bkg_list : list of ImageModels or filenames, default None
        exposures to include as part of median background model used to find bad pixels,
        and that are flagged and returned as background exposures

    Returns
    -------
    input : JWST data model or association
        input science data to be corrected

    selfcal_list : list of ImageModels or filenames to use for selfcal
    """
    if selfcal_list is None:
        selfcal_list = []
    selfcal_list = [dm.open(selfcal) for selfcal in selfcal_list]
    if bkg_list is None:
        bkg_list = []
    bkg_list = [dm.open(bkg) for bkg in bkg_list]
    selfcal_list = selfcal_list + bkg_list

    with dm.open(input) as input_data:

        # find science and background exposures in association file
        if isinstance(input_data, dm.ModelContainer):

            sci_models, bkg_list_asn, selfcal_list_asn = split_container_by_asn_exptype(
                input_data, exptypes=['science', 'background', 'selfcal'])

            selfcal_list = selfcal_list + list(bkg_list_asn) + list(selfcal_list_asn)
            bkg_list += list(bkg_list_asn)

            if len(sci_models) > 1:
                raise ValueError("Input data contains multiple science exposures. "
                                 "This is not supported in calwebb_spec2 steps.")
            input_sci = sci_models[0]

        elif isinstance(input_data, dm.IFUImageModel) or isinstance(input_data, dm.ImageModel):

            input_sci = input_data

        else:
            raise TypeError("Input data is not a ModelContainer, ImageModel, or IFUImageModel. "
                            "Cannot continue.")

    return input_sci, selfcal_list, bkg_list


def split_container_by_asn_exptype(container: dm.ModelContainer, exptypes: list) -> list:
    """
    Split a ModelContainer into a list of ImageModels based on exposure type.

    Parameters
    ----------
    container : ModelContainer
        input data

    exptype : str
        exposure type to split on

    Returns
    -------
    split_list : list of lists
        lists of ImageModels
    """
    split_list = []
    for exptype in exptypes:
        good_models = [container[j] for j in range(len(container)) if container[j].meta.asn.exptype == exptype]
        split_list.append(good_models)

    return split_list

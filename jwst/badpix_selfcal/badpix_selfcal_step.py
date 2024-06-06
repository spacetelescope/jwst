

from ..stpipe import Step
from . import badpix_selfcal
import numpy as np
from jwst import datamodels as dm
from jwst.master_background.master_background_step import split_container

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["BadpixSelfcalStep"]


class BadpixSelfcalStep(Step):
    """
    BadpixSelfcalStep: Flags residual artifacts as bad pixels in the DQ array
    of a JWST exposure using a median filter and percentile cutoffs.

    All input exposures in the association file (or manually-provided bkg_list) are combined
    into a single background model using a MIN operation. The bad pixels are then identified
    using a median filter and percentile cutoffs, and applied to the science data by setting
    the flagged pixels to NaN and the DQ flag to DO_NOT_USE + OTHER_BAD_PIXEL.
    """

    class_alias = "badpix_selfcal"
    bkg_suffix = "badpix_selfcal"

    spec = """
    flagfrac = float(default=0.001)  #fraction of pixels to flag on each of low and high end
    kernel_size = integer(default=15)  #size of kernel for median filter
    force_single = boolean(default=False)  #force single input exposure
    skip = boolean(default=True)
    """

    def process(self, input, selfcal_list=[]):
        """
        Flag residual artifacts as bad pixels in the DQ array of a JWST exposure

        Parameters
        ----------
        input: JWST data model or association
            input science data to be corrected

        selfcal_list: list of ImageModels or filenames to use for selfcal

        Returns
        -------
        output: JWST data model or association
            data model with CRs flagged

        Notes
        -----
        If selfcal_list is specified manually, it overrides any non-science exposures
        in the association file.
        If selfcal_list is set to None and an association file is read in, all exposures in the
        association file, including science, background, and selfcal exposures,
        are included in the MIN frame from which outliers are detected.
        If selfcal_list is set to None and input is a single science exposure, the step will
        be skipped with a warning unless the force_single parameter is set True.
        In that case, the input exposure will be used as the sole background exposure,
        i.e., true self-calibration.
        """
        input_sci, selfcal_list, bkg_list_asn = _parse_inputs(input, selfcal_list)

        # ensure that there are background exposures to use, otherwise skip the step
        # unless forced
        if (len(selfcal_list) == 0) and (not self.force_single):
            log.warning("No background exposures provided for self-calibration. Skipping step.")
            self.record_step_status(input_sci, "badpix_selfcal", success=False)
            return input_sci, bkg_list_asn

        # get the dispersion axis
        try:
            dispaxis = input_sci.meta.wcsinfo.dispersion_direction
        except AttributeError:
            log.warning("Dispersion axis not found in input science image metadata.\
                        Kernel for self-calibration will be two-dimensional.")
            dispaxis = None

        # open all selfcal exposures
        # note that selfcal_list includes the science exposure. This is expected.
        # all exposures are combined into a single background model using a MIN operation.
        selfcal_list = [input_sci] + [dm.open(k) for k in selfcal_list]

        # collapse background dithers into a single background model
        selfcal_3d = []
        for i, selfcal_model in enumerate(selfcal_list):
            selfcal_3d.append(selfcal_model.data)
        minimg = np.nanmin(np.asarray(selfcal_3d), axis=0)
        bad_indices = badpix_selfcal.badpix_selfcal(minimg, self.flagfrac, self.kernel_size, dispaxis)

        # apply the flags to the science data
        input_sci = badpix_selfcal.apply_flags(input_sci, bad_indices)

        # apply the flags to the background data to be passed to background sub step
        if len(bkg_list_asn) > 0:
            for i, background_model in enumerate(bkg_list_asn):
                bkg_list_asn[i] = badpix_selfcal.apply_flags(background_model, bad_indices)

        self.record_step_status(input_sci, "badpix_selfcal", success=True)

        return input_sci, *bkg_list_asn


def _parse_inputs(input, selfcal_list=[]):
    """
    Parse the input to the step. This is a helper function that is used in the
    command line interface to the step.

    Parameters
    ----------
    input: str
        input file or association

    selfcal_list: list
        ImageModels or filenames to use for selfcal

    Returns
    -------
    input: JWST data model or association
        input science data to be corrected

    selfcal_list: list of ImageModels or filenames to use for selfcal
    """
    with dm.open(input) as input_data:

        # find science and background exposures.
        if isinstance(input_data, dm.ModelContainer):

            sci_models, bkg_list_asn, selfcal_list_asn = split_container(input_data,
                                                                         fields=['science', 'background', 'selfcal'])

            if len(selfcal_list) == 0:
                selfcal_list = list(bkg_list_asn) + list(selfcal_list_asn)
            else:
                log.warning("selfcal_list provided directly as input, ignoring background \
                            and selfcal exposure types in the association file")

            # in calwebb_spec2 there should be only a single science exposure in an association
            input_sci = sci_models[0]

        elif isinstance(input_data, dm.IFUImageModel) or isinstance(input_data, dm.ImageModel):

            input_sci = input_data
            bkg_list_asn = []

        else:
            raise TypeError("Input data is not a ModelContainer, ImageModel, or IFUImageModel.\
                            Cannot continue.")

    return input_sci, selfcal_list, bkg_list_asn

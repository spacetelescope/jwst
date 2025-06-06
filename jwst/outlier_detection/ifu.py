"""
Submodule defined for performing outlier detection on IFU data.

This is the controlling routine for the outlier detection process.
It loads and sets the various input data and parameters needed to flag
outliers.  Pixel are flagged as outliers based on the MINIMUM difference
a pixel has with its neighbor across all the input cal files.

Notes
-----
This routine performs the following operations::

  1. Extracts parameter settings from input ModelContainer and merges
     them with any user-provided values
  2. Loop over cal files
     a. read in science data
     b. Store computed neighbor differences for all the pixels.
        The neighbor pixel  differences are defined by the dispersion axis.
        For MIRI, with the dispersion axis along the y axis, the neighbors that are used to
        to find the differences are to the left and right of each pixel being examined.
        For NIRSpec, with the dispersion along the x axis, the neighbors that are used to
        find the differences are above and below the pixel being examined.
  3. For each input file store the  minimum of the pixel neighbor differences
  4. Comparing all the differences from all the input data find the minimum neighbor difference
  5. Normalize minimum difference to local median of difference array
  6. Select outliers by flagging those normalized minimum values > threshold_percent
  7. Updates input ImageModel DQ arrays with mask of detected outliers.
"""

import logging
import warnings

import numpy as np

from jwst.datamodels import ModelContainer
from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.stpipe.utilities import record_step_status
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from stcal.outlier_detection.utils import medfilt

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["detect_outliers"]


def detect_outliers(
    input_models,
    save_intermediate_results,
    kernel_size,
    ifu_second_check,
    threshold_percent,
    make_output_path,
):
    """
    Flag outliers in ifu data.

    Parameters
    ----------
    input_models : ModelContainer or str
        A container of data models or an association file readable into a ModelContainer.
    save_intermediate_results : bool
        If True, save intermediate results.
    kernel_size : str
        The size of the kernel to use to normalize the pixel differences.
        Must only contain odd values. Valid values are a pair of ints in a single string
        (for example “7 7”, the step default).
    ifu_second_check : bool
        If True, perform a secondary check for outliers. This will set outliers
        wherever the difference array of adjacent pixels is a NaN.
    threshold_percent : float
        The threshold (in percent) of the normalized minimum pixel difference
        used to identify bad pixels. Pixels with a normalized minimum difference
        above this percentage are flagged as outliers.

    Returns
    -------
    input_models : ModelContainer
        The input data with DQ flags set for detected outliers.
    """
    if not isinstance(input_models, ModelContainer):
        input_models = ModelContainer(input_models)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models

    sizex, sizey = [int(val) for val in kernel_size.split()]
    kern_size = np.zeros(2, dtype=int)
    kern_size[0] = sizex
    kern_size[1] = sizey

    # check if kernel size is an odd value
    if kern_size[0] % 2 == 0:
        log.info(
            "X kernel size is given as an even number. This value must be an odd number. "
            "Increasing number by 1"
        )
        kern_size[0] = kern_size[0] + 1
        log.info(f"New x kernel size is {kern_size[0]}: ")
    if kern_size[1] % 2 == 0:
        log.info(
            "Y kernel size is given as an even number. This value must be an odd number. "
            "Increasing number by 1"
        )
        kern_size[1] = kern_size[1] + 1
        log.info(f"New y kernel size is {kern_size[1]}: ")

    (diffaxis, ny, nx) = _find_detector_parameters(input_models)

    nfiles = len(input_models)
    detector = np.empty(nfiles, dtype="<U15")
    for i, model in enumerate(input_models):
        detector[i] = model.meta.instrument.detector.lower()

    exptype = input_models[0].meta.exposure.type
    log.info(f"Performing IFU outlier_detection for exptype {exptype}")
    # How many unique values of detector?
    uq_det = np.unique(detector)
    ndet = len(uq_det)
    for idet in range(ndet):
        indx = (np.where(detector == uq_det[idet]))[0]
        ndet_files = int(len(indx))
        flag_outliers(
            input_models,
            idet,
            uq_det,
            ndet_files,
            diffaxis,
            nx,
            ny,
            kern_size,
            threshold_percent,
            save_intermediate_results,
            ifu_second_check,
            make_output_path,
        )
    return input_models


def flag_outliers(
    input_models,
    idet,
    uq_det,
    ndet_files,
    diffaxis,
    nx,
    ny,
    kern_size,
    threshold_percent,
    save_intermediate_results,
    ifu_second_check,
    make_output_path,
):
    """
    Flag outlier pixels on IFU.

    In general we are searching for pixels that
    are a form of a bad pixel but not in bad pixel mask, because the bad pixels vary with
    time. This program will flag the DQ of input images as DO_NOT_USE and OUTLIER and set
    the associated science pixel to a Nan. This routine only works on data from one detector.

    Parameters
    ----------
    idet : int
        Integer indicating which detector we are working with
    uq_det : np.array[str]
        Array of (unique) detector names found input data
    ndet_files : int
        Number of files for the detector we are working on
    diffaxis : int
        The axis to form the adjacent pixel differences
    nx : int
        Size of input data on x axis
    ny : int
        Since of input data on y axis
    kern_size : tuple
        Size of the kernel to use for median filtering
    threshold_percent : float
        Percent for flagging outliers. Flags pixels where the minimum difference between
        adjacent pixels for all the input data for a detector is above this percentage. The
        percentage is based on using all the pixels except a 4 X 4 row and column region around
        the detector that is often noisy.
    save_intermediate_results : bool
        If True then save intermediate output data
    ifu_second_check : bool
        If True then perform a secondary check searching for outliers. This will set outliers
        where ever the difference array of adjacent pixels is a Nan.
    make_output_path : function
        The functools.partial instance to pass to save_median. Has no effect if
        save_intermediate_results is False.
    """
    # set up array to hold group differences
    diffarr = np.zeros([ndet_files, ny, nx])
    j = 0
    for model in input_models:
        detector = model.meta.instrument.detector.lower()
        # only use data from the same detector
        if detector == uq_det[idet]:
            sci = model.data
            dq = model.dq
            bad = np.bitwise_and(dq, dqflags.pixel["DO_NOT_USE"]).astype(bool)
            # set all science data that have DO_NOT_USE to NAN
            sci[bad] = np.nan

            # Compute left and right differences (MIRI dispersion axis = 1 along y axis)
            # For NIRSpec dispersion axis = 0 (along the x axis), these differences are top, bottom
            # prepend = 0 has the effect of keeping the same shape as sci and
            # for MIRI data (disp axis = 1) the first column = sci data
            # OR
            # for NIRSpec data (disp axis = 0) the first row = sci data

            leftdiff = np.diff(sci, axis=diffaxis, prepend=0)
            flip = np.flip(sci, axis=diffaxis)
            rightdiff = np.diff(flip, axis=diffaxis, prepend=0)
            rightdiff = np.flip(rightdiff, axis=diffaxis)

            # Combine left and right differences with minimum of the abs value
            # to avoid artifacts from bright edges
            comb = np.zeros([2, ny, nx])
            comb[0, :, :] = np.abs(leftdiff)
            comb[1, :, :] = np.abs(rightdiff)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
                combdiff = np.nanmin(comb, axis=0)
            diffarr[j, :, :] = combdiff
            j = j + 1

    # minarr final minimum combined differences, size: ny X nx
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        minarr = np.nanmin(diffarr, axis=0)

    # Normalise the differences to a local median image to deal with ultra-bright sources
    normarr = medfilt(minarr, kern_size)
    nfloor = np.nanmedian(minarr) / 3
    normarr[normarr < nfloor] = nfloor  # Ensure we never divide by a tiny number
    minarr_norm = minarr / normarr
    # Percentile cut of the central region (cutting out weird detector edge effects)
    pctmin = np.nanpercentile(minarr_norm[4 : ny - 4, 4 : nx - 4], threshold_percent)
    log.debug(f"Flag pixels with values above {threshold_percent} {pctmin}: ")
    # Flag everything above this percentile value. Using np.where here because we count
    # the number of pixels flagged using len(indx[0])
    indx = minarr_norm > pctmin
    num_above = indx.sum()

    if save_intermediate_results:
        detector_name = uq_det[idet]
        opt_info = (
            kern_size[0],
            kern_size[1],
            threshold_percent,
            diffarr,
            minarr,
            normarr,
            minarr_norm,
        )
        opt_model = create_optional_results_model(opt_info)
        opt_model.meta.filename = make_output_path(
            basepath=input_models.meta.asn_table.products[0].name,
            suffix=detector_name + "_outlier_output",
        )
        log.info(f"Writing out intermediate outlier file {opt_model.meta.filename}")
        opt_model.save(opt_model.meta.filename)

    del diffarr

    # store some information if the second flagging step is to be done.
    if ifu_second_check:
        # store where the minarr is nan (neighbor pixels have nan so differences produces a nan)
        nanminarr = np.isnan(minarr)
        nanindx = np.where(nanminarr)

    # Update DQ flag
    for i in range(len(input_models)):
        detector = input_models[i].meta.instrument.detector.lower()
        # only use data from the same detector
        if detector == uq_det[idet]:
            model = input_models[i]
            sci = model.data
            dq = model.dq

            # There could be a large number of pixels with a sci value of NaN
            # but the dq flag of DO_NOT_USE has not been set.
            # This can occur in Non-science regions of the detector.
            check = np.where(
                np.logical_and(
                    ~np.bitwise_and(dq, dqflags.pixel["DO_NOT_USE"]).astype(bool), np.isnan(sci)
                )
            )
            log.debug(
                "Number of pixels for which DQ was not set to DO_NOT_USE "
                f"and SCI array was NaN: {len(check[0])} "
            )
            # set all pixels with dq = DO_NOT_USE to have sci values of Nan
            bad = np.bitwise_and(dq, dqflags.pixel["DO_NOT_USE"]).astype(bool)
            sci[bad] = np.nan

            # Basic setting outliers: flagging those at are found in from Percentage cut
            sci[indx] = np.nan
            dq[indx] = np.bitwise_or(dq[indx], dqflags.pixel["DO_NOT_USE"])
            dq[indx] = np.bitwise_or(dq[indx], dqflags.pixel["OUTLIER"])

            nadditional = 0
            # Second level of setting outliers: flagging pixels were minarr was a Nan
            # This will also catch pixels that have a sci of Nan but the DQ flags did
            # not have DO_NOT_USE set
            if ifu_second_check:
                # For counting purposes, count number of science values that were valid (not Nan)
                # after basic flagging in the nanminarr region that will now  be flagged as a Nan.
                nadditional = (~np.isnan(sci[nanindx])).sum()
                sci[nanindx] = np.nan
                dq[nanindx] = np.bitwise_or(dq[nanindx], dqflags.pixel["DO_NOT_USE"])
                dq[nanindx] = np.bitwise_or(dq[nanindx], dqflags.pixel["OUTLIER"])
                log.info(
                    "Number of outlier pixels flagged main ifu outlier flagging: "
                    f"{len(indx[0])} on detector {uq_det[idet]} "
                )
                log.info(
                    "Number of outlier pixels flagged in second check: "
                    f"{nadditional} on detector {uq_det[idet]} "
                )

            total_bad = num_above + nadditional
            percent_cr = total_bad / (model.data.shape[0] * model.data.shape[1]) * 100
            log.info(f"Total #  pixels flagged as outliers: {total_bad} ({percent_cr:.2f}%)")

            # Make sure all error and variance arrays also have matching
            # NaNs and DQ flags
            match_nans_and_flags(model)

            # update model
            input_models[i] = model


def _find_detector_parameters(input_models):
    """
    Find the size of data and the axis to form the differences (perpendicular to dispaxis).

    Parameters
    ----------
    input_models : ModelContainer
        The input data models.

    Returns
    -------
    diffaxis : int
        The axis perpendicular to dispaxis.
    ny : int
        The size of the data on the y axis.
    nx : int
        The size of the data on the x axis.
    """
    if input_models[0].meta.instrument.name.upper() == "MIRI":
        diffaxis = 1
    elif input_models[0].meta.instrument.name.upper() == "NIRSPEC":
        diffaxis = 0
    ny, nx = input_models[0].data.shape
    return (diffaxis, ny, nx)


def create_optional_results_model(opt_info):
    """
    Create an OutlierOutputModel from the computed arrays from outlier detection on IFU data.

    Parameters
    ----------
    opt_info : tuple
        The output arrays needed for the OutlierOutputModel.

    Returns
    -------
    opt_model : OutlierIFUOutputModel
        The optional OutlierIFUOutputModel to be returned from the outlier_detection_ifu step.
    """
    (kernsize_x, kernsize_y, threshold_percent, diffarr, minarr, normarr, minnorm) = opt_info
    opt_model = datamodels.OutlierIFUOutputModel(
        diffarr=diffarr, minarr=minarr, normarr=normarr, minnorm=minnorm
    )
    opt_model.meta.kernel_xsize = kernsize_x
    opt_model.meta.kernel_ysize = kernsize_y
    opt_model.meta.threshold_percent = threshold_percent
    return opt_model

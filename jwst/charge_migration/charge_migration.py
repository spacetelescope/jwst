import logging
import warnings

import numpy as np
from astropy.stats import SigmaClippedStats
from numpy.polynomial.polynomial import Polynomial, polyfit
from scipy.ndimage import binary_dilation
from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)

GOOD = dqflags.group["GOOD"]
DNU = dqflags.group["DO_NOT_USE"]
CHLO = dqflags.group["CHARGELOSS"]

CHLO_DNU = CHLO + DNU

__all__ = ["flag_charge_migration", "correct_charge_migration"]


def flag_charge_migration(output_model, signal_threshold):
    """
    Flag charge migration.

    Each group in each ramp exceeding the signal threshold is flagged as
    CHARGELOSS and DO_NOT_USE, skipping groups already flagged as DO_NOT_USE.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        The data model to flag. Updated in place.
    signal_threshold : float
        Science value above which a group will be flagged as CHARGELOSS.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        Data model with flagging applied. CHARGELOSS and DO_NOT_USE flags
        are added to groups exceeding ``signal_threshold``.
    """
    data = output_model.data
    gdq = output_model.groupdq

    log.info("Using signal_threshold: %.2f", signal_threshold)

    n_ints, n_grps, n_rows, n_cols = gdq.shape
    chargeloss_pix = np.where((data > signal_threshold) & (gdq != DNU))

    for k in range(len(chargeloss_pix[0])):
        integ, group = chargeloss_pix[0][k], chargeloss_pix[1][k]
        row, col = chargeloss_pix[2][k], chargeloss_pix[3][k]
        gdq[integ, group:, row, col] |= CHLO_DNU

        # North
        if row > 0:
            gdq[integ, group:, row - 1, col] |= CHLO_DNU

        # South
        if row < (n_rows - 1):
            gdq[integ, group:, row + 1, col] |= CHLO_DNU

        # East
        if col < (n_cols - 1):
            gdq[integ, group:, row, col + 1] |= CHLO_DNU

        # West
        if col > 0:
            gdq[integ, group:, row, col - 1] |= CHLO_DNU

    return output_model


def correct_charge_migration(
    output_model, readnoise_model, signal_threshold=10000.0, mask_pad=2, fit_info=False
):
    """
    Correct for charge migration.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        The data model to correct. Updated in place.
    readnoise_model : `~stdatamodels.jwst.datamodels.ReadnoiseModel`
        Read noise reference data associated with ``output_model``.
    signal_threshold : float, optional
        Science value above which the correction is triggered.
    mask_pad : int, optional
        Number of pixels to mask on either side of a pixel with charge loss.
    fit_info : bool, optional
        If True, return a diagnostic table of the ramp fits performed.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        Data model with correction applied and CHARGELOSS flags added to
        groups exceeding ``signal_threshold``.
    """
    # Input science and DQ data
    sci = output_model.data
    gdq = output_model.groupdq
    nints, ngroups, ny, nx = sci.shape

    # Readnoise reference file data
    rnvals = readnoise_model.data

    # Identify anywhere that the DO_NOT_USE or SATURATED bits are set
    bad = dqflags.pixel["DO_NOT_USE"] | dqflags.pixel["SATURATED"]
    badindx = (gdq & bad) != 0

    # Group by group differences in the ramp, masking bad pixels
    masked_ramp_diff = np.diff(sci, axis=1)
    masked_ramp_diff[badindx[:, 1:, :, :]] = np.nan
    masked_ramp_diff[badindx[:, :-1, :, :]] = np.nan
    del badindx

    # Median the diffs across integrations, with loop over first dimension to save memory
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        med_masked_ramp_diff = np.empty(masked_ramp_diff.shape[1:], dtype=masked_ramp_diff.dtype)
        for i in range(med_masked_ramp_diff.shape[0]):
            np.nanmedian(
                masked_ramp_diff[:, i, ...],
                axis=0,
                overwrite_input=True,
                out=med_masked_ramp_diff[i, ...],
            )
    del masked_ramp_diff

    # Find the maximum DN value reached anywhere on the ramps across groups and ints
    # (Generally this is the last pixel, but not always given linearity
    # corrections, saturation flagging, etc)
    maxvals = np.nanmax(sci, axis=(0, 1))

    # Identify pixels where max is greater than threshold
    mask = maxvals > signal_threshold

    # Grow the mask with the specified padding
    box_size = 2 * mask_pad + 1
    structure = np.full((box_size, box_size), True)
    mask = binary_dilation(mask, structure)

    # Get clipped stats for the median diff array
    valid_for_fit = np.isfinite(med_masked_ramp_diff) & mask
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        warnings.filterwarnings("ignore", "Degrees of freedom", RuntimeWarning)
        stats = SigmaClippedStats(med_masked_ramp_diff, axis=0, mask=~valid_for_fit)
        clipped_med = stats.median()
        sigma = stats.std()

    # Clip values where n_diffs > 5
    n_diffs = np.sum(valid_for_fit, axis=0)
    sigma[n_diffs <= 5] = 0.0
    clip_mask = np.abs(med_masked_ramp_diff - clipped_med) < 5 * sigma
    valid_for_fit[:, n_diffs > 5] &= clip_mask[:, n_diffs > 5]

    # Fit clipped pixels in the thresholded mask
    yfit, xfit = np.where(mask)
    nfit = len(xfit)
    log.info("Identified %s pixels to fit for BFE", str(nfit))

    ramp_diff = np.diff(sci, axis=1)
    gvec = np.arange(ngroups - 1)  # group diff index
    fit_info_table = []
    for ii in range(nfit):
        xx = xfit[ii]
        yy = yfit[ii]

        one_diff = med_masked_ramp_diff[:, yy, xx]
        indx = valid_for_fit[:, yy, xx]

        # Only proceed if at least 2 good differences
        if n_diffs[yy, xx] < 2:
            continue

        # Fit a polynomial to the group differences; only a constant term and a slope
        coeff = polyfit(gvec[indx], one_diff[indx], 1)
        full_model = Polynomial(coeff)
        slope = coeff[1]

        # First quality cuts on the polynomial fit: require that value
        # change due to slope over the range of fit values is greater than 3x the read noise
        # and that the full_model value at 0 group difference is positive
        good_slope = (np.abs(n_diffs[yy, xx] * slope) > np.abs(3 * rnvals[yy, xx])) & (
            full_model(0) > 0
        )
        if not good_slope:
            continue

        # Otherwise, good fit: make a model that uses only the slope term
        model = Polynomial([0, slope])

        # Second quality cut on the polynomial fit; if ndiffs > 5,
        # require that rms improves significantly
        # This deals with some strange cases with poorly conditioned ramps,
        # e.g. pixel (291,489) in jw04829001001_07101_00001_mirifushort
        if n_diffs[yy, xx] > 5:
            sigma_orig = SigmaClippedStats(one_diff[indx]).std()
            sigma_new = SigmaClippedStats(one_diff[indx] - model(gvec[indx])).std()

            # RMS of corrected data must be < 45% the rms of the original
            # (dialed based on visual inspection of thousands of MRS ramps)
            if sigma_new > 0.45 * sigma_orig:
                continue

        # Save fit data for inspection
        # xx, yy, gvec, one_diff[indx], full_model(gvec), sigma, sigma_new
        # value change: ndiffs * slope / rnvals[yy, xx]
        # pass/fail: model(10) != 0
        if fit_info:
            one_diff[~indx] = np.nan
            value_change = n_diffs[yy, xx] * slope / rnvals[yy, xx]
            fit_info_table.append(
                (xx, yy, one_diff, full_model(gvec), sigma[yy, xx], value_change, True)
            )

        # Subtract off the slope term from the original group differences
        ramp_diff[:, :, yy, xx] -= model(gvec)

    # Now reconstruct ramps from the modified differences:
    # each group builds on the previous corrected one.
    for gg in range(1, ngroups):
        sci[:, gg, :, :] = sci[:, gg - 1, :, :] + ramp_diff[:, gg - 1, :, :]

    output_model.data = sci
    return output_model

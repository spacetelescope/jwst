import logging
import warnings

import numpy as np
from astropy.stats import SigmaClippedStats
from astropy.table import Table
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
    output_model, readnoise_model, correction_threshold=10000.0, mask_pad=2, fit_info=False
):
    """
    Correct for the brighter-fatter effect (BFE).

    The process is:

    1. Identify pixels with maximum counts over the signal threshold, along with
       neighbors within the specified padding window.
    2. Compute diffs of all groups in the ramp, masking any bad pixels, and
       median-combine across all integrations.
    3. If more than 5 diffs are available, reject outliers with sigma clipping.
    4. For each identified pixel, fit a line to the clipped diffs.
    5. Reject the fit if it does not meet these criteria:

       a. The value change due to the slope over the range of the fit values is
          greater than 3 times the read noise.
       b. The fit value for the first group difference is positive.
       c. If there are more than 5 diffs, the RMS of the corrected data is less
          than 45% of the original RMS.

    6. For pixels with valid corrections, subtract the fit slope from the group
       differences.
    7. Reconstruct the ramps from the modified differences, adding each new diff
       into the previous corrected group.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        The data model to correct. Updated in place.
    readnoise_model : `~stdatamodels.jwst.datamodels.ReadnoiseModel`
        Read noise reference data associated with ``output_model``.
    correction_threshold : float, optional
        Science value above which the correction is triggered.
    mask_pad : int, optional
        Number of pixels to correct on either side of a pixel exceeding the threshold.
    fit_info : bool, optional
        If True, return a diagnostic table of the ramp fits performed.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        Data model with correction applied.
    fit_info_table : `~astropy.table.Table` or None, optional
        Table containing information about all fits attempted. Columns are:
        ``x, y, fit_success, value_change, sigma_orig, sigma_new, diff, model``.
        Returned only if ``fit_info`` is True. If there were no pixels to correct,
        the return value is None.
    """
    # Input science and DQ data
    sci = output_model.data
    gdq = output_model.groupdq
    nints, ngroups, ny, nx = sci.shape

    # Readnoise reference file data
    rnvals = readnoise_model.data

    # Find the maximum DN value reached anywhere on the ramps across groups and ints
    # (Generally this is the last pixel, but not always given linearity
    # corrections, saturation flagging, etc)
    maxvals = np.nanmax(sci, axis=(0, 1))

    # Identify pixels where max is greater than threshold
    mask = maxvals > correction_threshold

    # Check if there is anything to do
    if not np.any(mask):
        log.info("No pixels to correct.")
        if fit_info:
            return output_model, None
        else:
            return output_model

    # Grow the mask with the specified padding
    box_size = 2 * mask_pad + 1
    structure = np.full((box_size, box_size), True)
    mask = binary_dilation(mask, structure)

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

    # Get clipped stats for the median diff array
    valid_for_fit = np.isfinite(med_masked_ramp_diff) & mask
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        warnings.filterwarnings("ignore", "Degrees of freedom", RuntimeWarning)
        stats = SigmaClippedStats(med_masked_ramp_diff, axis=0, mask=~valid_for_fit)
        clipped_med = stats.median()
        clipped_sigma = stats.std()

    # Clip values where n_diffs > 5
    n_diffs = np.sum(valid_for_fit, axis=0)
    clip_mask = np.abs(med_masked_ramp_diff - clipped_med) < 5 * clipped_sigma
    valid_for_fit[:, n_diffs > 5] &= clip_mask[:, n_diffs > 5]
    del clipped_med, clip_mask, clipped_sigma

    # Clip again for a good measure of sigma
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        warnings.filterwarnings("ignore", "Degrees of freedom", RuntimeWarning)
        stats = SigmaClippedStats(med_masked_ramp_diff, axis=0, mask=~valid_for_fit)
        sigma = stats.std()

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
        if good_slope:
            # Good fit: make a model that uses only the slope term
            model = Polynomial([0, slope])
        else:
            model = None

        # Second quality cut on the polynomial fit; if ndiffs > 5,
        # require that rms improves significantly
        # This deals with some strange cases with poorly conditioned ramps,
        # e.g. pixel (291,489) in jw04829001001_07101_00001_mirifushort
        sigma_new = np.nan
        if n_diffs[yy, xx] > 5 and model is not None:
            sigma_new = SigmaClippedStats(one_diff[indx] - model(gvec[indx])).std()

            # RMS of corrected data must be < 45% the rms of the original
            # (dialed based on visual inspection of thousands of MRS ramps)
            if sigma_new > 0.45 * sigma[yy, xx]:
                model = None

        # Save fit data for inspection
        if fit_info:
            # median diff for this pixel is not needed anymore: fill clipped values with NaN
            one_diff[~indx] = np.nan
            value_change = n_diffs[yy, xx] * slope / rnvals[yy, xx]
            success = model is not None
            row = (
                xx,
                yy,
                success,
                value_change,
                sigma[yy, xx],
                sigma_new,
                one_diff,
                full_model(gvec),
            )
            fit_info_table.append(row)

        if model is not None:
            # Subtract off the slope term from the original group differences
            ramp_diff[:, :, yy, xx] -= model(gvec)

    # Now reconstruct ramps from the modified differences:
    # each group builds on the previous corrected one.
    for gg in range(1, ngroups):
        sci[:, gg, :, :] = sci[:, gg - 1, :, :] + ramp_diff[:, gg - 1, :, :]
    output_model.data = sci

    if fit_info:
        # Return a diagnostic table along with the corrected output
        fit_info_table = Table(
            rows=fit_info_table,
            names=[
                "x",
                "y",
                "fit_success",
                "value_change",
                "sigma_orig",
                "sigma_new",
                "diff",
                "model",
            ],
        )
        return output_model, fit_info_table
    else:
        # Just return the model
        return output_model

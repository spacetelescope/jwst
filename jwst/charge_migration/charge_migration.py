import logging
import warnings

import numpy as np
from astropy.stats import sigma_clipped_stats
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


def correct_charge_migration(output_model, readnoise_model, signal_threshold=10000.0, mask_pad=2):
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

    # Make a copy of the ramp where we NaN everything that's bad
    masked_ramp = sci.copy()
    masked_ramp[badindx] = np.nan

    # Group by group differences in the original and working ramps
    ramp_diff = np.diff(sci, axis=1)
    masked_ramp_diff = np.diff(masked_ramp, axis=1)

    # Median across integrations, with loop over first dimension to save memory
    with warnings.catch_warnings():
        med_masked_ramp_diff = np.empty(masked_ramp_diff.shape[1:], dtype=masked_ramp_diff.dtype)
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        for i in range(med_masked_ramp_diff.shape[0]):
            np.nanmedian(
                masked_ramp_diff[:, i, ...],
                axis=0,
                overwrite_input=True,
                out=med_masked_ramp_diff[i, ...],
            )

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

    # Fit all pixels in this mask
    yfit, xfit = np.where(mask)
    nfit = len(xfit)
    log.info("Identified %s pixels to fit for BFE", str(nfit))

    gvec = np.arange(ngroups - 1)  # group diff index
    for ii in range(nfit):
        xx = xfit[ii]
        yy = yfit[ii]

        thisdiff = med_masked_ramp_diff[:, yy, xx]
        indx = np.where(np.isfinite(thisdiff))

        # Only proceed if at least 2 good differences
        ndiffs = len(indx[0])
        if ndiffs >= 2:
            # If fewer than 5 good differences don't do rejection
            if ndiffs < 5:
                indx = np.where(np.isfinite(thisdiff))

            # If more than 5 good differences sigma clip bad values
            else:
                _, themed, thesigma = sigma_clipped_stats(thisdiff[indx])
                indx = np.where(
                    (np.isfinite(thisdiff)) & (np.abs(thisdiff - themed) < 5 * thesigma)
                )

            # Fit a polynomial to the group differences; only a constant term and a slope
            coeff = np.polyfit(gvec[indx], thisdiff[indx], 1, cov=False)
            fullmodel = np.poly1d(coeff)

            # First quality cuts on the polynomial fit; require that value
            # change due to slope over the range of fit values is greater than 3x the read noise
            # and that the fullmodel value at 0 group difference is positive
            if (np.abs(ndiffs * coeff[0]) > np.abs(3 * rnvals[yy, xx])) & (fullmodel(0) > 0):
                # Make a model that uses only the slope term
                model = np.poly1d([coeff[0], 0])
            else:
                model = np.poly1d([0, 0])  # Otherwise use a zero model

            # Second quality cut on the polynomial fit; if ndiffs > 5 and the
            # model isn't zero everywhere already, require that rms improves significantly
            # This deals with some strange cases with poorly conditioned ramps,
            # e.g. pixel (291,481) in jw04829001001_07101_00001_mirifushort
            if (ndiffs > 5) & (model(10) != 0):
                # model(10) is arbitrary, it's just evaluating the model
                # somewhere other than the origin
                _, _, sigmaorig = sigma_clipped_stats(thisdiff[indx])
                _, _, sigmanew = sigma_clipped_stats((thisdiff - model(gvec))[indx])

                # RMS of corrected data must be < 45% the rms of the original
                # (dialed based on visual inspection of thousands of MRS ramps)
                # Otherwise set the model to zero
                if sigmanew > 0.45 * sigmaorig:
                    # Otherwise use a zero model
                    model = np.poly1d([0, 0])

            # Save fit data for inspection
            # xx, yy, gvec, thisdiff[indx], fullmodel(gvec), sigmaorig, sigmanew
            # value change: ndiffs * coeff[0] / rnvals[yy, xx]
            # pass/fail: model(10) != 0

            # Subtract off the slope term from the original group differences
            ramp_diff[:, :, yy, xx] -= model(gvec)

    # Now reconstruct ramps from the modified differences:
    # each group builds on the previous corrected one.
    for gg in range(1, ngroups):
        sci[:, gg, :, :] = sci[:, gg - 1, :, :] + ramp_diff[:, gg - 1, :, :]

    output_model.data = sci
    return output_model

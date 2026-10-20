import logging
import warnings

import numpy as np
from astropy.stats import sigma_clipped_stats
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
    from matplotlib import pyplot as plt

    plotfits = False  # Plot all cases where we tried to fit a model and succeeded
    plotfails = False  # Plot all cases where we tried to fit a model but rejected the fit

    sci = output_model.data
    gdq = output_model.groupdq

    # What's the xy array size?
    nx, ny = sci.shape[3], sci.shape[2]

    # Readnoise ref file
    rnvals = readnoise_model.data

    nints, ngroups, _, _ = sci.shape
    bit_dnu = dqflags.pixel["DO_NOT_USE"]
    bit_sat = dqflags.pixel["SATURATED"]

    # Identify anywhere that the DO_NOT_USE or SATURATED bits are set
    badindx = np.where(((gdq & bit_dnu) != 0) | ((gdq & bit_sat) != 0))
    # Make a copy of the ramp where we NaN everything that's bad
    workramp = sci.copy()
    workramp[badindx] = np.nan

    # Group by group differences in the original and working ramps
    rampdiff = np.diff(sci, axis=1)
    workrampdiff = np.diff(workramp, axis=1)

    # Median across integrations
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
        medrampdiff = np.nanmedian(rampdiff, axis=0)
        medworkrampdiff = np.nanmedian(workrampdiff, axis=0)

    # Vector of group numbers
    gvec = np.arange(medrampdiff.shape[0])
    # Array of corrected group differences
    newdiff = rampdiff.copy()

    # Find the maximum DN value reached anywhere on the ramps across groups and ints
    # (Generally this is the last pixel, but not always given linearity
    # corrections, saturation flagging, etc)
    maxvals = np.nanmax(np.nanmax(sci, axis=0), axis=0)
    triggerpix = np.where(maxvals > signal_threshold)
    # We'll make a mask of such pixels
    mask = np.zeros_like(maxvals, dtype=int)
    mask[triggerpix] = 1
    # Loop over all triggering pixels growing the mask
    triggerx = triggerpix[1]
    triggery = triggerpix[0]
    for ii in range(0, len(triggerx)):
        x1 = np.nanmax([triggerx[ii] - mask_pad, 0])
        x2 = np.nanmin([triggerx[ii] + mask_pad + 1, nx])
        y1 = np.nanmax([triggery[ii] - mask_pad, 0])
        y2 = np.nanmin([triggery[ii] + mask_pad + 1, ny])
        mask[y1:y2, x1:x2] = 1

    # Fit all pixels in this mask
    tofit = np.where(mask == 1)
    yfit = tofit[0]
    xfit = tofit[1]
    nfit = len(xfit)
    log.info("Identified " + str(nfit) + " pixels to fit for BFE")

    for ii in range(0, nfit):
        xx = xfit[ii]
        yy = yfit[ii]

        thisdiff = medworkrampdiff[:, yy, xx]
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

            # Default values for sigma if we don't need to calculate them
            sigmaorig, sigmanew = 0, 0

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
                    model = np.poly1d([0, 0])  # Otherwise use a zero model

            # Plot cases where a model was fit successfully
            if plotfits & (model(10) != 0):
                log.info(
                    "(xx,yy) = ",
                    xx,
                    yy,
                )
                log.info("Value change in multiples of RN: ", ndiffs * coeff[0] / rnvals[yy, xx])
                log.info("Sigma_orig, sigma_new = ", sigmaorig, sigmanew)
                plt.plot(gvec[indx], thisdiff[indx], "x")
                plt.plot(gvec, fullmodel(gvec))
                plt.show()

            # Plot cases where a model was fit but rejected
            if plotfails & (model(10) == 0):
                log.info(
                    "(xx,yy) = ",
                    xx,
                    yy,
                )
                log.info("Value change in multiples of RN: ", ndiffs * coeff[0] / rnvals[yy, xx])
                log.info("Sigma_orig, sigma_new = ", sigmaorig, sigmanew)
                plt.plot(gvec[indx], thisdiff[indx], "x")
                plt.plot(gvec, fullmodel(gvec))
                plt.show()

            # Subtract off the slope term from the original group differences
            for ii in range(0, nints):
                newdiff[ii, :, yy, xx] = rampdiff[ii, :, yy, xx] - model(gvec)

                # Now reconstruct ramps from the modified differences
    for gg in range(1, ngroups):
        for ii in range(0, nints):
            sci[ii, gg, :, :] = sci[ii, gg - 1, :, :] + newdiff[ii, gg - 1, :, :]

    # Clean up
    output_model.data = sci

    return output_model

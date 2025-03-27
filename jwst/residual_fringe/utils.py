import logging

import numpy as np
import math
import numpy.polynomial.polynomial as poly

from scipy.interpolate import pchip
from astropy.timeseries import LombScargle
from BayesicFitting import SineModel
from BayesicFitting import LevenbergMarquardtFitter
from BayesicFitting import RobustShell
from BayesicFitting import ConstantModel
from BayesicFitting import Fitter

from jwst.residual_fringe.fitter import spline_fitter


log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

# Number of knots for bkg model if no other info provided
# Hard coded parameter, has been selected based on testing but can be changed
NUM_KNOTS = 80


def slice_info(slice_map, channel):
    """
    Identify pixels by slice.

    Parameters
    ----------
    slice_map : ndarray of int
        2D image containing slice identification values by pixel.
        Slice ID values are integers with the value 100 * channel number
        + slice number.  Pixels not included in a slice have value 0.
    channel : int
        Channel number.

    Returns
    -------
    slices_in_channel : ndarray of int
        1D array of slice IDs included in the channel.
    xrange_channel : ndarray of int
        1D array with two elements: minimum and maximum x indices
        for the channel.
    slice_x_ranges : ndarray of int
        N x 3 array for N slices, where the first column is the slice ID,
        second column is the minimum x index for the slice,
        and the third column is the maximum x index for the slice.
    all_slice_masks : ndarray of int
        N x nx x ny for N slices, matching the x and y shape of the
        input slice_map.  Values are 1 for pixels included in the slice,
        0 otherwise.
    """
    slice_inventory = np.unique(slice_map)
    slices_in_channel = slice_inventory[
        np.where((slice_inventory >= 100 * channel) & (slice_inventory < 100 * (channel + 1)))
    ]

    log.info(f"Number of slices in channel {slices_in_channel.shape[0]} ")
    slice_x_ranges = np.zeros((slices_in_channel.shape[0], 3), dtype=int)
    all_slice_masks = np.zeros((slices_in_channel.shape[0], slice_map.shape[0], slice_map.shape[1]))
    for n, s in enumerate(slices_in_channel):
        # create a mask of the slice
        pixels = np.where(slice_map == s)
        slice_mask = np.zeros(slice_map.shape)

        slice_mask[pixels] = 1

        # add this to the all_slice_mask array
        all_slice_masks[n] = slice_mask

        # get the indices at the start and end of the slice
        collapsed_slice = np.sum(slice_mask, axis=0)
        indices = np.where(collapsed_slice[:-1] != collapsed_slice[1:])[0]
        slice_x_ranges[n, 0], slice_x_ranges[n, 1], slice_x_ranges[n, 2] = (
            int(s),
            int(np.amin(indices)),
            int(np.amax(indices) + 1),
        )

        log.debug(
            f"For slice {slice_x_ranges[n, 0]} x ranges of slices "
            f"region {slice_x_ranges[n, 1]}, {slice_x_ranges[n, 2]}"
        )

    log.debug(
        "Min and max x pixel values of all slices "
        f"in channel {np.amin(slice_x_ranges[:, 1])} {np.amax(slice_x_ranges[:, 2])}"
    )

    xrange_channel = np.zeros(2)
    xrange_channel[0] = np.amin(slice_x_ranges[:, 1])
    xrange_channel[1] = np.amax(slice_x_ranges[:, 2])
    return slices_in_channel, xrange_channel, slice_x_ranges, all_slice_masks


def fill_wavenumbers(wnums):
    """
    Fill in missing wavenumber values.

    Given a wavenumber array with missing values (e.g., columns with
    on-slice and off-slice pixels), fit the good points using a
    polynomial, then use the coefficients to estimate wavenumbers
    on the off-slice pixels.

    Note that these new values are physically meaningless but having them
    in the wavenum array stops the BayesicFitting package from crashing with
    a LinAlgErr.

    Parameters
    ----------
    wnums : ndarray
        The wavenumber array.

    Returns
    -------
    wnums_filled : ndarray
        The wavenumber array with off-slice pixels filled
    """
    # set the off-slice pixels to nans and get their indices
    wnums[wnums == 0] = np.nan
    idx = np.isfinite(wnums)

    # fit the on-slice wavenumbers
    coefs = poly.polyfit(np.arange(wnums.shape[0])[idx], wnums[idx], 3)
    wnums_filled = poly.polyval(np.arange(wnums.shape[0]), coefs)

    # keep the original wnums for the on-slice pixels
    wnums_filled[idx] = wnums[idx]

    # clean
    del idx, coefs

    return wnums_filled


def multi_sine(n_sines):
    """
    Create a multi-sine model.

    Parameters
    ----------
    n_sines : int
        Number of sines to include.

    Returns
    -------
    model : BayesicFitting.SineModel
        The model composed of n sines.
    """
    # make the first sine
    mdl = SineModel()

    # make a copy
    model = mdl.copy()

    # add the copy n - 1 times
    for _ in range(1, n_sines):
        mdl.addModel(model.copy())

    # clean
    del model

    return mdl


def fit_envelope(wavenum, signal):
    """
    Fit the upper and lower envelope of signal using a univariate spline.

    Parameters
    ----------
    wavenum : ndarray
        Wavenumber values.
    signal : ndarray
        Signal values

    Returns
    -------
    lower_fit : ndarray
        Fit to the lower envelope.
    l_x : list
        Input lower wavenum values.
    l_y : list
        Input lower signal values.
    upper_fit : ndarray
        Fit to the upper envelope.
    u_x : list
        Input upper wavenum values.
    u_y : list
        Input lower wavenum values.
    """
    # Detect troughs and mark their location. Define endpoints
    l_x = [wavenum[0]]
    l_y = [signal[0]]
    u_x = [wavenum[0]]
    u_y = [signal[0]]

    for k in np.arange(1, len(signal) - 1):
        if (np.sign(signal[k] - signal[k - 1]) == -1) and (
            (np.sign(signal[k] - signal[k + 1])) == -1
        ):
            l_x.append(wavenum[k])
            l_y.append(signal[k])
        if (np.sign(signal[k] - signal[k - 1]) == 1) and (
            (np.sign(signal[k] - signal[k + 1])) == 1
        ):
            u_x.append(wavenum[k])
            u_y.append(signal[k])

    # Append the last value of (s) to the interpolating values.
    # This forces the model to use the same ending point
    l_x.append(wavenum[-1])
    l_y.append(signal[-1])
    u_x.append(wavenum[-1])
    u_y.append(signal[-1])

    # fit a model
    pcl = pchip(l_x, l_y)
    pcu = pchip(u_x, u_y)

    return pcl(wavenum), l_x, l_y, pcu(wavenum), u_x, u_y


def find_lines(signal, max_amp):
    """
    Determine the location of large spectral features.

    Parameters
    ----------
    signal : ndarray
        Signal data.
    max_amp : ndarray
        Maximum amplitude, by column.  Features larger than
        this value are flagged.

    Returns
    -------
    weights : ndarray
        1D array matching signal dimensions, containing 0 values
        for large features and 1 values where no features were
        detected.
    """
    r_x = np.arange(signal.shape[0] - 1)

    # setup the output arrays
    signal_check = signal.copy()
    weights_factors = np.ones(signal.shape[0])

    # Detect peaks
    u_y, u_x, l_y, l_x = [], [], [], []

    for x in r_x:
        if (np.sign(signal_check[x] - signal_check[x - 1]) == 1) and (
            np.sign(signal_check[x] - signal_check[x + 1]) == 1
        ):
            u_y.append(signal_check[x])
            u_x.append(x)

        if (np.sign(signal_check[x] - signal_check[x - 1]) == -1) and (
            np.sign(signal_check[x] - signal_check[x + 1]) == -1
        ):
            l_y.append(signal[x])
            l_x.append(x)

    for n, amp in enumerate(u_y):
        max_amp_val = max_amp[u_x[n]]
        log.debug("find_lines: check if peak above max amp")
        if amp > max_amp_val:
            # peak in x
            # log.debug("find_lines: flagging neighbours")
            xpeaks = [u_x[n] - 1, u_x[n], u_x[n] + 1]

            # log.debug("find_lines:  find neareast troughs")
            # find nearest troughs

            for xp in xpeaks:
                log.debug(f"find_lines:  checking ind {xp}")

                try:
                    x1 = l_x[np.argsort(np.abs(l_x - xp))[0]]

                    try:
                        x2 = l_x[np.argsort(np.abs(l_x - xp))[1]]

                        if x1 < x2:
                            xlow = x1
                            xhigh = x2
                        if x1 > x2:
                            xhigh = x1
                            xlow = x2

                    except IndexError:
                        # raised if x1 is at the edge
                        xlow = x1
                        xhigh = x1

                    # set the weights to 0 and signal 1
                    log.debug("find_lines: setting weights between troughs to 0")
                    signal_check[xlow:xhigh] = 0
                    weights_factors[xlow:xhigh] = 0

                except IndexError:
                    pass

    log.debug(f"find_lines: Found {len(u_x)} peaks   {len(l_x)} troughs")
    weights_factors[signal_check > max_amp * 2] = 0  # catch any remaining
    # weights_factors[signal_check > np.amax(max_amp)] = 0

    return weights_factors


def check_res_fringes(res_fringe_fit, max_amp):
    """
    Check for regions with bad fringe fits.

    Set the beat where this happens to 0 to avoid making the fringes worse.

    Parameters
    ----------
    res_fringe_fit : ndarray
        The residual fringe fit.
    max_amp : ndarray
        The maximum amplitude array.

    Returns
    -------
    res_fringe_fit: ndarray
        The residual fringe fit with bad fit regions removed and replaced
        with 0.
    flags: ndarray
        1D flag array indicating where the fit was altered, matching
        the size of the first dimension of `res_fringe_fit`.
        1 indicates a bad fit region; 0 indicates a good region, left
        unchanged.
    """
    flags = np.zeros(res_fringe_fit.shape[0])

    # get fit envelope
    npix = np.arange(res_fringe_fit.shape[0])
    lenv_fit, _, _, uenv_fit, _, _ = fit_envelope(npix, res_fringe_fit)

    # get the indices of the nodes (where uenv slope goes from
    # negative to positive), add 0 and 1023
    node_ind = [0]
    for k in np.arange(1, len(uenv_fit) - 1):
        if (np.sign(uenv_fit[k] - uenv_fit[k - 1]) == -1) and (
            (np.sign(uenv_fit[k] - uenv_fit[k + 1])) == -1
        ):
            node_ind.append(k)
    node_ind.append(res_fringe_fit.shape[0] - 1)
    node_ind = np.asarray(node_ind)
    log.debug(f"check_res_fringes: found {len(node_ind)} nodes")

    # find where res_fringes goes above max_amp
    runaway_rfc = np.argwhere((np.abs(lenv_fit) + np.abs(uenv_fit)) > (max_amp * 2))

    # check which signal env the blow ups are located in and set to 1, and set a flag array
    if len(runaway_rfc) > 0:
        log.debug(f"check_res_fringes: {len(runaway_rfc)} data points exceed threshold")
        log.debug("check_res_fringes: resetting fits to related beats")
        for i in runaway_rfc:
            # find where the index is compared to the nodes
            node_loc = np.searchsorted(node_ind, i)

            # set the res_fringes between the nodes to 0
            lind = node_ind[node_loc - 1]
            uind = node_ind[node_loc]
            res_fringe_fit[lind[0] : uind[0]] = 0
            flags[lind[0] : uind[0]] = 1  # set flag to 1 for reject fit region

    return res_fringe_fit, flags


def interp_helper(mask):
    """
    Create a convenience function for indexing low-weight values.

    Low-weight is defined to be a value < 1e-5.

    Parameters
    ----------
    mask : ndarray
        The 1D mask array (weights).

    Returns
    -------
    index_array : ndarray of bool
        Boolean index array for low weight pixels.
    index_function : callable
        A function, with signature indices = index_function(index_array),
        to convert logical indices to equivalent direct index values.
    """
    return mask < 1e-05, lambda z: z.nonzero()[0]


def fit_1d_background_complex(flux, weights, wavenum, ffreq=None, channel=1):
    """
    Fit the background signal using a piecewise spline.

    Note that this will also try to identify obvious emission lines
    and flag them, so they aren't considered in the fitting.

    Parameters
    ----------
    flux : ndarray
        1D array of fluxes.
    weights : ndarray
        1D array of weights.
    wavenum : ndarray
        1D array of wavenumbers.
    ffreq : float, optional
        The expected fringe frequency, used to determine number of knots.
        If None, defaults to NUM_KNOTS constant
    channel : int, optional
        The channel to process. Used to determine if other arrays
        need to be reversed given the direction of increasing
        wavelength down the detector in MIRIFULONG.

    Returns
    -------
    bg_fit : ndarray
        The fitted background.
    bgindx: ndarray
        The location of the knots.
    """
    # first get the weighted pixel fraction
    weighted_pix_frac = (weights > 1e-05).sum() / flux.shape[0]

    # define number of knots using fringe freq, want 1 knot per period
    if ffreq is not None:
        log.debug(f"fit_1d_background_complex: knot positions for {ffreq} cm-1")
        nknots = int((np.amax(wavenum) - np.amin(wavenum)) / (ffreq))

    else:
        log.debug(f"fit_1d_background_complex: using num_knots={NUM_KNOTS}")
        nknots = int((flux.shape[0] / 1024) * NUM_KNOTS)

    log.debug(f"fit_1d_background_complex: number of knots = {nknots}")

    # recale wavenums to around 1 for bayesicfitting
    factor = np.amin(wavenum)
    wavenum_scaled = wavenum.copy() / factor

    # get number of fringe periods in array
    nper = (np.amax(wavenum) - np.amin(wavenum)) // ffreq
    log.debug(f"fit_1d_background_complex: column is {nper} fringe periods")

    # now reduce by the weighted pixel fraction to see how many can be fitted
    nper_cor = int(nper * weighted_pix_frac)
    log.debug(f"fit_1d_background_complex: column has {nper_cor} weighted fringe periods")

    # require at least 5 sine periods to fit
    if nper < 5:
        log.info(" not enough weighted data, no fit performed")
        return flux.copy(), np.zeros(flux.shape[0]), None

    bgindx = make_knots(flux.copy(), int(nknots), weights=weights.copy())
    bgknots = wavenum_scaled[bgindx].astype(float)

    # Reverse (and clip) the fit data as scipy/astropy need monotone increasing data for SW detector
    if channel == 3 or channel == 4:
        t = bgknots[1:-1]
        x = wavenum_scaled
        y = flux
        w = weights
    elif channel == 1 or channel == 2:
        t = bgknots[::-1][1:-1]
        x = wavenum_scaled[::-1]
        y = flux[::-1]
        w = weights[::-1]
    else:
        raise ValueError("channel not in 1-4")

    # Fit the spline
    if ffreq > 1.5:
        bg_model = spline_fitter(x, y, w, t, 2, reject_outliers=True)
    else:
        # robust fitting causing problems for fringe 2 in channels 3 and 4,
        # just use the fitter class
        bg_model = spline_fitter(x, y, w, t, 1, reject_outliers=False)

    # fit the background
    bg_fit = bg_model(wavenum_scaled)
    bg_fit *= np.where(weights.copy() > 1e-07, 1, 1e-08)

    # linearly interpolate over the feature gaps if possible, stops issues later
    try:
        nz, z = interp_helper(weights)
        bg_fit[nz] = np.interp(z(nz), z(~nz), bg_fit[~nz])
        del nz, z
    except ValueError:
        pass

    return bg_fit, bgindx


def fit_1d_fringes_bayes_evidence(
    res_fringes, weights, wavenum, ffreq, dffreq, max_nfringes, pgram_res, col_snr2
):
    """
    Fit the residual fringe signal.

    Takes an input 1D array of residual fringes and fits using the
    supplied mode in the BayesicFitting package.

    Parameters
    ----------
    res_fringes : ndarray
        The 1D array with residual fringes.
    weights : ndarray
        The 1D array of weights
    wavenum : ndarray
        The 1D array of wavenum.
    ffreq : float
        The central scan frequency
    dffreq : float
        The one-sided interval of scan frequencies.
    max_nfringes : int
        The maximum number of fringes to check.
    pgram_res : float
        Resolution of the periodogram scan in cm-1.
    col_snr2 : ndarray
        Location of pixels with sufficient SNR to fit.

    Returns
    -------
    res_fringe_fit : ndarray
        The residual fringe fit data.
    """
    # initialize output to none
    res_fringe_fit = None

    # get the number of weighted pixels
    weighted_pix_num = (weights > 1e-05).sum()

    # get scan res
    res = np.around((2 * dffreq) / pgram_res).astype(int)
    log.debug(f"fit_1d_fringes_bayes: scan res = {res}")

    factor = np.amin(wavenum)
    wavenum = wavenum.copy() / factor
    ffreq = ffreq / factor
    dffreq = dffreq / factor

    # setup frequencies to scan
    freq = np.linspace(ffreq - dffreq, ffreq + dffreq, res)

    # handle out of slice pixels
    res_fringes = np.nan_to_num(res_fringes)
    res_fringes[res_fringes == np.inf] = 0
    res_fringes[res_fringes == -np.inf] = 0

    # initialise some parameters
    res_fringes_proc = res_fringes.copy()
    nfringes = 0
    keep_dict = {}
    fitted_frequencies = []

    # get the initial evidence from ConstantModel
    sdml = ConstantModel(values=1.0)
    sftr = Fitter(wavenum, sdml)
    _ = sftr.fit(res_fringes, weights=weights)
    evidence1 = sftr.getEvidence(limits=[-3, 10], noiseLimits=[0.001, 10])
    log.debug(f"fit_1d_fringes_bayes_evidence: Initial Evidence: {evidence1}")

    for f in np.arange(max_nfringes):
        log.debug(f"Starting fringe {f + 1}")

        # get the scan arrays
        weights *= col_snr2
        res_fringe_scan = res_fringes_proc[np.where(weights > 1e-05)]
        wavenum_scan = wavenum[np.where(weights > 1e-05)]

        # use a Lomb-Scargle periodogram to get PSD and identify the strongest frequency
        log.debug("fit_1d_fringes_bayes_evidence: get the periodogram")
        pgram = LombScargle(wavenum_scan[::-1], res_fringe_scan[::-1]).power(1 / freq)

        log.debug(
            "fit_1d_fringes_bayes_evidence: get the most significant frequency in the periodogram"
        )
        peak = np.argmax(pgram)
        freqs = 1.0 / freq[peak]

        # fix the most significant frequency in the fixed dict that is passed to fitter
        keep_ind = nfringes * 3
        keep_dict[keep_ind] = freqs

        log.debug(
            f"fit_1d_fringes_bayes_evidence: creating multisine model of {nfringes + 1} freqs"
        )
        mdl = multi_sine(nfringes + 1)

        # fit the multi-sine model and get evidence
        if ffreq * factor > 1.5:
            fitter = LevenbergMarquardtFitter(wavenum, mdl, verbose=0, keep=keep_dict)
            ftr = RobustShell(fitter, domain=10)
            try:
                pars = ftr.fit(res_fringes, weights=weights)

                # free the parameters and refit
                mdl = multi_sine(nfringes + 1)
                mdl.parameters = pars
                fitter = LevenbergMarquardtFitter(wavenum, mdl, verbose=0)
                ftr = RobustShell(fitter, domain=10)
                pars = ftr.fit(res_fringes, weights=weights)

                # try to get evidence (may fail for large component
                # fits to noisy data, set to very negative value)
                try:
                    evidence2 = fitter.getEvidence(limits=[-3, 10], noiseLimits=[0.001, 10])
                except ValueError:
                    evidence2 = -1e9
            except Exception:
                evidence2 = -1e9

        else:
            fitter = LevenbergMarquardtFitter(wavenum, mdl, verbose=0, keep=keep_dict)
            try:
                pars = fitter.fit(res_fringes, weights=weights)

                # free the parameters and refit
                mdl = multi_sine(nfringes + 1)
                mdl.parameters = pars
                fitter = LevenbergMarquardtFitter(wavenum, mdl, verbose=0)
                pars = fitter.fit(res_fringes, weights=weights)

                # try to get evidence (may fail for large component
                # fits to noisy data, set to very negative value)
                try:
                    evidence2 = fitter.getEvidence(limits=[-3, 10], noiseLimits=[0.001, 10])
                except ValueError:
                    evidence2 = -1e9
            except Exception:
                evidence2 = -1e9

        log.debug(
            f"fit_1d_fringes_bayes_evidence: nfringe={nfringes + 1} "
            f"ev={evidence2} chi={fitter.chisq}"
        )

        bayes_factor = evidence2 - evidence1
        log.debug(f"fit_1d_fringes_bayes_evidence: bayes factor={bayes_factor}")
        if bayes_factor > 1:
            # strong evidence threshold (log(bayes factor)>1, Kass and Raftery 1995)
            evidence1 = evidence2
            best_mdl = mdl.copy()
            fitted_frequencies.append(freqs)
            log.debug(
                f"fit_1d_fringes_bayes_evidence: strong evidence for nfringes={nfringes + 1} "
            )
        else:
            log.debug(f"fit_1d_fringes_bayes_evidence: no evidence for nfringes={nfringes + 1}")
            break

        # subtract the fringes for this frequency
        res_fringe_fit = best_mdl(wavenum)
        res_fringes_proc = res_fringes.copy() - res_fringe_fit
        nfringes += 1

    log.debug(f"fit_1d_fringes_bayes_evidence: optimal={nfringes} fringes")

    # create outputs to return
    fitted_frequencies = (1 / np.asarray(fitted_frequencies)) * factor
    peak_freq = fitted_frequencies[0]
    freq_min = np.amin(fitted_frequencies)
    freq_max = np.amax(fitted_frequencies)

    return res_fringe_fit, weighted_pix_num, nfringes, peak_freq, freq_min, freq_max


def make_knots(flux, nknots=20, weights=None):
    """
    Define knot positions for piecewise models.

    This function simply splits the array into sections. It does
    NOT take into account the shape of the data.

    Parameters
    ----------
    flux : ndarray
        The flux array or any array of the same dimension.
    nknots : int, optional
        The number of knots to create (excluding 0 and 1023).
    weights : ndarray or None, optional
        Optionally supply a weights array. This will be used to
        add knots at the edge of bad pixels or features.

    Returns
    -------
    knot_idx : ndarray
        The indices of the knots.
    """
    log.debug(f"make_knots: creating {nknots} knots on flux array")

    # handle nans or infs that may exist
    flux = np.nan_to_num(flux, posinf=1e-08, neginf=1e-08)
    flux[flux < 0] = 1e-08

    if weights is not None:
        weights = np.nan_to_num(weights, posinf=1e-08, neginf=1e-08)
        weights[weights < 0] = 1e-08

    # create an array of indices
    npoints = flux.shape[0]

    # kstep is the number of points / number of knots
    knot_step = npoints / nknots

    # define an initial knot index array
    init_knot_idx = np.zeros(nknots + 1)
    for n in range(nknots):
        init_knot_idx[n] = round(n * knot_step)
    init_knot_idx[-1] = npoints - 1

    # get difference between the indices
    knot_split = np.ediff1d(init_knot_idx)

    # the last diff will sometimes be different than the others
    last_split = knot_split[-1]

    # create new knot array with knots at 0 and 1023, then initial and final splits = last_split/2
    knot_idx = np.ones(nknots + 2)
    for n in range(nknots):
        knot_idx[n + 1] = (n * knot_step) + (last_split / 2)
    knot_idx[0] = 0
    knot_idx[-1] = npoints - 1

    # if the weights array is supplied, determine the edges of good data and set knots there
    if weights is not None:
        log.debug("make_knots: adding knots at edges of bad pixels in weights array")

        # if there are bad pixels in the flux array with flux~0,
        # add these to weights array if not already there
        weights *= (flux > 1e-03).astype(int)

        # use a two-point difference method
        weights_diff = np.ediff1d(weights)

        # set edges where diff should be almost equal to the largest of the two datapoints used
        # iterate over the diffs and compare to the datapoints
        edges_idx_list = []
        for n, wd in enumerate(weights_diff):
            # get the data points used for the diff
            datapoints = np.array([weights[n], weights[n + 1]])

            # get the value and index of the larges
            largest = np.amax(datapoints)
            largest_idx = np.argmax(datapoints)

            # we don't need knots in the bad pixels so ignore these
            if largest > 1e-03:
                # check if the absolute values are almost equal
                if math.isclose(largest, np.abs(wd), rel_tol=1e-01):
                    # if so, set the index and adjust depending on whether the
                    # first or second datapoint is the largest
                    idx = n + largest_idx

                    # check if this is right next to another index already defined
                    # causes problems in fitting, minimal difference
                    if (idx - 1 in knot_idx) | (idx + 1 in knot_idx):
                        pass
                    else:
                        # append to the index list
                        edges_idx_list.append(idx)

                else:
                    pass

        # convert the list to array, add to the knot_idx array, remove duplicates and sort
        edges_idx = np.asarray(edges_idx_list)
        knot_idx = np.sort(np.concatenate((knot_idx, edges_idx), axis=0), axis=0)
        knot_idx = np.unique(knot_idx.astype(int))

    return knot_idx.astype(int)


# The below functions were added to enable residual fringe correction
# in 1D extracted data.

# Define some constants describing the two central fringe frequencies
# (primary fringe, and dichroic fringe) and a range around them to search for residual fringes,
# along with the maximum number of fringes to fit and the max amplitude allowed for those fringes.
FFREQ_1d = [2.9, 0.4]
DFFREQ_1d = [1.5, 0.15]
MAX_NFRINGES_1d = [10, 15]
MAXAMP_1d = 0.2


def fit_1d_background_complex_1d(flux, weights, wavenum, ffreq=None):
    """
    Fit the background signal using a piecewise spline of n knots.

    Note that this will also try to identify obvious emission lines and
    flag them so they aren't considered in the fitting.

    Parameters
    ----------
    flux : ndarray
        The 1D array of fluxes.
    weights : ndarray
        The 1D array of weights.
    wavenum : ndarray
        The 1D array of wavenum.
    ffreq : float or None, optional
        The expected fringe frequency, used to determine number of knots.
        If None, defaults to NUM_KNOTS constant.

    Returns
    -------
    bg_fit : ndarray
        The fitted background.
    bgindx : ndarray
        The location of the knots.
    """
    # first get the weighted pixel fraction
    weighted_pix_frac = (weights > 1e-05).sum() / flux.shape[0]

    # define number of knots using fringe freq, want 1 knot per period
    if ffreq is not None:
        log.debug(f"fit_1d_background_complex: knot positions for {ffreq} cm-1")
        nknots = int((np.amax(wavenum) - np.amin(wavenum)) / (ffreq))

    else:
        log.debug(f"fit_1d_background_complex: using num_knots={NUM_KNOTS}")
        nknots = int((flux.shape[0] / 1024) * NUM_KNOTS)

    log.debug(f"fit_1d_background_complex: number of knots = {nknots}")

    # recale wavenums to around 1 for bayesicfitting
    factor = np.amin(wavenum)
    wavenum_scaled = wavenum.copy() / factor

    # get number of fringe periods in array
    nper = (np.amax(wavenum) - np.amin(wavenum)) // ffreq
    log.debug(f"fit_1d_background_complex: column is {nper} fringe periods")

    # now reduce by the weighted pixel fraction to see how many can be fitted
    nper_cor = int(nper * weighted_pix_frac)
    log.debug(f"fit_1d_background_complex: column has {nper_cor} weighted fringe periods")

    # require at least 5 sine periods to fit
    if nper < 5:
        log.info(" not enough weighted data, no fit performed")
        return flux.copy(), np.zeros(flux.shape[0]), None

    bgindx = make_knots(flux.copy(), int(nknots), weights=weights.copy())
    bgknots = wavenum_scaled[bgindx].astype(float)

    # Reverse (and clip) the fit data as scipy/astropy need monotone increasing data for SW detector
    t = bgknots[::-1][1:-1]
    x = wavenum_scaled[::-1]
    y = flux[::-1]
    w = weights[::-1]

    # Fit the spline
    if ffreq > 1.5:
        bg_model = spline_fitter(x, y, w, t, 2, reject_outliers=True)
    else:
        # robust fitting causing problems for fringe 2, change to just using fitter there
        bg_model = spline_fitter(x, y, w, t, 1, reject_outliers=False)

    # fit the background
    bg_fit = bg_model(wavenum_scaled)
    bg_fit *= np.where(weights.copy() > 1e-07, 1, 1e-08)

    # linearly interpolate over the feature gaps if possible, stops issues later
    try:
        nz, z = interp_helper(weights)
        bg_fit[nz] = np.interp(z(nz), z(~nz), bg_fit[~nz])
        del nz, z
    except ValueError:
        pass

    return bg_fit, bgindx


def fit_1d_fringes_bayes_evidence_1d(
    res_fringes, weights, wavenum, ffreq, dffreq, max_nfringes, pgram_res
):
    """
    Fit the residual fringe signal in 1D.

    Takes an input 1D array of residual fringes and fits them using
    the supplied mode in the BayesicFitting package.

    Parameters
    ----------
    res_fringes : ndarray
        The 1D array with residual fringes.
    weights : ndarray
        The 1D array of weights.
    wavenum : ndarray
        The 1D array of wavenum.
    ffreq : float
        The central scan frequency.
    dffreq :  float
        The one-sided interval of scan frequencies.
    max_nfringes : int
        The maximum number of fringes to check.
    pgram_res : float
        Resolution of the periodogram scan in cm-1.

    Returns
    -------
    res_fringe_fit : ndarray
        The residual fringe fit data.
    """
    # initialize output to none
    res_fringe_fit = None

    # get the number of weighted pixels
    weighted_pix_num = (weights > 1e-05).sum()

    # get scan res
    res = np.around((2 * dffreq) / pgram_res).astype(int)

    factor = np.amin(wavenum)
    wavenum = wavenum.copy() / factor
    ffreq = ffreq / factor
    dffreq = dffreq / factor

    # setup frequencies to scan
    freq = np.linspace(ffreq - dffreq, ffreq + dffreq, res)

    # handle out of slice pixels
    res_fringes = np.nan_to_num(res_fringes)

    # initialise some parameters
    res_fringes_proc = res_fringes.copy()
    nfringes = 0
    keep_dict = {}
    fitted_frequencies = []

    # get the initial evidence from ConstantModel
    sdml = ConstantModel(values=1.0)
    sftr = Fitter(wavenum, sdml)
    _ = sftr.fit(res_fringes, weights=weights)
    evidence1 = sftr.getEvidence(limits=[-2, 1000], noiseLimits=[0.001, 1])

    for _ in range(max_nfringes):
        # get the scan arrays
        res_fringe_scan = res_fringes_proc[np.where(weights > 1e-05)]
        wavenum_scan = wavenum[np.where(weights > 1e-05)]

        # use a Lomb-Scargle periodogram to get PSD and identify the strongest frequency
        pgram = LombScargle(wavenum_scan[::-1], res_fringe_scan[::-1]).power(1 / freq)

        peak = np.argmax(pgram)
        freqs = 1.0 / freq[peak]

        # fix the most significant frequency in the fixed dict that is passed to fitter
        keep_ind = nfringes * 3
        keep_dict[keep_ind] = freqs

        mdl = multi_sine(nfringes + 1)

        # fit the multi-sine model and get evidence
        fitter = LevenbergMarquardtFitter(wavenum, mdl, verbose=0, keep=keep_dict)
        ftr = RobustShell(fitter, domain=10)
        try:
            pars = ftr.fit(res_fringes, weights=weights)

            # free the parameters and refit
            mdl = multi_sine(nfringes + 1)
            mdl.parameters = pars
            fitter = LevenbergMarquardtFitter(wavenum, mdl, verbose=0)
            ftr = RobustShell(fitter, domain=10)
            pars = ftr.fit(res_fringes, weights=weights)

            # try to get evidence (may fail for large component
            # fits to noisy data, set to very negative value)
            try:
                evidence2 = fitter.getEvidence(limits=[-2, 1000], noiseLimits=[0.001, 1])
            except ValueError:
                evidence2 = -1e9
        except RuntimeError:
            evidence2 = -1e9

        bayes_factor = evidence2 - evidence1
        if bayes_factor > 1:  # strong evidence thresh (log(bayes factor)>1, Kass and Raftery 1995)
            evidence1 = evidence2
            best_mdl = mdl.copy()
            fitted_frequencies.append(freqs)
        else:
            break

        # subtract the fringes for this frequency
        res_fringe_fit = best_mdl(wavenum)
        res_fringes_proc = res_fringes.copy() - res_fringe_fit
        nfringes += 1

    # create outputs to return
    fitted_frequencies = (1 / np.asarray(fitted_frequencies)) * factor
    peak_freq = fitted_frequencies[0]
    freq_min = np.amin(fitted_frequencies)
    freq_max = np.amax(fitted_frequencies)

    return res_fringe_fit, weighted_pix_num, nfringes, peak_freq, freq_min, freq_max


def fit_residual_fringes_1d(flux, wavelength, channel=1, dichroic_only=False, max_amp=None):
    """
    Fit residual fringes in 1D.

    Parameters
    ----------
    flux : ndarray
        The 1D array of fluxes.
    wavelength : ndarray
        The 1D array of wavelengths.
    channel : int, optional
        The MRS spectral channel.
    dichroic_only : bool, optional
        Fit only dichroic fringes.
    max_amp : ndarray, optional
        The maximum amplitude array.

    Returns
    -------
    output : ndarray
        Modified version of input flux array.
    """
    # Restrict to just the non-zero positive fluxes
    indx = np.where(flux > 0)
    useflux = flux[indx]
    usewave = wavelength[indx]

    wavenum = 10000.0 / usewave

    weights = useflux / np.nanmedian(useflux)
    weights[weights == np.inf] = 0
    weights[np.isnan(weights)] = 0

    # Zero out any weights longward of 27.6 microns as the calibration is too uncertain
    # and can bias the fringe finding
    weights[usewave > 27.6] = 0

    # get the maxamp of the fringes
    if max_amp is None:
        max_amp = np.ones(useflux.shape) * MAXAMP_1d

    # find spectral features (env is spline fit of troughs and peaks)
    # smooth the data slightly first to avoid noisy broad lines being missed
    env, l_x, l_y, _, _, _ = fit_envelope(np.arange(useflux.shape[0]), useflux)
    mod = np.abs(useflux / env) - 1

    # given signal in mod find location of lines > col_max_amp * 2
    weight_factors = find_lines(mod, max_amp * 2)
    weights_feat = weights * weight_factors

    if dichroic_only is True:
        if channel not in [3, 4]:
            raise ValueError(
                "Dichroic fringe should only be removed from channels 3 and 4, stopping!"
            )

        ffreq_vals = [FFREQ_1d[1]]
        dffreq_vals = [DFFREQ_1d[1]]
        max_nfringes_vals = [MAX_NFRINGES_1d[1]]

    else:
        # check the channel and remove second fringe for channels 1 and 2
        if channel == 1 or channel == 2:
            ffreq_vals = [FFREQ_1d[0]]
            dffreq_vals = [DFFREQ_1d[0]]
            max_nfringes_vals = [MAX_NFRINGES_1d[0]]

        else:
            ffreq_vals = FFREQ_1d
            dffreq_vals = DFFREQ_1d
            max_nfringes_vals = MAX_NFRINGES_1d

    # BayesicFitting doesn't like 0s at data or weight array edges so set to small value
    useflux[useflux <= 0] = 1e-08
    weights_feat[weights_feat <= 0] = 1e-08

    # check for off-slice pixels and send to be filled with interpolated/extrapolated wnums
    # to stop BayesicFitting crashing, will not be fitted anyway
    found_bad = np.logical_or(np.isnan(wavenum), np.isinf(wavenum))
    num_bad = len(np.where(found_bad)[0])

    if num_bad > 0:
        wavenum[found_bad] = 0
        wavenum = fill_wavenumbers(wavenum)

    # do the processing
    proc_arr = [useflux.copy()]

    for m, proc_data in enumerate(proc_arr):
        for n, ffreq in enumerate(ffreq_vals):
            bg_fit, bgindx = fit_1d_background_complex_1d(
                proc_data, weights_feat, wavenum, ffreq=ffreq
            )

            # get the residual fringes as fraction of signal
            res_fringes = np.divide(
                proc_data, bg_fit, out=np.zeros_like(proc_data), where=bg_fit != 0
            )
            res_fringes = np.subtract(res_fringes, 1, where=res_fringes != 0)
            res_fringes *= np.where(weights > 1e-07, 1, 1e-08)

            # fit the residual fringes
            res_fringe_fit, wpix_num, opt_nfringes, peak_freq, freq_min, freq_max = (
                fit_1d_fringes_bayes_evidence_1d(
                    res_fringes,
                    weights_feat,
                    wavenum,
                    ffreq,
                    dffreq_vals[n],
                    max_nfringes_vals[n],
                    0.001,
                )
            )

            # check for fit blowing up, reset rfc fit to 0, raise a flag
            res_fringe_fit, res_fringe_fit_flag = check_res_fringes(res_fringe_fit, max_amp)

            # correct for residual fringes
            _, _, _, env, u_x, u_y = fit_envelope(
                np.arange(res_fringe_fit.shape[0]), res_fringe_fit
            )
            rfc_factors = 1 / (res_fringe_fit * (weights > 1e-05).astype(int) + 1)
            proc_data *= rfc_factors

            # handle nans or infs that may exist
            proc_data[proc_data == np.inf] = 0
            proc_data = np.nan_to_num(proc_data)
            proc_arr[m] = proc_data

    # Embed output back in a full-size array
    output = flux.copy()
    output[indx] = proc_arr[0]

    return output

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

from astropy.modeling.models import Spline1D
from astropy.modeling.fitting import SplineExactKnotsFitter

from .fitter import ChiSqOutlierRejectionFitter

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
#

# hard coded parameters, have been selected based on testing but can be changed
NUM_KNOTS = 80  # number of knots for bkg model if no other info provided


def find_nearest(array, value):
    """ Utility function to find the index of pixel with value in 'array' nearest to 'value'.

    Used by det_pixel_trace function
    """
    idx = (np.abs(array - value)).argmin()
    return idx


def slice_info(slice_map, c):
    """ Function to take slice map and channel and find pixels in the slice and xrange of each slice
    """

    slice_inventory = np.unique(slice_map)
    slices_in_band = slice_inventory[np.where((slice_inventory >= 100 * c)
                                              & (slice_inventory < 100 * (c + 1)))]

    log.info("Number of slices in band {} ".format(slices_in_band.shape[0]))
    slice_x_ranges = np.zeros((slices_in_band.shape[0], 3), dtype=int)
    all_slice_masks = np.zeros((slices_in_band.shape[0], slice_map.shape[0], slice_map.shape[1]))
    for n, s in enumerate(slices_in_band):
        # create a mask of the slice
        pixels = np.where(slice_map == s)
        slice = np.zeros(slice_map.shape)

        slice[pixels] = 1

        # add this to the all_slice_mask array
        all_slice_masks[n] = slice

        # get the indices at the start and end of the slice
        collapsed_slice = np.sum(slice, axis=0)
        indices = np.where(collapsed_slice[:-1] != collapsed_slice[1:])[0]
        slice_x_ranges[n, 0], slice_x_ranges[n, 1], slice_x_ranges[n, 2] = int(s), \
            int(np.amin(indices)), int(np.amax(indices) + 1)

        log.debug("For slice {} x ranges of slices region {}, {}".
                  format(slice_x_ranges[n, 0], slice_x_ranges[n, 1], slice_x_ranges[n, 2]))

    log.debug("Min and max x pixel values of all slices in channel {} {}".
              format(np.amin(slice_x_ranges[:, 1]), np.amax(slice_x_ranges[:, 2])))

    xrange_channel = np.zeros(2)
    xrange_channel[0] = np.amin(slice_x_ranges[:, 1])
    xrange_channel[1] = np.amax(slice_x_ranges[:, 2])
    result = (slices_in_band, xrange_channel, slice_x_ranges, all_slice_masks)
    return result


def fill_wavenumbers(wnums):
    """
    Function to take a wavenumber array with missing values (e.g., columns with on-slice and
    off-slice pixels), fit the good points using a polynomial, then run use the coefficients
    to estimate wavenumbers on the off-slice pixels.
    Note that these new values are physically meaningless but having them in the wavenum array
    stops the BayesicFitting package from crashing with a LinAlgErr.

    :Parameters:

    wnums: numpy array, required
        the wavenumber array

    :Returns:

    wnums_filled: numpy array
        the wavenumber array with off-slice pixels filled

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


def fit_quality_stats(stats):
    """Get simple statistics for the fits

    :Parameters:

    fit_stats: np.array
        the fringe contrast per fit

    :Returns:

   median, stddev, max of stat numpy array

    """
    return np.mean(stats), np.median(stats), np.std(stats),  np.amax(stats)


def multi_sine(n):
    """
    Return a model composed of n sines

    :Parameters:

    n: int, required
        number of sines

    :Returns:

    mdl, BayesFitting model
        the model composed of n sines
    """

    # make the first sine
    mdl = SineModel()

    # make a copy
    model = mdl.copy()

    # add the copy n-1 times
    for i in range(1, n):
        mdl.addModel(model.copy())

    # clean
    del model

    return mdl


def fit_envelope(wavenum, signal):
    """ Fit the upper and lower envelope of signal using a univariate spline

    :param wavenum:
    :param signal:
    :return:
    """

    # Detect troughs and mark their location. Define endpoints
    l_x = [wavenum[0]]
    l_y = [signal[0]]
    u_x = [wavenum[0]]
    u_y = [signal[0]]

    for k in np.arange(1, len(signal) - 1):
        if ((np.sign(signal[k] - signal[k - 1]) == -1) and
                ((np.sign(signal[k] - signal[k + 1])) == -1)):
            l_x.append(wavenum[k])
            l_y.append(signal[k])
        if ((np.sign(signal[k] - signal[k - 1]) == 1) and
                ((np.sign(signal[k] - signal[k + 1])) == 1)):
            u_x.append(wavenum[k])
            u_y.append(signal[k])

    # Append the last value of (s) to the interpolating values. This forces the model to use the same ending point
    l_x.append(wavenum[-1])
    l_y.append(signal[-1])
    u_x.append(wavenum[-1])
    u_y.append(signal[-1])

    # fit a model
    pcl = pchip(l_x, l_y)
    pcu = pchip(u_x, u_y)

    return pcl(wavenum), l_x, l_y, pcu(wavenum), u_x, u_y


def find_lines_resfringe(signal, max_amp):
    """
    *** Replaced with find_lines below. This version does not include some of the
    feature finding functionality***

    Take signal and max amp array, determine location of spectral
    features with amplitudes greater than max amp

    :param signal:
    :param max_amp:
    :return:
    """

    r_x = np.arange(signal.shape[0] - 1)

    # setup the output arrays
    signal_check = signal.copy()
    weights_factors = np.ones(signal.shape[0])

    # Detect peaks

    u_y, u_x, l_y, l_x = [], [], [], []

    for x in r_x:
        if ((np.sign(signal_check[x] - signal_check[x - 1]) == 1) and
                (np.sign(signal_check[x] - signal_check[x + 1]) == 1)):
            u_y.append(signal_check[x])
            u_x.append(x)

        if ((np.sign(signal_check[x] - signal_check[x - 1]) == -1) and
                (np.sign(signal_check[x] - signal_check[x + 1]) == -1)):
            l_y.append(signal[x])
            l_x.append(x)

    log.debug("find_lines: Found {} peaks   {} troughs".format(len(u_x), len(l_x)))
    weights_factors[signal_check > max_amp] = 0
    return weights_factors


def find_lines(signal, max_amp):
    """
    Take signal and max amp array, determine location of spectral
    features with amplitudes greater than max amp

    :param signal:
    :param max_amp:
    :return:
    """

    r_x = np.arange(signal.shape[0] - 1)

    # setup the output arrays
    signal_check = signal.copy()
    weights_factors = np.ones(signal.shape[0])

    # Detect peaks
    u_y, u_x, l_y, l_x = [], [], [], []

    for x in r_x:
        if (np.sign(signal_check[x] - signal_check[x - 1]) == 1) and \
                (np.sign(signal_check[x] - signal_check[x + 1]) == 1):
            u_y.append(signal_check[x])
            u_x.append(x)

        if (np.sign(signal_check[x] - signal_check[x - 1]) == -1) and \
                (np.sign(signal_check[x] - signal_check[x + 1]) == -1):
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
                log.debug("find_lines:  checking ind {}".format(xp))

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

    log.debug("find_lines: Found {} peaks   {} troughs".format(len(u_x), len(l_x)))
    weights_factors[signal_check > max_amp * 2] = 0  # catch any remaining
    # weights_factors[signal_check > np.amax(max_amp)] = 0

    return weights_factors


def check_res_fringes(res_fringe_fit, max_amp):
    """
    Check for regions where res fringe fit runs away (greater than max amp),
    set the beat where this happens to 0 to avoid  making the fringes worse

    :Parameters:

    res_fringe_fit:  numpy array, required
        the residual fringe fit

    max_amp:  numpy array, required
        the maximum amplitude array

    :Returns:

    res_fringe_fit: numpy array
        the residual fringe fit with exploding fit regions removed

    flats: numpy array
        flags where the fit was rejected

    """

    flags = np.zeros(res_fringe_fit.shape[0])

    # get fit envelope
    npix = np.arange(res_fringe_fit.shape[0])
    lenv_fit, _, _, uenv_fit, _, _ = fit_envelope(npix, res_fringe_fit)

    # get the indices of the nodes (where uenv slope goes from negative to positive), add 0 and 1023
    node_ind = [0]
    for k in np.arange(1, len(uenv_fit) - 1):
        if (np.sign(uenv_fit[k] - uenv_fit[k - 1]) == -1) and ((np.sign(uenv_fit[k] - uenv_fit[k + 1])) == -1):
            node_ind.append(k)
    node_ind.append(res_fringe_fit.shape[0] - 1)
    node_ind = np.asarray(node_ind)
    log.debug("check_res_fringes: found {} nodes".format(len(node_ind)))

    # find where res_fringes goes above max_amp
    runaway_rfc = np.argwhere((np.abs(lenv_fit) + np.abs(uenv_fit)) > (max_amp * 2))

    # check which signal env the blow ups are located in and set to 1, and set a flag array
    if len(runaway_rfc) > 0:
        log.debug("check_res_fringes: {} data points exceed threhold".format(len(runaway_rfc)))
        log.debug("check_res_fringes: resetting fits to related beats")
        for i in runaway_rfc:
            # find where the index is compared to the nodes
            node_loc = np.searchsorted(node_ind, i)

            # set the res_fringes between the nodes to 1
            lind = node_ind[node_loc - 1]
            uind = node_ind[node_loc]
            res_fringe_fit[lind[0]:uind[0]] = 0
            flags[lind[0]:uind[0]] = 1  # set flag to 1 for reject fit region

    return res_fringe_fit, flags


def interp_helper(mask):
    """Helper function to for interpolating in feature gaps.

    :Parameters:

    mask:  numpy array, required
        the 1D mask array (weights)

    :Returns:

        - logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices to 'equivalent' indices

    """
    return mask < 1e-05, lambda z: z.nonzero()[0]


def fit_1d_background_complex(flux, weights, wavenum, order=2, ffreq=None, channel=1, test=False):
    """Fit the background signal using a pieceweise spline of n knots. Note that this will also try to identify
    obvious emission lines and flag them so they aren't considered in the fitting.

    :Parameters:

    flux:  numpy array, required
        the 1D array of fluxes

    weights: numpy array, required
        the 1D array of weights

    wavenum: numpy array, required
        the 1D array of wavenum

    order: int, optional, default=2
        the order of the Splines model

    ffreq: float, optional, default=None
        the expected fringe frequency, used to determine number of knots. If None,
        defaults to NUM_KNOTS constant

    channel: int, optional, default=1
        the channel processed. used to determine if other arrays need to be reversed given the direction of increasing
        wavelength down the detector in MIRIFULONG

    :Returns:

    bg_fit: numpy array
        the fitted background

    bgindx: numpy array
        the location of the knots

    fitter: BayesicFitting object
        fitter object, mainly used for testing

    """

    # first get the weighted pixel fraction
    weighted_pix_frac = (weights > 1e-05).sum() / flux.shape[0]

    # define number of knots using fringe freq, want 1 knot per period
    if ffreq is not None:
        log.debug("fit_1d_background_complex: knot positions for {} cm-1".format(ffreq))
        nknots = int((np.amax(wavenum) - np.amin(wavenum)) / (ffreq))

    else:
        log.debug("fit_1d_background_complex: using num_knots={}".format(NUM_KNOTS))
        nknots = int((flux.shape[0] / 1024) * NUM_KNOTS)

    log.debug("fit_1d_background_complex: number of knots = {}".format(nknots))

    # recale wavenums to around 1 for bayesicfitting
    factor = np.amin(wavenum)
    wavenum_scaled = wavenum.copy() / factor

    # get number of fringe periods in array
    nper = (np.amax(wavenum) - np.amin(wavenum)) // ffreq
    log.debug("fit_1d_background_complex: column is {} fringe periods".format(nper))

    # now reduce by the weighted pixel fraction to see how many can be fitted
    nper_cor = int(nper * weighted_pix_frac)
    log.debug("fit_1d_background_complex: column has {} weighted fringe periods".format(nper_cor))

    # require at least 5 sine periods to fit
    if nper < 5:
        log.info(" not enough weighted data, no fit performed")
        return flux.copy(), np.zeros(flux.shape[0]), None

    bgindx = new_make_knots(flux.copy(), int(nknots), weights=weights.copy())
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
        raise ValueError('channel not in 1-4')

    # Fit the spline
    # robust fitting causing problems for fringe 2 in channels 3 and 4, just use the fitter class
    if ffreq > 1.5:
        spline_model = Spline1D(knots=t, degree=2, bounds=[x[0], x[-1]])
        fitter = SplineExactKnotsFitter()
        robust_fitter = ChiSqOutlierRejectionFitter(fitter)
        bg_model = robust_fitter(spline_model, x, y, weights=w)
    else:
        spline_model = Spline1D(knots=t, degree=1, bounds=[x[0], x[-1]])
        fitter = SplineExactKnotsFitter()
        bg_model = fitter(spline_model, x, y, weights=w)

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


def fit_quality(wavenum, res_fringes, weights, ffreq, dffreq, save_results=False):
    """Determine the post correction fringe residual

    Fit a single sine model to the corrected array to get the post correction fringe residual

    :Parameters:

    wavenum: numpy array
        the wavenum array

    res_fringes: numpy array
        the residual fringe fit data

    weights: numpy array
        the weights array

    ffreq: float, required
        the central scan frequency

    dffreq:  float, required
        the one-sided interval of scan frequencies

    :Returns:

    fringe_res_amp: numpy array
        the post correction fringe residual amplitude

    """
    ffreq, dffreq = 2.8, 0.2

    # fit the residual with a single sine model
    # use a Lomb-Scargle periodogram to get PSD and identify the strongest frequency
    freq = np.linspace(ffreq - dffreq, ffreq + dffreq, 100)

    # handle out of slice pixels
    res_fringes = np.nan_to_num(res_fringes)
    res_fringe_scan = res_fringes[np.where(weights > 1e-05)]
    wavenum_scan = wavenum[np.where(weights > 1e-05)]
    pgram = LombScargle(wavenum_scan[::-1], res_fringe_scan[::-1]).power(1 / freq)
    peak = np.argmax(pgram)
    peak_freq = freq[peak]
    log.debug("fit_quality: strongest frequency is {}".format(peak_freq))

    # create the model
    mdl = SineModel(pars=[0.1, 0.1], fixed={0: 1 / peak_freq})

    fitter = LevenbergMarquardtFitter(wavenum[10:-10], mdl)
    ftr = RobustShell(fitter, domain=10)

    fr_par = ftr.fit(res_fringes[10:-10], weights=weights[10:-10])
    log.debug("fit_quality: best fit pars: {}".format(fr_par))

    if np.abs(fr_par[0]) > np.abs(fr_par[1]):
        contrast = np.abs(round(fr_par[0] * 2, 3))
    else:
        contrast = np.abs(round(fr_par[1] * 2, 3))

    # make data to return for fit quality
    quality = None
    if save_results:
        best_mdl = SineModel(fixed={0: 1 / peak_freq, 1: fr_par[0], 2: fr_par[1]})
        fit = best_mdl.result(wavenum)
        quality = np.array([(10000.0 / wavenum), res_fringes, fit])

    return contrast, quality


def new_fit_1d_fringes_bayes_evidence(res_fringes, weights, wavenum, ffreq, dffreq, min_nfringes, max_nfringes,
                                      pgram_res, col_snr2):

    """Fit the residual fringe signal.- Improved method
    Takes an input 1D array of residual fringes and fits using the supplied mode in the BayesicFitting package:
    :Parameters:
    res_fringes:  numpy array, required
        the 1D array with residual fringes
    weights: numpy array, required
        the 1D array of weights
    ffreq: float, required
        the central scan frequency
    dffreq:  float, required
        the one-sided interval of scan frequencies
    min_nfringes: int, required
        the minimum number of fringes to check
    max_nfringes: int, required
        the maximum number of fringes to check
    pgram_res: float, optional
        resolution of the periodogram scan in cm-1
    wavenum: numpy array, required
        the 1D array of wavenum
    :Returns:
    res_fringe_fit: numpy array
        the residual fringe fit data
    """
    # initialize output to none
    res_fringe_fit = None
    weighted_pix_num = None
    peak_freq = None
    freq_min = None
    freq_max = None

    # get the number of weighted pixels
    weighted_pix_num = (weights > 1e-05).sum()
    # set the maximum array size, always 1024

    # get scan res
    res = np.around((2 * dffreq) / pgram_res).astype(int)
    log.debug("fit_1d_fringes_bayes: scan res = {}".format(res))

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
    best_mdl = None
    fitted_frequencies = []

    # get the initial evidence from ConstantModel
    sdml = ConstantModel(values=1.0)
    sftr = Fitter(wavenum, sdml)
    _ = sftr.fit(res_fringes, weights=weights)
    evidence1 = sftr.getEvidence(limits=[-3, 10], noiseLimits=[0.001, 10])
    log.debug(
        "fit_1d_fringes_bayes_evidence: Initial Evidence: {}".format(evidence1))

    for f in np.arange(max_nfringes):
        log.debug(
            "Starting fringe {}".format(f + 1))

        # get the scan arrays
        weights *= col_snr2
        res_fringe_scan = res_fringes_proc[np.where(weights > 1e-05)]
        wavenum_scan = wavenum[np.where(weights > 1e-05)]

        # use a Lomb-Scargle periodogram to get PSD and identify the strongest frequency
        log.debug("fit_1d_fringes_bayes_evidence: get the periodogram")
        pgram = LombScargle(wavenum_scan[::-1], res_fringe_scan[::-1]).power(1 / freq)

        log.debug("fit_1d_fringes_bayes_evidence: get the most significant frequency in the periodogram")
        peak = np.argmax(pgram)
        freqs = 1. / freq[peak]

        # fix the most significant frequency in the fixed dict that is passed to fitter
        keep_ind = nfringes * 3
        keep_dict[keep_ind] = freqs

        log.debug("fit_1d_fringes_bayes_evidence: creating multisine model of {} freqs".format(nfringes + 1))
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

                # try get evidence (may fail for large component fits to noisy data, set to very negative value
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

                # try get evidence (may fail for large component fits to noisy data, set to very negative value
                try:
                    evidence2 = fitter.getEvidence(limits=[-3, 10], noiseLimits=[0.001, 10])
                except ValueError:
                    evidence2 = -1e9
            except Exception:
                evidence2 = -1e9

        log.debug("fit_1d_fringes_bayes_evidence: nfringe={} ev={} chi={}".format(nfringes, evidence2, fitter.chisq))

        bayes_factor = evidence2 - evidence1
        log.debug(
            "fit_1d_fringes_bayes_evidence: bayes factor={}".format(bayes_factor))
        if bayes_factor > 1:  # strong evidence thresh (log(bayes factor)>1, Kass and Raftery 1995)
            evidence1 = evidence2
            best_mdl = mdl.copy()
            fitted_frequencies.append(freqs)
            log.debug(
                "fit_1d_fringes_bayes_evidence: strong evidence for nfringes={} ".format(nfringes + 1))
        else:
            log.debug(
                "fit_1d_fringes_bayes_evidence: no evidence for nfringes={}".format(nfringes + 1))
            break

        # subtract the fringes for this frequency
        res_fringe_fit = best_mdl(wavenum)
        res_fringes_proc = res_fringes.copy() - res_fringe_fit
        nfringes += 1

    log.debug("fit_1d_fringes_bayes_evidence: optimal={} fringes".format(nfringes))

    # create outputs to return
    fitted_frequencies = (1 / np.asarray(fitted_frequencies)) * factor
    peak_freq = fitted_frequencies[0]
    freq_min = np.amin(fitted_frequencies)
    freq_max = np.amax(fitted_frequencies)

    return res_fringe_fit, weighted_pix_num, nfringes, peak_freq, freq_min, freq_max


def new_make_knots(flux, nknots=20, weights=None):
    """Defines knot positions for piecewise models. This simply splits the array into sections. It does
    NOT take into account the shape of the data.

    :Parameters:

    flux: numpy array, required
        the flux array or any array of the same dimension

    nknots: int, optional, default=20
        the number of knots to create (excluding 0 and 1023)

    weights: numpy array, optional, default=None
        optionally supply a weights array. This will be used to add knots at the edge of bad pixels or features

    :Returns:

    knot_idx, numpy array
        the indices of the knots

    """
    log.debug("new_make_knots: creating {} knots on flux array".format(nknots))

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

    # create new knot array with knots at 0 and 1023, then inital and final splits = last_split/2
    knot_idx = np.ones(nknots + 2)
    for n in range(nknots):
        knot_idx[n + 1] = (n * knot_step) + (last_split / 2)
    knot_idx[0] = 0
    knot_idx[-1] = npoints - 1

    # if the weights array is supplied, determine the edges of good data and set knots there
    if weights is not None:

        log.debug("new_make_knots: adding knots at edges of bad pixels in weights array")

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

# RFC1D additions ======================================

# Define some constants describing the two central fringe frequencies
# (primary fringe, and dichroic fringe) and a range around them to search for residual fringes,
# along with the maximum number of fringes to fit and the max amplitude allowed for those fringes.
FFREQ_1d = [2.9, 0.4]
DFFREQ_1d = [1.5, 0.15]
MAX_NFRINGES_1d = [10, 15]
MAXAMP_1d = 0.2

# functions
def fit_1d_background_complex_1d(flux, weights, wavenum, order=2, ffreq=None, channel=1, test=False):
    """Fit the background signal using a pieceweise spline of n knots. Note that this will also try to identify
    obvious emission lines and flag them so they aren't considered in the fitting.

    Parameters
    ----------
    flux :  numpy array, required
        the 1D array of fluxes

    weights : numpy array, required
        the 1D array of weights

    wavenum : numpy array, required
        the 1D array of wavenum

    order : int, optional, default=2
        the order of the Splines model

    ffreq : float, optional, default=None
        the expected fringe frequency, used to determine number of knots. If None,
        defaults to NUM_KNOTS constant

    channel : int, optional, default=1
        the channel processed. used to determine if other arrays need to be reversed given the direction of increasing
        wavelength down the detector in MIRIFULONG

    Returns
    -------

    bg_fit : numpy array
        the fitted background

    bgindx : numpy array
        the location of the knots

    fitter : BayesicFitting object
        fitter object, mainly used for testing

    """

    # first get the weighted pixel fraction
    weighted_pix_frac = (weights > 1e-05).sum() / flux.shape[0]

    # define number of knots using fringe freq, want 1 knot per period
    if ffreq is not None:
        log.debug("fit_1d_background_complex: knot positions for {} cm-1".format(ffreq))
        nknots = int((np.amax(wavenum) - np.amin(wavenum)) / (ffreq))

    else:
        log.debug("fit_1d_background_complex: using num_knots={}".format(NUM_KNOTS))
        nknots = int((flux.shape[0] / 1024) * NUM_KNOTS)

    log.debug("fit_1d_background_complex: number of knots = {}".format(nknots))

    # recale wavenums to around 1 for bayesicfitting
    factor = np.amin(wavenum)
    wavenum_scaled = wavenum.copy() / factor

    # get number of fringe periods in array
    nper = (np.amax(wavenum) - np.amin(wavenum)) // ffreq
    log.debug("fit_1d_background_complex: column is {} fringe periods".format(nper))

    # now reduce by the weighted pixel fraction to see how many can be fitted
    nper_cor = int(nper * weighted_pix_frac)
    log.debug("fit_1d_background_complex: column has {} weighted fringe periods".format(nper_cor))

    # require at least 5 sine periods to fit
    if nper < 5:
        log.info(" not enough weighted data, no fit performed")
        return flux.copy(), np.zeros(flux.shape[0]), None

    bgindx = new_make_knots(flux.copy(), int(nknots), weights=weights.copy())
    bgknots = wavenum_scaled[bgindx].astype(float)

    # Reverse (and clip) the fit data as scipy/astropy need monotone increasing data for SW detector

    t = bgknots[::-1][1:-1]
    x = wavenum_scaled[::-1]
    y = flux[::-1]
    w = weights[::-1]

    # Fit the spline
    # TODO: robust fitting causing problems for fringe 2, change to just using fitter there
    if ffreq > 1.5:
        spline_model = Spline1D(knots=t, degree=2, bounds=[x[0], x[-1]])
        fitter = SplineExactKnotsFitter()
        robust_fitter = ChiSqOutlierRejectionFitter(fitter)
        bg_model = robust_fitter(spline_model, x, y, weights=w)
    else:
        spline_model = Spline1D(knots=t, degree=1, bounds=[x[0], x[-1]])
        fitter = SplineExactKnotsFitter()
        bg_model = fitter(spline_model, x, y, weights=w)

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


def new_fit_1d_fringes_bayes_evidence_1d(res_fringes, weights, wavenum, ffreq, dffreq, min_nfringes, max_nfringes,
                                         pgram_res):
    """Fit the residual fringe signal.- 1d version
    Takes an input 1D array of residual fringes and fits using the supplied mode in the BayesicFitting package.

    Parameters
    ----------
    res_fringes :  numpy array, required
        the 1D array with residual fringes
    weights : numpy array, required
        the 1D array of weights
    ffreq : float, required
        the central scan frequency
    dffreq :  float, required
        the one-sided interval of scan frequencies
    min_nfringes : int, required
        the minimum number of fringes to check
    max_nfringes : int, required
        the maximum number of fringes to check
    pgram_res : float, optional
        resolution of the periodogram scan in cm-1
    wavenum : numpy array, required
        the 1D array of wavenum

    Returns
    -------
    res_fringe_fit : numpy array
        the residual fringe fit data
    """
    # initialize output to none
    res_fringe_fit = None
    weighted_pix_num = None
    peak_freq = None
    freq_min = None
    freq_max = None

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
    best_mdl = None
    fitted_frequencies = []

    # get the initial evidence from ConstantModel
    sdml = ConstantModel(values=1.0)
    sftr = Fitter(wavenum, sdml)
    _ = sftr.fit(res_fringes, weights=weights)
    evidence1 = sftr.getEvidence(limits=[-2, 1000], noiseLimits=[0.001, 1])

    for n in np.arange(max_nfringes):
        # get the scan arrays
        res_fringe_scan = res_fringes_proc[np.where(weights > 1e-05)]
        wavenum_scan = wavenum[np.where(weights > 1e-05)]

        # use a Lomb-Scargle periodogram to get PSD and identify the strongest frequency
        pgram = LombScargle(wavenum_scan[::-1], res_fringe_scan[::-1]).power(1 / freq)

        peak = np.argmax(pgram)
        freqs = 1. / freq[peak]

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

            # try get evidence (may fail for large component fits to noisy data, set to very negative value
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
    """This is the wrapper function for 1d residual fringe correction.

    Parameters
    ----------
    flux : numpy array, required
        The 1D array of fluxes
    wavelength : numpy array, required
        The 1D array of wavelengths
    channel : integer, optional
        The MRS spectral channel
    dichroic_only : boolean, optional
        Fit only dichroic fringes
    max_amp : numpy array, optional
        The maximum amplitude array

    Returns
    -------
    output : numpy array
        Modified version of input flux array
    """

    # Restrict to just the non-zero fluxes
    indx = np.where(flux != 0)
    useflux = flux[indx]

    wavenum = 10000.0 / wavelength[indx]

    weights = useflux / np.nanmedian(useflux)
    weights[weights == np.inf] = 0
    weights[np.isnan(weights)] = 0

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
            raise ValueError('Dichroic fringe should only be removed from channels 3 and 4, stopping!')

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

            # fit background
            try:
                bg_fit, bgindx = fit_1d_background_complex_1d(proc_data, weights_feat,
                                                              wavenum, ffreq=ffreq, channel=1)
            except Exception as e:
                raise e

            # get the residual fringes as fraction of signal
            res_fringes = np.divide(proc_data, bg_fit, out=np.zeros_like(proc_data),
                                    where=bg_fit != 0)
            res_fringes = np.subtract(res_fringes, 1, where=res_fringes != 0)
            res_fringes *= np.where(weights > 1e-07, 1, 1e-08)

            # fit the residual fringes
            try:
                res_fringe_fit, wpix_num, opt_nfringes, peak_freq, freq_min, freq_max = new_fit_1d_fringes_bayes_evidence_1d(
                    res_fringes, weights_feat,
                    wavenum, ffreq, dffreq_vals[n], min_nfringes=0,
                    max_nfringes=max_nfringes_vals[n], pgram_res=0.001)

            except Exception as e:
                raise e

            # check for fit blowing up, reset rfc fit to 0, raise a flag
            res_fringe_fit, res_fringe_fit_flag = check_res_fringes(res_fringe_fit, max_amp)

            # correct for residual fringes
            _, _, _, env, u_x, u_y = fit_envelope(np.arange(res_fringe_fit.shape[0]),
                                                  res_fringe_fit)
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

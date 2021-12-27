import numpy as np
import math
import numpy.polynomial.polynomial as poly

from scipy.interpolate import pchip
from astropy.timeseries import LombScargle
from BayesicFitting import SplinesModel
from BayesicFitting import Fitter
from BayesicFitting import SineModel
from BayesicFitting import LevenbergMarquardtFitter
from BayesicFitting import RobustShell

from numpy.linalg.linalg import LinAlgError

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
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
    slice_x_ranges = np.zeros((slices_in_band.shape[0], 3),dtype=int)
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
            int(np.amin(indices)),int(np.amax(indices) + 1)

        log.debug("For slice {} x ranges of slices region {}, {}".
                  format(slice_x_ranges[n,0],slice_x_ranges[n,1],slice_x_ranges[n,2]))

    log.debug("Min and max x pixel values of all slices in channel {} {}".
              format(np.amin(slice_x_ranges[:,1]), np.amax(slice_x_ranges[:,2])))

    xrange_channel = np.zeros(2)
    xrange_channel[0] = np.amin(slice_x_ranges[:,1])
    xrange_channel[1] = np.amax(slice_x_ranges[:,2])
    result = (slices_in_band, xrange_channel, slice_x_ranges,all_slice_masks)
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
    coeffs = poly.polyfit(np.arange(wnums.shape[0])[idx], wnums[idx], 3)
    wnums_filled = poly.polyval(np.arange(wnums.shape[0]), coeffs)

    # keep the original wnums for the on-slice pixels
    wnums_filled[idx] = wnums[idx]

    # clean
    del idx, coeffs

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
        if (np.sign(signal[k] - signal[k - 1]) == -1) and ((np.sign(signal[k] - signal[k + 1])) == -1):
            l_x.append(wavenum[k])
            l_y.append(signal[k])
        if (np.sign(signal[k] - signal[k - 1]) == 1) and ((np.sign(signal[k] - signal[k + 1])) == 1):
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
        if (np.sign(signal_check[x] - signal_check[x - 1]) == 1) and (np.sign(signal_check[x] - signal_check[x + 1]) == 1):
            u_y.append(signal_check[x])
            u_x.append(x)

        if (np.sign(signal_check[x] - signal_check[x - 1]) == -1) and (np.sign(signal_check[x] - signal_check[x + 1]) == -1):
            l_y.append(signal[x])
            l_x.append(x)

    log.debug("find_lines: Found {} peaks   {} troughs".format(len(u_x), len(l_x)))
    weights_factors[signal_check > np.amax(max_amp)] = 0
    return weights_factors


def check_res_fringes(res_fringe_fit, max_amp):
    """
    Check for regions where res fringe fit runs away, set to 0

    :param res_fringes:
    :return:
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
    runaway_rfc = np.argwhere((np.abs(lenv_fit) + np.abs(uenv_fit)) > max_amp * 2)

    # check which signal env the blow ups are located in and set to 1, and set a flag array
    if len(runaway_rfc) > 0:
        log.debug("check_res_fringes: {} data points exceed threshold".format(len(runaway_rfc)))
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


def make_knots(flux, nknots=20, weights=None):
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
    log.debug("make_knots: creating {} knots on flux array".format(nknots))

    # create an array of indices
    npoints = flux.shape[0]

    # kstep is the number of points / number of knots
    knot_step = npoints // nknots  # + 1

    # define an initial knot index array
    init_knot_idx = np.zeros(nknots + 1)
    for n in range(nknots):
        init_knot_idx[n] = n * knot_step
    init_knot_idx[-1] = npoints - 1

    # get difference between the indices
    knot_split = np.ediff1d(init_knot_idx)

    # the last diff will always be different than the others
    last_split = knot_split[-1]

    # create new knot array with knots at 0 and 1023, then inital and final splits = last_split/2
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


def fit_1d_background_complex(flux, weights, wavenum, order=2, ffreq=None):
    """Fit the background signal using a piecewise spline of n knots. Note that this will also try to identify
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
    if nper_cor >= 5:
        bgindx = make_knots(flux.copy(), int(nknots), weights=weights.copy())
        bgknots = wavenum_scaled[bgindx].astype(float)
        bg_model = SplinesModel(bgknots, order=order)
        fitter = Fitter(wavenum_scaled, bg_model)
    else:
        log.info(" not enough weighted data, no fit performed")
        return flux.copy(), np.zeros(flux.shape[0]), None
    ftr = RobustShell(fitter, domain=10)
    # sometimes the fits will fail for small segments because the model is too complex for the data,
    # in this case return the original array
    # NOTE: since iso-alpha no longer supported this shouldn't be an issue, but will leave here for now
    try:
        ftr.fit(flux.copy(), weights=weights.copy())
    except LinAlgError:
        return flux.copy(), np.zeros(flux.shape[0]), None

    # fit the background
    bg_fit = bg_model.result(wavenum_scaled)
    bg_fit *= np.where(weights.copy() > 1e-07, 1, 1e-08)

    # linearly interpolate over the feature gaps if possible, stops issues later
    try:
        nz, z = interp_helper(weights)
        bg_fit[nz] = np.interp(z(nz), z(~nz), bg_fit[~nz])
        del nz, z
    except ValueError:
        pass

    return bg_fit, bgindx, fitter


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

    fitter = LevenbergMarquardtFitter(wavenum[50:-50], mdl)
    ftr = RobustShell(fitter, domain=10)
    fr_par = ftr.fit(res_fringes[50:-50], weights=weights[50:-50])
    log.debug("fit_quality: best fit pars: {}".format(fr_par))

    if np.abs(fr_par[0]) > np.abs(fr_par[1]):
        contrast = np.abs(round(fr_par[0] * 2, 3))
    else:
        contrast = np.abs(round(fr_par[1] * 2, 3))

    # make data to return for fit quality
    quality = None
    if save_results:
        best_mdl = SineModel(fixed={0: 1 / peak_freq, 1:fr_par[0], 2:fr_par[1]})
        fit = best_mdl.result(wavenum)
        quality = np.array([(10000.0 / wavenum), res_fringes, fit])

    return contrast, quality


def fit_1d_fringes_bayes_evidence(res_fringes, weights, wavenum, ffreq, dffreq, min_nfringes, max_nfringes, pgram_res):

    """Fit the residual fringe signal.

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
    opt_nfringes = None
    peak_freq = None
    freq_min = None
    freq_max = None

    # get the number of weighted pixels
    weighted_pix_num = (weights > 1e-05).sum()
    # set the maximum array size, always 1024
    max_arr_size = int(res_fringes.shape[0])

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
    res_fringe_scan = res_fringes[np.where(weights > 1e-05)]
    wavenum_scan = wavenum[np.where(weights > 1e-05)]

    # use a Lomb-Scargle periodogram to get PSD and identify the strongest frequencies
    log.debug("fit_1d_fringes_bayes: get periodogram")
    pgram = LombScargle(wavenum_scan[::-1], res_fringe_scan[::-1]).power(1 / freq)

    res_fringe_thresh = np.mean(np.abs(pgram))
    peaks = np.argwhere(pgram > res_fringe_thresh)[:, 0]
    if len(peaks) == 0:
        log.debug("fit_1d_fringes_bayes: no significant frequencies found")
        return res_fringes, weighted_pix_num, opt_nfringes, peak_freq, freq_min, freq_max

    else:
        log.debug("fit_1d_fringes_bayes: set significant frequencies")
        # get the peaks
        peaks_freq = freq[peaks]
        # peaks_power = pgram[peaks]

        # limit to max number of fringes to MAX_NFRINGES / (weighted_pix_num / 1024), lower limit of 2
        log.debug("fit_1d_fringes_bayes: max_array_size = {} ".format(max_arr_size))
        log.debug("fit_1d_fringes_bayes: max_nfringes = {} ".format(max_nfringes))
        log.debug("fit_1d_fringes_bayes: fractional weighted pixels: {}/{} ".format(weighted_pix_num, max_arr_size))
        max_fringe_num = int(max_nfringes * (weighted_pix_num / max_arr_size))

        if max_fringe_num < 2:
            max_fringe_num = 2

        log.debug("fit_1d_fringes_bayes: using max fringe number = {} ".format(max_fringe_num))

        if peaks_freq.shape[0] > max_fringe_num:
            log.debug("fit_1d_fringes_bayes: limit scan to {} strongest freqs".format(max_fringe_num))
            peak_ind = pgram.argsort()[-max_fringe_num:][::-1]
            freqs = 1. / freq[peak_ind]
        else:
            log.debug("fit_1d_fringes_bayes: limit scan to {} strongest freqs".format(peaks_freq.shape[0]))
            peak_ind = pgram.argsort()[-peaks_freq.shape[0]:][::-1]
            freqs = 1. / freq[peak_ind]

        # start from 2 frequency, add 1 each loop and check evidence
        log.debug("fit_1d_fringes_bayes: fit {} freqs incrementally, check bayes evidence".format(freqs.shape[0]))

        evidence1 = 1e-5  # arbitrarily small
        opt_nfringes = min_nfringes  # initialise
        if min_nfringes > freqs.shape[0]:
            warning = 'Number of fringes found is less then minimum number of required fringes for column'
            raise Exception(warning)

        for nfringes in np.arange(min_nfringes, freqs.shape[0]):
            # use the significant frequencies setup a multi-sine model of that number of sines
            mdl_fit = multi_sine(nfringes)

            # initialise some variables used in the fitting
            pars = []
            keep_dict = {}

            # fill the variables with parameters to be frozen (freqs) and initialised (amps)
            for n in np.arange(nfringes):
                pars.append(freqs[n])
                pars.append(1.0)
                pars.append(1.0)
                keep_index = n * 3
                keep_dict[keep_index] = freqs[n]

            # fit the multi-sine model and get evidence
            fitter = LevenbergMarquardtFitter(wavenum, mdl_fit, verbose=0, keep=keep_dict)
            ftr = RobustShell(fitter, domain=10)
            try:
                ftr.fit(res_fringes, weights=weights)

                # try get evidence (may fail for large component fits to noisy data, set to very negative value
                try:
                    evidence2 = fitter.getEvidence(limits=[-1, 1], noiseLimits=[0.001, 1])
                except ValueError:
                    evidence2 = -1e9
            except RuntimeError:
                evidence2 = -1e9

            log.debug("fit_1d_fringes_bayes_evidence: nfringe={} ev={} chi={}".format(nfringes, evidence2, fitter.chisq))

            bayes_factor = evidence2 - evidence1
            log.debug(
                "fit_1d_fringes_bayes_evidence: bayes factor={}".format(bayes_factor))
            if bayes_factor > 1:  # strong evidence thresh (log(bayes factor)>1, Kass and Raftery 1995)
                evidence1 = evidence2
                opt_nfringes = nfringes
                log.debug(
                    "fit_1d_fringes_bayes_evidence: strong evidence for nfringes={} ".format(nfringes))
            else:
                log.debug(
                    "fit_1d_fringes_bayes_evidence: no evidence for nfringes={}".format(nfringes))
                break

        mdl_fit = multi_sine(opt_nfringes)
        log.debug("fit_1d_fringes_bayes_evidence: optimal={} fringes".format(opt_nfringes))

        # initialise some variables used in the fitting
        pars = []
        keep_dict = {}

        # fill the variables with parameters to be frozen (freqs) and initialised (amps)
        for n in range(opt_nfringes):
            pars.append(freqs[n])
            pars.append(1.0)
            pars.append(1.0)
            keep_index = n * 3
            keep_dict[keep_index] = freqs[n]

        # fit the optimal multi-sine model
        fitter = LevenbergMarquardtFitter(wavenum, mdl_fit, verbose=0, keep=keep_dict)
        ftr = RobustShell(fitter, domain=10)
        fr_par = ftr.fit(res_fringes, weights=weights)
        best_fringe_model = mdl_fit.copy()
        best_fringe_model.parameters = fr_par
        res_fringe_fit = best_fringe_model(wavenum)

        # create outputs to return
        fitted_frequencies = (1 / freqs[:opt_nfringes + 1]) * factor
        peak_freq = fitted_frequencies[0]
        freq_min = np.amin(fitted_frequencies)
        freq_max = np.amax(fitted_frequencies)

        return res_fringe_fit, weighted_pix_num, opt_nfringes, peak_freq, freq_min, freq_max

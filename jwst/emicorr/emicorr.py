#
#  Module for applying MIRI EMI correction
#

import numpy as np
from scipy import interpolate
import logging
from astropy.stats import sigma_clipped_stats as scs
from stdatamodels.jwst import datamodels
import matplotlib.pyplot as plt
plt.style.use('/Users/tbrandt/tim.mplstyle')

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

subarray_clocks = {

    "SLITLESSPRISM": {
        "rowclocks": 28,
        "frameclocks": 15904},

    "MASKLYOT": {
        "rowclocks": 90,
        "frameclocks": 32400},

    "SUB64": {
        "rowclocks": 28,
        "frameclocks": 8512},

    "SUB128": {
        "rowclocks": 44,
        "frameclocks": 11904},

    "MASK1140": {
        "rowclocks": 82,
        "frameclocks": 23968},

    "MASK1550": {
        "rowclocks": 82,
        "frameclocks": 23968},

    "MASK1065": {
        "rowclocks": 82,
        "frameclocks": 23968},

    "FULL_FAST": {
        "rowclocks": 271,
        "frameclocks": 277504},

    "FULL_SLOW": {
        "rowclocks": 2333,
        "frameclocks": 2388992},

    "BRIGHTSKY": {
        "rowclocks": 162,
        "frameclocks": 86528},

    "SUB256": {
        "rowclocks": 96,
        "frameclocks": 29952},
}


def do_correction(input_model, emicorr_model, save_onthefly_reffile, **pars):
    """
    EMI-correct a JWST data model using an emicorr model

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        Input science data model to be emi-corrected

    emicorr_model : `~jwst.datamodels.EmiModel`
        Data model containing emi correction

    save_onthefly_reffile : str or None
        Full path and root name of on-the-fly reference file

    pars : dict
        Optional user-specified parameters to modify how emicorr
        will operate.  Valid parameters include:
            save_intermediate_results - saves the output into a file and the
                                        reference file (if created on-the-fly)

    Returns
    -------
    output_model : JWST data model
        emi-corrected science data model

    """
    save_intermediate_results = pars['save_intermediate_results']
    nints_to_phase = pars['nints_to_phase']
    nbins = pars['nbins']
    scale_reference = pars['scale_reference']
    onthefly_corr_freq = pars['onthefly_corr_freq']
    use_n_cycles = pars['use_n_cycles']

    output_model = apply_emicorr(input_model, emicorr_model,
                        onthefly_corr_freq, save_onthefly_reffile,
                        save_intermediate_results=save_intermediate_results,
                        nints_to_phase=nints_to_phase,
                        nbins_all=nbins,
                        scale_reference=scale_reference,
                        use_n_cycles=use_n_cycles
                        )

    return output_model


def apply_emicorr(output_model, emicorr_model,
        onthefly_corr_freq, save_onthefly_reffile,
        save_intermediate_results=False, nints_to_phase=None,
        nbins_all=None, scale_reference=True, use_n_cycles=3):
    """
    -> NOTE: This is translated from IDL code fix_miri_emi.pro

    EMI-corrects data by the following procedure:
         [repeat recursively for each discrete EMI frequency desired]
    1) Read image data.
    2) Make very crude slope image and fixed pattern "super"bias for each
       integration, ignoring everything (nonlin, saturation, badpix, etc).
    3) Subtract scaled slope image and bias from each frame of each integration.
    4) Calculate phase of every pixel in the image at the desired emi frequency
       (e.g. 390 Hz) relative to the first pixel in the image.
    5) Make a binned, phased amplitude (pa) wave from the cleaned data (plot
       cleaned vs phase, then bin by phase).
    6)* Measure the phase shift between the binned pa wave and the input
       reference wave
    7)* Use look-up table to get the aligned reference wave value for each pixel
       (the noise array corresponding to the input image).

     6a-7a) Alternately, use the binned pa wave instead of the reference wave to
       "self-correct"

    8) Subtract the noise array from the input image and return the cleaned result.

         [repeat for next frequency using cleaned output as input]
      [removing highest amplitude emi first will probably give best results]

    Parameters
    ----------
    output_model : `~jwst.datamodels.JwstDataModel`
        input science data model to be emi-corrected

    emicorr_model : `~jwst.datamodels.EmiModel`
         ImageModel of emi

    onthefly_corr_freq : float
        Frequency number to do correction with on-the-fly reference file

    save_onthefly_reffile : str or None
        Full path and root name of on-the-fly reference file

    save_intermediate_results : bool
        Saves the output into a file and the reference file (if created on-the-fly)

    nints_to_phase : int
        Number of integrations to phase

    nbins_all : int
        Number of bins in one phased wave (this number will be used for all
        frequencies to be corrected)

    scale_reference : bool
        If True, the reference wavelength will be scaled to the data's phase amplitude

    use_n_cycles : int
        Only use N cycles to calculate the phase to reduce code running time

    Returns
    -------
    output_model : JWST data model
        input science data model which has been emi-corrected
    """
    # get the subarray case and other info
    detector = output_model.meta.instrument.detector
    subarray = output_model.meta.subarray.name
    readpatt = output_model.meta.exposure.readpatt
    xsize = output_model.meta.subarray.xsize   # SUBSIZE1 keyword
    xstart = output_model.meta.subarray.xstart   # SUBSTRT1 keyword
    # get the number of samples, 10us sample times per pixel (1 for fastmode, 9 for slowmode)
    nsamples = output_model.meta.exposure.nsamples

    # get the subarray case from either the ref file or set default values
    freqs_numbers = []
    if emicorr_model is not None:
        log.info('Using reference file to get subarray case.')
        subname, rowclocks, frameclocks, freqs2correct = get_subarcase(emicorr_model, subarray, readpatt, detector)
        reference_wave_list = []

        log.info('With configuration: Subarray={}, Read_pattern={}, Detector={}'.format(subarray, readpatt, detector))
        log.info('Will correct data for the following {} frequencies: '.format(len(freqs2correct)))
        log.info('   {}'.format(freqs2correct))

        if freqs2correct is not None:
            for fnme in freqs2correct:
                freq, ref_wave = get_frequency_info(emicorr_model, fnme)
                freqs_numbers.append(freq)
                reference_wave_list.append(ref_wave)

                if readpatt.upper() == 'FASTR1' or readpatt.upper() == 'SLOWR1':
                    _frameclocks = frameclocks
                else:
                    _frameclocks = 0

                period_in_pixels = (1./freq) / 10.0e-6 #* (1 + 4.3e-6)

                # New routine gets called here.
                d = emicorr_refwave(output_model.data, output_model.pixeldq,
                                    np.array(ref_wave),
                                    nsamples, rowclocks, _frameclocks,
                                    period_in_pixels,
                                    fit_ints_separately=True,
                                    makeplots=True, plotfile='chisq.pdf')
                output_model.data[:] = d

            return output_model
    else:
        # if we got here, the user requested to do correction with on-the-fly reference file
        subname = subarray
        if subname == 'FULL':
            if 'FAST' in readpatt.upper():
                subname += '_FAST'
            elif 'SLOW' in readpatt.upper():
                subname += '_SLOW'
        if subname in subarray_clocks:
            rowclocks = subarray_clocks[subname]['rowclocks']
            frameclocks = subarray_clocks[subname]['frameclocks']
            freqs_numbers = onthefly_corr_freq
            freqs2correct = [repr(freq) for freq in onthefly_corr_freq]
        else:
            subname, rowclocks, frameclocks, freqs2correct = None, None, None, None

    log.info('With configuration: Subarray={}, Read_pattern={}, Detector={}'.format(subarray, readpatt, detector))
    if rowclocks is None or len(freqs_numbers) == 0:
        # no subarray or read pattern match found, print to log and skip correction
        return subname

    # Initialize the output model as a copy of the input
    nints, ngroups, ny, nx = np.shape(output_model.data)

    # create the dictionary to store the frequencies and corresponding phase amplitudes
    if save_intermediate_results and save_onthefly_reffile is not None:
        freq_pa_dict = {'frequencies': {}, 'subarray_cases': {}}

    # Loop over the frequencies to correct
    log.info('Will correct data for the following {} frequencies: '.format(len(freqs2correct)))
    log.info('   {}'.format(freqs2correct))
    for fi, frequency_name in enumerate(freqs2correct):
        frequency = freqs_numbers[fi]
        log.info('Correcting for frequency: {} Hz  ({} out of {})'.format(frequency, fi+1, len(freqs2correct)))

        # Set up some variables

        # Correspondence of array order in IDL
        # sz[0] = 4 in idl
        # sz[1] = nx
        # sz[2] = ny
        # sz[3] = ngroups
        # sz[4] = nints
        nx4 = int(nx/4)

        dd_all = np.zeros((nints, ngroups, ny, nx4))
        log.info('Subtracting self-superbias from each group of each integration and')

        # Calculate times of all pixels in the input integration, then use that to calculate
        # phase of all pixels. Times here is in integer numbers of 10us pixels starting from
        # the first data pixel in the input image. Times can be a very large integer, so use
        # a big datatype. Phaseall (below) is just 0-1.0.

        # A safer option is to calculate times_per_integration and calculate the phase at each
        # int separately. That way times array will have a smaller number of elements at each
        # int, with less risk of datatype overflow. Still, use the largest datatype available
        # for the time_this_int array.

        times_this_int = np.zeros((ngroups, ny, nx4), dtype='ulonglong')
        phaseall = np.zeros((nints, ngroups, ny, nx4))

        # non-roi rowclocks between subarray frames (this will be 0 for fullframe)
        extra_rowclocks = (1024. - ny) * (4 + 3.)
        # e.g. ((1./390.625) / 10e-6) = 256.0 pix and ((1./218.52055) / 10e-6) = 457.62287 pix
        period_in_pixels = (1./frequency) / 10.0e-6

        if nints_to_phase is None and use_n_cycles is None:  # user wats to use all integrations
            # use all integrations
            nints_to_phase = nints
        elif nints_to_phase is None and use_n_cycles is not None: # user wants to use nints_to_phase
            # Calculate how many integrations you need to get that many cycles for
            # a given frequency (rounding up to the nearest Nintegrations)
            nints_to_phase = (use_n_cycles * period_in_pixels) / (frameclocks * ngroups)
            nints_to_phase = int(np.ceil(nints_to_phase))
        elif nints_to_phase is not None and use_n_cycles == 3:
            # user wants to use nints_to_phase because default value for use_n_cycles is 3
            # make sure to never use more integrations than data has
            if nints_to_phase > nints:
                nints_to_phase = nints

        start_time, ref_pix_sample = 0, 3

        # Need colstop for phase calculation in case of last refpixel in a row. Technically,
        # this number comes from the subarray definition (see subarray_cases dict above), but
        # calculate it from the input image header here just in case the subarray definitions
        # are not available to this routine.
        colstop = int(xsize/4 + xstart - 1)
        log.info('doing phase calculation per integration')

        for ninti in range(nints):
            log.debug('  Working on integration: {}'.format(ninti+1))
            # Read in this integration
            data = output_model.data[ninti].copy()

            # Remove source signal and fixed bias from each integration ramp
            # (linear is good enough for phase finding)

            # do linear fit for source + sky
            s0, mm0 = sloper(data[1:ngroups-1, :, :])

            # subtract source+sky from each frame of this ramp
            for ngroupi in range(ngroups):
                data[ngroupi, ...] = output_model.data[ninti, ngroupi, ...] - (s0 * ngroupi)

            # make a self-superbias
            m0 = minmed(data[1:ngroups-1, :, :])

            # subtract self-superbias from each frame of this ramp
            for ngroupi in range(ngroups):
                data[ngroupi, ...] = data[ngroupi, ...] - m0

                # de-interleave each frame into the 4 separate output channels and
                # average (or median) them together for S/N
                d0 = data[ngroupi, :, 0:nx:4]
                d1 = data[ngroupi, :, 1:nx:4]
                d2 = data[ngroupi, :, 2:nx:4]
                d3 = data[ngroupi, :, 3:nx:4]
                dd = (d0 + d1 + d2 + d3)/4.

                # fix a bad ref col
                dd[:, 1] = (dd[:, 0] + dd[:, 3])/2
                dd[:, 2] = (dd[:, 0] + dd[:, 3])/2
                # This is the quad-averaged, cleaned, input image data for the exposure
                dd_all[ninti, ngroupi, ...] = dd - np.median(dd)

            for k in range(ngroups):   # frames
                for j in range(ny):    # rows
                    # nsamples= 1 for fast, 9 for slow (from metadata)
                    times_this_int[k, j, :] = np.arange(nx4, dtype='ulonglong') * nsamples + start_time

                    # If the last pixel in a row is a reference pixel, need to push it out
                    # by ref_pix_sample sample times. The same thing happens for the first
                    # ref pix in each row, but that gets absorbed into the inter-row pad and
                    # can be ignored here. Since none of the current subarrays hit the
                    # right-hand reference pixel, this correction is not in play, but for
                    # fast and slow fullframe (e.g. 10Hz) it should be applied. And even
                    # then, leaving this out adds just a *tiny* phase error on the last ref
                    # pix in a row (only) - it does not affect the phase of the other pixels.

                    if colstop == 258:
                        ulonglong_ref_pix_sample = ref_pix_sample + 2**32
                        times_this_int[k, j, nx4-1] = times_this_int[k, j, nx4-1] + ulonglong_ref_pix_sample

                    # point to the first pixel of the next row (add "end-of-row" pad)
                    start_time += rowclocks

                # point to the first pixel of the next frame (add "end-of-frame" pad)
                start_time += extra_rowclocks

            # Convert "times" to phase each integration. Note that times has units of
            # number of 10us from the first data pixel in this integration, so to
            # convert to phase, divide by the waveform *period* in float pixels
            phase_this_int = times_this_int / period_in_pixels
            phaseall[ninti, ...] = phase_this_int - phase_this_int.astype('ulonglong')

            # add a frame time to account for the extra frame reset between MIRI integrations
            if readpatt.upper() == 'FASTR1' or readpatt.upper() == 'SLOWR1':
                start_time += frameclocks

        # use phaseall vs dd_all

        # Define the sizew of 1 wave of the phased waveform vector, then bin the whole
        # dataset at this interval. This is essentially the desired number of bins along
        # the waveform. This can be any number that is at least 1 less than the period in
        # pixels. If larger, some bins could end up sparsely sampled. Fewer bins results
        # in a smoother waveform but lower resolution. Can always be smoothed and/or
        # stretched to a different dimension later. By default, use period_in_pixels/2.0

        # number of bins in one phased wave (e.g. 255, 218, etc)
        if nbins_all is None:
            # the IDL code sets nbins as ulong type (ulonglong in python)
            nbins = int(period_in_pixels/2.0)
        else:
            nbins = nbins_all
        if nbins > 501:
            nbins = 500

        # bin the whole set
        log.info('Calculating the phase amplitude for {} bins'.format(nbins))
        # Define the binned waveform amplitude (pa = phase amplitude)
        pa = np.arange(nbins, dtype=float)
        # keep track of n per bin to check for low n
        nb_over_nbins = [nb/nbins for nb in range(nbins)]
        nbp1_over_nbins = [(nb + 1)/nbins for nb in range(nbins)]
        # Construct a phase map and dd map for only the nints_to_phase
        phase_temp = phaseall[0: nints_to_phase, :, :, :]
        dd_temp = dd_all[0: nints_to_phase, :, :, :]
        for nb in range(nbins):
            u = (phase_temp > nb_over_nbins[nb]) & (phase_temp <= nbp1_over_nbins[nb])
            # calculate the sigma-clipped mean
            dmean, _, _ = scs(dd_temp[u])
            pa[nb] = dmean   # amplitude in this bin

        pa -= np.median(pa)

        # pa_phase is the phase corresponding to each bin in pa. The +0.5 here is to center in
        # the bin. Together, (pa, pa_phase) define the waveform to phase against the
        # reference_wave, or to use directly for a self correction (likely noiser than the
        # reference_wave, but perfectly matched in phase by definition so no need for a
        # crosscor or other phase-finder against a reference_wave). This may not be
        # necessary. Can just assume pa and the reference_wave are exactly 1.0 waves long.
        #pa_phase = (np.arange(nbins, dtype=np.float32) + 0.5) / nbins

        # Correction using direct look-up from either pa (self-correction) or
        # the reference wave (if provided)
        # Make a vector to hold a version of pa that is close to an integer number of
        # pixels long. This will be the noise amplitude look-up table.
        lut = rebin(pa, [period_in_pixels])   # shift and resample reference_wave at pa's phase

        # If a reference_wave_filename was provided, use it to measure phase shift between
        # pa (binned phased wave) from the input image vs. the reference_wave using xcorr.
        # These two methods give the same integer-reference_wave-element resolution results.
        if emicorr_model is not None:
            log.info('Using reference file to measure phase shift')
            reference_wave = np.array(reference_wave_list[fi])
            reference_wave_size = np.size(reference_wave)
            rebinned_pa = rebin(pa, [reference_wave_size])
            cc = np.zeros(reference_wave_size)
            for i in range(reference_wave_size):
                shifted_ref_wave = np.roll(reference_wave, i)
                pears_coeff = np.corrcoef(shifted_ref_wave, rebinned_pa)
                cc[i] = pears_coeff[0, 1]

            # Use the reference file (shifted to match the phase of pa,
            # and optionally amplitude scaled)
            # shift and resample reference_wave at pa's phase
            # u[0] is the phase shift of reference_wave *to* pa
            u = np.argmax(cc)
            lut_reference = rebin(np.roll(reference_wave, u), [period_in_pixels])

            # Scale reference wave amplitude to match the pa amplitude from this dataset by
            # fitting a line rather than taking a mean ratio so that any DC offset drops out
            # in the intercept (and to avoid potential divide-by-zeros)
            m, b = np.polyfit(lut_reference, lut, 1)   # the slope is the amplitude ratio

            # check pa-to-reference_wave phase match and optional scale
            if scale_reference:
                lut_reference = lut_reference * m

            lut = lut_reference

        if save_intermediate_results and save_onthefly_reffile is not None:
            freq_pa_dict['frequencies'][frequency_name] = {'frequency': frequency,
                                                           'phase_amplitudes': pa}

        log.info('Creating phased-matched noise model to subtract from data')
        # This is the phase matched noise model to subtract from each pixel of the input image
        dd_noise = lut[(phaseall * period_in_pixels).astype(int)]

        # Safety catch; anywhere the noise value is not finite, set it to zero
        dd_noise[~np.isfinite(dd_noise)] = 0.0

        # Subtract EMI noise from the input data
        log.info('Subtracting EMI noise from data')

        # Interleave (straight copy) into 4 amps
        for k in range(4):
            output_model.data[..., k::4] -= dd_noise

        # clean up
        del dd_all
        del times_this_int
        del phaseall
        del dd_noise

    if save_intermediate_results and save_onthefly_reffile is not None:
        if 'FAST' in readpatt:
            freqs_dict = {readpatt: freqs2correct}
        else:
            freqs_dict = {readpatt: {detector: freqs2correct} }
        on_the_fly_subarr_case = {}
        on_the_fly_subarr_case[subarray] = {
            'rowclocks': rowclocks,
            'frameclocks': frameclocks,
            'freqs': freqs_dict
        }
        freq_pa_dict['subarray_cases'] = on_the_fly_subarr_case
        mk_reffile(freq_pa_dict, save_onthefly_reffile)

    return output_model


def sloper(data):
    """ Fit slopes to all pix of a ramp, using numerical recipes plane-adding
     returning intercept image.

    Parameters
    ----------
    data : numpy array
        3-D integration data array

    Returns
    -------
    outarray : numpy array
        3-D integration data array

    intercept: numpy array
        slope intercept values
    """
    ngroups, ny, nx = np.shape(data)
    frametime = 1.0   # use 1.0 for per-frame slopes, otherwise put a real time here
    grouptimes = np.arange(ngroups) * frametime

    sxy = np.zeros((ny, nx))
    sx = np.zeros((ny, nx))
    sy = np.zeros((ny, nx))
    sxx = np.zeros((ny, nx))

    for groupcount in range(ngroups):
        # do the incremental slope calculation operations
        sxy = sxy + grouptimes[groupcount] * data[groupcount, ...]
        sx = sx + grouptimes[groupcount]
        sy = sy + data[groupcount, ...]
        sxx = sxx + grouptimes[groupcount]**2

    # calculate the final per-pixel slope values
    outarray = (ngroups * sxy - sx * sy)/(ngroups * sxx - sx * sx)

    # calculate the final per-pixel intercept values
    intercept = (sxx * sy - sx * sxy)/(ngroups * sxx - sx * sx)

    return outarray, intercept


def minmed(data, minval=False, avgval=False, maxval=False):
    """ Returns the median image of a stack of images, or if there are 2 or less
    non-zero frames, returns the minimum (or this minimum can be forced).

    Parameters
    ----------
    data : numpy array
        3-D integration data array

    minval : bool
        Returned image will be the minimum

    avgval : bool
        Returned image will be the mean

    maxval : bool
        Returned image will be the maximum

    Returns
    -------
    medimg : numpy array
        Median image of a stack of images
    """

    n = data.size
    if n > 0:
        if n <= 2 or minval:
            medimg = np.nanmin(data, axis=0)
        if maxval:
            medimg = np.nanmax(data, axis=0)
        if not minval and not maxval and not avgval:
            medimg = np.nanmedian(data, axis=0)
        if avgval:
            medimg = np.nanmean(data, axis=0)
    return medimg


def get_subarcase(subarray_cases, subarray, readpatt, detector):
    """ Get the rowclocks and frameclocks values for the given configuration.

    Parameters
    ----------
    subarray_cases : dict or model object
        Either default corrections dictionary or datamodel

    subarray : str
        Keyword value

    readpatt : str
        Keyword value

    detector : str
        Keyword value

    Returns
    -------
    subname : str
        Modified subarray name

    rowclocks : int

    frameclocks : int

    frequencies : list
        List of frequencies to correct according to subarray name
    """
    subname, rowclocks, frameclocks, frequencies = None, None, None, None

    # make sure the readpattern is defined as expected to read data from reference file
    readpatt, no_slow_match_readpatt, no_fast_match_readpatt = readpatt.upper(), False, False
    if "SLOW" not in readpatt:
        no_slow_match_readpatt = True
    if "FAST" not in readpatt:
        no_fast_match_readpatt = True
    if no_slow_match_readpatt and no_fast_match_readpatt:
        log.info('Read pattern {} does not include expected string FAST or SLOW'.format(readpatt))
        return subname, rowclocks, frameclocks, frequencies

    # search and return the specific values for the configuration
    frequencies = []
    mdl_dict = subarray_cases.to_flat_dict()
    for subname in subarray_cases.subarray_cases:
        subname = subname.split(sep='.')[0]
        dataconfig = subarray
        if 'FULL' in dataconfig:
            dataconfig = subarray + '_' + readpatt.replace('R1', '')
        if dataconfig != subname:
            continue
        log.debug('Found subarray case {}!'.format(subname))
        for item, val in mdl_dict.items():
            if subname in item:
                if "rowclocks" in item:
                    rowclocks = val
                elif "frameclocks" in item:
                    frameclocks = val
                else:
                    if "SLOW" in readpatt and "SLOW" in item and detector in item:
                        frequencies.append(val)
                    elif "FAST" in readpatt and "FAST" in item:
                        frequencies.append(val)
            if subname is not None and rowclocks is not None and frameclocks is not None and frequencies is not None:
                break
    return subname, rowclocks, frameclocks, frequencies


def get_frequency_info(freqs_names_vals, frequency_name):
    """Get the frequency number from the given dictionary
    Parameters
    ----------
    freqs_names_vals : dict or model object
        Either default corrections dictionary or datamodel

    frequency_name : str
        Frequency of interest

    Returns
    -------
    frequency_number : float
        Frequency

    phase_amplitudes : array
        1-D array of the corresponding phase amplidues for this frequency
    """
    if isinstance(freqs_names_vals, dict):
        for freq_nme, val in freqs_names_vals.items():
            if freq_nme == frequency_name:
                return val
    else:
        freq_number, phase_amplitudes = None, None
        mdl_dict = freqs_names_vals.to_flat_dict()
        for item, val in mdl_dict.items():
            if frequency_name in item:
                if 'frequency' in item:
                    freq_number = val
                if 'phase_amplitudes' in item:
                    phase_amplitudes = val
            if freq_number is not None and phase_amplitudes is not None:
                break
        return freq_number, phase_amplitudes


def rebin(arr, newshape):
    """Rebin an array to a new shape.

    Parameters
    ----------
    arr : numpy array
        array to rebin

    newshape : list
        New shape for the input array

    Returns
    -------
    arr : numpy array
        rebinned array
    """
    assert len(arr.shape) == len(newshape)

    slices = [ slice(0,old, float(old)/new) for old,new in zip(arr.shape,newshape) ]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')   # choose the biggest smaller integer index

    return arr[tuple(indices)]


def mk_reffile(freq_pa_dict, emicorr_ref_filename):
    """ Create the reference file hdulist object.

    Parameters
    ----------
    freq_pa_dict : dictionary
        Dictionary containing the phase amplitude (pa) array
        corresponding to the appropriate frequency

    emicorr_ref_filename : str
        Full path and root name of on-the-fly reference file

    Returns
    -------
    Nothing
    """
    # save the reference file if save_intermediate_results=True and no CRDS file exists
    emicorr_model = datamodels.EmiModel(freq_pa_dict)
    if 'fits' in emicorr_ref_filename:
        emicorr_ref_filename = emicorr_ref_filename.replace('fits', 'asdf')
    emicorr_model.save(emicorr_ref_filename)
    emicorr_model.close()
    log.info('On-the-fly reference file written as: %s', emicorr_ref_filename)





def emicorr_refwave(data, pdq, model, nsamples, rowclocks, frameclocks,
                    period_in_pixels, fit_ints_separately=False,
                    nphases_opt=500, makeplots=False, plotfile='chisq.pdf'):

    """
    Derive the best amplitude and phase for the EMI waveform, subtract it off.

    """

    nints, ngroups, ny, nx = data.shape
    nx4 = nx//4
    extra_rowclocks = (1024. - ny) * (4 + 3.)

    frametime = ny*rowclocks + extra_rowclocks

    # Times of the individual reads in pixel clock times
    t0_arr = np.zeros((ny, nx4))
    for i in range(ny):
        t0_arr[i] = i*rowclocks + np.arange(nx4)*nsamples

    phase = (t0_arr/period_in_pixels)%1

    # Phase gap between groups
    dphase = (frametime/period_in_pixels)%1

    # Phase gap between integrations
    dphase_frame = dphase*ngroups + frameclocks/period_in_pixels

    # Time-like index for group number
    grouptimes = np.arange(ngroups)

    # Phase(model) for interpolation.  Wrap the phases to avoid
    # extrapolation, and use a cubic spline.

    model_extended = np.array([model[-1]] + list(model) + [model[0]])
    phase_extended = (np.arange(len(model) + 2) - 0.5)/len(model)
    phasefunc = interpolate.interp1d(phase_extended, model_extended, kind='cubic')

    phases_template = (phase_extended[1:-1, np.newaxis] + grouptimes*dphase)%1
    nphases = phases_template.shape[0]

    # These arrays hold the sum of the values of good pixels at a
    # given phase and the total number of good pixels at a given
    # phase, respectively.  The phase refers to the first group; other
    # groups will have the appropriate phase delay added.

    all_y = np.zeros((nints, nphases, ngroups))
    all_N = np.zeros((nints, nphases))

    # "Good" pixel here has no more than twice the median standard
    # deviation among group differences and is not flagged in the pdq
    # array.  This should discard most bad and high-flux pixels.

    pixel_std = np.std(data, axis=1)
    pixel_ok = (pixel_std < 2*np.median(pixel_std))&(pdq == 0)

    for j in range(ny):
        for k in range(nx4):
            # Bad reference column?  I am not sure whether this is needed.
            # I tried to carry it from the current routine.
            if k == 1 or k == 2:
                continue
            # Choose the index corresponding to our phase.
            indx = int(phase[j, k]*nphases)
            # All four output channels.
            for l in range(4):
                pixok = pixel_ok[:, j, k*4 + l]
                all_y[:, indx] += data[:, :, j, k*4 + l]*pixok[:, np.newaxis]
                all_N[:, indx] += pixok

    # We'll compute chi2 at nphases_opt evenly spaced phases.
    phaselist = np.arange(nphases_opt)*1./nphases_opt

    # Class that computes and holds all of the intermediate
    # information needed for the fits.

    emifitter = EMIfitter(all_y, all_N, nints, ngroups, phasefunc,
                          phases_template, phaselist, dphase_frame)

    # diagnostic
    if makeplots:
        fig, (ax) = plt.subplots(1, 1)

    # The case where each integration gets its own phase and amplitude
    if fit_ints_separately:

        amplitudes_to_correct = []
        phases_to_correct = []
        for i in range(nints):

            # Compute chi squared at all phase offsets, then find the
            # phase offset that gives the best answer, then get the
            # corresponding amplitude.

            chisq, _ = calc_chisq_amplitudes(emifitter, ints=[i])
            best_phase = get_best_phase(phaselist, chisq)
            phases_to_correct += [best_phase + i * dphase_frame]

            _, c = calc_chisq_amplitudes(emifitter, phases=[best_phase], ints=[i])
            amplitudes_to_correct += c

            if makeplots:
                cval = plt.cm.viridis_r([i/(nints - 1.)])
                ax.plot(phaselist, chisq, color=cval)

    # One phase and amplitude fit to all integrations
    else:

        chisq, _ = calc_chisq_amplitudes(emifitter)
        best_phase = get_best_phase(phaselist, chisq)
        _, c = calc_chisq_amplitudes(emifitter, phases=[best_phase])

        # Phases appropriate for each integration
        phases_to_correct = (best_phase + np.arange(nints)*dphase_frame)%1

        # Output from calc_chisq_amplitudes here is a one-element list.
        # Repeat it for each integration.

        amplitudes_to_correct = c*nints

        if makeplots:
            ax.plot(phaselist, chisq)

    # Diagnostic
    if makeplots:
        ax.set_xlabel("Phase")
        ax.set_ylabel(r"$\Delta \chi^2$")

        if fit_ints_separately:
            sm = plt.cm.ScalarMappable(cmap='viridis_r')
            cb = plt.colorbar(sm, ax=ax, label="Integration Number", cmap='viridis_r')
            locs = np.arange(nints) + 1
            cb.set_ticks([], minor=True)
            cb.set_ticks((locs - 1)/(locs[-1] - 1), labels=locs)
        plt.tight_layout()
        plt.savefig(plotfile)

    for i in range(nints):
        for j in range(ngroups):

            # Place the reference waveform at the appropriate phase,
            # scale, and subtract from each output channel.

            phased_emi = phasefunc((phase + dphase*j + phases_to_correct[i])%1)

            for k in range(4):
                data[i, j, :, k::4] -= amplitudes_to_correct[i]*phased_emi

    return data



def get_best_phase(phases, chisq):
    """
    Fit a parabola to get the phase corresponding to the best chi2.

    Phases input should be periodic: phases[0] comes next after phases[-1]
    """

    ibest = np.argmin(chisq)

    # If ibest is zero, ibest_m1 will give the final element.
    ibest_m1 = ibest - 1
    ibest_p1 = ibest + 1

    # If ibest is the last element, we have to explicitly wrap.
    if ibest_p1 > len(chisq) - 1:
        ibest_p1 = 0

    # Calculate the vertex of the parabola.
    chisq_opt = [chisq[ibest_m1], chisq[ibest], chisq[ibest_p1]]
    x = [phases[ibest_m1], phases[ibest], phases[ibest_p1]]

    # Just in case all of these chi squared values are equal...
    if np.std(chisq_opt) == 0:
        log.info('Chi squared values as a function of phase offset are '
                 'unexpectedly equal in EMIcorr')
        return phases[ibest]

    a, b, _ = np.polyfit(x, chisq_opt, 2)
    bestphase = -b/(2*a)

    return bestphase


def calc_chisq_amplitudes(emifitter, ints=None, phases=None):
    """
    Compute chi2 and amplitude of EMI waveform at phase(s)

    This has the math from the writeup

    """

    E = emifitter

    # By default, compute a single phase and amplitudes over all
    # integrations.

    if ints is None:
        ints = np.arange(E.nints)

    # By default, calculate chi squared and the best-fit amplitude
    # for every phase in the input EMIfitter's phaselist.

    if phases is None:
        phases = E.phaselist

    chisq = []
    amplitudes = []

    # Compute the best chi squared and the best amplitude at every
    # requested phase using the math in the writeup.

    for phase in phases:
        A = 0
        B = 0

        for i in ints:

            y = E.all_y[i]
            N = E.all_N[i]
            Sy = E.all_Sy[i]
            Sty = E.all_Sty[i]

            phase_diff = E.dphase_frame*i
            k = np.argmin(np.abs((E.phaselist - phase - phase_diff)%1))

            z = E.zlist[k]
            Stz = E.Stzlist[k]
            Sz = E.Szlist[k]
            Szz = E.Szzlist[k]
            Syz = np.sum(z * y, axis=1)

            A += np.sum(N / E.Delta*(-E.Stt * Sz**2 + 2 * E.St * Sz * Stz
                                     - E.ngroups * Stz**2 + Szz * E.Delta))
            B += (2 / E.Delta)*np.sum(E.Stt * Sz * Sy - E.St * Stz * Sy
                                      - E.St * Sz * Sty + E.ngroups * Stz * Sty)
            B -= 2 * np.sum(Syz)

        c = -B / (2 * A)
        chisq += [A * c**2 + B * c]
        amplitudes += [c]

    return chisq, amplitudes



class EMIfitter:
    """
    Compute and store quantities needed for chi2, amplitude calculation.

    Some of the math from the writeup goes in here.
    """
    def __init__(self, all_y, all_N, nints, ngroups, phasefunc,
                 phases_template, phaselist, dphase_frame):

        self.all_y = all_y
        self.all_N = all_N
        self.nints = nints
        self.ngroups = ngroups

        self.grouptimes = np.arange(self.ngroups)
        self.Stt = np.sum(self.grouptimes**2)
        self.St = np.sum(self.grouptimes)
        self.Delta = self.ngroups*self.Stt - self.St**2

        self.all_Sy = [np.sum(self.all_y[i], axis=1)
                       for i in range(self.nints)]
        self.all_Sty = [np.sum(self.all_y[i]*self.grouptimes, axis=1)
                        for i in range(self.nints)]

        self.phaselist = phaselist
        self.phases_template = phases_template

        self.phasefunc = phasefunc
        self.dphase_frame = dphase_frame
        self.nphases = len(phaselist)

        self.zlist = []
        self.Szlist = []
        self.Stzlist = []
        self.Szzlist = []

        # Waveform for each pixel when the first one is at dphaseval.
        for dphaseval in phaselist:
            z = self.phasefunc((self.phases_template + dphaseval)%1)
            self.zlist += [z]
            self.Stzlist += [np.sum(self.grouptimes*z, axis=1)]
            self.Szlist += [np.sum(z, axis=1)]
            self.Szzlist += [np.sum(z**2, axis=1)]

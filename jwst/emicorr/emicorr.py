#
#  Module for applying MIRI EMI correction
#

import numpy as np
import logging
from astropy.stats import sigma_clipped_stats as scs
from stdatamodels.jwst import datamodels

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
            user_supplied_reffile - reference file supplied by the user

    Returns
    -------
    output_model : JWST data model
        emi-corrected science data model

    """
    save_intermediate_results = pars['save_intermediate_results']
    user_supplied_reffile = pars['user_supplied_reffile']
    nints_to_phase = pars['nints_to_phase']
    nbins = pars['nbins']
    scale_reference = pars['scale_reference']
    onthefly_corr_freq = pars['onthefly_corr_freq']
    use_n_cycles = pars['use_n_cycles']

    output_model = apply_emicorr(input_model, emicorr_model,
                        onthefly_corr_freq, save_onthefly_reffile,
                        save_intermediate_results=save_intermediate_results,
                        user_supplied_reffile=user_supplied_reffile,
                        nints_to_phase=nints_to_phase,
                        nbins_all=nbins,
                        scale_reference=scale_reference,
                        use_n_cycles=use_n_cycles
                        )

    return output_model


def apply_emicorr(input_model, emicorr_model,
        onthefly_corr_freq, save_onthefly_reffile,
        save_intermediate_results=False, user_supplied_reffile=None,
        nints_to_phase=None, nbins_all=None, scale_reference=True,
        use_n_cycles=3):
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
    input_model : `~jwst.datamodels.JwstDataModel`
        input science data model to be emi-corrected

    emicorr_model : `~jwst.datamodels.EmiModel`
         ImageModel of emi

    onthefly_corr_freq : float
        Frequency number to do correction with on-the-fly reference file

    save_onthefly_reffile : str or None
        Full path and root name of on-the-fly reference file

    save_intermediate_results : bool
        Saves the output into a file and the reference file (if created on-the-fly)

    user_supplied_reffile : str
        Reference file supplied by the user

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
    detector = input_model.meta.instrument.detector
    subarray = input_model.meta.subarray.name
    readpatt = input_model.meta.exposure.readpatt
    xsize = input_model.meta.subarray.xsize   # SUBSIZE1 keyword
    xstart = input_model.meta.subarray.xstart   # SUBSTRT1 keyword

    # get the subarray case from either the ref file or set default values
    freqs_numbers = []
    if emicorr_model is not None:
        log.info('Using reference file to get subarray case.')
        subname, rowclocks, frameclocks, freqs2correct = get_subarcase(emicorr_model, subarray, readpatt, detector)
        reference_wave_list = []
        if freqs2correct is not None:
            for fnme in freqs2correct:
                freq, ref_wave = get_frequency_info(emicorr_model, fnme)
                freqs_numbers.append(freq)
                reference_wave_list.append(ref_wave)
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

    # get the number of samples, 10us sample times per pixel (1 for fastmode, 9 for slowmode)
    nsamples = input_model.meta.exposure.nsamples

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()
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

        # Read image data and set up some variables
        orig_data = output_model.data
        data = orig_data.copy()

        # Correspondance of array order in IDL
        # sz[0] = 4 in idl
        # sz[1] = nx
        # sz[2] = ny
        # sz[3] = ngroups
        # sz[4] = nints
        nx4 = int(nx/4)

        dd_all = np.ones((nints, ngroups, ny, nx4))
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
        colstop = int( xsize/4 + xstart - 1 )
        log.info('doing phase calculation per integration')

        for ninti in range(nints):
            log.debug('  Working on integration: {}'.format(ninti+1))

            # Remove source signal and fixed bias from each integration ramp
            # (linear is good enough for phase finding)

            # do linear fit for source + sky
            s0, mm0 = sloper(data[ninti, 1:ngroups-1, :, :])

            # subtract source+sky from each frame of this ramp
            for ngroupi in range(ngroups):
                data[ninti, ngroupi, ...] = orig_data[ninti, ngroupi, ...] - (s0 * ngroupi)

            # make a self-superbias
            m0 = minmed(data[ninti, 1:ngroups-1, :, :])

            # subtract self-superbias from each frame of this ramp
            for ngroupi in range(ngroups):
                data[ninti, ngroupi, ...] = data[ninti, ngroupi, ...] - m0

                # de-interleave each frame into the 4 separate output channels and
                # average (or median) them together for S/N
                d0 = data[ninti, ngroupi, :, 0:nx:4]
                d1 = data[ninti, ngroupi, :, 1:nx:4]
                d2 = data[ninti, ngroupi, :, 2:nx:4]
                d3 = data[ninti, ngroupi, :, 3:nx:4]
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
        # streched to a different dimension later. By default, use period_in_pixels/2.0

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
        phase_temp = phaseall[0:nints_to_phase,:,:,:]
        dd_temp = dd_all[0:nints_to_phase,:,:,:]
        for nb in range(nbins):
            u = np.where((phase_temp > nb_over_nbins[nb]) & (phase_temp <= nbp1_over_nbins[nb]))
            # calculate the sigma-clipped mean
            dmean,_,_ = scs(dd_temp[u])
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
            freq_pa_dict['frequencies'][frequency_name] = {'frequency' : frequency,
                                                        'phase_amplitudes' : pa}

        log.info('Creating phased-matched noise model to subtract from data')
        # This is the phase matched noise model to subtract from each pixel of the input image
        dd_noise = lut[(phaseall * period_in_pixels).astype(int)]

        # Interleave (straight copy) into 4 amps
        noise = np.ones((nints, ngroups, ny, nx))   # same size as input data
        noise_x = np.arange(nx4) * 4
        for k in range(4):
            noise[:, :, :, noise_x + k] = dd_noise

        # Subtract EMI noise from the input data
        log.info('Subtracting EMI noise from data')
        corr_data = orig_data - noise
        output_model.data = corr_data

        # clean up
        del data
        del dd_all
        del times_this_int
        del phaseall

    if save_intermediate_results and save_onthefly_reffile is not None:
        if 'FAST' in readpatt:
            freqs_dict = {readpatt: freqs2correct}
        else:
            freqs_dict = {readpatt: {detector: freqs2correct} }
        on_the_fly_subarr_case = {}
        on_the_fly_subarr_case[subarray] = {
            'rowclocks': rowclocks,
            'frameclocks': frameclocks,
            'freqs': freqs_dict }
        freq_pa_dict['subarray_cases'] = on_the_fly_subarr_case
        mk_reffile(freq_pa_dict, save_onthefly_reffile)

    return output_model


def sloper(data):
    """ Fit slopes to all pix of a ramp, using numerical recipies plane-adding
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

    ngroups, ny, nx = np.shape(data)
    medimg = np.zeros((ny, nx))
    # use a mask to ignore nans for calculations
    vec = np.ma.array(data, mask=np.isnan(data))
    n = vec.size
    if n > 0:
        if n <= 2 or minval:
            medimg = np.ma.min(vec, axis=0)
        if maxval:
            medimg = np.ma.max(vec, axis=0)
        if not minval and not maxval and not avgval:
            medimg = np.ma.median(vec, axis=0)
        if avgval:
            medimg = np.ma.mean(vec, axis=0)
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
                    elif "FAST" in readpatt  and "FAST" in item:
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



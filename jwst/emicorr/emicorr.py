#
#  Module for applying MIRI EMI correction
#

import numpy as np
import logging

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


default_emi_freqs = {
    "Hz390": 390.625,
    "Hz390_sub128": 390.625,
    "Hz10": 10.039216,
    "Hz10_slow_MIRIMAGE": 10.039216,
    "Hz10_slow_MIRIFULONG": 10.039216,
    "Hz10_slow_MIRIFUSHORT": 10.039216
}

default_subarray_cases = {

    # 390Hz out-of-phase - these all need 390hz correction

    "SLITLESSPRISM": {
        "rowclocks": 28,
        "frameclocks": 15904,
        "freqs":  {"FAST": ["Hz390", "Hz10"],
                   "SLOW":  {"MIRIMAGE" : ["Hz390", "Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz390", "Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz390", "Hz10_slow_MIRIFUSHORT"]}}},

    "MASKLYOT": {
        "rowclocks": 90,
        "frameclocks": 32400,
        "freqs":  {"FAST": ["Hz390", "Hz10"],
                   "SLOW":  {"MIRIMAGE" : ["Hz390", "Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz390", "Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz390", "Hz10_slow_MIRIFUSHORT"]}}},

    "SUB64": {
        "rowclocks": 28,
        "frameclocks": 8512,
        "freqs":  {"FAST": ["Hz390", "Hz10"],
                   "SLOW":  {"MIRIMAGE" : ["Hz390", "Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz390", "Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz390", "Hz10_slow_MIRIFUSHORT"]}}},

    "SUB128": {
        "rowclocks": 44,
        "frameclocks": 11904,
        "freqs":  {"FAST": ["Hz390_sub128", "Hz10"],
                   "SLOW":  {"MIRIMAGE" : ["Hz390", "Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz390", "Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz390", "Hz10_slow_MIRIFUSHORT"]}}},

    "MASK1140": {
        "rowclocks": 82,
        "frameclocks": 23968,
        "freqs":  {"FAST": ["Hz390", "Hz10"],
                   "SLOW":  {"MIRIMAGE" : ["Hz390", "Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz390", "Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz390", "Hz10_slow_MIRIFUSHORT"]}}},

    "MASK1550": {
        "rowclocks": 82,
        "frameclocks": 23968,
        "freqs":  {"FAST": ["Hz390", "Hz10"],
                   "SLOW":  {"MIRIMAGE" : ["Hz390", "Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz390", "Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz390", "Hz10_slow_MIRIFUSHORT"]}}},

    # 390Hz already in-phase for these, but may need corr for other
    # frequencies (e.g. 10Hz heater noise)

    "FULL_FAST": {
        "rowclocks": 271,
        "frameclocks": 277504,
        "freqs": {"FAST" : ["Hz10"]}},

    "FULL_SLOW": {
        "rowclocks": 2333,
        "frameclocks": 2388992,
        "freqs": {"SLOW":  {"MIRIMAGE" : ["Hz10_slow_MIRIMAGE"],
                              "MIRIFULONG" : ["Hz10_slow_MIRIFULONG"],
                              "MIRIFUSHORT" : ["Hz10_slow_MIRIFUSHORT"]}}},

    "BRIGHTSKY": {
        "rowclocks": 162,
        "frameclocks": 86528,
        "freqs": {"FAST" : ["Hz10"],
                  "SLOW": {"MIRIMAGE" : ["Hz10_slow_MIRIMAGE"],
                            "MIRIFULONG" : ["Hz10_slow_MIRIFULONG"],
                            "MIRIFUSHORT" : ["Hz10_slow_MIRIFUSHORT"]}}},

    "SUB256": {
        "rowclocks": 96,
        "frameclocks": 29952,
        "freqs": {"FAST" : ["Hz10"],
                  "SLOW": {"MIRIMAGE" : ["Hz10_slow_MIRIMAGE"],
                            "MIRIFULONG" : ["Hz10_slow_MIRIFULONG"],
                            "MIRIFUSHORT" : ["Hz10_slow_MIRIFUSHORT"]}}},
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

    output_model = apply_emicorr(input_model, emicorr_model, save_onthefly_reffile,
                        save_intermediate_results=save_intermediate_results,
                        user_supplied_reffile=user_supplied_reffile,
                        nints_to_phase=nints_to_phase,
                        nbins_all=nbins,
                        scale_reference=scale_reference
                        )

    return output_model


def apply_emicorr(input_model, emicorr_model, save_onthefly_reffile,
        save_intermediate_results=False, user_supplied_reffile=None,
        nints_to_phase=None, nbins_all=None, scale_reference=True):
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
        for fnme in freqs2correct:
            freq, ref_wave = get_frequency_info(emicorr_model, fnme)
            freqs_numbers.append(freq)
            reference_wave_list.append(ref_wave)
    else:
        log.info('Using default subarray case corrections.')
        subname, rowclocks, frameclocks, freqs2correct = get_subarcase(default_subarray_cases, subarray, readpatt, detector)
        for fnme in freqs2correct:
            freq = get_frequency_info(default_emi_freqs, fnme)
            freqs_numbers.append(freq)

    log.info('With configuration: Subarray={}, Read_pattern={}, Detector={}'.format(subarray, readpatt, detector))
    if rowclocks is None or len(freqs_numbers) == 0:
        # no subarray or read pattern match found, print to log and skip correction
        return subname

    # get the number of samples, 10us sample times per pixel (1 for fastmode, 9 for slowmode)
    nsamples = input_model.meta.exposure.nsamples

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

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
        nints, ngroups, ny, nx = np.shape(data)
        if nints_to_phase is None:
            nints_to_phase = nints

        # Correspondance of array order in IDL
        # sz[0] = 4 in idl
        # sz[1] = nx
        # sz[2] = ny
        # sz[3] = ngroups
        # sz[4] = nints
        nx4 = int(nx/4)

        dd_all = np.ones((nints, ngroups, ny, nx4))
        log.info('Subtracting self-superbias from each group of each integration')
        for ninti in range(nints_to_phase):
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

        start_time, ref_pix_sample = 0, 3

        # Need colstop for phase calculation in case of last refpixel in a row. Technically,
        # this number comes from the subarray definition (see subarray_cases dict above), but
        # calculate it from the input image header here just in case the subarray definitions
        # are not available to this routine.
        colstop = int( xsize/4 + xstart - 1 )
        log.info('Phase calculation per integration')
        for l in range(nints_to_phase):
            log.debug('  Working on integration: {}'.format(l+1))
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

            # Convert "times" to phase each integration. Nothe times has units of
            # number of 10us from the first data pixel in this integration, so to
            # convert to phase, divide by the waveform *period* in float pixels
            phase_this_int = times_this_int / period_in_pixels
            phaseall[l, ...] = phase_this_int - phase_this_int.astype('ulonglong')

            # add a frame time to account for the extra frame reset between MIRI integrations
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
        nu = np.arange(nbins)
        for nb in range(nbins):
            u = np.where((phaseall > nb/nbins) & (phaseall <= (nb + 1)/nbins) & (np.isfinite(dd_all)))
            nu[nb] = phaseall[u].size
            # calculate the sigma-clipped mean
            dmean, _, _, _ = iter_stat_sig_clip(dd_all[u])
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
            reference_wave = reference_wave_list[fi]
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
            u = np.where(cc >= max(cc))
            lut_reference = rebin(np.roll(reference_wave, u[0]), [period_in_pixels])

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
        for k in range(4):
            noise_x = np.arange(nx4)*4 + k
            noise[..., noise_x] = dd_noise

        # Subtract EMI noise from the input data
        log.info('Subtracting EMI noise from data')
        corr_data = orig_data - noise
        output_model.data = corr_data

    if save_intermediate_results and save_onthefly_reffile is not None:
        if readpatt == 'FAST':
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
    masked_data = np.ma.array(data, mask=np.isnan(data))

    for i in range(nx):
        for j in range(ny):
            vec = masked_data[:, j, i]
            u = np.where(vec != 0)
            n = vec[u].size
            if n > 0:
                if n <= 2 or minval:
                    medimg[j, i] = np.ma.min(vec[u])
                if maxval:
                    medimg[j, i] = np.ma.max(vec[u])
                if not minval and not maxval and not avgval:
                    medimg[j, i] = np.ma.median(vec[u])
                if avgval:
                    dmean , _, _, _ = iter_stat_sig_clip(vec[u])
                    medimg[j, i] = dmean
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

    freqs_names : list
        List of frequency names to use
    """
    subname, rowclocks, frameclocks, frequencies = None, None, None, None

    # make sure the readpattern is defined as expected
    readpatt = readpatt.upper()
    if "SLOW" in readpatt:
        readpatt = "SLOW"
    elif "FAST" in readpatt:
        readpatt = "FAST"
    else:
        return subname, rowclocks, frameclocks, frequencies

    # search and return the specific values for the configuration
    if isinstance(subarray_cases, dict):
        for subname in subarray_cases:
            if subarray not in subname:
                continue
            if subname == 'FULL':
                subname = subname + '_' + readpatt
            rowclocks = subarray_cases[subname]["rowclocks"]
            frameclocks = subarray_cases[subname]["frameclocks"]
            if readpatt == "SLOW":
                frequencies = subarray_cases[subname]["freqs"]["SLOW"][detector]
            else:
                frequencies = subarray_cases[subname]["freqs"]["FAST"]
            break
    else:
        frequencies = []
        mdl_dict = subarray_cases.to_flat_dict()
        for subname in subarray_cases.subarray_cases:
            subname = subname.split(sep='.')[0]
            if subarray not in subname:
                continue
            log.debug('Found subarray case {}!'.format(subname))
            if 'FULL' in subname:
                subname = subname + '_' + readpatt
            for item, val in mdl_dict.items():
                if subname in item:
                    if "rowclocks" in item:
                        rowclocks = val
                    elif "frameclocks" in item:
                        frameclocks = val
                    else:
                        if readpatt == "SLOW" and "SLOW" in item and detector in item:
                            frequencies.append(val)
                        elif readpatt == "FAST"  and "FAST" in item:
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


def iter_stat_sig_clip(data, sigrej=3.0, maxiter=10):
    """ Compute the mean, mediand and/or sigma of data with iterative sigma clipping.
    This funtion is based on djs_iterstat.pro (authors therein)

    Parameters
    ----------
    data : numpy array

    sigrej : float
        Sigma for rejection

    maxiter: int
        Maximum number of sigma rejection iterations

    Returns
    -------
    dmean : float
        Computed mean

    dmedian : float
        Computed median

    dsigma : float
        Computed sigma

    dmask : numpy array
        Mask set to 1 for good points, and 0 for rejected points
    """
    # special cases of 0 or 1 data points
    ngood = np.size(data)
    dmean, dmedian, dsigma, dmask = 0.0, 0.0, 0.0, np.zeros(np.shape(data))
    if ngood == 0:
        log.debug('No data points for sigma clipping')
        return dmean, dmedian, dsigma, dmask
    elif ngood == 1:
        log.debug('Only 1 data point for sigma clipping')
        return dmean, dmedian, dsigma, dmask

    # Compute the mean + standard deviation of the entire data array,
    # these values will be returned if there are fewer than 2 good points.
    dmask = np.ones(ngood, dtype='b') + 1
    dmean = sum(data * dmask) / ngood
    dsigma = np.sqrt(sum((data - dmean)**2) / (ngood - 1))
    dsigma = dsigma
    iiter = 1

    # Iteratively compute the mean + stdev, updating the sigma-rejection thresholds
    # each iteration
    nlast = -1
    while iiter < maxiter and nlast != ngood and ngood >= 2:
        loval = dmean - sigrej * dsigma
        hival = dmean + sigrej * dsigma
        nlast = ngood

        dmask[data < loval] = 0
        dmask[data > hival] = 0
        ngood = sum(dmask)

        if ngood >= 2:
            dmean = sum(data*dmask) / ngood
            dsigma = np.sqrt( sum((data - dmean)**2 * dmask) / (ngood - 1) )
            dsigma = dsigma

        iiter += 1

    return dmean, dmedian, dsigma, dmask


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



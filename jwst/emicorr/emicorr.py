import datetime
import logging
import warnings

import numpy as np
from scipy import interpolate
from astropy.stats import sigma_clipped_stats as scs
from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

subarray_clocks = {
    "SLITLESSPRISM": {"rowclocks": 28, "frameclocks": 15904},
    "MASKLYOT": {"rowclocks": 90, "frameclocks": 32400},
    "SUB64": {"rowclocks": 28, "frameclocks": 8512},
    "SUB128": {"rowclocks": 44, "frameclocks": 11904},
    "MASK1140": {"rowclocks": 82, "frameclocks": 23968},
    "MASK1550": {"rowclocks": 82, "frameclocks": 23968},
    "MASK1065": {"rowclocks": 82, "frameclocks": 23968},
    "FULL_FAST": {"rowclocks": 271, "frameclocks": 277504},
    "FULL_SLOW": {"rowclocks": 2333, "frameclocks": 2388992},
    "BRIGHTSKY": {"rowclocks": 162, "frameclocks": 86528},
    "SUB256": {"rowclocks": 96, "frameclocks": 29952},
}


def apply_emicorr(
    input_model,
    emicorr_model,
    save_onthefly_reffile=None,
    algorithm="sequential",
    nints_to_phase=None,
    nbins=None,
    scale_reference=True,
    onthefly_corr_freq=None,
    use_n_cycles=3,
    fit_ints_separately=False,
):
    """
    Apply an EMI correction to MIRI ramps.

    The sequential algorithm corrects data for EMI by the following procedure:
         [repeat iteratively for each discrete EMI frequency desired]

    1) Make very crude slope image and fixed pattern superbias for each
       integration, ignoring everything (nonlin, saturation, badpix, etc).
    2) Subtract scaled slope image and bias from each frame of each integration.
    3) Calculate phase of every pixel in the image at the desired EMI frequency
       (e.g. 390 Hz) relative to the first pixel in the image.
    4) Make a binned, phased amplitude (PA) wave from the cleaned data.
    5) Either:
       -  Measure the phase shift between the binned pa wave and the input
          reference wave. Use look-up table to get the aligned reference wave
          value for each pixel (the noise array corresponding to the input image).
       Or:
       - Use the binned PA wave instead of the reference wave to "self-correct".
    6) Subtract the noise array from the input image and return the cleaned result.

         [repeat for next frequency using cleaned output as input]

    The joint algorithm proceeds similarly, except that the linear ramps and
    the EMI noise are fit simultaneously. It works by choosing pixels with modest
    scatter among the reads, and then finding the amplitude and phase of a supplied
    reference waveform that, when subtracted, makes these pixels' ramps as straight
    as possible. The straightness of the ramps is measured by chi squared after
    fitting lines to each one.  As for the sequential algorithm, the EMI signal at
    each frequency is fit and removed iteratively, for each desired frequency.

    For either algorithm, removing the highest amplitude EMI first will probably
    give best results.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        Input science data model to be EMI-corrected.
    emicorr_model : `~jwst.datamodels.EmiModel` or None
        Data model containing EMI correction.
    save_onthefly_reffile : str or None, optional
        Full path and root name to save an on-the-fly reference file to.
    algorithm : {'sequential', 'joint'}, optional
        Algorithm for fitting EMI noise in ramps.  If 'sequential', ramps are fit
        and then EMI noise is fit to the residuals.  If 'joint', ramps and noise
        are fit simultaneously.  The 'sequential' algorithm can be used without
        a reference waveform, generating a new reference file on-the-fly, but it
        requires 10 or more groups for a reliable fit.  The 'joint' algorithm
        requires a reference waveform but can successfully fit EMI in ramps with
        3 or more groups.
    nints_to_phase : int or None, optional
        Number of integrations to phase, when `algorithm` is 'sequential'.
    nbins : int or None, optional
        Number of bins to use in one phased wave, when `algorithm` is 'sequential'.
    scale_reference : bool, optional
        If True and `algorithm` is 'sequential', the reference wavelength will be scaled
        to the data's phase amplitude.
    onthefly_corr_freq : list of float or None, optional
        Frequency values to use to create a correction on-the-fly.  If provided,
        any input `emicorr_model` is ignored and the `algorithm` is set to
        'sequential'.
    use_n_cycles : int, optional
        Only use N cycles to calculate the phase to reduce code running time,
        when `algorithm` is 'sequential'.
    fit_ints_separately : bool, optional
        If True, fit each integration separately, when `algorithm` is 'joint'.

    Returns
    -------
    output_model : JWST data model
        Input science data model, corrected for EMI.
    """
    # Get the subarray case and other info
    detector = input_model.meta.instrument.detector
    subarray = input_model.meta.subarray.name
    readpatt = input_model.meta.exposure.readpatt

    # Get the number of samples, 10us sample times per pixel (1 for fastmode, 9 for slowmode)
    nsamples = input_model.meta.exposure.nsamples

    # Ignore emicorr model if user frequencies provided
    if onthefly_corr_freq is not None:
        emicorr_model = None

    # Check algorithm against input data shape
    ngroups = input_model.data.shape[1]
    if ngroups < 3:
        log.warning(f"EMI correction cannot be performed for ngroups={ngroups}")
        return None
    if ngroups < 10 and algorithm == "sequential":
        log.warning(f"The 'sequential' algorithm is selected and ngroups={ngroups}.")
        log.warning("The 'joint' algorithm is recommended for ngroups < 10.")
    if algorithm == "joint" and emicorr_model is None:
        log.warning(
            "The 'joint' algorithm cannot be used without an EMICORR "
            "reference file.  Setting the algorithm to 'sequential'."
        )
        algorithm = "sequential"

    # Get the subarray case from the ref file
    freq_numbers = []
    reference_wave_list = []
    subname, rowclocks, frameclocks, freqs2correct = None, None, None, None
    if emicorr_model is not None:
        log.info("Using reference file to get subarray case.")
        subname, rowclocks, frameclocks, freqs2correct = get_subarcase(
            emicorr_model, subarray, readpatt, detector
        )

        log.info(
            f"With configuration: Subarray={subarray}, Read_pattern={readpatt}, Detector={detector}"
        )
        if freqs2correct is not None:
            log.info(f"Will correct data for the following {len(freqs2correct)} frequencies: ")
            log.info(f"   {freqs2correct}")

            for freq_val in freqs2correct:
                # Collect the frequency numbers and reference wave list
                freq, ref_wave = get_frequency_info(emicorr_model, freq_val)
                freq_numbers.append(freq)
                reference_wave_list.append(ref_wave)
    else:
        # if we got here, the user requested to do correction with on-the-fly reference file
        subname = subarray
        if subname == "FULL":
            if "FAST" in readpatt.upper():
                subname += "_FAST"
            elif "SLOW" in readpatt.upper():
                subname += "_SLOW"
        if subname in subarray_clocks:
            rowclocks = subarray_clocks[subname]["rowclocks"]
            frameclocks = subarray_clocks[subname]["frameclocks"]
            freq_numbers = onthefly_corr_freq
            freqs2correct = [repr(freq) for freq in onthefly_corr_freq]

            log.info(f"Will correct data for the following {len(freqs2correct)} frequencies: ")
            log.info(f"   {freqs2correct}")

    # No subarray or read pattern match found, print to log and skip correction
    if rowclocks is None or len(freq_numbers) == 0:
        log.warning("No correction match for this configuration")
        return None

    # Run either the joint or sequential fitting algorithm
    # to fit and correct for EMI data
    log.info(f"Running EMI fit with algorithm = '{algorithm}'.")
    if algorithm == "joint":
        output_model = _run_joint_algorithm(
            input_model,
            freq_numbers,
            reference_wave_list,
            nsamples,
            rowclocks,
            frameclocks,
            fit_ints_separately=fit_ints_separately,
        )
    else:
        output_model = _run_sequential_algorithm(
            input_model,
            freqs2correct,
            freq_numbers,
            reference_wave_list,
            nsamples,
            rowclocks,
            frameclocks,
            save_onthefly_reffile=save_onthefly_reffile,
            nints_to_phase=nints_to_phase,
            nbins=nbins,
            scale_reference=scale_reference,
            use_n_cycles=use_n_cycles,
        )

    return output_model


def _run_joint_algorithm(
    input_model,
    freq_numbers,
    reference_wave_list,
    nsamples,
    rowclocks,
    frameclocks,
    fit_ints_separately=False,
):
    """
    Remove EMI noise with a joint fit to ramps and EMI signal.

    Parameters
    ----------
    input_model : DataModel
        Input datamodel.  Will be modified in place.
    freq_numbers : list of float
        Frequency values to correct.
    reference_wave_list : list of ndarray
        Reference waveforms.
    nsamples : int
        Number of samples of each pixel in each group:
        1 for fast, 9 for slow
    rowclocks : int
        Extra pixel times in each row before reading out the following row
    frameclocks : int
        Pixel clock cycles in a reset frame
    fit_ints_separately : bool, optional
        Fit the integrations separately? If True, fit amplitude and phase
        for refwave independently for each integration.  If False, fit
        for a single amplitude and phase across all integrations.

    Returns
    -------
    DataModel
        The input datamodel, with noise fit and subtracted.
    """
    # Additional frame time to account for the extra frame reset between
    # MIRI integrations
    readpatt = str(input_model.meta.exposure.readpatt).upper()
    if readpatt == "FASTR1" or readpatt == "SLOWR1":
        _frameclocks = frameclocks
    else:
        _frameclocks = 0

    for freq, ref_wave in zip(freq_numbers, reference_wave_list, strict=True):
        period_in_pixels = (1.0 / freq) / 10.0e-6

        corrected_data = emicorr_refwave(
            input_model.data,
            input_model.pixeldq,
            np.array(ref_wave),
            nsamples,
            rowclocks,
            _frameclocks,
            period_in_pixels,
            fit_ints_separately=fit_ints_separately,
        )

        # Data is updated in place, so it is corrected iteratively
        input_model.data[:] = corrected_data

    # Return the corrected model.
    return input_model


def _run_sequential_algorithm(
    input_model,
    freqs2correct,
    freq_numbers,
    reference_wave_list,
    nsamples,
    rowclocks,
    frameclocks,
    save_onthefly_reffile=None,
    nints_to_phase=None,
    nbins=None,
    scale_reference=True,
    use_n_cycles=3,
):
    """
    Remove EMI noise with a sequential fit to ramps and EMI signal.

    Parameters
    ----------
    input_model : DataModel
        Input datamodel.  Will be modified in place.
    freqs2correct : list of str
        Names of the frequencies to correct.
    freq_numbers : list of float
        Frequency values to correct.
    reference_wave_list : list of ndarray
        Reference waveforms.
    nsamples : int
        Number of samples of each pixel in each group:
        1 for fast, 9 for slow
    rowclocks : int
        Extra pixel times in each row before reading out the following row
    frameclocks : int
        Pixel clock cycles in a reset frame
    save_onthefly_reffile : str or None, optional
        Full path and root name to save an on-the-fly reference file to.
    nints_to_phase : int or None, optional
        Number of integrations to phase, when `algorithm` is 'sequential'.
    nbins : int or None, optional
        Number of bins to use in one phased wave, when `algorithm` is 'sequential'.
    scale_reference : bool, optional
        If True and `algorithm` is 'sequential', the reference wavelength will be scaled
        to the data's phase amplitude.
    use_n_cycles : int, optional
        Only use N cycles to calculate the phase to reduce code running time,
        when `algorithm` is 'sequential'.

    Returns
    -------
    DataModel
        The output datamodel, with noise fit and subtracted.

    Notes
    -----
    The 'sequential' algorithm was originally translated from the IDL
    procedure 'fix_miri_emi.pro', written by E. Bergeron.
    """
    # Get the shape and metadata for the input data
    nints, ngroups, ny, nx = np.shape(input_model.data)
    xsize = input_model.meta.subarray.xsize  # SUBSIZE1 keyword
    xstart = input_model.meta.subarray.xstart  # SUBSTRT1 keyword
    detector = input_model.meta.instrument.detector
    subarray = input_model.meta.subarray.name
    readpatt = input_model.meta.exposure.readpatt

    # create the dictionary to store the frequencies and corresponding phase amplitudes
    if save_onthefly_reffile is not None:
        freq_pa_dict = {"frequencies": {}, "subarray_cases": {}}

    # Copy the input nbins setting, so it can be set to different values for different frequencies
    nbins_all = nbins

    # Loop over the frequencies to correct
    for fi, frequency_name in enumerate(freqs2correct):
        frequency = freq_numbers[fi]
        log.info(
            f"Correcting for frequency: {frequency} Hz  ({fi + 1} out of {len(freqs2correct)})"
        )

        # Set up some variables

        # Correspondence of array order in IDL
        # sz[0] = 4 in idl
        # sz[1] = nx
        # sz[2] = ny
        # sz[3] = ngroups
        # sz[4] = nints
        nx4 = int(nx / 4)

        dd_all = np.zeros((nints, ngroups, ny, nx4))
        log.info("Subtracting self-superbias from each group of each integration")

        # Calculate times of all pixels in the input integration, then use that to calculate
        # phase of all pixels. Times here is in integer numbers of 10us pixels starting from
        # the first data pixel in the input image. Times can be a very large integer, so use
        # a big datatype. Phaseall (below) is just 0-1.0.

        # A safer option is to calculate times_per_integration and calculate the phase at each
        # int separately. That way times array will have a smaller number of elements at each
        # int, with less risk of datatype overflow. Still, use the largest datatype available
        # for the time_this_int array.

        times_this_int = np.zeros((ngroups, ny, nx4), dtype="ulonglong")
        phaseall = np.zeros((nints, ngroups, ny, nx4))

        # non-roi rowclocks between subarray frames (this will be 0 for fullframe)
        extra_rowclocks = (1024.0 - ny) * (4 + 3.0)
        # e.g. ((1./390.625) / 10e-6) = 256.0 pix and ((1./218.52055) / 10e-6) = 457.62287 pix
        period_in_pixels = (1.0 / frequency) / 10.0e-6

        if nints_to_phase is None and use_n_cycles is None:  # user wants to use all integrations
            # use all integrations
            nints_to_phase = nints
        elif (
            nints_to_phase is None and use_n_cycles is not None
        ):  # user wants to use nints_to_phase
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
        colstop = int(xsize / 4 + xstart - 1)
        log.info("Doing phase calculation per integration")

        for ninti in range(nints):
            log.debug(f"  Working on integration: {ninti + 1}")
            # Read in this integration
            data = input_model.data[ninti].copy()

            # Remove source signal and fixed bias from each integration ramp
            # (linear is good enough for phase finding)

            # do linear fit for source + sky
            s0, mm0 = sloper(data[1 : ngroups - 1, :, :])

            # subtract source+sky from each frame of this ramp
            for ngroupi in range(ngroups):
                data[ngroupi, ...] = input_model.data[ninti, ngroupi, ...] - (s0 * ngroupi)

            # make a self-superbias
            m0 = minmed(data[1 : ngroups - 1, :, :])

            # subtract self-superbias from each frame of this ramp
            for ngroupi in range(ngroups):
                data[ngroupi, ...] = data[ngroupi, ...] - m0

                # de-interleave each frame into the 4 separate output channels and
                # average (or median) them together for S/N
                d0 = data[ngroupi, :, 0:nx:4]
                d1 = data[ngroupi, :, 1:nx:4]
                d2 = data[ngroupi, :, 2:nx:4]
                d3 = data[ngroupi, :, 3:nx:4]
                dd = (d0 + d1 + d2 + d3) / 4.0

                # fix a bad ref col
                dd[:, 1] = (dd[:, 0] + dd[:, 3]) / 2
                dd[:, 2] = (dd[:, 0] + dd[:, 3]) / 2
                # This is the quad-averaged, cleaned, input image data for the exposure
                dd_all[ninti, ngroupi, ...] = dd - np.median(dd)

            for k in range(ngroups):  # frames
                for j in range(ny):  # rows
                    # nsamples= 1 for fast, 9 for slow (from metadata)
                    times_this_int[k, j, :] = (
                        np.arange(nx4, dtype="ulonglong") * nsamples + start_time
                    )

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
                        times_this_int[k, j, nx4 - 1] = (
                            times_this_int[k, j, nx4 - 1] + ulonglong_ref_pix_sample
                        )

                    # point to the first pixel of the next row (add "end-of-row" pad)
                    start_time += rowclocks

                # point to the first pixel of the next frame (add "end-of-frame" pad)
                start_time += extra_rowclocks

            # Convert "times" to phase each integration. Note that times has units of
            # number of 10us from the first data pixel in this integration, so to
            # convert to phase, divide by the waveform *period* in float pixels
            phase_this_int = times_this_int / period_in_pixels
            phaseall[ninti, ...] = phase_this_int - phase_this_int.astype("ulonglong")

            # add a frame time to account for the extra frame reset between MIRI integrations
            if readpatt.upper() == "FASTR1" or readpatt.upper() == "SLOWR1":
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
            nbins = int(period_in_pixels / 2.0)
        else:
            nbins = nbins_all
        if nbins > 501:
            nbins = 500

        # bin the whole set
        log.info(f"Calculating the phase amplitude for {nbins} bins")
        # Define the binned waveform amplitude (pa = phase amplitude)
        pa = np.arange(nbins, dtype=float)
        # keep track of n per bin to check for low n
        nb_over_nbins = [nb / nbins for nb in range(nbins)]
        nbp1_over_nbins = [(nb + 1) / nbins for nb in range(nbins)]
        # Construct a phase map and dd map for only the nints_to_phase
        phase_temp = phaseall[0:nints_to_phase, :, :, :]
        dd_temp = dd_all[0:nints_to_phase, :, :, :]
        for nb in range(nbins):
            u = (phase_temp > nb_over_nbins[nb]) & (phase_temp <= nbp1_over_nbins[nb])
            # calculate the sigma-clipped mean
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                dmean, _, _ = scs(dd_temp[u])
            pa[nb] = dmean  # amplitude in this bin

        pa -= np.median(pa)

        # pa_phase is the phase corresponding to each bin in pa. The +0.5 here is to center in
        # the bin. Together, (pa, pa_phase) define the waveform to phase against the
        # reference_wave, or to use directly for a self correction (likely noiser than the
        # reference_wave, but perfectly matched in phase by definition so no need for a
        # crosscor or other phase-finder against a reference_wave). This may not be
        # necessary. Can just assume pa and the reference_wave are exactly 1.0 waves long.
        # pa_phase = (np.arange(nbins, dtype=np.float32) + 0.5) / nbins

        # Correction using direct look-up from either pa (self-correction) or
        # the reference wave (if provided)
        # Make a vector to hold a version of pa that is close to an integer number of
        # pixels long. This will be the noise amplitude look-up table.
        lut = rebin(pa, [period_in_pixels])  # shift and resample reference_wave at pa's phase

        # If a reference waveform was provided, use it to measure phase shift between
        # pa (binned phased wave) from the input image vs. the reference_wave using xcorr.
        # These two methods give the same integer-reference_wave-element resolution results.
        if len(reference_wave_list) > 0:
            log.info("Using reference file to measure phase shift")
            reference_wave = np.array(reference_wave_list[fi])
            reference_wave_size = np.size(reference_wave)
            rebinned_pa = rebin(pa, [reference_wave_size])
            cc = np.zeros(reference_wave_size)
            for i in range(reference_wave_size):
                shifted_ref_wave = np.roll(reference_wave, i)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
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
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                m, b = np.polyfit(lut_reference, lut, 1)  # the slope is the amplitude ratio

            # check pa-to-reference_wave phase match and optional scale
            if scale_reference:
                lut_reference = lut_reference * m

            lut = lut_reference

        if save_onthefly_reffile is not None:
            freq_pa_dict["frequencies"][frequency_name] = {
                "frequency": frequency,
                "phase_amplitudes": pa,
            }

        log.info("Creating phased-matched noise model to subtract from data")
        # This is the phase matched noise model to subtract from each pixel of the input image
        dd_noise = lut[(phaseall * period_in_pixels).astype(int)]

        # Safety catch; anywhere the noise value is not finite, set it to zero
        dd_noise[~np.isfinite(dd_noise)] = 0.0

        # Subtract EMI noise from the input data
        log.info("Subtracting EMI noise from data")

        # Interleave (straight copy) into 4 amps
        for k in range(4):
            input_model.data[..., k::4] -= dd_noise

        # clean up
        del dd_all
        del times_this_int
        del phaseall
        del dd_noise

    if save_onthefly_reffile is not None:
        if "FAST" in readpatt:
            freqs_dict = {readpatt: freqs2correct}
        else:
            freqs_dict = {readpatt: {detector: freqs2correct}}
        on_the_fly_subarr_case = {}
        on_the_fly_subarr_case[subarray] = {
            "rowclocks": rowclocks,
            "frameclocks": frameclocks,
            "freqs": freqs_dict,
        }
        freq_pa_dict["subarray_cases"] = on_the_fly_subarr_case
        mk_reffile(freq_pa_dict, save_onthefly_reffile)

    return input_model


def sloper(data):
    """
    Fit slopes to all pix of a ramp.

    Parameters
    ----------
    data : ndarray
        3-D integration data array

    Returns
    -------
    outarray : ndarray
        3-D integration data array
    intercept: ndarray
        Slope intercept values
    """
    ngroups, ny, nx = np.shape(data)
    frametime = 1.0  # use 1.0 for per-frame slopes, otherwise put a real time here
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
        sxx = sxx + grouptimes[groupcount] ** 2

    # calculate the final per-pixel slope values
    outarray = (ngroups * sxy - sx * sy) / (ngroups * sxx - sx * sx)

    # calculate the final per-pixel intercept values
    intercept = (sxx * sy - sx * sxy) / (ngroups * sxx - sx * sx)

    return outarray, intercept


def minmed(data):
    """
    Calculate the median image from a stack of images.

    If there are 2 or fewer frames, the minimum image is computed
    instead of the median.

    Parameters
    ----------
    data : ndarray
        3D data array, to be stacked along the first axis.

    Returns
    -------
    medimg : ndarray
        2D median or minimum image.
    """
    if data.shape[0] <= 2:
        medimg = np.nanmin(data, axis=0)
    else:
        medimg = np.nanmedian(data, axis=0)
    return medimg


def get_subarcase(emi_model, subarray, readpatt, detector):
    """
    Get the rowclocks and frameclocks values for the given configuration.

    If no match is found for the configuration, None will be returned
    for all values.

    Parameters
    ----------
    emi_model : EmiModel
        EMI datamodel containing subarray_cases.
    subarray : str
        Keyword value
    readpatt : str
        Keyword value
    detector : str
        Keyword value

    Returns
    -------
    subname : str or None
        Modified subarray name.
    rowclocks : int or None
        Row clock value.
    frameclocks : int or None
        Frame clock value.
    frequencies : list of str or None
        List of frequencies to correct according to subarray name.
    """
    subname, rowclocks, frameclocks, frequencies = None, None, None, None

    # make sure the read pattern is defined as expected to read data from reference file
    readpatt = readpatt.upper()
    if "SLOW" in readpatt:
        readout_speed = "SLOW"
    elif "FAST" in readpatt:
        readout_speed = "FAST"
    else:
        log.warning(f"Read pattern {readpatt} does not include expected string FAST or SLOW")
        return subname, rowclocks, frameclocks, frequencies

    # modify the subarray name for matching if needed,
    # appending FAST or SLOW for full frame data
    subname = subarray
    if "FULL" in subarray:
        subname = f"{subarray}_{readout_speed}"

    try:
        subarray_case = getattr(emi_model.subarray_cases, subname)
        rowclocks = subarray_case.rowclocks
        frameclocks = subarray_case.frameclocks
        freqs = subarray_case.freqs
        rp_freqs = getattr(freqs, readout_speed)
        if readout_speed == "FAST":
            frequencies = rp_freqs
        else:
            frequencies = getattr(rp_freqs, detector)
    except AttributeError:
        # match not found, will return None for all values
        log.debug(f"Subarray {subname} not found")
        pass

    return subname, rowclocks, frameclocks, frequencies


def get_frequency_info(emi_model, frequency_name):
    """
    Get the frequency number from the given EMI model.

    Parameters
    ----------
    emi_model : EmiModel
        EMI reference datamodel.
    frequency_name : str
        Frequency of interest.

    Returns
    -------
    frequency_number : float
        Frequency
    phase_amplitudes : array
        1-D array of the corresponding phase amplidues for this frequency

    Raises
    ------
    AttributeError
        If the `frequency_name` was not found in the `emi_model`.
    """
    freq_set = getattr(emi_model.frequencies, frequency_name)
    freq_number = freq_set.frequency
    phase_amplitudes = freq_set.phase_amplitudes
    return freq_number, phase_amplitudes


def rebin(arr, newshape):
    """
    Rebin an array to a new shape.

    Rebinned values are taken from the input at nearest-neighbor indices.
    No interpolation is performed.  If a value falls between two input
    pixels, the smaller index value is used.

    Parameters
    ----------
    arr : ndarray
        Array to rebin.
    newshape : list or tuple
        New shape for the input array.

    Returns
    -------
    arr : ndarray
        Rebinned array.

    Notes
    -----
    The intent for this function is to replicate the behavior of the
    IDL function `congrid`, without the /INTERP flag.
    """
    if len(arr.shape) != len(newshape):
        raise ValueError("New shape dimensions must match old shape")

    # get coordinates for old values in the new array
    slices = [
        slice(0, old, float(old) / new) for old, new in zip(arr.shape, newshape, strict=False)
    ]
    coordinates = np.mgrid[slices]

    # take the nearest smaller integer value (truncate to an integer index)
    indices = coordinates.astype("i")

    # take elements from the array at the specified indices
    return arr[tuple(indices)]


def mk_reffile(freq_pa_dict, emicorr_ref_filename):
    """
    Create the reference file hdulist object.

    Parameters
    ----------
    freq_pa_dict : dict
        Dictionary containing the phase amplitude (pa) array
        corresponding to the appropriate frequency
    emicorr_ref_filename : str
        Full path and root name of on-the-fly reference file
    """
    # Save the reference file if desired
    emicorr_model = datamodels.EmiModel()
    emicorr_model.frequencies = freq_pa_dict["frequencies"]
    emicorr_model.subarray_cases = freq_pa_dict["subarray_cases"]

    # Add some placeholder values required for validation
    emicorr_model.meta.reftype = "emicorr"
    emicorr_model.meta.author = "JWST calibration pipeline"
    emicorr_model.meta.description = "EMI correction file"
    emicorr_model.meta.pedigree = "GROUND"
    emicorr_model.meta.useafter = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%S")

    if emicorr_ref_filename.endswith(".fits"):
        emicorr_ref_filename = emicorr_ref_filename.replace(".fits", ".asdf")
    log.info(f"Writing on-the-fly reference file to: {emicorr_ref_filename}")
    emicorr_model.save(emicorr_ref_filename)
    emicorr_model.close()


def emicorr_refwave(
    data,
    pdq,
    refwave,
    nsamples,
    rowclocks,
    frameclocks,
    period_in_pixels,
    fit_ints_separately=False,
    nphases_opt=500,
):
    """
    Derive the best amplitude and phase for the EMI waveform, subtract it off.

    This routine works by choosing pixels with modest scatter among
    the reads, and then finding the amplitude and phase of a supplied
    reference waveform that, when subtracted, makes these pixels'
    ramps as straight as possible. The straightness of the ramps is
    measured by chi squared after fitting lines to each one.

    Parameters
    ----------
    data : ndarray
        4D data to correct, shape (nints, ngroups, nrows, ncols).
        Modified in-place.
    pdq : ndarray
        2D pixeldq array, shape (nrows, ncols)
    refwave : ndarray
        Reference waveform, 1D array
    nsamples : int
        Number of samples of each pixel in each group:
        1 for fast, 9 for slow
    rowclocks : int
        Extra pixel times in each row before reading out the following row
    frameclocks : int
        Pixel clock cycles in a reset frame
    period_in_pixels : float
        Period in pixel clock times of the waveform given by refwave
    fit_ints_separately : bool, optional
        Fit the integrations separately? If True, fit amplitude and phase
        for refwave independently for each integration.  If False, fit
        for a single amplitude and phase across all integrations.
    nphases_opt : int, optional
        Number of phases to sample chi squared as a function of phase

    Returns
    -------
    data : ndarray
        The input 4D data array with the model waveform phased, scaled,
        and subtracted.
    """
    nints, ngroups, ny, nx = data.shape

    # divide ncols by the four output channels (which are read simultaneously)
    nx4 = nx // 4

    # non-roi rowclocks between subarray frames (this will be 0 for fullframe)
    extra_rowclocks = (1024.0 - ny) * (4 + 3.0)

    frametime = ny * rowclocks + extra_rowclocks

    # Times of the individual reads in pixel clock times
    t0_arr = np.zeros((ny, nx4))
    for i in range(ny):
        t0_arr[i] = i * rowclocks + np.arange(nx4) * nsamples
    phase = (t0_arr / period_in_pixels) % 1

    # Phase gap between groups
    dphase = (frametime / period_in_pixels) % 1

    # Phase gap between integrations
    dphase_frame = dphase * ngroups + frameclocks / period_in_pixels

    # Time-like index for group number
    grouptimes = np.arange(ngroups)

    # Phase(model) for interpolation.  Wrap the phases to avoid
    # extrapolation, and use a cubic spline.

    refwave_extended = np.array([refwave[-1]] + list(refwave) + [refwave[0]])
    phase_extended = (np.arange(len(refwave) + 2) - 0.5) / len(refwave)
    phasefunc = interpolate.interp1d(phase_extended, refwave_extended, kind="cubic")

    phases_template = (phase_extended[1:-1, np.newaxis] + grouptimes * dphase) % 1
    nphases = phases_template.shape[0]

    # These arrays hold the sum of the values of good pixels at a
    # given phase and the total number of good pixels at a given
    # phase, respectively.  The phase refers to the first group; other
    # groups will have the appropriate phase delay added.

    all_y = np.zeros((nints, nphases, ngroups))
    all_n = np.zeros((nints, nphases))

    # "Good" pixel here has no more than twice the median standard
    # deviation among group values and is not flagged in the pdq
    # array.  This should discard most bad and high-flux pixels.

    pixel_std = np.std(data, axis=1)
    pixel_ok = (pixel_std < 2 * np.median(pixel_std)) & (pdq == 0)

    for j in range(ny):
        for k in range(nx4):
            # Bad reference column
            if k == 1 or k == 2:
                continue

            # Choose the index corresponding to our phase.
            # phases_template has the midpoints of the intervals,
            # so rounding down here is appropriate.

            indx = int(phase[j, k] * nphases)
            # All four output channels.
            for l in range(4):
                pixok = pixel_ok[:, j, k * 4 + l]
                all_y[:, indx] += data[:, :, j, k * 4 + l] * pixok[:, np.newaxis]
                all_n[:, indx] += pixok

    # We'll compute chi2 at nphases_opt evenly spaced phases.
    phaselist = np.arange(nphases_opt) * 1.0 / nphases_opt

    # Class that computes and holds all of the intermediate
    # information needed for the fits.

    emifitter = EMIfitter(all_y, all_n, phasefunc, phases_template, phaselist, dphase_frame)

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

    # One phase and amplitude fit to all integrations
    else:
        chisq, _ = calc_chisq_amplitudes(emifitter)
        best_phase = get_best_phase(phaselist, chisq)
        _, c = calc_chisq_amplitudes(emifitter, phases=[best_phase])

        # Phases appropriate for each integration
        phases_to_correct = (best_phase + np.arange(nints) * dphase_frame) % 1

        # Output from calc_chisq_amplitudes here is a one-element list.
        # Repeat it for each integration.

        amplitudes_to_correct = c * nints

    for i in range(nints):
        for j in range(ngroups):
            # Place the reference waveform at the appropriate phase,
            # scale, and subtract from each output channel.

            phased_emi = phasefunc((phase + dphase * j + phases_to_correct[i]) % 1)
            for k in range(4):
                data[i, j, :, k::4] -= amplitudes_to_correct[i] * phased_emi

    return data


def get_best_phase(phases, chisq):
    """
    Fit a parabola to get the phase corresponding to the best chi squared.

    Parameters
    ----------
    phases : ndarray
        1D array of phases. Assumed to be periodic: phases[0] comes
        after phases[-1]
    chisq : ndarray
        1D array of chi squared values.

    Returns
    -------
    bestphase : float
        Phase where chisq is minimized, computed by fitting a parabola
        to the three points centered on the input phase with the
        lowest chisq.
    """
    if np.all(~np.isfinite(chisq)):
        log.warning("Chi squared values are all NaN")
        return phases[0]

    ibest = np.nanargmin(chisq)

    # If ibest is zero, ibest_m1 will give the final element.
    ibest_m1 = ibest - 1
    ibest_p1 = ibest + 1

    # If ibest is the last element, we have to explicitly wrap.
    if ibest_p1 > len(chisq) - 1:
        ibest_p1 = 0

    # Calculate the vertex of the parabola.
    chisq_opt = [chisq[ibest_m1], chisq[ibest], chisq[ibest_p1]]
    if np.any(~np.isfinite(chisq_opt)):
        log.warning("Chi squared values in final optimization are unexpectedly NaN")
        return phases[ibest]

    x = [phases[ibest_m1], phases[ibest], phases[ibest_p1]]

    # Just in case all of these chi squared values are equal
    # (in which case fitting a parabola would give a=0)
    if np.nanstd(chisq_opt) == 0:
        log.warning(
            "Chi squared values as a function of phase offset are unexpectedly equal in EMIcorr"
        )
        return phases[ibest]

    a, b, _ = np.polyfit(x, chisq_opt, 2)
    bestphase = -b / (2 * a)

    return bestphase


def calc_chisq_amplitudes(emifitter, ints=None, phases=None):
    """
    Compute chi2 and amplitude of EMI waveform at phase(s).

    This routine has the math from the writeup in it; see that for the
    definitions and derivations.

    Parameters
    ----------
    emifitter : EMIFitter
        Convenience object containing all the parameters and intermediate
        arrays needed to compute the best amplitude and the corresponding
        chi squared value at one or more input phases.
    ints : list of int or None
        Integration(s) to combine for computing amplitude and chi
        squared. If None, use all integrations.
    phases : list of float or None
        Phase(s) at which to compute the amplitude and best chi squared.
        If None, use all phases in emifitter.phaselist.

    Returns
    -------
    chisq : list of float
        Best chi squared value for each phase used. Will be the same
        length as phases (or emifitter.phaselist, if phases is None).
    amplitudes : list of float
        Best-fit amplitudes for each phase used. Will be the same
        length as phases (or emifitter.phaselist, if phases is None).
    """
    ef = emifitter  # Alias to simplify subsequent references

    # By default, compute a single phase and amplitudes over all
    # integrations.

    if ints is None:
        ints = np.arange(ef.nints)

    # By default, calculate chi squared and the best-fit amplitude
    # for every phase in the input EMIfitter's phaselist.

    if phases is None:
        phases = ef.phaselist

    chisq = []
    amplitudes = []

    # Compute the best chi squared and the best amplitude at every
    # requested phase using the math in the writeup.

    for phase in phases:
        a_ = 0
        b_ = 0

        for i in ints:
            y = ef.all_y[i]
            n_ = ef.all_n[i]
            s_y = ef.all_sy[i]
            s_ty = ef.all_sty[i]

            # Phase difference between the start of this integration
            # and the start of the first integration

            phase_diff = ef.dphase_frame * i

            # Choose the closest phase in emifitter's phaselist

            k = np.argmin(np.abs((ef.phaselist - phase - phase_diff) % 1))

            z = ef.zlist[k]
            s_tz = ef.stzlist[k]
            s_z = ef.szlist[k]
            s_zz = ef.szzlist[k]
            s_yz = np.sum(z * y, axis=1)

            a_ += np.sum(
                n_
                / ef.delta
                * (
                    -ef.s_tt * s_z**2
                    + 2 * ef.s_t * s_z * s_tz
                    - ef.ngroups * s_tz**2
                    + s_zz * ef.delta
                )
            )
            b_ += (2 / ef.delta) * np.sum(
                ef.s_tt * s_z * s_y
                - ef.s_t * s_tz * s_y
                - ef.s_t * s_z * s_ty
                + ef.ngroups * s_tz * s_ty
            )
            b_ -= 2 * np.sum(s_yz)

        if a_ == 0 or ~np.isfinite(a_) or ~np.isfinite(b_):
            chisq.append(np.nan)
            amplitudes.append(0.0)
        else:
            c = -b_ / (2 * a_)
            chisq.append(a_ * c**2 + b_ * c)
            amplitudes.append(c)

    return chisq, amplitudes


class EMIfitter:
    """Convenience class to facilitate chi2, amplitude calculation."""

    def __init__(self, all_y, all_n, phasefunc, phases_template, phaselist, dphase_frame):
        """
        Compute and store quantities needed for chi2, amplitude calculation.

        This class precomputes a series of arrays so that the sums in
        `calc_chisq_amplitudes` do not need to be fully recomputed at every
        possible phase of the EMI signal.

        Sums of pixel read values are precomputed over pixels that share the
        same EMI phase to avoid double sums over pixels and reads later. Sums over
        the EMI waveform itself are also precomputed at each trial phase offset for
        the same reason.

        Parameters
        ----------
        all_y : ndarray
            3D array of phased, summed y values,
            shape (nints, nphases, ngroups)
        all_n : ndarray
            2D array of the number of pixels used for the calculation
            shape (nints, nphases): number of values summed to create
            all_y at each group
        phasefunc : callable
            Spline function to compute the waveform as a function of
            input phase (input phase should be between 0 and 1).
        phases_template : ndarray
            2D array of phases corresponding to all_y[0]
        phaselist : ndarray
            1D array of phases at which to pre-compute quantities
            needed for chi squared, amplitude calculation.  Typically
            uniformly spaced between 0 and 1.
        dphase_frame : float
            Phase difference between successive integrations
        """
        self.all_y = all_y
        self.all_n = all_n

        self.nints, self.nphases, self.ngroups = all_y.shape

        self.grouptimes = np.arange(self.ngroups)
        self.s_tt = np.sum(self.grouptimes**2)
        self.s_t = np.sum(self.grouptimes)
        self.delta = self.ngroups * self.s_tt - self.s_t**2

        self.all_sy = np.sum(self.all_y, axis=2)
        self.all_sty = np.sum(self.all_y * self.grouptimes, axis=2)

        self.phaselist = phaselist
        self.phases_template = phases_template

        self.phasefunc = phasefunc
        self.dphase_frame = dphase_frame

        self.zlist = []
        self.szlist = []
        self.stzlist = []
        self.szzlist = []

        # Waveform for each pixel when the first one is at dphaseval.
        for dphaseval in phaselist:
            # The transposes help ensure that similar phases are evaluated
            # consecutively, which significantly improves runtime when
            # there is a very large number of groups.

            z = self.phasefunc((self.phases_template.T + dphaseval) % 1).T

            self.zlist += [z]
            self.stzlist += [np.sum(self.grouptimes * z, axis=1)]
            self.szlist += [np.sum(z, axis=1)]
            self.szzlist += [np.sum(z**2, axis=1)]

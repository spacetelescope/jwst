import logging
import warnings

import numpy as np

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.datamodels.utils.flat_multispec import expand_flat_spec
from jwst.extract_1d.spec_wcs import create_spectral_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# attributes that we want to copy unmodified from input to output spectra
SPECMETA_ATTRIBUTES = [
    "source_id",
    "dispersion_direction",
    "source_type",
    "source_ra",
    "source_dec",
]


class InputSpectrumModel:
    """
    Model an input spectrum.

    Attributes
    ----------
    wavelength : ndarray
        Input wavelength.
    flux : ndarray
        Input flux.
    flux_error : ndarray
        Input error on the flux.
    surf_bright : ndarray
        Input surface brightness.
    sb_error : ndarray
        Input error on the surface brightness.
    dq : ndarray
        Input DQ array.
    nelem : int
        Number of spectral elements.
    weight : ndarray
        Weight value for each spectral element.
    unit_weight : bool
        Flag to indicate uniform weights are used.
    right_ascension : ndarray
        RA value for each spectral element.
    declination : ndarray
        Dec value for each spectral element.
    name : str
        Slit name for the spectrum.
    source_id : int
        Source ID for the spectrum.
    source_type : str
        Source type for the spectrum.
    source_ra : float
        Right ascension of the source.
    source_dec : float
        Declination of the source.
    dispersion_direction : str
        Dispersion direction for the spectrum.
    flux_unit : str
        Unit for the flux values.
    sb_unit : str
        Unit for the surface brightness values.
    """

    def __init__(self, ms, spec, exptime_key):
        """
        Create an InputSpectrumModel object.

        Parameters
        ----------
        ms : `~jwst.datamodels.JwstDataModel`, MultiSpecModel or SpecModel
            This is used to get the integration time.

        spec : `~jwst.datamodels.JwstDataModel`, SpecModel table
            The table containing columns "wavelength" and "flux".
            The `ms` object may contain more than one spectrum, but `spec`
            should be just one of those.

        exptime_key : str
            A string identifying which keyword to use to get the exposure
            time, which is used as a weight; or "unit_weight", which means
            to use weight = 1.
        """
        self.wavelength = spec.spec_table.field("wavelength").copy()

        self.flux = spec.spec_table.field("flux").copy()
        self.flux_error = spec.spec_table.field("flux_error").copy()
        self.surf_bright = spec.spec_table.field("surf_bright").copy()
        self.sb_error = spec.spec_table.field("sb_error").copy()
        self.dq = spec.spec_table.field("dq").copy()
        self.nelem = self.wavelength.shape[0]
        self.unit_weight = False  # may be reset below
        self.right_ascension = np.zeros_like(self.wavelength)
        self.declination = np.zeros_like(self.wavelength)
        self.name = spec.name
        for attr in SPECMETA_ATTRIBUTES:
            setattr(self, attr, getattr(spec, attr))
        self.flux_unit = spec.spec_table.columns["flux"].unit
        self.sb_unit = spec.spec_table.columns["surf_bright"].unit

        self.weight = np.ones_like(self.wavelength)
        if exptime_key == "integration_time":
            self.weight *= ms.meta.exposure.integration_time
        elif exptime_key == "exposure_time":
            self.weight *= ms.meta.exposure.exposure_time
        elif exptime_key == "unit_weight":
            self.unit_weight = True
        else:
            raise RuntimeError(f"Don't understand exptime_key = '{exptime_key}'")

        try:
            self.right_ascension[:], self.declination[:], _ = spec.meta.wcs(0.0)
        except AttributeError:
            self.right_ascension[:] = ms.meta.target.ra
            self.declination[:] = ms.meta.target.dec
            # This exception is hit for NIRISS and NIRCam WFSS data,
            # for which it doesn't matter anyway, since the RA and Dec are not used
            # in any meaningful way to combine the spectra. A future refactor should
            # make it so the WCS is not expected in the input spectra for those modes.
            log.debug("There is no WCS in the input. Getting RA, Dec from target metadata.")

    def close(self):
        """Set data attributes to null values."""
        self.wavelength = None
        self.flux = None
        self.flux_error = None
        self.surf_bright = None
        self.sb_error = None
        self.dq = None
        self.nelem = 0
        self.weight = 1.0
        self.unit_weight = False
        self.right_ascension = None
        self.declination = None
        self.source_id = None


class OutputSpectrumModel:
    """
    Model an output spectrum.

    Attributes
    ----------
    wavelength : ndarray
        Output wavelength.
    flux : ndarray
        Output flux.
    flux_error : ndarray
        Output error on the flux.
    surf_bright : ndarray
        Output surface brightness.
    sb_error : ndarray
        Output error on the surface brightness.
    dq : ndarray
        Output DQ array.
    weight : ndarray
        Weight value for each spectral element.
    count : ndarray
        Input value count for each output spectral element.
    wcs : gwcs.WCS
        Output spectral WCS.
    normalized : bool
        Flag to indicate data has been combined (sums are normalized).
    """

    def __init__(self):
        self.wavelength = None
        self.flux = None
        self.flux_error = None
        self.surf_bright = None
        self.sb_error = None
        self.dq = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.normalized = False
        self.source_id = None
        self.flux_unit = None
        self.sb_unit = None

    def assign_wavelengths(self, input_spectra):
        """
        Create an array of wavelengths to use for the output spectrum.

        Take the union of all input wavelengths, then call compute_output_wl
        to bin wavelengths in groups of the number of overlapping spectra.

        Parameters
        ----------
        input_spectra : list of InputSpectrumModel objects
            List of input spectra.
        """
        (wl, n_input_spectra) = count_input(input_spectra)

        self.wavelength = np.sort(compute_output_wl(wl, n_input_spectra))

        self.wcs = create_spectral_wcs(
            input_spectra[0].right_ascension[0], input_spectra[0].declination[0], self.wavelength
        )

    def accumulate_sums(self, input_spectra, sigma_clip=None):
        """
        Compute a weighted sum of all the input spectra.

        Each pixel of each input spectrum will be added to one pixel of
        the output spectrum.  The wavelength spacing of the input and
        output should be comparable, so each input pixel will usually
        correspond to one output pixel (i.e. have the same wavelength,
        to within half a pixel).  However, for any given input spectrum,
        there can be output pixels that are incremented by more than one
        input pixel, and there can be output pixels that are not incremented.
        If there are several input spectra, such gaps will hopefully be
        filled in.

        Pixels in an input spectrum that are flagged (via the DQ array)
        as DO_NOT_USE will not be used, i.e. they will simply be ignored
        when incrementing output arrays from an input spectrum.  If an
        input pixel is flagged with a non-zero value other than DO_NOT_USE,
        the value will be propagated to the output DQ array via bitwise OR.

        Parameters
        ----------
        input_spectra : list of InputSpectrumModel objects
            List of input spectra.
        sigma_clip : float, optional
            Factor for clipping outliers in spectral combination.
        """
        # This is the data type for the output spectrum.  We'll use double
        # precision for accumulating sums for most columns, but for the DQ
        # array, use the correct output data type.
        cmb_dtype = datamodels.CombinedSpecModel().spec_table.dtype
        dq_dtype = cmb_dtype.fields["DQ"][0]

        nelem = self.wavelength.shape[0]
        nspec = len(input_spectra)

        dq = np.zeros(nelem, dtype=dq_dtype)

        flux = np.zeros((nspec, nelem), dtype=np.float64)
        flux_error = np.zeros((nspec, nelem), dtype=np.float64)
        surf_bright = np.zeros((nspec, nelem), dtype=np.float64)
        sb_error = np.zeros((nspec, nelem), dtype=np.float64)
        weight = np.zeros((nspec, nelem), dtype=np.float64)
        count = np.zeros((nspec, nelem), dtype=np.float64)

        self.flux_unit = input_spectra[0].flux_unit
        self.sb_unit = input_spectra[0].sb_unit

        n_nan = 0  # initial value
        ninputs = 0
        for s, in_spec in enumerate(input_spectra):
            ninputs += 1
            if in_spec.name is not None:
                slit_name = f"{ninputs}, slit {in_spec.name}"
            else:
                slit_name = ninputs
            log.info(f"Accumulating data from input spectrum {slit_name}")
            # Get the pixel numbers in the output corresponding to the
            # wavelengths of the current input spectrum.
            out_pixel = self.wcs.invert(
                in_spec.right_ascension, in_spec.declination, in_spec.wavelength
            )
            # i is a pixel number in the current input spectrum, and
            # k is the corresponding pixel number in the output spectrum.
            nan_flag = np.isnan(out_pixel)
            n_nan += nan_flag.sum()
            for i in range(len(out_pixel)):
                # Need to check on dq and nan flux because dq is not set for some x1d
                if (in_spec.dq[i] & datamodels.dqflags.pixel["DO_NOT_USE"] > 0) | np.isnan(
                    in_spec.flux[i]
                ):
                    continue
                # Round to the nearest pixel.
                if nan_flag[i]:  # skip if the pixel number is NaN
                    continue
                k = round(float(out_pixel[i]))
                if k < 0 or k >= nelem:
                    continue
                dq[k] |= in_spec.dq[i]
                flux[s, k] = in_spec.flux[i]
                flux_error[s, k] = in_spec.flux_error[i]
                surf_bright[s, k] = in_spec.surf_bright[i]
                sb_error[s, k] = in_spec.sb_error[i]
                weight[s, k] = in_spec.weight[i]
                count[s, k] = 1.0

        (flux, flux_error, surf_bright, sb_error, weight, count) = self.combine_spectra(
            flux, flux_error, surf_bright, sb_error, weight, count, sigma_clip=sigma_clip
        )

        if n_nan > 0:
            log.warning(f"{int(n_nan)} output pixel numbers were NaN")

        self.flux = flux
        self.flux_error = flux_error
        self.surf_bright = surf_bright
        self.sb_error = sb_error
        self.dq = dq
        self.weight = weight
        self.count = count

        # Since the output wavelengths will not usually be exactly the same
        # as the input wavelengths, it's possible that there will be output
        # pixels for which there is no corresponding pixel in any of the
        # input spectra.  Check for this case.
        index = np.where(self.count > 0.0)
        n_good = len(index[0])
        if nelem > n_good:
            log.warning(
                f"{int(nelem - n_good)} elements of output had no corresponding input data;"
            )
            log.warning("    these elements will be omitted.")
            self.wavelength = self.wavelength[index]
            self.flux = self.flux[index]
            self.flux_error = self.flux_error[index]
            self.surf_bright = self.surf_bright[index]
            self.sb_error = self.sb_error[index]
            self.dq = self.dq[index]
            self.weight = self.weight[index]
            self.count = self.count[index]
        del index

    def combine_spectra(
        self, flux, flux_error, surf_bright, sb_error, weight, count, sigma_clip=None
    ):
        """
        Combine accumulated spectra.

        Parameters
        ----------
        flux : ndarray, 2-D
            Tabulated fluxes for the input spectra in the format
            [N spectra, M wavelengths].
        flux_error : ndarray, 2-D
            Flux errors of input spectra.
        surf_bright : ndarray, 2-D
            Surface brightnesses of input spectra.
        sb_error : ndarray, 2-D
            Surface brightness errors for input spectra.
        weight : ndarray, 2-D
            Pixel weights for input spectra
        count : ndarray, 2-D
            Count of how many values at each index in the input arrays.
        sigma_clip : float, optional
            Factor for clipping outliers.  Compares input spectra to the
            median and medaian absolute devaition, by default None.

        Returns
        -------
        flux : ndarray, 1-D
            Combined 1-D fluxes.
        flux_error : ndarray, 1-D
            Combined 1-D flux errors.
        surf_bright : ndarray, 1-D
            Combined 1-D surface brightnesses.
        sb_error : ndarray, 1-D
            Combined 1-D surface brightness errors.
        weight : ndarray, 1-D
            Total, per wavelength weights.
        count : ndarray, 1-D
            Total count of spectra contributing to each wavelength.
        """
        # Catch warnings for all NaN slices in an array.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if sigma_clip is not None:
                # Copy the fluxes for modifying
                flux_2d = np.array(flux * weight)

                # Mask any missing pixels in the input spectra
                missing = (count < 1) | (flux_2d == 0.0)
                flux_2d[missing] = np.nan
                # Calculate median and median absolute deviation
                med_flux = np.nanmedian(flux_2d, axis=0)
                mad = np.nanmedian(np.abs(flux_2d - med_flux))

                # Clip any outlier pixels in the input spectra
                clipped = np.abs(flux * weight - med_flux) > sigma_clip * mad
                flux[clipped] = np.nan
                flux_error[clipped] = np.nan
                surf_bright[clipped] = np.nan
                sb_error[clipped] = np.nan
                count[clipped] = 0
                weight[clipped] = 0

            # Perform a weighted sum of the input spectra
            sum_weight = np.nansum(weight, axis=0)
            sum_weight_nonzero = np.where(sum_weight > 0.0, sum_weight, 1.0)

            flux = np.nansum(flux * weight, axis=0) / sum_weight_nonzero
            flux_error = np.sqrt(np.nansum((flux_error * weight) ** 2, axis=0)) / sum_weight_nonzero
            surf_bright = np.nansum(surf_bright * weight, axis=0) / sum_weight_nonzero
            sb_error = np.sqrt(np.nansum((sb_error * weight) ** 2, axis=0)) / sum_weight_nonzero
            count = np.nansum(count, axis=0)

        self.normalized = True

        return flux, flux_error, surf_bright, sb_error, sum_weight, count

    def create_output_data(self):
        """
        Create the output data.

        Returns
        -------
        output_model : `~jwst.datamodels.JwstDataModel`, CombinedSpecModel object
            A table of combined spectral data.
        """
        if not self.normalized:
            log.warning("Data have not been divided by the sum of the weights.")

        cmb_dtype = datamodels.CombinedSpecModel().spec_table.dtype

        # Note that these arrays have to be in the right order.
        data = np.array(
            list(
                zip(
                    self.wavelength,
                    self.flux,
                    self.flux_error,
                    self.surf_bright,
                    self.sb_error,
                    self.dq,
                    self.weight,
                    self.count,
                    strict=False,
                )
            ),
            dtype=cmb_dtype,
        )
        output_model = datamodels.CombinedSpecModel(spec_table=data)

        output_model.spec_table.columns["wavelength"].unit = "um"
        output_model.spec_table.columns["flux"].unit = self.flux_unit
        output_model.spec_table.columns["error"].unit = self.flux_unit
        output_model.spec_table.columns["surf_bright"].unit = self.sb_unit
        output_model.spec_table.columns["sb_error"].unit = self.sb_unit

        return output_model

    def close(self):
        """Set data attributes to null values."""
        self.wavelength = None
        self.flux = None
        self.flux_error = None
        self.surf_bright = None
        self.sb_error = None
        self.dq = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.normalized = False
        self.source_id = None


def count_input(input_spectra):
    """
    Determine the number of input spectra that cover each wavelength.

    For any given input spectrum, the array of wavelengths gives the
    wavelengths at the centers of the pixels.  In this context, the
    expression that an input spectrum "covers" some particular wavelength
    (typically not one of the elements in the input spectrum's array of
    wavelengths) means that this wavelength is within the interval between
    the left edge of the first pixel and the right edge of the last pixel
    of the input spectrum.

    Parameters
    ----------
    input_spectra : list of InputSpectrumModel objects
        List of input spectra.

    Returns
    -------
    wl : ndarray
        Sorted list of all the wavelengths in all the input spectra.

    n_input_spectra : ndarray, 1-D
        For each element of `wl`, the corresponding element of
        `n_input_spectra` is the number of input spectra that cover the
        wavelength in `wl`.
    """
    # Create an array with all the input wavelengths (i.e. the union
    # of the input wavelengths).
    wl = None
    for in_spec in input_spectra:
        input_wl = in_spec.wavelength
        # only include spectra that have more than 1 data point
        if len(input_wl) > 1:
            if wl is None:
                wl = input_wl.copy()
            else:
                wl = np.hstack((input_wl, wl))
    wl.sort()
    nwl = len(wl)

    # n_input_spectra will be the number of input spectra that cover the
    # corresponding wavelength in wl.
    n_input_spectra = np.zeros(nwl, dtype=np.int64)
    for in_spec in input_spectra:
        input_wl = in_spec.wavelength

        # Check for degenerate spectrum. Skip with log.
        if len(input_wl) < 2:
            log.warning(f"Spectrum {in_spec} is degenerate with length {len(input_wl)}")
            log.warning("Skipping...")
            continue

        # wl0 and wl1 will be about a half pixel wider on either side
        # of the wavelength range for the current input spectrum.
        if input_wl[1] > input_wl[0]:  # wavelengths are increasing
            wl0 = input_wl[0] - 0.5 * (input_wl[1] - input_wl[0])
            wl1 = input_wl[-1] + 0.5 * (input_wl[-1] - input_wl[-2])
        elif input_wl[1] < input_wl[0]:  # wavelengths are decreasing
            wl0 = input_wl[-1] - 0.5 * (input_wl[-2] - input_wl[-1])
            wl1 = input_wl[0] + 0.5 * (input_wl[0] - input_wl[1])
        else:
            log.warning(f"Spectrum {in_spec} has a monotonic wavelength solution.")
            log.warning("Skipping...")
            continue
        temp = np.where(wl >= wl0, 1, 0)
        temp = np.where(wl >= wl1, 0, temp)
        n_input_spectra += temp
        del temp

    # This shouldn't happen.
    if np.any(n_input_spectra <= 0.0):
        raise RuntimeError("Problem with input wavelengths.")

    return wl, n_input_spectra


def compute_output_wl(wl, n_input_spectra):
    """
    Compute output wavelengths.

    In summary, the output wavelengths are computed by binning the
    input wavelengths in groups of the number of overlapping spectra.
    The merged and sorted input wavelengths are in `wl`.

    If the wavelength arrays of the input spectra are nearly all the
    same, the values in wl will be in clumps of N nearly identical
    values (for N input files).  We want to average the wavelengths
    in those clumps to get the output wavelengths.  However, if the
    wavelength scales of the input files differ slightly, the values
    in wl may clump in some regions but not in others.  The algorithm
    below tries to find clumps; a clump is identified by having a
    small standard deviation of all the wavelengths in a slice of wl,
    compared with neighboring slices.  It is intended that this should
    not find clumps where they're not present, i.e. in regions where
    the wavelengths of the input spectra are not well aligned with each
    other.  These regions will be gaps in the initial determination
    of output wavelengths based on clumps.  In such regions, output
    wavelengths will then be computed by just averaging groups of
    input wavelengths to fill in the gaps.

    Parameters
    ----------
    wl : 1-D array
        An array containing all the wavelengths from all the input
        spectra, sorted in increasing order.

    n_input_spectra : 1-D array
        An integer array of the same length as `wl`.  For a given
        array index k (for example), n_input_spectra[k] is the number of
        input spectra that cover wavelength wl[k].

    Returns
    -------
    wavelength :  1-D array
        Array of wavelengths for the output spectrum.
    """
    nwl = len(wl)

    # sigma is an array of the standard deviation at each element
    # of wl, over n_input_spectra elements.  A small value implies that
    # there's a clump, i.e. several elements of wl with nearly the
    # same wavelength.
    sigma = np.zeros(nwl, dtype=np.float64) + 9999.0

    # mean_wl is the mean wavelength over the same slice of wl that we
    # used to compute sigma.  If sigma is small enough that it looks
    # as if there's a clump, we'll copy the mean_wl value to temp_wl
    # to be one element of the output wavelengths.
    mean_wl = np.zeros(nwl, dtype=np.float64) - 99.0

    # temp_wl has the same number of elements as wl, but we expect the
    # array of output wavelengths to be significantly smaller, so
    # temp_wl is initialized to a negative value as a flag.  Positive
    # elements will be copied to the array of output wavelengths.
    temp_wl = np.zeros(nwl, dtype=np.float64) - 99.0

    for k in range(nwl):
        n = n_input_spectra[k]
        if n == 1:
            sigma[k] = 0.0
            mean_wl[k] = wl[k]
            temp_wl[k] = mean_wl[k]
        else:
            k1 = k + n // 2 + 1
            k0 = k1 - n
            if k0 >= 0 and k1 <= nwl:
                sigma[k] = wl[k0:k1].std()
                mean_wl[k] = wl[k0:k1].mean()
                if sigma[k] == 0.0:
                    temp_wl[k] = mean_wl[k]

    cutoff = 0.8
    for k in range(nwl):
        # If sigma[k] equals 0, temp_wl has already been assigned.
        if sigma[k] > 0.0:
            if k == 0:
                if sigma[k] < cutoff * sigma[1]:
                    temp_wl[k] = mean_wl[k]
            elif k == nwl - 1:
                if sigma[k] < cutoff * sigma[nwl - 2]:
                    temp_wl[k] = mean_wl[k]
            else:
                if sigma[k] < cutoff * (sigma[k - 1] + sigma[k + 1]) / 2.0:
                    temp_wl[k] = mean_wl[k]

    # Fill gaps in the output wavelengths by taking averages of the
    # input wavelengths.  If there are n overlapping input spectra,
    # average a block of n elements of `wl`.
    done = False
    i = 0
    while not done:
        if i >= nwl:
            done = True
        else:
            n = n_input_spectra[i]
            middle = n // 2
            if i + middle < nwl:
                n = n_input_spectra[i + middle]
                if i + n < nwl:
                    # The nominal range is i - 1 to i + n (inclusive) for
                    # checking whether an output wavelength has already
                    # been assigned.
                    low = max(i - 1, 0)
                    high = min(i + n + 1, nwl - 1)
                    if temp_wl[low:high].max() <= 0.0:
                        temp_wl[i] = wl[i : i + n].mean()
            i += n

    return temp_wl[np.where(temp_wl > 0.0)].copy()


def check_exptime(exptime_key):
    """
    Check exptime_key for validity.

    This function checks exptime_key.  If it is valid, the corresponding
    value used by the metadata interface will be returned.  This will be
    either "integration_time" or "exposure_time".  If it is not valid,
    "unit weight" will be returned (meaning that a weight of 1 will be
    used when averaging spectra), and a warning will be logged.

    Parameters
    ----------
    exptime_key : str
        A keyword or string indicating what value (integration time or
        exposure time) should be used as a weight when combing spectra.

    Returns
    -------
    exptime_key : str
        The value will be either "integration_time", "exposure_time",
        or "unit_weight".
    """
    exptime_lwr = exptime_key.lower()
    if exptime_lwr.startswith("integration") or exptime_lwr == "effinttm":
        exptime_key = "integration_time"
        log.info("Using integration time as the weight.")
    elif exptime_lwr.startswith("exposure") or exptime_lwr == "effexptm":
        exptime_key = "exposure_time"
        log.info("Using exposure time as the weight.")
    elif exptime_lwr == "unit_weight" or exptime_lwr == "unit weight":
        exptime_key = "unit_weight"
        log.info("Using weight = 1.")
    else:
        log.warning(f"Don't understand exptime_key = '{exptime_key}'; using unit weight.")
        log.info("The options for exptime_key are:")
        log.info("  integration_time, effinttm, exposure_time, effexptm, unit_weight, unit weight")
        exptime_key = "unit_weight"

    return exptime_key


def _read_input_spectra(input_model, exptime_key, input_spectra):
    """
    Read input spectra from a datamodel.

    Parameters
    ----------
    input_model : MultiSpecModel, TSOMultiSpecModel, MRSSpecModel
        A datamodel with a ``spec`` attribute, containing spectra.
        If TSOMultiSpecModel, integrations in the spectral table rows
        are expanded into separate spectra.
    exptime_key : str
        Exposure time key to use for weighting.
    input_spectra : dict
        Dictionary to hold input spectra, keyed by spectral order.
        Updated in place.

    Returns
    -------
    input_spectra : dict
        The updated dictionary, holding all spectra in the input model.

    Raises
    ------
    TypeError
        If the input datamodel does not have a ``spec`` attribute.
    """
    if not hasattr(input_model, "spec"):
        raise TypeError(f"Invalid input datamodel: {type(input_model)}")
    if isinstance(input_model, datamodels.TSOMultiSpecModel):
        spectra = expand_flat_spec(input_model).spec
    else:
        spectra = input_model.spec
    for in_spec in spectra:
        if not np.any(np.isfinite(in_spec.spec_table.field("flux"))):
            log.warning(
                f"Input spectrum {in_spec.source_id} order {in_spec.spectral_order} "
                f"from group_id {in_spec.meta.group_id} has no valid flux values; skipping."
            )
            continue
        spectral_order = in_spec.spectral_order
        if spectral_order not in input_spectra:
            input_spectra[spectral_order] = []
        input_spectra[spectral_order].append(InputSpectrumModel(input_model, in_spec, exptime_key))
    return input_spectra


def combine_1d_spectra(input_model, exptime_key, sigma_clip=None):
    """
    Combine the input spectra.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        The input spectra.  This will likely be a ModelContainer object,
        but may also be a multi-spectrum model, such as MultiSpecModel or
        TSOMultiSpecModel.  Input spectra may have different spectral orders
        or wavelengths but should all share the same target.
    exptime_key : str
        A string identifying which keyword to use to get the exposure time,
        which is used as a weight when combining spectra.  The value should
        be one of:  "exposure_time" (the default), "integration_time",
        or "unit_weight".

    Returns
    -------
    output_model : `~jwst.datamodels.MultiCombinedSpecModel`
        A combined spectra datamodel.
    """
    log.debug(f"Using exptime_key = {exptime_key}.")

    exptime_key = check_exptime(exptime_key)

    input_spectra = {}
    output_spectra = {}
    if isinstance(input_model, ModelContainer):
        for ms in input_model:
            _read_input_spectra(ms, exptime_key, input_spectra)
    else:
        _read_input_spectra(input_model, exptime_key, input_spectra)

    if len(input_spectra) == 0:
        log.error("No valid input spectra found for source. Skipping.")
        result = input_model.copy()
        result.meta.cal_step.combine_1d = "SKIPPED"
        return result

    for order in input_spectra:
        output_spectra[order] = OutputSpectrumModel()
        output_spectra[order].assign_wavelengths(input_spectra[order])
        output_spectra[order].accumulate_sums(input_spectra[order], sigma_clip=sigma_clip)

    output_model = datamodels.MultiCombinedSpecModel()

    for order in output_spectra:
        output_order = output_spectra[order].create_output_data()
        output_order.spectral_order = order
        for attr in SPECMETA_ATTRIBUTES:
            setattr(output_order, attr, getattr(input_spectra[order][0], attr))
        output_model.spec.append(output_order)

    # Copy one of the input headers to output.
    if isinstance(input_model, ModelContainer):
        output_model.update(input_model[0], only="PRIMARY")
    else:
        output_model.update(input_model, only="PRIMARY")

    # Looks clunky, but need an output_spec instance to copy wcs
    output_model.meta.wcs = output_spectra[list(output_spectra)[0]].wcs
    output_model.meta.cal_step.combine_1d = "COMPLETE"

    for order in input_spectra:
        for in_spec in input_spectra[order]:
            in_spec.close()

    for order in output_spectra:
        output_spectra[order].close()

    return output_model

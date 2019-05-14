#
# Module for combining 1-D spectra
#

import logging

import numpy as np

from .. import datamodels
from ..extract_1d.spec_wcs import create_spectral_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class InputSpectrumModel:
    """Attributes:
        wavelength
        flux
        error
        surf_bright
        sb_error
        dq
        nelem
        weight
        unit_weight
        right_ascension
        declination
    """

    def __init__(self, ms, spec, exptime_key):
        """Create an InputSpectrumModel object.

        Parameters
        ----------
        ms : `~jwst.datamodels.DataModel`, MultiSpecModel or SpecModel
            This is used to get the integration time.

        spec : `~jwst.datamodels.DataModel`, SpecModel table
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
        self.error = spec.spec_table.field("error").copy()
        try:
            self.surf_bright = spec.spec_table.field("surf_bright").copy()
            self.sb_error = spec.spec_table.field("sb_error").copy()
        except KeyError:
            self.surf_bright = np.zeros_like(self.flux)
            self.sb_error = np.zeros_like(self.flux)
            log.warning("There is no SURF_BRIGHT column in the input.")
        self.dq = spec.spec_table.field("dq").copy()
        self.nelem = self.wavelength.shape[0]
        self.unit_weight = False        # may be reset below
        self.right_ascension = np.zeros_like(self.wavelength)
        self.declination = np.zeros_like(self.wavelength)

        self.weight = np.ones_like(self.wavelength)
        if exptime_key == "integration_time":
            self.weight *= ms.meta.exposure.integration_time
        elif exptime_key == "exposure_time":
            self.weight *= ms.meta.exposure.exposure_time
        elif exptime_key == "unit_weight":
            self.unit_weight = True
        else:
            raise RuntimeError("Don't understand exptime_key = '%s'" %
                               exptime_key)

        got_wcs = True
        try:
            wcs = spec.meta.wcs
        except AttributeError:
            got_wcs = False
        if got_wcs:
            self.right_ascension[:], self.declination[:], _ = wcs(0.)
        else:
            self.right_ascension[:] = ms.meta.target.ra
            self.declination[:] = ms.meta.target.dec
            log.warning("There is no WCS in the input.")

    def close(self):
        self.wavelength = None
        self.flux = None
        self.error = None
        self.surf_bright = None
        self.sb_error = None
        self.dq = None
        self.nelem = 0
        self.weight = 1.
        self.unit_weight = False
        self.right_ascension = None
        self.declination = None


class OutputSpectrumModel:
    """Attributes:
        wavelength
        flux
        error
        surf_bright
        sb_error
        dq
        weight
        count
        wcs
        wavelength_dtype
        flux_dtype
        dq_dtype
        normalized
    """

    def __init__(self):

        self.wavelength = None
        self.flux = None
        self.error = None
        self.surf_bright = None
        self.sb_error = None
        self.dq = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.wavelength_dtype = None
        self.flux_dtype = None
        self.dq_dtype = None
        self.normalized = False

    def assign_wavelengths(self, input_spectra):
        """Create an array of wavelengths to use for the output spectrum.

        Take the union of all input wavelengths, then call method
        compute_output_wl to bin wavelengths in groups of the number of
        overlapping spectra.

        Parameters
        ----------
        input_spectra : list of InputSpectrumModel objects
            List of input spectra.
        """

        # Save these, so we'll know what data type to use for the output.
        # The types used for accumulating sums and taking averages may not
        # be the same as these types.
        self.wavelength_dtype = input_spectra[0].wavelength.dtype
        self.flux_dtype = input_spectra[0].flux.dtype
        self.dq_dtype = input_spectra[0].dq.dtype

        nwl = 0
        for in_spec in input_spectra:
            nwl += in_spec.nelem

        # Create an array with all the input wavelengths (i.e. the union
        # of the input wavelengths).
        wl = np.zeros(nwl, dtype=np.float64)
        i = 0
        for in_spec in input_spectra:
            nelem = in_spec.nelem
            # Concatenate current input wavelengths to wl array.
            wl[i:i + nelem] = in_spec.wavelength.copy()
            i += nelem
        wl.sort()

        # count_input will be the number of input spectra that cover the
        # corresponding wavelength in wl.
        count_input = np.zeros(nwl, dtype=np.int64)
        for in_spec in input_spectra:
            input_wl = in_spec.wavelength
            # wl0 and wl1 will be about a half pixel wider on either side
            # of the wavelength range for the current input spectrum.
            if input_wl[1] > input_wl[0]:       # wavelengths are increasing
                wl0 = input_wl[0] - 0.5 * (input_wl[1] - input_wl[0])
                wl1 = input_wl[-1] + 0.5 * (input_wl[-1] - input_wl[-2])
            elif input_wl[1] < input_wl[0]:     # wavelengths are decreasing
                wl0 = input_wl[-1] - 0.5 * (input_wl[-2] - input_wl[-1])
                wl1 = input_wl[0] + 0.5 * (input_wl[0] - input_wl[1])
            else:
                raise RuntimeError("Wavelength increment must not be zero.")
            temp = np.where(wl >= wl0, 1, 0)
            temp = np.where(wl >= wl1, 0, temp)
            count_input += temp
            del temp
        # This shouldn't happen.
        if np.any(count_input <= 0.):
            raise RuntimeError("Problem with input wavelengths.")

        self.wavelength = self.compute_output_wl(wl, count_input)

        self.wcs = create_spectral_wcs(input_spectra[0].right_ascension[0],
                                       input_spectra[0].declination[0],
                                       self.wavelength)

    def compute_output_wl(self, wl, count_input):
        """Compute output wavelengths.

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

        count_input : 1-D array
            An integer array of the same length as `wl`.  For a given
            array index k (for example), count_input[k] is the number of
            input spectra that cover wavelength wl[k].

        Returns
        -------
        wavelength :  1-D array
            Array of wavelengths for the output spectrum.
        """

        nwl = len(wl)

        # sigma is an array of the standard deviation at each element
        # of wl, over count_input elements.  A small value implies that
        # there's a clump, i.e. several elements of wl with nearly the
        # same wavelength.
        sigma = np.zeros(nwl, dtype=np.float64) + 9999.

        # mean_wl is the mean wavelength over the same slice of wl that we
        # used to compute sigma.  If sigma is small enough that it looks
        # as if there's a clump, we'll copy the mean_wl value to temp_wl
        # to be one element of the output wavelengths.
        mean_wl = np.zeros(nwl, dtype=np.float64) - 99.

        # temp_wl has the same number of elements as wl, but we expect the
        # array of output wavelengths to be significantly smaller, so
        # temp_wl is initialized to a negative value as a flag.  Positive
        # elements will be copied to the array of output wavelengths.
        temp_wl = np.zeros(nwl, dtype=np.float64) - 99.

        for k in range(nwl):
            n = count_input[k]
            if n == 1:
                sigma[k] = 0.
                mean_wl[k] = wl[k]
                temp_wl[k] = mean_wl[k]
            else:
                k1 = k + n // 2 + 1
                k0 = k1 - n
                if k0 >= 0 and k1 <= nwl:
                    sigma[k] = wl[k0:k1].std()
                    mean_wl[k] = wl[k0:k1].mean()
                    if sigma[k] == 0.:
                        temp_wl[k] = mean_wl[k]

        cutoff = 0.8
        for k in range(nwl):
            # If sigma[k] equals 0, temp_wl has already been assigned.
            if sigma[k] > 0.:
                if k == 0:
                    if sigma[k] < cutoff * sigma[1]:
                        temp_wl[k] = mean_wl[k]
                elif k == nwl - 1:
                    if sigma[k] < cutoff * sigma[nwl - 2]:
                        temp_wl[k] = mean_wl[k]
                else:
                    if sigma[k] < cutoff * (sigma[k - 1] + sigma[k + 1]) / 2.:
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
                n = count_input[i]
                middle = n // 2
                if i + middle < nwl:
                    n = count_input[i + middle]
                    if i + n < nwl:
                        # The nominal range is i - 1 to i + n (inclusive) for
                        # checking whether an output wavelength has already
                        # been assigned.
                        low = max(i - 1, 0)
                        high = min(i + n + 1, nwl - 1)
                        if temp_wl[low:high].max() <= 0.:
                            temp_wl[i] = wl[i:i + n].mean()
                i += n

        return temp_wl[np.where(temp_wl > 0.)].copy()

    def accumulate_sums(self, input_spectra):
        """Compute a weighted sum of all the input spectra.

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
        """

        nelem = self.wavelength.shape[0]

        self.flux = np.zeros(nelem, dtype=np.float64)
        self.error = np.zeros(nelem, dtype=np.float64)
        self.surf_bright = np.zeros(nelem, dtype=np.float64)
        self.sb_error = np.zeros(nelem, dtype=np.float64)
        self.dq = np.zeros(nelem, dtype=self.dq_dtype)
        self.weight = np.zeros(nelem, dtype=np.float64)
        self.count = np.zeros(nelem, dtype=np.float64)

        for in_spec in input_spectra:
            # Get the pixel numbers in the output corresponding to the
            # wavelengths of the current input spectrum.
            out_pixel = self.wcs.invert(in_spec.right_ascension,
                                        in_spec.declination,
                                        in_spec.wavelength)
            # i is a pixel number in the current input spectrum, and
            # k is the corresponding pixel number in the output spectrum.
            for i in range(len(out_pixel)):
                if in_spec.dq[i] & datamodels.dqflags.pixel['DO_NOT_USE'] > 0:
                    continue
                # Round to the nearest pixel.
                k = round(float(out_pixel[i]))
                weight = in_spec.weight[i]
                self.dq[k] |= in_spec.dq[i]
                self.flux[k] += in_spec.flux[i] * weight
                self.error[k] += (in_spec.error[i] * weight)**2
                self.surf_bright[k] += (in_spec.surf_bright[i] * weight)
                self.sb_error[k] += (in_spec.sb_error[i] * weight)**2
                self.weight[k] += weight
                self.count[k] += 1.

        # Since the output wavelengths will not usually be exactly the same
        # as the input wavelengths, it's possible that there will be output
        # pixels for which there is no corresponding pixel in any of the
        # input spectra.  Check for this case.
        index = np.where(self.count > 0.)
        n_good = len(index[0])
        if nelem > n_good:
            log.warning("%d elements of output had no corresponding"
                        " input data;" % (nelem - n_good,))
            log.warning("    these elements will be omitted.")
            self.wavelength = self.wavelength[index]
            self.flux = self.flux[index]
            self.error = self.error[index]
            self.surf_bright = self.surf_bright[index]
            self.sb_error = self.error[index]
            self.dq = self.dq[index]
            self.weight = self.weight[index]
            self.count = self.count[index]
        del index

        self.normalized = False

    def compute_combination(self):
        """Compute the combined values."""

        if not self.normalized:
            sum_weight = np.where(self.weight > 0., self.weight, 1.)
            self.surf_bright /= sum_weight
            self.flux /= sum_weight
            self.error = np.sqrt(self.error / sum_weight)
            self.sb_error = np.sqrt(self.sb_error / sum_weight)
            self.normalized = True

    def create_output(self):
        """Create the output model.

        Returns
        -------
        output_model : `~jwst.datamodels.DataModel`, CombinedSpecModel object
            A table of combined spectral data.
        """

        if not self.normalized:
            log.warning("Data have not been divided by"
                        " the sum of the weights.")

        dtype = [('WAVELENGTH', self.wavelength_dtype),
                 ('FLUX', self.flux_dtype),
                 ('ERROR', self.flux_dtype),
                 ('SURF_BRIGHT', self.flux_dtype),
                 ('SB_ERROR', self.flux_dtype),
                 ('DQ', self.dq_dtype),
                 ('WEIGHT', self.wavelength_dtype),
                 ('N_INPUT', np.float64)]

        data = np.array(list(zip(self.wavelength,
                                 self.flux,
                                 self.error,
                                 self.surf_bright,
                                 self.sb_error,
                                 self.dq,
                                 self.weight,
                                 self.count)), dtype=dtype)
        output_model = datamodels.CombinedSpecModel(spec_table=data)

        return output_model

    def close(self):
        self.wavelength = None
        self.flux = None
        self.error = None
        self.surf_bright = None
        self.sb_error = None
        self.dq = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.wavelength_dtype = None
        self.flux_dtype = None
        self.dq_dtype = None
        self.normalized = False


def check_exptime(exptime_key):
    """Check exptime_key for validity.

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
    if exptime_lwr.startswith("integration") or \
       exptime_lwr == "effinttm":
        exptime_key = "integration_time"
        log.info("Using integration time as the weight.")
    elif exptime_lwr.startswith("exposure") or \
         exptime_lwr == "effexptm":
        exptime_key = "exposure_time"
        log.info("Using exposure time as the weight.")
    elif exptime_lwr == "unit_weight" or exptime_lwr == "unit weight":
        exptime_key = "unit_weight"
        log.info("Using weight = 1.")
    else:
        log.warning("Don't understand exptime_key = '%s';"
                    " using unit weight." % exptime_key)
        log.info("The options for exptime_key are:")
        log.info("  integration_time, effinttm, exposure_time, effexptm,"
                 " unit_weight, unit weight")
        exptime_key = "unit_weight"

    return exptime_key


def combine_1d_spectra(input_model, exptime_key):
    """Combine the input spectra.

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        The input spectra.  This will likely be a ModelContainer object.

    exptime_key : str
        A string identifying which keyword to use to get the exposure time,
        which is used as a weight when combining spectra.  The value should
        be one of:  "exposure_time" (the default), "integration_time",
        or "unit_weight".

    Returns
    -------
    output_model : `~jwst.datamodels.DataModel`
        A datamodels.CombinedSpecModel object.
    """

    log.debug("Using exptime_key = {}.".format(exptime_key))

    exptime_key = check_exptime(exptime_key)

    input_spectra = []
    if isinstance(input_model, datamodels.ModelContainer):
        for ms in input_model:
            for in_spec in ms.spec:
                input_spectra.append(InputSpectrumModel(
                                ms, in_spec, exptime_key))
    else:
        for in_spec in input_model.spec:
            input_spectra.append(InputSpectrumModel(
                                input_model, in_spec, exptime_key))

    output_spec = OutputSpectrumModel()
    output_spec.assign_wavelengths(input_spectra)
    output_spec.accumulate_sums(input_spectra)
    output_spec.compute_combination()

    for in_spec in input_spectra:
        in_spec.close()

    output_model = output_spec.create_output()

    # Copy one of the input headers to output.
    if isinstance(input_model, datamodels.ModelContainer):
        output_model.update(input_model[0], only="PRIMARY")
    else:
        output_model.update(input_model, only="PRIMARY")

    output_model.meta.wcs = output_spec.wcs
    output_model.meta.cal_step.combine_1d = 'COMPLETE'

    output_spec.close()

    return output_model

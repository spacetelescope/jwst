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
        net
        dq
        nelem
        weight
        unit_weight
        right_ascension
        declination
    """

    def __init__(self, ms, spec, exptime_key, background):
        """Create an InputSpectrumModel object.

        Parameters
        ----------
        ms : MultiSpecModel or SpecModel object
            This is used to get the integration time.

        spec : SpecModel table
            The table containing columns "wavelength" and "net".
            The `ms` object may contain more than one spectrum, but `spec`
            should be just one of those.

        exptime_key : str
            A string identifying which keyword to use to get the exposure
            time, which is used as a weight; or "unit_weight", which means
            to use weight = 1.

        background : bool
            If the flux data are actually background rather than a target
            spectrum, `background` should be set to True.  In this case, the
            values read from the flux column of each input spectrum will be
            divided by the npixels column (if that column exists).  This is
            to convert the values to background per pixel.
        """

        self.wavelength = spec.spec_table.field("wavelength").copy()

        if background:
            try:
                npixels = spec.spec_table.field("npixels").copy()
            except KeyError:
                npixels = np.ones_like(self.wavelength)

        if background:
            # Convert to background value per pixel.
            self.flux = spec.spec_table.field("flux") / npixels
            self.error = spec.spec_table.field("error") / npixels
            self.net = spec.spec_table.field("net") / npixels
        else:
            self.flux = spec.spec_table.field("flux").copy()
            self.error = spec.spec_table.field("error").copy()
            self.net = spec.spec_table.field("net").copy()
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
        if background and exptime_key != "unit_weight":
            self.weight *= npixels

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
        self.net = None
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
        net
        dq
        flux_weight
        weight
        count
        wcs
        wavelength_dtype
        net_dtype
        dq_dtype
        normalized
    """

    def __init__(self):

        self.wavelength = None
        self.flux = None
        self.error = None
        self.net = None
        self.dq = None
        self.flux_weight = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.wavelength_dtype = None
        self.net_dtype = None
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
        self.net_dtype = input_spectra[0].net.dtype
        self.dq_dtype = input_spectra[0].dq.dtype

        nwl = 0
        for in_spec in input_spectra:
            nwl += in_spec.nelem

        # Create an array with all the input wavelengths (i.e. the union
        # of the input wavelengths).
        wl = np.zeros(nwl, dtype=np.float)
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
        sigma = np.zeros(nwl, dtype=np.float) + 9999.

        # mean_wl is the mean wavelength over the same slice of wl that we
        # used to compute sigma.  If sigma is small enough that it looks
        # as if there's a clump, we'll copy the mean_wl value to temp_wl
        # to be one element of the output wavelengths.
        mean_wl = np.zeros(nwl, dtype=np.float) - 99.

        # temp_wl has the same number of elements as wl, but we expect the
        # array of output wavelengths to be significantly smaller, so
        # temp_wl is initialized to a negative value as a flag.  Positive
        # elements will be copied to the array of output wavelengths.
        temp_wl = np.zeros(nwl, dtype=np.float) - 99.

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

        The net will be weighted by the exposure (or integration) time.
        The flux will also be weighted by the sensitivity.  The ratio of
        net to flux is the sensitivity.  For both net and flux, the idea
        is that it is counts that should be added up when accumulating
        sums.  If unit weight was specified, however, both net and flux
        will be weighted by one.

        Parameters
        ----------
        input_spectra : list of InputSpectrumModel objects
            List of input spectra.
        """

        nelem = self.wavelength.shape[0]

        self.flux = np.zeros(nelem, dtype=np.float)
        self.error = np.zeros(nelem, dtype=np.float)
        self.flux_weight = np.zeros(nelem, dtype=np.float)
        self.dq = np.zeros(nelem, dtype=self.dq_dtype)
        self.net = np.zeros(nelem, dtype=np.float)
        self.weight = np.zeros(nelem, dtype=np.float)
        self.count = np.zeros(nelem, dtype=np.float)

        # The flux should be weighted by sensitivity (as well as exposure
        # time), but if the input net columns are not populated, we can't
        # compute the sensitivity.
        weight_flux_by_sensitivity = True
        for in_spec in input_spectra:
            if in_spec.net.min() == 0. and in_spec.net.max() == 0.:
                weight_flux_by_sensitivity = False
                log.warning("The NET column is all zero in one or more "
                            "input tables, so FLUX will not be weighted by "
                            "sensitivity.")
                break

        for in_spec in input_spectra:
            if weight_flux_by_sensitivity:
                # Replace zeros so we can divide by the flux.
                temp_flux = np.where(in_spec.flux == 0., 1., in_spec.flux)
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
                self.net[k] += (in_spec.net[i] * in_spec.weight[i])
                self.weight[k] += in_spec.weight[i]
                self.dq[k] |= in_spec.dq[i]
                if in_spec.unit_weight:
                    flux_wgt = 1.
                elif weight_flux_by_sensitivity:
                    # net / flux is the sensitivity
                    flux_wgt = (in_spec.weight[i] *
                                in_spec.net[i] / temp_flux[i])
                    flux_wgt = max(flux_wgt, 0.)
                else:
                    flux_wgt = in_spec.weight[i]
                self.flux[k] += in_spec.flux[i] * flux_wgt
                self.error[k] += (in_spec.error[i] * flux_wgt)**2
                self.flux_weight[k] += flux_wgt
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
            self.net = self.net[index]
            self.weight = self.weight[index]
            self.flux_weight = self.flux_weight[index]
            self.error = self.error[index]
            self.count = self.count[index]
        del index

        self.normalized = False

    def compute_combination(self):
        """Compute the combined values."""

        if not self.normalized:
            sum_weight = np.where(self.weight > 0., self.weight, 1.)
            sum_flux_wgt = np.where(self.flux_weight > 0.,
                                    self.flux_weight, 1.)
            self.net /= sum_weight
            self.flux /= sum_flux_wgt
            self.error = np.sqrt(self.error / sum_flux_wgt)
            self.normalized = True


    def create_output(self):
        """Create the output model.

        Returns
        -------
        output_model : CombinedSpecModel object
            A table of combined spectral data.
        """

        if not self.normalized:
            log.warning("Data have not been divided by"
                        " the sum of the weights.")

        dtype = [('wavelength', self.wavelength_dtype),
                 ('flux', self.net_dtype),
                 ('error', self.net_dtype),
                 ('net', self.net_dtype),
                 ('dq', self.dq_dtype),
                 ('weight', self.wavelength_dtype),
                 ('n_input', np.float)]

        data = np.array(list(zip(self.wavelength,
                                 self.flux,
                                 self.error,
                                 self.net,
                                 self.dq,
                                 self.weight,
                                 self.count)), dtype=dtype)
        output_model = datamodels.CombinedSpecModel(spec_table=data)

        return output_model


    def close(self):
        self.wavelength = None
        self.flux = None
        self.error = None
        self.net = None
        self.dq = None
        self.flux_weight = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.wavelength_dtype = None
        self.net_dtype = None
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


def combine_1d_spectra(input_model, exptime_key, background=False):
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

    background : bool, default=False
        If the flux data are actually background rather than a target
        spectrum, `background` should be set to True.  In this case, the
        values read from the flux column of each input spectrum will be
        divided by the npixels column (if that column exists).  This is
        to convert the values to background per pixel.

    Returns
    -------
    output_model : `~jwst.datamodels.DataModel`
        A datamodels.CombinedSpecModel object.
    """

    log.debug("Using exptime_key = {}.".format(exptime_key))
    if background:
        log.debug("The FLUX data will be treated as background data.")

    exptime_key = check_exptime(exptime_key)

    input_spectra = []
    if isinstance(input_model, datamodels.ModelContainer):
        for ms in input_model:
            for in_spec in ms.spec:
                input_spectra.append(InputSpectrumModel(
                                ms, in_spec, exptime_key, background))
    else:
        for in_spec in input_model.spec:
            input_spectra.append(InputSpectrumModel(
                                input_model, in_spec, exptime_key, background))

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

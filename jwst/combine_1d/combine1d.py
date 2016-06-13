#
# Module for combining 1-D spectra
#

import logging
import math

import numpy as np
import json
from .. import datamodels

from . import temp_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class InputSpectrumModel(object):
    """Attributes:
        wavelength
        countrate
        nelem
        wcs
        weight
    """

    def __init__(self, ms, spec, exptime_key):
        """Create an InputSpectrumModel object.

        Parameters
        ----------
        ms: MultiSpecModel or SpecModel object
            This is used to get the integration time.

        spec: SpecModel table
            The table containing columns "wavelength" and "countrate".
            The `ms` object may contain more than one spectrum, but `spec`
            should be just one of those.

        exptime_key: str
            A string identifying which keyword to use to get the exposure
            time, which is used as a weight; or "unit_weight", which means
            to use weight = 1.
        """

        self.wavelength = spec.spec_table.field("wavelength").copy()
        self.countrate = spec.spec_table.field("countrate").copy()
        self.nelem = self.wavelength.shape[0]
        self.wcs = temp_wcs.WCS(self.wavelength)

        if exptime_key == "integration_time":
            self.weight = ms.meta.exposure.integration_time
        elif exptime_key == "exposure_time":
            self.weight = ms.meta.exposure.exposure_time
        elif exptime_key == "unit_weight":
            self.weight = 1.
        else:
            raise RuntimeError("Don't understand exptime_key = '%s'" %
                               exptime_key)

    def close(self):
        if self.wavelength is not None:
            del self.wavelength
            self.wavelength = None
        if self.countrate is not None:
            del self.countrate
            self.countrate = None
        if self.wcs:
            self.wcs.close()
            del self.wcs
            self.wcs = None
        self.nelem = 0
        self.weight = 1.

class OutputSpectrumModel(object):
    """Attributes:
        wavelength
        countrate
        weight
        count
        wcs
        wavelength_dtype
        countrate_dtype
        normalized
    """

    def __init__(self):

        self.wavelength = None
        self.countrate = None
        self.weight = None
        self.count = None
        self.wcs = None
        self.wavelength_dtype = None
        self.countrate_dtype = None
        self.normalized = False

    def assign_wavelengths(self, input_spectra):
        """Create an array of wavelengths to use for the output spectrum.

        Take the union of all input wavelengths, then call method
        compute_output_wl to bin wavelengths in groups of the number of
        overlapping spectra.

        Parameters
        ----------
        input_spectra: list of InputSpectrumModel objects
            List of input spectra.
        """

        # Save these, so we'll know what data type to use for the output.
        # The types used for accumulating sums and taking averages may not
        # be the same as these types.
        self.wavelength_dtype = input_spectra[0].wavelength.dtype
        self.countrate_dtype = input_spectra[0].countrate.dtype

        nspectra = len(input_spectra)
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

        self.compute_output_wl(wl, count_input)

        self.wcs = temp_wcs.WCS(self.wavelength)

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

        The output wavelengths will be assigned to the `wavelengths`
        attribute.

        Parameters
        ----------
        wl: 1-D array
            An array containing all the wavelengths from all the input
            spectra, sorted in increasing order.

        count_input: 1-D array
            An integer array of the same length as `wl`.  For a given
            array index k (for example), count_input[k] is the number of
            input spectra that cover wavelength wl[k].
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

        self.wavelength = temp_wl[np.where(temp_wl > 0.)].copy()

    def accumulate_sums(self, input_spectra, interpolation):
        """Compute a weighted sum of all the input spectra.

        Each pixel of each input spectrum will be added to one pixel of
        the output spectrum.  The wavelength spacing of the input and
        output should be comparable, so each input pixel will usually
        correspond to one output pixel (i.e. have the same wavelength,
        to within half a pixel).  However, if the input and output
        wavelength spacings are not identical, then for any given input
        spectrum, there can be output pixels that are incremented by more
        than one input pixel, and in that case there will be at least one
        output pixel that is not incremented.  If there are several input
        spectra, such gaps will hopefully be filled in.

        Parameters
        ----------
        input_spectra: list of InputSpectrumModel objects
            List of input spectra.

        interpolation: str
            The type of interpolation to use.
            This is currently ignored.
        """

        nelem = self.wavelength.shape[0]

        self.countrate = np.zeros(nelem, dtype=np.float64)
        self.weight = np.zeros(nelem, dtype=np.float64)
        self.count = np.zeros(nelem, dtype=np.float32)

        for in_spec in input_spectra:
            # Get the pixel numbers in the output arrays corresponding to
            # the wavelengths of the current input spectrum.
            out_pixel = self.wcs.to_pixel(in_spec.wavelength)
            for i in range(len(out_pixel)):
                # Nearest pixel "interpolation."
                k = int(round(out_pixel[i]))
                self.countrate[k] += (in_spec.countrate[i] * in_spec.weight)
                self.weight[k] += in_spec.weight
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
            self.countrate = self.countrate[index]
            self.weight = self.weight[index]
            self.count = self.count[index]
        del index

        self.normalized = False

    def compute_combination(self):
        """Compute the combined count rate."""

        if not self.normalized:
            weight = np.where(self.weight > 0., self.weight, 1.)
            self.countrate /= weight
            self.normalized = True

    def create_output(self):
        """Create the output model.

        Returns
        -------
        out_model: CombinedSpecModel object
            A table of combined spectral data.
        """

        if not self.normalized:
            log.warning("Data have not been divided by"
                        " the sum of the weights.")

        dtype = [('wavelength', self.wavelength_dtype),
                 ('countrate', self.countrate_dtype),
                 ('weight', self.wavelength_dtype),
                 ('n_input', np.float32)]

        data = np.array(list(zip(self.wavelength,
                            self.countrate,
                            self.weight,
                            self.count)), dtype=dtype)
        out_model = datamodels.CombinedSpecModel(spec_table=data)

        return out_model

    def close(self):
        if self.wavelength is not None:
            del self.wavelength
            self.wavelength = None
        if self.countrate is not None:
            del self.countrate
            self.countrate = None
        if self.weight is not None:
            del self.weight
            self.weight = None
        if self.count is not None:
            del self.count
            self.count = None
        if self.wcs:
            self.wcs.close()
            del self.wcs
            self.wcs = None
        self.wavelength_dtype = None
        self.countrate_dtype = None
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
    exptime_key: str
        A keyword or string indicating what value (integration time or
        exposure time) should be used as a weight when combing spectra.

    Returns
    -------
    exptime_key: str
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
        exptime_key = "unit_weight"

    return exptime_key

def do_combine1d(asn_file, exptime_key, interpolation):
    """Combine the input spectra.

    Parameters
    ----------
    asn_file: str
        The association file name.

    exptime_key: str
        A string identifying which keyword to use to get the exposure time,
        which is used as a weight when combining spectra.

    interpolation: str
        The type of interpolation between input pixels.
    """

    asn = json.load(open(asn_file))
    # Get the name to use for the output file.
    output_name = asn["products"][0]["name"]
    # Get the dictionary that contains the input file names.
    input_files = asn["products"][0]["members"]
    del asn

    exptime_key = check_exptime(exptime_key)

    input_spectra = []
    for file_dict in input_files:
        file = file_dict["expname"]
        ms = datamodels.MultiSpecModel(file)
        nspec = len(ms.spec)
        for in_spec in ms.spec:
            input_spectra.append(InputSpectrumModel(ms, in_spec, exptime_key))
        ms.close()

    output_spec = OutputSpectrumModel()
    output_spec.assign_wavelengths(input_spectra)
    output_spec.accumulate_sums(input_spectra, interpolation)
    output_spec.compute_combination()

    for in_spec in input_spectra:
        in_spec.close()

    out_model = output_spec.create_output()
    output_spec.close()

    # Copy one of the input headers to output.
    ms = datamodels.MultiSpecModel(input_files[0]["expname"])
    out_model.update(ms)# , primary_only=True)
    ms.close()

    out_model.meta.filename = output_name
    out_model.meta.cal_step.combine_1d = 'COMPLETE'

    out_model.save(output_name)
    out_model.close()

def correct_model(input_file, exptime_key, interpolation):
    """Combine 1-D spectra."""

    do_combine1d(input_file, exptime_key, interpolation)

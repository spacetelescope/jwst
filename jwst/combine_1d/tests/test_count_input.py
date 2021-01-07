"""
Test for combine1d.count_input
"""
import numpy as np

# from jwst import datamodels
from jwst.combine_1d import combine1d

class DummySpectra:
    def __init__(self, wavelength):
        # The argument to count_input is a list of InputSpectrumModel
        # objects, but the only attribute we need is wavelength, so this
        # DummySpectra class will suffice.
        self.wavelength = wavelength.copy()


def test_1a():
    input_spectra = create_input_spectra_1()

    (wl, n_input_spectra) = combine1d.count_input(input_spectra)

    truth_1a = truth_array_1a()

    # Compare wl with the correct values.
    assert np.allclose(wl, truth_1a, rtol=1.e-10)


def test_1b():
    input_spectra = create_input_spectra_1()
    (wl, n_input_spectra) = combine1d.count_input(input_spectra)

    truth_1b = truth_array_1b()

    # Compare n_input_spectra with the correct values.
    assert np.allclose(n_input_spectra, truth_1b, rtol=1.e-10)


def test_1c():
    input_spectra = create_input_spectra_1()
    (wl, n_input_spectra) = combine1d.count_input(input_spectra)

    output_wl = combine1d.compute_output_wl(wl, n_input_spectra)

    truth_1c = truth_array_1c()

    assert np.allclose(output_wl, truth_1c, rtol=1.e-10)


def create_input_spectra_1():

    input_spectra = []

    # This is for the case that there are four input spectra, each one
    # six pixels long.  The wavelengths are roughly the same in each
    # spectrum.  The differences in wavelength between spectra is small
    # compared to the spacing between wavelengths in any one spectrum.
    # That is, the union of the wavelengths in all the input spectra
    # will be clumped.  For this test, first specify the wavelengths in
    # these clumps.  Then assign the wavelengths for each input spectrum
    # by picking one element from each clump.
    clumps = np.array([[1.372, 1.385, 1.391, 1.402],
                       [1.701, 1.718, 1.719, 1.723],
                       [1.979, 1.984, 1.991, 2.013],
                       [2.212, 2.223, 2.228, 2.310],
                       [2.475, 2.479, 2.487, 2.496],
                       [2.716, 2.723, 2.729, 2.741]], dtype=np.float64)

    spec = [clumps[0, 3], clumps[1, 0], clumps[2, 2], clumps[3, 1], clumps[4, 2], clumps[5, 1]]
    input_spectra.append(DummySpectra(np.array(spec)))
    spec = [clumps[0, 0], clumps[1, 3], clumps[2, 0], clumps[3, 3], clumps[4, 0], clumps[5, 2]]
    input_spectra.append(DummySpectra(np.array(spec)))
    spec = [clumps[0, 2], clumps[1, 1], clumps[2, 1], clumps[3, 2], clumps[4, 3], clumps[5, 0]]
    input_spectra.append(DummySpectra(np.array(spec)))
    spec = [clumps[0, 1], clumps[1, 2], clumps[2, 3], clumps[3, 0], clumps[4, 1], clumps[5, 3]]
    input_spectra.append(DummySpectra(np.array(spec)))

    return input_spectra


def truth_array_1a():
    """Merged and sorted wavelengths from all input files."""

    # Wavelengths (merged and sorted) in input files.
    return np.array([1.372, 1.385, 1.391, 1.402, 1.701, 1.718, 1.719, 1.723,
                     1.979, 1.984, 1.991, 2.013, 2.212, 2.223, 2.228, 2.31,
                     2.475, 2.479, 2.487, 2.496, 2.716, 2.723, 2.729, 2.741],
                    dtype=np.float64)


def truth_array_1b():
    """Number of input models that cover each of the wavelengths in the
       combined array.

    The first element in the wl array (see function truth_array_1a) has
    value 1.372, and that's not actually within the array of wavelengths
    for the first input model, which are (see function
    create_input_spectra_1):
        [1.402 1.701 1.991 2.223 2.487 2.723]
    But the array that we'll return has 4 (the number of input models)
    for every element.  How can that be correct?  It's because we assume
    the pixels in the input models have finite width, but the wavelengths
    in the input models are the values at the centers of the pixels.  So
    the range of wavelengths that are actually covered by an input model
    is not just the interval between the first and last wavelengths in
    the input model's wavelength array, it's the range from the wavelength
    at the left edge of the first pixel to the wavelength at the right edge
    of the last pixel.
    """

    # n_input_spectra
    return np.array([4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], dtype=np.int64)


def truth_array_1c():
    """Computed wavelengths for the output.

    These are averages of the input wavelengths within clusters, sets
    of input wavelengths that are closer together than the spacing
    between pixels.
    """

    return np.array([1.3875, 1.71525, 1.99175, 2.24325, 2.48425, 2.72725],
                    dtype=np.float64)

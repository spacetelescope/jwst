"""
Create a miniature, slightly simplified model of the SOSS detector/optics.

Use the model detector to instantiate an extraction engine, use the engine to create mock data
from a known input spectrum, and then check if the engine can retrieve that spectrum
from the data.

The model has the following features:
- Factor-of-10 smaller along each dimension
- Similar wavelengths for each of the two orders
- Partially overlapping traces for the two orders
- Randomly-selected bad pixels in the data
- Wave grid of size ~100 with varying resolution
- Triangle function throughput for each spectral order
- Kernel is also a triangle function peaking at the center, or else unity for certain tests
- (partial) Mock of the Pastasoss reference model
"""

import pytest
import numpy as np
from scipy.signal import savgol_filter
from functools import partial
from jwst.extract_1d.soss_extract import atoca
from jwst.extract_1d.soss_extract import atoca_utils as au
from stdatamodels.jwst.datamodels import PastasossModel

PWCPOS = 245.85932900002442
DATA_SHAPE = (25, 200)
WAVE_BNDS_O1 = [2.8, 0.8]
WAVE_BNDS_O2 = [1.4, 0.5]
WAVE_BNDS_GRID = [0.7, 2.7]
ORDER1_SCALING = 20.0
ORDER2_SCALING = 2.0
TRACE_END_IDX = [DATA_SHAPE[1], 180]
SPECTRAL_SLOPE = 2


@pytest.fixture(scope="package")
def wave_map():
    """
    Make a mock wave map.

    Returns
    -------
    list of array[float]
        2D wavelength arrays for each order
    """
    wave_ord1 = np.linspace(WAVE_BNDS_O1[0], WAVE_BNDS_O1[1], DATA_SHAPE[1])
    wave_ord1 = np.ones(DATA_SHAPE) * wave_ord1[np.newaxis, :]

    wave_ord2 = np.linspace(WAVE_BNDS_O2[0], WAVE_BNDS_O2[1], DATA_SHAPE[1])
    wave_ord2 = np.ones(DATA_SHAPE) * wave_ord2[np.newaxis, :]
    # add a small region of zeros to mimic what is input to the step from ref files
    wave_ord2[:, TRACE_END_IDX[1] :] = 0.0

    return [wave_ord1, wave_ord2]


@pytest.fixture(scope="package")
def trace_profile():
    """
    Make a mock trace profile.

    Order 2 is partially on top of, and partially not on top of, order 1.
    Give order 2 some slope here to simulate that; Order 1 is flat.

    Returns
    -------
    list of array[float]
        2D trace profiles for each order
    """
    # order 1
    ord1 = np.zeros(DATA_SHAPE[0])
    ord1[3:9] = 0.2
    ord1[2] = 0.1
    ord1[9] = 0.1
    profile_ord1 = np.ones(DATA_SHAPE) * ord1[:, np.newaxis]

    # order 2
    yy, xx = np.meshgrid(np.arange(DATA_SHAPE[0]), np.arange(DATA_SHAPE[1]))
    yy = yy.astype(np.float32) - xx.astype(np.float32) * 0.08
    yy = yy.T

    profile_ord2 = np.zeros_like(yy)
    full = (yy >= 3) & (yy < 9)
    half0 = (yy >= 9) & (yy < 11)
    half1 = (yy >= 1) & (yy < 3)
    profile_ord2[full] = 0.2
    profile_ord2[half0] = 0.1
    profile_ord2[half1] = 0.1

    return [profile_ord1, profile_ord2]


@pytest.fixture(scope="package")
def trace1d(wave_map, trace_profile):
    """
    Build 1D traces from the 2D traces.

    Returns
    -------
    list of tuple of array[float]
        (x, y, wavelength) for each order
    """
    trace_list = []
    for order in [0, 1]:
        profile = trace_profile[order]
        wave2d = wave_map[order].copy()  # avoid modifying wave_map, it's needed elsewhere!
        end_idx = TRACE_END_IDX[order]

        # find mean y-index at each x where trace_profile is nonzero
        shp = profile.shape
        xx, yy = np.mgrid[: shp[0], : shp[1]]
        xx = xx.astype("float")
        xx[profile == 0] = np.nan
        mean_trace = np.nanmean(xx, axis=0)
        # same strategy for wavelength
        wave2d[profile == 0] = np.nan
        mean_wave = np.nanmean(wave2d, axis=0)

        # smooth it. we know it should be linear so use 1st order poly and large box size
        ytrace = savgol_filter(mean_trace, mean_trace.size - 1, 1)
        wavetrace = savgol_filter(mean_wave, mean_wave.size - 1, 1)

        # apply cutoff
        wavetrace = wavetrace[:end_idx]
        ytrace = ytrace[:end_idx]
        xtrace = np.arange(0, end_idx)

        trace_list.append((xtrace, ytrace, wavetrace))

    return trace_list


@pytest.fixture(scope="package")
def wave_grid():
    """
    Create an irregularly-spaced 1-D wavelength grid.

    wave_grid has smaller spacings in some places than others
    and is not in long-to-short order like the wave map.
    Two wavelength values are duplicated here on purpose, for testing.

    Returns
    -------
    array[float]
        1-D wavelength grid
    """
    lo0 = np.linspace(WAVE_BNDS_GRID[0], 1.2, 16)
    hi = np.linspace(1.2, 1.7, 46)
    lo2 = np.linspace(1.7, WAVE_BNDS_GRID[1], 31)
    return np.concatenate([lo0, hi, lo2])


@pytest.fixture(scope="package")
def throughput():
    """
    Make a triangle function for each order but with different peak wavelength.

    Returns
    -------
    list of function
        Throughput functions for each order
    """

    def filter_function(wl, wl_max):
        # Set free parameters to roughly mimic throughput functions on main.
        maxthru = 0.4
        thresh = 0.01
        scaling = 0.3
        dist = np.abs(wl - wl_max)
        thru = maxthru - dist * scaling
        thru[thru < thresh] = thresh
        return thru

    thru_o1 = partial(filter_function, wl_max=1.7)
    thru_o2 = partial(filter_function, wl_max=0.7)

    return [thru_o1, thru_o2]


@pytest.fixture(scope="package")
def kernels_unity():
    """
    Make a trivial kernel.

    Returns
    -------
    list of array[float]
        Kernel for each order, both set to unity
    """
    return [
        np.array(
            [
                1.0,
            ]
        ),
        np.array(
            [
                1.0,
            ]
        ),
    ]


@pytest.fixture(scope="package")
def webb_kernels(wave_map):
    """
    Toy model of the JWST kernels.

    Let the kernel be a triangle function along the pixel position axis,
    peaking at the center, and independent of wavelength.

    Same for both orders except the wavelengths and wave_trace are different.

    Returns
    -------
    list of WebbKernel
        WebbKernel for each order
    """
    n_os = 5
    n_pix = 15  # full pixel width of kernel
    n_wave = 10
    peakiness = 4
    min_val = 1
    kernel_width = n_os * n_pix - (n_os - 1)
    ctr_idx = kernel_width // 2

    kernels = []
    for order in [0, 1]:
        # set up wavelength grid over which kernel is defined
        wave_trace = wave_map[order][0]
        wave_range = (np.min(wave_trace), np.max(wave_trace))
        wavelengths = np.linspace(*wave_range, n_wave)
        wave_kernel = np.ones((kernel_width, wavelengths.size), dtype=float) * wavelengths[None, :]

        # model kernel as simply a triangle function that peaks at the center
        triangle_function = ctr_idx - peakiness * np.abs(ctr_idx - np.arange(0, kernel_width))
        triangle_function[triangle_function <= min_val] = min_val
        kernel = np.ones((kernel_width, wavelengths.size), dtype=float) * triangle_function[:, None]
        kernel /= np.sum(kernel)

        kernels.append(au.WebbKernel(wave_kernel, kernel, wave_trace, n_pix))

    return kernels


@pytest.fixture(scope="package")
def mask_trace_profile(trace_profile):
    """
    Make masks corresponding to the trace profiles, set to unity where outside the profile.

    Returns
    -------
    list of array[bool]
        Mask for each order
    """

    def mask_from_trace(trace_in, cut_low=None, cut_hi=None):
        trace = trace_in.copy()
        trace[trace_in <= 0] = 1
        trace[trace_in > 0] = 0
        if cut_low is not None:
            trace[:, :cut_low] = 1
        if cut_hi is not None:
            trace[:, cut_hi:] = 1
        return trace.astype(bool)

    trace_o1 = mask_from_trace(trace_profile[0], cut_low=0, cut_hi=199)
    trace_o2 = mask_from_trace(trace_profile[1], cut_low=0, cut_hi=175)
    return [trace_o1, trace_o2]


@pytest.fixture(scope="package")
def detector_mask():
    """
    Make a mock detector mask made up of a few random bad pixels.

    Returns
    -------
    array[bool]
        Detector mask
    """
    rng = np.random.default_rng(42)
    mask = np.zeros(DATA_SHAPE, dtype=bool)
    bad = rng.choice(mask.size, 100)
    bad = np.unravel_index(bad, DATA_SHAPE)
    mask[bad] = 1
    return mask


@pytest.fixture(scope="package")
def engine(
    wave_map,
    trace_profile,
    throughput,
    kernels_unity,
    wave_grid,
    mask_trace_profile,
    detector_mask,
):  # numpydoc ignore=RT01
    """Instantiate the extraction engine from the mock detector model."""
    return atoca.ExtractionEngine(
        wave_map,
        trace_profile,
        throughput,
        kernels_unity,
        wave_grid,
        mask_trace_profile,
        global_mask=detector_mask,
    )


def f_lam(wl, m=SPECTRAL_SLOPE, b=0):
    """
    Return a linear model of flux as function of wavelength.

    Parameters
    ----------
    wl : array[float]
        Wavelength
    m : float
        Slope of linear model
    b : float
        Intercept of linear model

    Returns
    -------
    array[float]
        Flux
    """
    return m * wl + b


@pytest.fixture(scope="package")
def imagemodel(engine, detector_mask):
    """
    Use engine.rebuild to make an image model from an expected f(lambda).

    Returns
    -------
    array[float], array[float]
        Data and error arrays

    Notes
    -----
    The data coming out of this fixture display some unphysical striping.
    This is likely due to inadequate wavelength sampling in the mock data.
    It has no real effect for the unit testing suite, but may be worth
    in the future if this test suite is expanded.
    """
    rng = np.random.default_rng(seed=42)
    shp = engine.trace_profile[0].shape

    # make the detector bad values NaN, but leave the trace masks alone
    # in reality, of course, this is backward: detector bad values
    # would be determined from data
    data = engine.rebuild(f_lam, fill_value=0.0)
    data[detector_mask] = np.nan

    # add noise
    noise_scaling = 3e-5
    data += noise_scaling * rng.standard_normal(shp)

    # error random, all positive, but also always larger than a certain value
    # to avoid very large values of data / error
    error = noise_scaling * (rng.standard_normal(shp) ** 2 + 0.5)

    return data, error


@pytest.fixture(scope="module")
def refmodel(trace1d):
    """
    Make a mock Pastasoss reference model.

    Spatial dimensions are scaled down by a factor of 10.
    Since the traces are just linear, the polynomials
    also have coefficients equal to 0 except for the constant and linear terms.

    Returns
    -------
    PastasossModel
        Mock reference model
    """
    model = PastasossModel()
    model.meta.pwcpos_cmd = 245.76

    trace0 = {
        "pivot_x": 189.0,
        "pivot_y": 5.0,
        "spectral_order": 1,
        "trace": np.array([trace1d[0][0], trace1d[0][1]], dtype=np.float64).T,
        "padding": 0,
    }
    trace1 = {
        "pivot_x": 168.0,
        "pivot_y": 20.0,
        "spectral_order": 2,
        "trace": np.array([trace1d[1][0], trace1d[1][1]], dtype=np.float64).T,
    }
    model.traces = [trace0, trace1]

    wavecal0 = {
        "coefficients": [
            WAVE_BNDS_O1[0],
            -2.0,
        ]
        + [0.0 for i in range(19)],
        "polynomial_degree": 5,
        "scale_extents": [[0, -1.03552000e-01], [DATA_SHAPE[1], 1.62882080e-01]],
    }
    wavecal1 = {
        "coefficients": [
            WAVE_BNDS_O2[0],
            -1.0,
        ]
        + [0.0 for i in range(8)],
        "polynomial_degree": 3,
        "scale_extents": [[0, 245.5929], [DATA_SHAPE[1], 245.9271]],
    }
    model.wavecal_models = [wavecal0, wavecal1]

    thru0 = {
        "spectral_order": 1,
        "wavelength": np.linspace(0.5, 5.5, 501),
        "throughput": np.ones((501,)),
    }  # peaks around 1.22 at value 0.37
    thru1 = {
        "spectral_order": 2,
        "wavelength": np.linspace(0.5, 5.5, 501),
        "throughput": np.ones((501,)),
    }  # peaks around 0.7 at value of 0.16
    model.throughputs = [thru0, thru1]

    return model


@pytest.fixture
def ref_files(refmodel):
    """
    Make a dictionary of reference files.

    Returns
    -------
    dict
        Dictionary of reference files
    """
    ref_files = {"pastasoss": refmodel}
    ref_files["subarray"] = "SUBSTRIP256"
    ref_files["pwcpos"] = PWCPOS
    return ref_files

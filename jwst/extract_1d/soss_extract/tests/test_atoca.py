import pytest
import numpy as np
from functools import partial
from scipy.sparse import csr_matrix
from jwst.extract_1d.soss_extract import atoca
from jwst.extract_1d.soss_extract.tests.conftest import (
    SPECTRAL_SLOPE,
    f_lam,
    DATA_SHAPE,
    WAVE_BNDS_O1,
    WAVE_BNDS_O2,
)

"""Tests for the ATOCA extraction engine, taking advantage of the miniature
model set up by conftest.py.
The test_call() function ensures that the engine can retrieve the spectrum
with SPECTRAL_SLOPE that we put into the data, which implicitly checks
a lot of the matrix math."""


def test_extraction_engine_init(
    wave_map,
    trace_profile,
    throughput,
    wave_grid,
    mask_trace_profile,
    detector_mask,
    engine,
):
    """Test the init of the engine with default/good inputs"""

    # test wave_grid became unique
    assert engine.wave_grid.dtype == np.float64
    unq = np.unique(wave_grid)
    assert np.allclose(engine.wave_grid, unq)
    assert engine.n_wavepoints == unq.size

    for order in [0, 1]:
        # test assignment of attributes and conversion to expected float64 dtype
        assert engine.wave_map[order].dtype == np.float64
        assert engine.trace_profile[order].dtype == np.float64
        assert engine.kernels[order].dtype == np.float64
        assert engine.mask_trace_profile[order].dtype == np.bool_

        assert np.allclose(engine.wave_map[order], wave_map[order])
        assert np.allclose(engine.trace_profile[order], trace_profile[order])
        assert np.allclose(engine.mask_trace_profile[order], mask_trace_profile[order])

        # test derived attributes
        assert engine.data_shape == DATA_SHAPE
        assert engine.n_orders == 2

    # test wave_p and wave_m. separate unit test for their calculation
    for att in ["wave_p", "wave_m"]:
        wave = getattr(engine, att)
        assert wave.dtype == np.float64
        assert wave.shape == (2,) + DATA_SHAPE

    # test _get_i_bounds
    assert len(engine.i_bounds) == 2
    for order in [0, 1]:
        assert len(engine.i_bounds[order]) == 2
        assert engine.i_bounds[order][0] >= 0
        assert engine.i_bounds[order][1] < DATA_SHAPE[1]

    # test to ensure that wave_map is considered for bounds
    # in order 1 the wave_map is more restrictive on the shortwave end
    # in order 2 the wave_grid is more restrictive on the longwave end
    # in order 1 no restriction on the longwave end so get the full extent
    # in order 2 no restriction on the shortwave end so get the full extent
    assert engine.i_bounds[0][0] > 0
    assert engine.i_bounds[1][1] < engine.n_wavepoints

    # test _get_masks
    # ensure they all include the input detector bad pixel mask
    for mask in [engine.mask_ord[0], engine.mask_ord[1], engine.mask, engine.general_mask]:
        assert mask.dtype == np.bool_
        assert mask.shape == DATA_SHAPE
        assert np.all(mask[detector_mask == 1])

    # general_mask should be a copy of mask
    assert np.allclose(engine.mask, engine.general_mask)

    wave_bnds = [WAVE_BNDS_O1, WAVE_BNDS_O2]
    for order in [0, 1]:
        # ensure wavelength bounds from wave_grid are respected in mask_ord
        mask = engine.mask_ord[order]
        wls = wave_map[order]
        lo, hi = wave_bnds[order]
        outside = (wls > lo) | (wls < hi)
        assert np.all(mask[outside])

        # a bit paradoxically, engine.mask_ord does not contain a single order's mask_trace_profile
        # instead, it's mask_trace_profile[0] AND mask_trace_profile[1], i.e.,
        # the trace profiles of BOTH orders are UNmasked in the masks of each order
        # the general mask and wavelength bounds are then applied, so the
        # only difference between mask_ord[0] and mask_ord[1] are the wavelength bounds
        # test that NOT all locations masked by mask_trace_profile are masked in mask_ord
        assert not np.all(mask[mask_trace_profile[order]])
        # test that all locations masked by both profiles are masked in mask_ord
        combined_profile = mask_trace_profile[0] & mask_trace_profile[1]
        assert np.all(mask[combined_profile])

    # test throughput function conversion to array
    for order in [0, 1]:
        thru = engine.throughput[order]
        assert thru.size == engine.n_wavepoints
        assert np.all(thru >= 0)
        assert thru.dtype == np.float64

    # test kernel is cast to proper shape for input (trivial) kernel
    for order in [0, 1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        expected_shape = (n_valid, engine.n_wavepoints)
        assert engine.kernels[order].shape == expected_shape
        # for trivial kernel only one element per row is nonzero
        assert engine.kernels[order].count_nonzero() == expected_shape[0]

    # test weights. see separate unit tests to ensure the calculation is correct
    for order in [0, 1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        weights = engine.weights[order]
        k_idx = engine.weights_k_idx[order]
        assert weights.dtype == np.float64
        assert np.issubdtype(k_idx.dtype, np.integer)

        # why is weights the same size as ~engine.mask, and not ~engine.mask_ord?
        assert weights.shape == (np.sum(~engine.mask), n_valid)
        assert k_idx.shape[0] == weights.shape[0]

    # test assignment of empty attributes
    for att in ["w_t_wave_c", "tikho_mat", "_tikho_mat"]:
        assert hasattr(engine, att)
    assert getattr(engine, "pixel_mapping") == [None, None]


def test_extraction_engine_bad_inputs(
    wave_map,
    trace_profile,
    throughput,
    kernels_unity,
    wave_grid,
    mask_trace_profile,
    detector_mask,
):
    # not enough good pixels in order
    with pytest.raises(atoca.MaskOverlapError):
        detector_mask = np.ones_like(detector_mask)
        detector_mask[5:7, 50:55] = 0  # still a few good pixels but very few
        atoca.ExtractionEngine(
            wave_map,
            trace_profile,
            throughput,
            kernels_unity,
            wave_grid,
            mask_trace_profile,
            global_mask=detector_mask,
        )

    # wrong number of orders
    with pytest.raises(ValueError):
        atoca.ExtractionEngine(
            wave_map,
            trace_profile,
            throughput,
            kernels_unity,
            wave_grid,
            mask_trace_profile,
            global_mask=detector_mask,
            orders=[
                0,
            ],
        )


def test_get_attributes(engine):
    # test string input
    assert np.allclose(engine.get_attributes("wave_map"), engine.wave_map)

    # test list of strings input
    name_list = ["wave_map", "wave_grid"]
    att_list = engine.get_attributes(*name_list)
    expected = [engine.wave_map, engine.wave_grid]
    for i in range(len(expected)):
        for j in range(2):  # orders
            assert np.allclose(att_list[i][j], expected[i][j])

    # test i_order not None
    att_list = engine.get_attributes(*name_list, i_order=1)
    expected = [engine.wave_map[1], engine.wave_grid[1]]
    for i in range(len(expected)):
        assert np.allclose(att_list[i], expected[i])


def test_update_throughput(engine, throughput):
    old_thru = engine.throughput

    # test callable input
    new_thru = [throughput[1], throughput[1]]
    engine.update_throughput(new_thru)
    for i, thru in enumerate(engine.throughput):
        assert isinstance(thru, np.ndarray)
        assert thru.shape == engine.wave_grid.shape
    assert np.allclose(engine.throughput[0], engine.throughput[1])

    # test array input
    # first reset to old throughput
    engine.throughput = old_thru
    new_thru = [thru * 2 for thru in engine.throughput]
    engine.update_throughput(new_thru)
    for i, thru in enumerate(engine.throughput):
        assert np.allclose(thru, old_thru[i] * 2)

    # test fail on bad array shape
    new_thru = [thru[:-1] for thru in engine.throughput]
    with pytest.raises(ValueError):
        engine.update_throughput(new_thru)

    # test fail on callable that doesn't return correct array shape
    def new_thru_f(wl):
        return 1.0

    with pytest.raises(TypeError):
        engine.update_throughput([new_thru_f, new_thru_f])


def test_create_kernels(webb_kernels, engine):
    """test_atoca_utils.test_get_c_matrix already tests the creation
    of individual kernels for different input types, here just ensure
    the options get passed into that function properly"""

    kernels_0 = engine._create_kernels(webb_kernels)
    kernels_1 = engine._create_kernels([None, None])

    for kernel_list in [kernels_0, kernels_1]:
        assert len(kernel_list) == 2
        for order in [0, 1]:
            kern = kernel_list[order]
            assert isinstance(kern, csr_matrix)
            assert kern.dtype == np.float64


def test_wave_grid_c(engine):
    for order in [0, 1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        assert engine.wave_grid_c(order).size == n_valid


def test_set_w_t_wave_c(engine):
    """all this does is copy whatever is input"""
    product = np.zeros((1,))
    engine._set_w_t_wave_c(0, product)
    assert len(engine.w_t_wave_c) == engine.n_orders
    assert engine.w_t_wave_c[0] == product
    assert engine.w_t_wave_c[1] == []
    assert product is not engine.w_t_wave_c[0]


def test_get_pixel_mapping(engine):
    pixel_mapping_0 = engine.get_pixel_mapping(0)
    # check attribute is set and identical to output
    # check the second one is not set but there is space for it
    assert hasattr(engine, "pixel_mapping")
    assert len(engine.pixel_mapping) == engine.n_orders
    assert np.allclose(engine.pixel_mapping[0].data, pixel_mapping_0.data)
    assert engine.pixel_mapping[1] is None

    # set the second one so can check both at once
    engine.get_pixel_mapping(1)
    for order in [0, 1]:
        mapping = engine.pixel_mapping[order]
        assert mapping.dtype == np.float64
        # why is this the shape, instead of using mask_ord and only the valid wave_grid?
        expected_shape = (np.sum(~engine.mask), engine.wave_grid.size)
        assert mapping.shape == expected_shape

        # test that w_t_wave_c is getting saved
        w_t_wave_c = engine.w_t_wave_c[order]
        assert w_t_wave_c.dtype == np.float64
        assert w_t_wave_c.shape == expected_shape

        # check if quick=True works
        mapping_quick = engine.get_pixel_mapping(order, quick=True)
        assert np.allclose(mapping.data, mapping_quick.data)

    # check that quick=True does not work if w_t_wave_c unsaved
    engine.w_t_wave_c = None
    with pytest.raises(AttributeError):
        engine.get_pixel_mapping(1, quick=True)


def test_rebuild(engine):
    detector_model = engine.rebuild(f_lam)
    assert detector_model.dtype == np.float64
    assert detector_model.shape == engine.wave_map[0].shape

    # test that input spectrum is ok as either callable or array
    assert np.allclose(detector_model, engine.rebuild(f_lam(engine.wave_grid)))

    # test fill value
    detector_model_nans = engine.rebuild(f_lam, fill_value=np.nan)
    assert np.allclose(np.isnan(detector_model_nans), engine.general_mask)


def test_build_sys(imagemodel, engine):
    data, error = imagemodel
    matrix, result = engine.build_sys(data, error)
    assert result.size == engine.n_wavepoints
    assert matrix.shape == (result.size, result.size)


def test_get_detector_model(imagemodel, engine):
    data, error = imagemodel
    unmasked_size = np.sum(~engine.mask)
    b_matrix, data_matrix = engine.get_detector_model(data, error)

    assert data_matrix.shape == (1, unmasked_size)
    assert b_matrix.shape == (unmasked_size, engine.n_wavepoints)
    assert np.allclose(data_matrix.toarray()[0], (data / error)[~engine.mask])


def test_estimate_tikho_factors(engine):
    factor = engine.estimate_tikho_factors(f_lam)
    assert isinstance(factor, float)


@pytest.fixture(scope="module")
def tikho_tests(imagemodel, engine):
    data, error = imagemodel

    log_guess = np.log10(engine.estimate_tikho_factors(f_lam))
    factors = np.logspace(log_guess - 9, log_guess + 9, 19)
    return factors, engine.get_tikho_tests(factors, data, error)


def test_get_tikho_tests(tikho_tests, engine):
    factors, tests = tikho_tests
    unmasked_size = np.sum(~engine.mask)

    # test all the output shapes
    assert np.allclose(tests["factors"], factors)
    assert tests["solution"].shape == (len(factors), engine.n_wavepoints)
    assert tests["error"].shape == (len(factors), unmasked_size)
    assert tests["reg"].shape == (len(factors), engine.n_wavepoints - 1)
    assert tests["chi2"].shape == (len(factors),)
    assert tests["chi2_soft_l1"].shape == (len(factors),)
    assert tests["chi2_cauchy"].shape == (len(factors),)
    assert np.allclose(tests["grid"], engine.wave_grid)

    # test data type is preserved through solve
    for key in tests.keys():
        assert tests[key].dtype == np.dtype("float64")


def test_best_tikho_factor(engine, tikho_tests):
    input_factors, tests = tikho_tests
    fit_modes = ["all", "curvature", "chi2", "d_chi2"]
    best_factors = []
    for mode in fit_modes:
        factor = engine.best_tikho_factor(tests, mode)
        assert isinstance(factor, float)
        best_factors.append(factor)

    # ensure fit_mode=all found one of the three others
    assert best_factors[0] in best_factors[1:]

    # test the logic tree by manually changing the tests dict
    # is non-trivial because dchi2 and curvature
    # are both derived from other keys in the dictionary, not just the statistical metrics


def test_call(engine, tikho_tests, imagemodel):
    """
    Run the actual extract method.
    Ensure it can retrieve the input spectrum based on f_lam to within a few percent
    at all points on the wave_grid.

    Note this round-trip implicitly checks the math of the build_sys, get_detector_model,
    _solve, and _solve_tikho, at least at first blush.
    """
    data, error = imagemodel
    _, tests = tikho_tests
    best_factor = engine.best_tikho_factor(tests, "all")

    expected_spectrum = f_lam(engine.wave_grid)
    for tikhonov in [True, False]:
        spectrum = engine(data, error, tikhonov=tikhonov, factor=best_factor)
        diff = (spectrum - expected_spectrum) / expected_spectrum
        assert not np.all(np.isnan(diff))
        diff = diff[~np.isnan(diff)]
        assert np.all(np.abs(diff) < 0.05)

    # test bad input, failing to put factor in for Tikhonov solver
    with pytest.raises(ValueError):
        engine(data, error, tikhonov=True)


def test_compute_likelihood(engine, imagemodel):
    """Ensure log-likelihood is highest for the correct slope"""

    data, error = imagemodel
    test_slopes = np.arange(0, 5, 0.5)
    logl = []
    for slope in test_slopes:
        spectrum = partial(f_lam, m=slope)
        logl.append(engine.compute_likelihood(spectrum, data, error))

    assert np.argmax(logl) == np.argwhere(test_slopes == SPECTRAL_SLOPE)
    assert np.all(np.array(logl) < 0)

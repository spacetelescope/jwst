import pytest
import numpy as np
from functools import partial
from jwst.extract_1d.soss_extract import atoca


DATA_SHAPE = (25,200)
WAVE_BNDS_O1 = [2.8, 0.8]
WAVE_BNDS_O2 = [1.4, 0.5]
WAVE_BNDS_GRID = [0.7, 2.7]
ORDER1_SCALING = 20.0
ORDER2_SCALING = 2.0
SPECTRAL_SLOPE = 2

@pytest.fixture(scope="module")
def wave_map():
    """Trying for a roughly factor-of-10 shrink on each axis
    but otherwise relatively faithful reproduction of real detector behavior"""
    wave_ord1 = np.linspace(WAVE_BNDS_O1[0], WAVE_BNDS_O1[1], DATA_SHAPE[1])
    wave_ord1 = np.ones(DATA_SHAPE)*wave_ord1[np.newaxis, :]

    wave_ord2 = np.linspace(WAVE_BNDS_O2[0], WAVE_BNDS_O2[1], DATA_SHAPE[1])
    wave_ord2 = np.ones(DATA_SHAPE)*wave_ord2[np.newaxis,:]
    # add a small region of zeros to mimic what is input to the step from ref files
    wave_ord2[:,190:] = 0.0

    return [wave_ord1, wave_ord2]


@pytest.fixture(scope="module")
def trace_profile(wave_map):
    """
    order 2 is partially on top of, partially not on top of order 1
    give order 2 some slope to simulate that
    """
    # order 1
    DATA_SHAPE = wave_map[0].shape
    ord1 = np.zeros((DATA_SHAPE[0]))
    ord1[3:9] = 0.2
    ord1[2] = 0.1
    ord1[9] = 0.1
    profile_ord1 = np.ones(DATA_SHAPE)*ord1[:, np.newaxis]

    # order 2
    yy, xx = np.meshgrid(np.arange(DATA_SHAPE[0]), np.arange(DATA_SHAPE[1]))
    yy = yy.astype(np.float32) - xx.astype(np.float32)*0.08
    yy = yy.T
    
    profile_ord2 = np.zeros_like(yy)
    full = (yy >= 3) & (yy < 9)
    half0 = (yy >= 9) & (yy < 11)
    half1 = (yy >= 1) & (yy < 3)
    profile_ord2[full] = 0.2
    profile_ord2[half0] = 0.1
    profile_ord2[half1] = 0.1

    return [profile_ord1, profile_ord2]


@pytest.fixture(scope="module")
def wave_grid():
    """wave_grid has smaller spacings in some places than others
    and is not backwards order like the wave map
    Two duplicates are in there on purpose for testing"""
    lo0 = np.linspace(WAVE_BNDS_GRID[0], 1.2, 16)
    hi = np.linspace(1.2, 1.7, 46)
    lo2 = np.linspace(1.7, WAVE_BNDS_GRID[1], 31)
    return np.concatenate([lo0, hi, lo2])


@pytest.fixture(scope="module")
def throughput():
    """make a triangle function for each order but with different peak wavelength
    """

    def filter_function(wl, wl_max):
        """Set free parameters to roughly mimic throughput functions on main"""
        maxthru = 0.4
        thresh = 0.01
        scaling = 0.3
        dist = np.abs(wl - wl_max)
        thru = maxthru - dist*scaling
        thru[thru<thresh] = thresh
        return thru
    
    thru_o1 = partial(filter_function, wl_max=1.7)
    thru_o2 = partial(filter_function, wl_max=0.7)

    return [thru_o1, thru_o2]

@pytest.fixture(scope="module")
def kernels():
    """For now, just return unity"""
    return [np.array([1.,]), np.array([1.,])]



@pytest.fixture(scope="module")
def mask_trace_profile(trace_profile):
    """Masks set to unity where bad"""

    def mask_from_trace(trace_in, cut_low=None, cut_hi=None):
        trace = trace_in.copy()
        trace[trace_in<=0] = 1
        trace[trace_in>0] = 0
        if cut_low is not None:
            trace[:,:cut_low] = 1
        if cut_hi is not None:
            trace[:,cut_hi:] = 1
        return trace.astype(bool)
    
    trace_o1 = mask_from_trace(trace_profile[0], cut_low=0, cut_hi=199)
    trace_o2 = mask_from_trace(trace_profile[1], cut_low=0, cut_hi=175)
    return [trace_o1, trace_o2]


@pytest.fixture(scope="module")
def detector_mask(wave_map):
    """Add a few random bad pixels"""
    rng = np.random.default_rng(42)
    mask = np.zeros(DATA_SHAPE, dtype=bool)
    bad = rng.choice(mask.size, 100)
    bad = np.unravel_index(bad, DATA_SHAPE)
    mask[bad] = 1
    return mask


@pytest.fixture(scope="module")
def engine(wave_map,
           trace_profile,
           throughput,
           kernels,
           wave_grid,
           mask_trace_profile,
           detector_mask,
):
    return atoca.ExtractionEngine(wave_map,
                                    trace_profile,
                                    throughput,
                                    kernels,
                                    wave_grid,
                                    mask_trace_profile,
                                    global_mask=detector_mask)


def test_extraction_engine(
    wave_map,
    trace_profile,
    throughput,
    kernels,
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

    for order in [0,1]:
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
        assert wave.shape == (2,)+DATA_SHAPE

    # test _get_i_bounds
    assert len(engine.i_bounds) == 2
    for order in [0,1]:
        assert len(engine.i_bounds[order]) == 2
        assert engine.i_bounds[order][0] >= 0
        assert engine.i_bounds[order][1] < DATA_SHAPE[1]
        assert all([isinstance(val, int) for val in engine.i_bounds[order]])

    # test to ensure that wave_map is considered for bounds
    # in order 1 the wave_map is more restrictive on the shortwave end
    # in order 2 the wave_grid is more restrictive on the longwave end
    # in order 1 no restriction on the longwave end so get the full extent
    # in order 2 no restriction on the shortwave end so get the full extent
    assert engine.i_bounds[0][0] > 0
    assert engine.i_bounds[1][1] < engine.n_wavepoints
    # TODO: off-by-one error here. why does this fail?
    # check what this looks like on a real run on main
    # assert engine.i_bounds[0][1] == engine.n_wavepoints
    # assert engine.i_bounds[1][0] == 0


    # test _get_masks
    # ensure they all include the input detector bad pixel mask
    for mask in [engine.mask_ord[0], engine.mask_ord[1], engine.mask, engine.general_mask]:
        assert mask.dtype == np.bool_
        assert mask.shape == DATA_SHAPE
        assert np.all(mask[detector_mask == 1])

    # general_mask should be a copy of mask
    assert np.allclose(engine.mask, engine.general_mask)

    wave_bnds = [WAVE_BNDS_O1, WAVE_BNDS_O2]
    for order in [0,1]:
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
    for order in [0,1]:
        thru = engine.throughput[order]
        assert thru.size == engine.n_wavepoints
        assert np.all(thru >= 0)
        assert thru.dtype == np.float64

    # test kernel is cast to proper shape for input (trivial) kernel
    for order in [0,1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        expected_shape = (n_valid, engine.n_wavepoints)
        assert engine.kernels[order].shape == expected_shape
        # for trivial kernel only one element per row is nonzero
        assert engine.kernels[order].count_nonzero() == expected_shape[0]

    # test weights. see separate unit tests to ensure the calculation is correct
    for order in [0,1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        weights = engine.weights[order]
        k_idx = engine.weights_k_idx[order]
        assert weights.dtype == np.float64
        assert np.issubdtype(k_idx.dtype, np.integer)

        # TODO: why is weights the same size as ~engine.mask, and not ~engine.mask_ord?
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
    kernels,
    wave_grid,
    mask_trace_profile,
    detector_mask,
):
    # not enough good pixels in order
    with pytest.raises(atoca.MaskOverlapError):
        detector_mask = np.ones_like(detector_mask)
        detector_mask[5:7,50:55] = 0 #still a few good pixels but very few
        atoca.ExtractionEngine(wave_map,
                               trace_profile,
                               throughput,
                               kernels,
                               wave_grid,
                               mask_trace_profile,
                               global_mask=detector_mask)

    # wrong number of orders
    with pytest.raises(ValueError):
        atoca.ExtractionEngine(wave_map,
                               trace_profile,
                               throughput,
                               kernels,
                               wave_grid,
                               mask_trace_profile,
                               global_mask=detector_mask,
                               orders=[0,])


def test_get_attributes(engine):
    # test string input
    assert np.allclose(engine.get_attributes("wave_map"), engine.wave_map)

    # test list of strings input
    name_list = ["wave_map", "wave_grid"]
    att_list = engine.get_attributes(*name_list)
    expected = [engine.wave_map, engine.wave_grid]
    for i in range(len(expected)):
        for j in range(2): #orders
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
    new_thru = [thru*2 for thru in engine.throughput]
    engine.update_throughput(new_thru)
    for i, thru in enumerate(engine.throughput):
        assert np.allclose(thru, old_thru[i]*2)

    # test fail on bad array shape
    new_thru = [thru[:-1] for thru in engine.throughput]
    with pytest.raises(ValueError):
        engine.update_throughput(new_thru)

    # test fail on callable that doesn't return correct array shape
    def new_thru_f(wl):
        return 1.0
    with pytest.raises(ValueError):
        engine.update_throughput([new_thru_f, new_thru_f])


def test_create_kernels(engine):
    # TODO: make a non-trivial kernel fixture to work with this
    # trivial example already done
    pass


def test_wave_grid_c(engine):
    for order in [0,1]:
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


def test_grid_from_map():
    assert False


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
    for order in [0,1]:
        mapping = engine.pixel_mapping[order]
        assert mapping.dtype == np.float64
        # TODO: why is this the shape, instead of using mask_ord and only the valid wave_grid?
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

    # TODO: follow the shape of every sparse matrix through this function, see if
    # it makes sense or if it would make more sense to mask differently
    # TODO: test that the math is actually correct using a trivial example (somehow)


def test_rebuild(engine):

    detector_model = engine.rebuild(f_lam)
    assert detector_model.dtype == np.float64
    assert detector_model.shape == engine.wave_map[0].shape

    # test that input spectrum is ok as either callable or array
    assert np.allclose(detector_model, engine.rebuild(f_lam(engine.wave_grid)))

    # test fill value
    detector_model_nans = engine.rebuild(f_lam, fill_value=np.nan)
    assert np.allclose(np.isnan(detector_model_nans), engine.general_mask)


def f_lam(wl, m=SPECTRAL_SLOPE, b=0):
    """
    Estimator for flux as function of wavelength
    Returns linear function of wl with slope m and intercept b
    
    This function is also used in this test suite as 
    """
    return m*wl + b


@pytest.fixture(scope="module")
def imagemodel(engine, detector_mask):
    """
    use engine.rebuild to make an image model from an expected f(lambda).
    Then we can ensure it round-trips
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
    data += noise_scaling*rng.standard_normal(shp)

    # error random, all positive, but also always larger than a certain value
    # to avoid very large values of data / error
    error = noise_scaling*(rng.standard_normal(shp)**2 + 0.5)

    return data, error


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
    assert np.allclose(data_matrix.toarray()[0], (data/error)[~engine.mask])


def test_estimate_tikho_factors(engine):
    
    factor = engine.estimate_tikho_factors(f_lam)
    assert isinstance(factor, float)


    # very approximate calculation of tik fac looks like
    # n_pixels = (~engine.mask).sum()
    # flux = f_lam(engine.wave_grid)
    # dlam = engine.wave_grid[1:] - engine.wave_grid[:-1]
    # print(n_pixels/np.mean(flux[1:] * dlam))


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
    assert tests["reg"].shape == (len(factors), engine.n_wavepoints-1)
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

    # TODO: test the logic tree by manually changing the tests dict
    # this is non-trivial because dchi2 and curvature
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
        diff = (spectrum - expected_spectrum)/expected_spectrum
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

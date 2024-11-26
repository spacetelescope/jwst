import pytest
import numpy as np
from functools import partial
from jwst.extract_1d.soss_extract import atoca


DATA_SHAPE = (25,200)
WAVE_BNDS_O1 = [2.8, 0.8]
WAVE_BNDS_O2 = [1.4, 0.5]
WAVE_BNDS_GRID = [0.7, 2.7]

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
    return [np.array([1.,])]



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


def test_extraction_engine(wave_map,
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


    # test masks
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
    print(engine.kernels[0].shape)
    print(engine.kernels[0].dtype)
    print(engine.kernels[0].toarray())



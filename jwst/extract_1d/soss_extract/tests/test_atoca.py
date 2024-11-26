import pytest
import numpy as np
from functools import partial
from jwst.extract_1d.soss_extract import atoca


@pytest.fixture(scope="module")
def wave_map():
    """Trying for a roughly factor-of-10 shrink on each axis
    but otherwise relatively faithful reproduction of real detector behavior"""
    shp = (25, 200)
    wave_ord1 = np.linspace(2.8, 0.8, shp[1])
    wave_ord1 = np.ones(shp)*wave_ord1[np.newaxis, :]

    wave_ord2 = np.linspace(1.4, 0.5, shp[1])
    wave_ord2 = np.ones(shp)*wave_ord2[np.newaxis,:]
    wave_ord2[:,190:] = 0.0

    return [wave_ord1, wave_ord2]

@pytest.fixture(scope="module")
def trace_profile(wave_map):
    """
    order 2 is partially on top of, partially not on top of order 1
    give order 2 some slope to simulate that
    """
    # order 1
    shp = wave_map[0].shape
    ord1 = np.zeros((shp[0]))
    ord1[3:9] = 0.2
    ord1[2] = 0.1
    ord1[9] = 0.1
    profile_ord1 = np.ones(shp)*ord1[:, np.newaxis]

    # order 2
    yy, xx = np.meshgrid(np.arange(shp[0]), np.arange(shp[1]))
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
    lo0 = np.linspace(0.7, 1.2, 6)
    hi = np.linspace(1.2, 1.7, 19)
    lo2 = np.linspace(1.7, 2.8, 12)
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
    shp = wave_map[0].shape
    rng = np.random.default_rng(42)
    mask = np.zeros(shp, dtype=bool)
    bad = rng.choice(mask.size, 100)
    bad = np.unravel_index(bad, shp)
    mask[bad] = 1
    return mask


def test_extraction_engine(wave_map,
                           trace_profile,
                           throughput,
                           kernels,
                           wave_grid,
                           mask_trace_profile,
                           detector_mask):

    engine = atoca.ExtractionEngine(wave_map,
                                    trace_profile,
                                    throughput,
                                    kernels,
                                    wave_grid,
                                    mask_trace_profile,
                                    global_mask=detector_mask)
    shp = wave_map[0].shape
    
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
        assert engine.data_shape == shp
        assert engine.n_orders == 2
        
    # test wave_p and wave_m. separate unit test for their calculation
    for att in ["wave_p", "wave_m"]:
        wave = getattr(engine, att)
        assert wave.dtype == np.float64
        assert wave.shape == (2,)+shp

    # test masks
    # ensure they all include the input detector bad pixel mask
    for mask in [engine.mask_ord[0], engine.mask_ord[1], engine.mask, engine.general_mask]:
        assert mask.dtype == np.bool_
        assert mask.shape == shp
        assert np.all(mask[detector_mask == 1])
    # general_mask should be a copy of mask except it doesn't respect wavelength bounds
    #assert np.allclose(engine.mask, engine.general_mask)
    # order masks should include mask of trace
    for order in [0,1]:
        mask = engine.mask_ord[order]
        mask_profile = mask_trace_profile[order]
        import matplotlib.pyplot as plt
        fig, (ax0, ax1) = plt.subplots(2,1,figsize = (10,10))
        ax0.imshow(mask, origin='lower')
        ax1.imshow(mask_profile, origin='lower')
        plt.show()
        #assert np.all(mask[mask_profile == 1])





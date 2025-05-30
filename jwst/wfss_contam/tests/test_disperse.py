import pytest
import numpy as np
from jwst.wfss_contam.disperse import flux_interpolator_injector, determine_wl_spacing

'''
Note that main disperse.py call is tested in test_observations.py because 
it requires all the fixtures defined there.
'''

@pytest.mark.parametrize("lams, flxs, extrapolate_sed, expected_outside_bounds", 
                         [([1, 3], [1, 3], False, 0), 
                          ([2], [2], False, 2), 
                          ([1, 3], [1, 3], True, 4)])
def test_interpolate_fluxes(lams, flxs, extrapolate_sed, expected_outside_bounds):

    flux_interpf = flux_interpolator_injector(lams, flxs, extrapolate_sed)
    assert flux_interpf(2.0) == 2.0
    assert flux_interpf(4.0) == expected_outside_bounds


@pytest.mark.parametrize("lams, expected_dw", 
                         [([1, 1.2, 1.4], 0.05),
                         ([1, 1.02, 1.04], 0.01)
                         ])
def test_determine_wl_spacing(lams, expected_dw):
    
    dw = 0.1
    oversample_factor = 2
    dw_out = determine_wl_spacing(dw, np.array(lams), oversample_factor)

    assert np.isclose(dw_out, expected_dw, atol=1e-8)

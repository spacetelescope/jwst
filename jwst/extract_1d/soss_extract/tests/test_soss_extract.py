import pytest
import numpy as np
from numpy.polynomial import Polynomial

from jwst.extract_1d.soss_extract.soss_extract import (
    _estim_flux_first_order, 
)
from .conftest import SPECTRAL_SLOPE


def test_estim_flux_first_order(imagemodel,
                                detector_mask,
                                wave_map,
                                trace_profile,
                                throughput,
                                mask_trace_profile):

    data, err = imagemodel
    ref_file_args = [wave_map, trace_profile, throughput, None]

    func = _estim_flux_first_order(data,
                                   err,
                                   detector_mask,
                                   ref_file_args,
                                   mask_trace_profile,
                                   threshold=1e-4)
    
    # use the function to generate a test spectrum
    test_wl = np.linspace(1.0, 2.5, 100)
    test_flux = func(test_wl)

    # check slope against expected value to within a few percent
    p_fitted = Polynomial.fit(test_wl, test_flux, 1)
    b, m = p_fitted.coef
    assert np.isclose(m, SPECTRAL_SLOPE, rtol=0.05)



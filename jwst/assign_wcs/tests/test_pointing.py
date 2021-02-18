import numpy as np

from jwst.assign_wcs import pointing
from astropy.modeling.models import Identity


def test_dva_corr_noop_missing_meta_values():
    assert isinstance(pointing.dva_corr_model(va_scale=None, v2_ref=1, v3_ref=1), Identity)
    assert isinstance(pointing.dva_corr_model(va_scale=1, v2_ref=2, v3_ref=3), Identity)


def test_dva_corr_valid_match():
    m = pointing.dva_corr_model(va_scale=1.001, v2_ref=-380, v3_ref=-770)
    assert np.allclose(m(1, 10), (380 * 0.001 + 1.001, 770 * 0.001 + 10.01))


def test_dva_corr_valid_match_no_shifts():
    scale = 1.001
    m = pointing.dva_corr_model(va_scale=scale, v2_ref=0, v3_ref=0)
    assert np.allclose(m(1, 10), (scale, 10 * scale))

    m = pointing.dva_corr_model(va_scale=scale, v2_ref=None, v3_ref=None)
    assert np.allclose(m(1, 10), (scale, 10 * scale))


def test_dva_corr_inverse():
    v2_ref = -380
    v3_ref = -770
    va_scale = 1.001
    test_v2 = 2300
    test_v3 = 5600
    fm = pointing.dva_corr_model(va_scale=va_scale, v2_ref=v2_ref, v3_ref=v3_ref)
    im = pointing.dva_corr_model(va_scale=1 / va_scale, v2_ref=v2_ref, v3_ref=v3_ref)
    assert np.allclose(im(*fm(test_v2, test_v3)), (test_v2, test_v3))

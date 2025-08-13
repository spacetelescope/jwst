import numpy as np
from stdatamodels.jwst.datamodels import dqflags

from jwst.coron import median_replace_img as mri

DNU = dqflags.pixel["DO_NOT_USE"]
NONSCI = dqflags.pixel["NON_SCIENCE"]


def test_median_replace_img_miri(target_model_miri):
    # add bad science pixels, one flagged, one just NaN
    target_model_miri.dq[:, 100, 100] = DNU
    target_model_miri.data[:, 110, 110] = np.nan

    # add bad non-science pixels, one flagged, one just NaN
    target_model_miri.dq[:, 200, 200] = NONSCI | DNU
    target_model_miri.data[:, 210, 210] = np.nan
    target_model_miri.dq[:, 210, 210] = NONSCI

    result = mri.median_replace_img(target_model_miri, box_size=3, bad_bitvalue=DNU)

    # Updated in place
    assert result is target_model_miri

    # Science data is filled in from neighboring pixels
    assert not np.any(result.data[:, 100, 100] == 0)
    np.testing.assert_allclose(result.data[:, 100, 100], result.data[:, 101, 101], atol=0.1)
    assert not np.any(result.data[:, 110, 110] == 0)
    np.testing.assert_allclose(result.data[:, 110, 110], result.data[:, 111, 111], atol=0.1)

    # Non-science data is set to 0
    np.testing.assert_array_equal(result.data[:, 200, 200], 0.0)
    np.testing.assert_array_equal(result.data[:, 210, 210], 0.0)


def test_median_replace_img_nircam(target_model):
    # add bad science pixels, one flagged, one just NaN
    target_model.dq[:, 100, 100] = DNU
    target_model.data[:, 110, 110] = np.nan

    # add bad non-science pixels, one flagged, one just NaN
    target_model.dq[:, 200, 200] = NONSCI | DNU
    target_model.data[:, 210, 210] = np.nan
    target_model.dq[:, 210, 210] = NONSCI

    result = mri.median_replace_img(target_model, box_size=3, bad_bitvalue=DNU)

    # Updated in place
    assert result is target_model

    # All flagged data is filled in for NIRCam
    assert not np.any(result.data == 0)
    np.testing.assert_allclose(result.data[:, 100, 100], result.data[:, 101, 101], atol=0.1)
    np.testing.assert_allclose(result.data[:, 110, 110], result.data[:, 111, 111], atol=0.1)
    np.testing.assert_allclose(result.data[:, 200, 200], result.data[:, 201, 201], atol=0.1)
    np.testing.assert_allclose(result.data[:, 210, 210], result.data[:, 211, 211], atol=0.1)


def test_median_fill_all_flagged():
    data = np.full((10, 10), 1.0)
    dq = np.full((10, 10), 1)
    xc, yc = 5, 5
    result = mri.median_fill_value(data, dq, 3, 1, xc, yc)
    assert result == 0.0


def test_median_fill_no_good_data():
    data = np.full((10, 10), np.nan)
    dq = np.full((10, 10), 0)
    xc, yc = 5, 5
    result = mri.median_fill_value(data, dq, 3, 1, xc, yc)
    assert result == 0.0


def test_median_fill_out_of_range():
    data = np.full((10, 10), 1.0)
    dq = np.full((10, 10), 0)
    xc, yc = 12, 14
    result = mri.median_fill_value(data, dq, 3, 1, xc, yc)
    assert result == 0.0

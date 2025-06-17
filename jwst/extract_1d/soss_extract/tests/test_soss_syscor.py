import pytest
import numpy as np

from jwst.extract_1d.soss_extract.soss_syscor import (
    soss_background,
    make_background_mask,
)


def test_soss_background(imagemodel, detector_mask, mask_trace_profile):
    data, err = imagemodel
    bkg_mask = ~mask_trace_profile[0] | ~mask_trace_profile[1] | detector_mask

    data_bkg, col_bkg = soss_background(data, detector_mask, bkg_mask)
    assert data_bkg.shape == data.shape
    assert col_bkg.size == data.shape[1]

    # check background now has mean zero
    mean_bkg = np.mean(data_bkg[~bkg_mask])
    assert np.isclose(mean_bkg, 0.0)

    # check col_bkg are at least close to the non-sigma-clipped version which is much easier to calculate
    # For the test case, there are no outliers so this should be quite a close match
    data[bkg_mask] = np.nan
    col_bkg_unclipped = np.nanmean(data, axis=0)
    assert np.allclose(col_bkg, col_bkg_unclipped)


def test_make_background_mask():
    rng = np.random.default_rng(seed=42)
    for sub in [96, 256, 2048]:
        shape = (sub, 2048)
        width = int(sub / 4)
        data = rng.normal(0.0, 1.0, shape)

        mask = make_background_mask(data, width)

        if sub == 96:
            expected_bad_frac = 1 / 4
        else:
            expected_bad_frac = 1 / 2

        bad_frac = np.sum(mask) / mask.size
        # test that bad fraction is computed properly for all modes
        assert np.isclose(bad_frac, expected_bad_frac)

        # test that mask=True is the high-flux pixels
        assert np.mean(data[mask]) > np.mean(data)

    # test unrecognized shape
    with pytest.raises(ValueError):
        data = rng.normal(0.0, 1.0, (40, 2048))
        make_background_mask(data, width)

import numpy as np
import pytest
from photutils.aperture import CircularAnnulus, RectangularAperture

from jwst.extract_1d import ifu
from jwst.extract_1d.extract_1d_step import Extract1dStep


@pytest.fixture()
def setup_params(mock_miri_ifu):
    bkg_data = mock_miri_ifu.data[0, :, :]
    temp_weightmap = mock_miri_ifu.weightmap[0, :, :]
    temp_weightmap[temp_weightmap > 1] = 1
    shape = np.shape(mock_miri_ifu.data)
    bmask = np.zeros([shape[1], shape[2]], dtype=bool)
    bmask[:] = False
    bmask[np.where(temp_weightmap == 0)] = True
    bkg_sigma_clip = 3.0
    step = Extract1dStep()
    ref_file = step.get_reference_file(mock_miri_ifu, "extract1d")
    extract_params = ifu.get_extract_parameters(ref_file, bkg_sigma_clip)
    func_params = [bkg_data, temp_weightmap, shape, bmask, bkg_sigma_clip, extract_params]
    return func_params


def test_apply_bkg_sigma_clip_point(mock_miri_ifu, setup_params):
    mock_miri_ifu.meta.target.source_type = "POINT"
    (bkg_data, temp_weightmap, shape, bmask, bkg_sigma_clip, extract_params) = setup_params
    method = extract_params["method"]
    subpixels = extract_params["subpixels"]
    # set a mock inner and outer radius
    inner_bkg = 9.4
    outer_bkg = 12.7
    x_center = float(shape[-1]) / 2.0
    y_center = float(shape[-2]) / 2.0
    position = (x_center, y_center)
    annulus = CircularAnnulus(position, r_in=inner_bkg, r_out=outer_bkg)

    bkg_table, maskclip = ifu._apply_bkg_sigma_clip(
        bkg_data, temp_weightmap, bmask, bkg_sigma_clip, annulus, method, subpixels
    )

    expected_shape = (shape[1], shape[2])
    assert np.shape(maskclip) == expected_shape
    assert float(bkg_table["aperture_sum"][0]) == 292332.0


def setup_extended(extract_params, shape):
    method = extract_params["method"]
    subpixels = extract_params["subpixels"]
    width = float(shape[-1])
    height = float(shape[-2])
    x_center = width / 2.0 - 0.5
    y_center = height / 2.0 - 0.5
    theta = 0.0
    position = (x_center, y_center)
    aperture = RectangularAperture(position, width, height, theta)
    return aperture, method, subpixels


def test_apply_bkg_sigma_clip_extended(setup_params):
    (bkg_data, temp_weightmap, shape, bmask, bkg_sigma_clip, extract_params) = setup_params
    aperture, method, subpixels = setup_extended(extract_params, shape)

    bkg_table, maskclip = ifu._apply_bkg_sigma_clip(
        bkg_data, temp_weightmap, bmask, bkg_sigma_clip, aperture, method, subpixels
    )

    expected_shape = (shape[1], shape[2])
    assert np.shape(maskclip) == expected_shape
    assert float(bkg_table["aperture_sum"][0]) == 3123750.0


def test_apply_bkg_sigma_clip_no_bkg_stat_data(setup_params):
    (bkg_data, temp_weightmap, shape, bmask, bkg_sigma_clip, extract_params) = setup_params
    aperture, method, subpixels = setup_extended(extract_params, shape)
    temp_weightmap *= 0.0
    bmask[np.where(temp_weightmap == 0)] = True

    bkg_table, maskclip = ifu._apply_bkg_sigma_clip(
        bkg_data, temp_weightmap, bmask, bkg_sigma_clip, aperture, method, subpixels
    )

    assert maskclip is None
    assert float(bkg_table["aperture_sum"][0]) == 0.0

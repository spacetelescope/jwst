from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from jwst.cube_build.msa_cube import MSACubeData


@pytest.fixture
def basic_pars():
    """Default parameters fixture using plain dictionaries."""
    return {
        "source": {
            "slitname": ["slit1"],
            "number_slits": 1,
            "bunit": "Jy/sr",
            "bunit_err": "Jy/sr",
        },
        "linear_wave": True,
        "scalexy": 0.1,
        "scalew": 0.001,
        "slit_frac": 1.0,
        "ra_center": 100.0,
        "dec_center": 0.0,
        "cube_pa": 45.0,
        "nspax_x": 10,
        "nspax_y": 10,
        "debug_spaxel": "1 2 3",
        "coord_system": "skyalign",
        "wavemin": 1.0,
        "wavemax": 2.0,
        "weighting": "volume",
        "wavelength_table": [1.0, 1.5, 2.0],
        "suffix": "cube",
    }


def test_find_spaxel_flux_normalization(basic_pars):
    """Ensure flux and variance normalization correctly skip unpopulated spaxels."""
    cube_data = MSACubeData(**basic_pars)

    spaxel_iflux = np.array([0.0, 1.0, 2.0])
    spaxel_flux = np.array([0.0, 10.0, 20.0])
    spaxel_weight = np.array([0.0, 2.0, 4.0])
    spaxel_var = np.array([0.0, 4.0, 16.0])

    cube_data.find_spaxel_flux(spaxel_iflux, spaxel_flux, spaxel_weight, spaxel_var)

    # Unpopulated spaxel (iflux == 0) is skipped and remains unchanged
    assert spaxel_flux[0] == 0.0
    assert spaxel_var[0] == 0.0
    # Populated spaxel 1: 10.0 / 2.0 = 5.0
    assert spaxel_flux[1] == 5.0
    # Variance 1: 4.0 / (2.0^2) = 1.0
    assert spaxel_var[1] == 1.0
    # Populated spaxel 2: 20.0 / 4.0 = 5.0
    assert spaxel_flux[2] == 5.0


def test_set_slit_wcs_linear(basic_pars):
    """Test standard WCS grid setup for linear wavelength mode."""
    cube_data = MSACubeData(**basic_pars)

    # Initialize cdelt parameters set normally by find_footprint()
    cube_data.cdelt1 = cube_data.scalexy
    cube_data.cdelt2 = cube_data.scalexy
    cube_data.cdelt3 = cube_data.scalew

    corner_ra = [99.9, 100.1, 100.1, 99.9]
    corner_dec = [-0.1, -0.1, 0.1, 0.1]
    lambda_min = 1.0
    lambda_max = 2.0
    rot_angle = 0.0

    # Mock spatial coordinate conversion to standard tangent projection coordinates
    def mock_radec2std(crval1, crval2, ra, dec, rot):
        return (ra - crval1) * 3600.0, (dec - crval2) * 3600.0

    with patch("jwst.cube_build.msa_cube.coord.radec2std", side_effect=mock_radec2std):
        cube_data.set_slit_wcs(corner_ra, corner_dec, lambda_min, lambda_max, rot_angle)

    assert cube_data.crval1 == 100.0
    assert cube_data.crval2 == 0.0
    assert cube_data.naxis1 == 21  # 10*2 + 1 based on nspax_x=10
    assert cube_data.naxis2 == 21  # 10*2 + 1 based on nspax_y=10
    assert cube_data.naxis3 == 1000  # range_lambda (1.0) / scalew (0.001)
    assert len(cube_data.xcoord) == cube_data.naxis1
    assert len(cube_data.zcoord) == cube_data.naxis3
    assert len(cube_data.cdelt3_normal) == cube_data.naxis3


def test_set_slit_wcs_nonlinear(basic_pars):
    """Test WCS wavelength array generation for non-linear wavelength mode."""
    basic_pars["linear_wave"] = False
    basic_pars["wavelength_table"] = [1.0, 1.2, 1.5, 2.0]
    cube_data = MSACubeData(**basic_pars)

    # Initialize cdelt attributes normally populated by find_footprint()
    cube_data.cdelt1 = cube_data.scalexy
    cube_data.cdelt2 = cube_data.scalexy

    def mock_radec2std(crval1, crval2, ra, dec, rot):
        return (ra - crval1) * 3600.0, (dec - crval2) * 3600.0

    with patch("jwst.cube_build.msa_cube.coord.radec2std", side_effect=mock_radec2std):
        cube_data.set_slit_wcs([99.9, 100.1], [-0.1, 0.1], 1.0, 2.0, 0.0)

    assert cube_data.naxis3 == 4
    np.testing.assert_array_equal(cube_data.zcoord, [1.0, 1.2, 1.5, 2.0])
    assert len(cube_data.cdelt3_normal) == 4


def test_map_source_pixels_to_output_frame_no_good_data(basic_pars):
    """Verify map_source_pixels_to_output_frame handles empty/invalid pixel data gracefully."""
    cube_data = MSACubeData(**basic_pars)
    cube_data.crval1 = 100.0
    cube_data.crval2 = 0.0
    cube_data.rot_angle = 0.0

    # Mock a SlitModel with completely unusable pixel flags
    slit = MagicMock()
    slit.data = np.ones((5, 5))
    slit.err = np.ones((5, 5))
    slit.dq = np.full((5, 5), 1)  # Set bit 1 (DO_NOT_USE)
    slit.var_rnoise = np.ones((5, 5))

    # Generic mock transform that handles 2 or 3 input coordinate arrays
    def mock_transform(*args):
        x = args[0]
        return np.zeros_like(x), np.zeros_like(x), np.full_like(x, 1.5)

    mock_wcs = MagicMock()
    mock_wcs.bounding_box = ((0, 4), (0, 4))
    mock_wcs.side_effect = mock_transform
    mock_wcs.get_transform.return_value = mock_transform

    slit.meta.wcs = mock_wcs

    with patch(
        "jwst.cube_build.msa_cube.wcstools.grid_from_bounding_box", return_value=np.zeros((2, 5, 5))
    ):
        results = cube_data.map_source_pixels_to_output_frame(slit)

    x_array, y_array, flux_array = results[:3]
    assert len(x_array) == 0
    assert len(y_array) == 0
    assert len(flux_array) == 0


def test_setup_final_model(basic_pars):
    """Verify setup_final_model returns a correctly formatted IFUCubeModel."""
    cube_data = MSACubeData(**basic_pars)
    cube_wcs = (2, 2, 2, 0.1, 0.1, 0.001, 1.5, 1.5, 1.0, 100.0, 0.0, 1.0)

    # 2x2x2 cube flattening to size 8
    spaxel_flux = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    spaxel_iflux = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    spaxel_var = np.array([0.1] * 8)
    spaxel_dq = np.zeros(8, dtype=np.uint32)
    spaxel_counter = np.ones(8, dtype=np.uint32)
    zcoord = np.array([1.0, 1.001])

    # Create mock WCS with a valid ASDF tag to pass stdatamodels schema validation
    mock_wcs = MagicMock()
    mock_wcs._tag = "tag:stsci.edu:gwcs/wcs-1.0.0"

    with patch("jwst.cube_build.msa_cube.pointing.create_fitswcs", return_value=mock_wcs):
        cube_model = cube_data.setup_final_model(
            cube_wcs,
            0.0,
            spaxel_flux,
            spaxel_iflux,
            spaxel_var,
            spaxel_dq,
            spaxel_counter,
            zcoord,
            "Jy/sr",
            "Jy/sr",
        )

    assert cube_model.data.shape == (2, 2, 2)
    assert cube_model.meta.wcsinfo.crval1 == 100.0
    assert cube_model.meta.wcsinfo.crval2 == 0.0
    assert cube_model.meta.bunit_data == "Jy/sr"


def test_init_parses_debug_spaxel(basic_pars):
    """Verify initialization and string parsing of debug spaxel indices."""
    cube_data = MSACubeData(**basic_pars)
    assert cube_data.spaxel_x == 1
    assert cube_data.spaxel_y == 2
    assert cube_data.spaxel_z == 3
    assert cube_data.scalexy == 0.1
    assert cube_data.linear_wave is True

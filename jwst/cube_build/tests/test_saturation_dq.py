import gwcs
import numpy as np
import pytest
from astropy.modeling.models import Mapping
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.cube_build import ifu_cube

SHAPE = (1, 1)


# ---------------------------------------------------------------------------
# Pytest Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_wcs():
    """Simple GWCS pipeline instance matching test_wcs.py conventions."""
    input_frame = gwcs.Frame2D(name="detector")
    output_frame = gwcs.Frame2D(name="world")
    pipeline = [(input_frame, Mapping((0, 1, 1))), (output_frame, None)]
    return gwcs.WCS(pipeline)


@pytest.fixture
def mock_miri_image_model(mock_wcs):
    """Returns a real JWST IFUImageModel populated with MIRI metadata."""
    input_model = datamodels.IFUImageModel()
    input_model.meta.instrument.name = "MIRI"
    input_model.meta.instrument.detector = "MIRIFULONG"
    input_model.meta.instrument.channel = "34"
    input_model.meta.instrument.band = "SHORT"
    input_model.meta.filename = "test_drizzle_saturated.fits"

    input_model.data = np.zeros(SHAPE, dtype=np.float32)
    input_model.err = np.zeros(SHAPE, dtype=np.float32)
    input_model.dq = np.zeros(SHAPE, dtype=np.uint32)
    input_model.meta.wcs = mock_wcs
    return input_model


@pytest.fixture
def drizzle_cube_instance(mock_miri_image_model):
    """Fixture to create and configure an IFUCubeData instance."""
    cube = ifu_cube.IFUCubeData.__new__(ifu_cube.IFUCubeData)

    cube.output_name = "test_drizzle_cube.fits"
    cube.input_models = [mock_miri_image_model]

    # Instrument and coordinate metadata configuration
    cube.instrument = "MIRI"
    cube.coord_system = "sky"
    cube.list_par1 = ["1"]
    cube.list_par2 = ["SHORT"]

    # Initialize spatial/spectral dimensions
    cube.naxis1 = 1
    cube.naxis2 = 1
    cube.naxis3 = 1
    cube.linear_wave = True
    cube.interpolation = "drizzle"
    cube.weighting = "drizzle"
    cube.rois = 1.0
    cube.roiw = 1.0
    cube.weight_power = 1.0
    cube.soft_rad = 1.0
    cube.scalerad = 1.0
    cube.offsets = None

    # Reference coordinates
    cube.crval1 = 0.0
    cube.crval2 = 0.0
    cube.crval3 = 1.0
    cube.crpix1 = 1.0
    cube.crpix2 = 1.0
    cube.crpix3 = 1.0
    cube.cdelt1 = 0.1
    cube.cdelt2 = 0.1
    cube.cdelt3 = 0.1

    cube.rot_angle = 0.0
    cube.zcoord = np.array([1.0])
    cube.cdelt3_normal = np.array([0.1])

    # Pre-allocate output arrays
    total_num = cube.naxis1 * cube.naxis2 * cube.naxis3
    cube.spaxel_flux = np.zeros(total_num, dtype=np.float64)
    cube.spaxel_weight = np.zeros(total_num, dtype=np.float64)
    cube.spaxel_var = np.zeros(total_num, dtype=np.float64)
    cube.spaxel_iflux = np.zeros(total_num, dtype=np.float64)
    cube.spaxel_dq = np.zeros(total_num, dtype=np.uint32)

    return cube


# ---------------------------------------------------------------------------
# Unit Tests
# ---------------------------------------------------------------------------


def test_drizzle_saturated_pixel_dq_propagation(drizzle_cube_instance, mock_miri_image_model):
    """
    Test that detector pixels with DQ = SATURATED pass through properly
    during drizzle processing, bitwise-OR into spaxel_dq, and survive flux division.
    """
    cube = drizzle_cube_instance
    sat_flag = dqflags.pixel["SATURATED"]

    # Populate saturated input image using real datamodel
    input_model = mock_miri_image_model
    input_model.data[0, 0] = 500.0
    input_model.err[0, 0] = 2.0
    input_model.dq[0, 0] = sat_flag

    corner_coords = np.zeros((8, 1), dtype=np.float64)

    sky_result = (
        np.array([0]),  # x
        np.array([0]),  # y
        np.array([0.0]),  # ra
        np.array([0.0]),  # dec
        np.array([1.0]),  # wave
        np.array([1]),  # slice_no
        np.array([0.1]),  # dwave
        corner_coords,  # corner_coord_all
    )

    cube.map_miri_pixel_to_sky = lambda *args, **kwargs: sky_result

    res = cube.map_detector_to_outputframe("MEDIUM", input_model)

    dq_out = res[7]
    assert (dq_out[0] & sat_flag) == sat_flag

    drizzle_spaxel_dq = np.array([sat_flag], dtype=np.uint32)
    cube.spaxel_dq = np.bitwise_or(cube.spaxel_dq, drizzle_spaxel_dq)

    cube.spaxel_weight[0] = 0.85
    cube.spaxel_flux[0] = 425.0
    cube.spaxel_iflux[0] = 1.0

    cube.find_spaxel_flux()
    cube.set_final_dq_flags()

    assert cube.spaxel_flux[0] == pytest.approx(500.0)
    assert (cube.spaxel_dq[0] & sat_flag) == sat_flag


def test_drizzle_saturated_nan_flux_sets_dq(drizzle_cube_instance):
    """
    Test that when all input detector pixels are SATURATED (flux=NaN),
    the C routine's SATURATED DQ bit is preserved when set_final_dq_flags()
    appends NON_SCIENCE and DO_NOT_USE for weight == 0 spaxels.
    """

    cube = drizzle_cube_instance
    sat_flag = dqflags.pixel["SATURATED"]
    do_not_use_flag = dqflags.pixel["DO_NOT_USE"]
    non_science_flag = dqflags.pixel["NON_SCIENCE"]  # 512

    # Simulate C routine output (cube_match_sky_driz):
    # - Bitwise-OR detector DQ (SATURATED) into spaxel_dq
    # - npy_isnan(flux) causes spaxel_flux and spaxel_weight to stay 0.0
    cube.spaxel_dq[0] = sat_flag
    cube.spaxel_flux[0] = 0.0
    cube.spaxel_weight[0] = 0.0
    cube.spaxel_iflux[0] = 0.0

    # Run flux processing to catch NaN fluxes and append DO_NOT_USE
    cube.find_spaxel_flux()
    cube.set_final_dq_flags()

    # 3. Construct final output model
    dq_cube = cube.spaxel_dq.reshape((cube.naxis3, cube.naxis2, cube.naxis1))
    data_cube = cube.spaxel_flux.reshape((cube.naxis3, cube.naxis2, cube.naxis1))
    err_cube = np.sqrt(cube.spaxel_var).reshape((cube.naxis3, cube.naxis2, cube.naxis1))
    weight_cube = cube.spaxel_weight.reshape((cube.naxis3, cube.naxis2, cube.naxis1))

    final_model = datamodels.IFUCubeModel(
        data=data_cube,
        err=err_cube,
        dq=dq_cube,
        weight=weight_cube,
    )

    final_dq = final_model.dq[0, 0, 0]

    # Verify SATURATED (2), DO_NOT_USE (1), and NON_SCIENCE (512) are all set (total: 515)
    assert (final_dq & sat_flag) == sat_flag
    assert (final_dq & do_not_use_flag) == do_not_use_flag
    assert (final_dq & non_science_flag) == non_science_flag

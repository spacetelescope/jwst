from types import SimpleNamespace

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import dqflags

from jwst.cube_build.ifu_cube import IFUCubeData


class AutoNamespace(SimpleNamespace):
    """Dynamic namespace that auto-instantiates missing attributes as sub-namespaces."""

    def __getattr__(self, name):
        val = AutoNamespace()
        setattr(self, name, val)
        return val


class FakeWCS:
    """Fake WCS implementation returning static sky coordinates."""

    def __call__(self, x, y):
        return np.array([0.0]), np.array([0.0]), np.array([1.0])


class FakeMeta(AutoNamespace):
    """Fake meta object storing image attributes and nested metadata structures."""

    def __init__(self, filename="test_drizzle_saturated.fits"):
        super().__init__()
        self.filename = filename
        self.wcs = FakeWCS()
        self.instrument.name = "MIRI"


class FakeInputModel:
    """Fake input data model containing array attributes and metadata."""

    def __init__(self, data, err, dq):
        self.data = data
        self.err = err
        self.dq = dq
        self.meta = FakeMeta()


class FakeIFUCubeModel:
    """Fake output cube model recording initialization arguments and metadata updates."""

    def __init__(self, **kwargs):
        self.dq = kwargs.get("dq")
        self.data = kwargs.get("data")
        self.err = kwargs.get("err")
        self.weight = kwargs.get("weight")
        self.meta = kwargs.get("meta", FakeMeta())

    def update(self, reference_model, **kwargs):
        """Stub for datamodel update method to copy reference metadata."""
        pass

    def save(self, *args, **kwargs):
        pass


# ---------------------------------------------------------------------------
# Pytest Fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def drizzle_cube_instance(monkeypatch):
    """Fixture to create an IFUCubeData instance without executing __init__."""

    # Avoid real __init__ execution using monkeypatch
    monkeypatch.setattr(IFUCubeData, "__init__", lambda self: None)
    cube = IFUCubeData()

    # Initialize output and input metadata attributes
    cube.output_name = "test_drizzle_cube.fits"
    cube.input_models = [
        FakeInputModel(np.zeros((1, 1)), np.zeros((1, 1)), np.zeros((1, 1), dtype=np.uint32))
    ]

    # Instrument and coordinate metadata configuration
    cube.instrument = "MIRI"
    cube.coord_system = "sky"  # Coordinate frame: 'sky' or 'internal_cal'
    cube.list_par1 = ["1"]  # MIRI Channels (e.g. 1, 2, 3, 4)
    cube.list_par2 = ["SHORT"]  # MIRI Bands (e.g. SHORT, MEDIUM, LONG)

    # Initialize spatial/spectral dimensions and config for DRIZZLE
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

    # Reference coordinate values (CRVAL)
    cube.crval1 = 0.0
    cube.crval2 = 0.0
    cube.crval3 = 1.0

    # Reference pixel positions (CRPIX)
    cube.crpix1 = 1.0
    cube.crpix2 = 1.0
    cube.crpix3 = 1.0

    # Coordinate step sizes (CDELT)
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


def test_drizzle_saturated_pixel_dq_propagation(drizzle_cube_instance, monkeypatch):
    """
    Test that detector pixels with DQ = SATURATED pass through properly
    during drizzle processing, bitwise-OR into spaxel_dq, and survive flux division.
    """
    cube = drizzle_cube_instance
    sat_flag = dqflags.pixel["SATURATED"]

    # 1. Create stub input data model with a SATURATED detector pixel
    input_model = FakeInputModel(
        data=np.array([[500.0]]), err=np.array([[2.0]]), dq=np.array([[sat_flag]], dtype=np.uint32)
    )

    # 8 corner coordinates: ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4
    corner_coords = np.zeros((8, 1), dtype=np.float64)

    # Fake mapping result from map_miri_pixel_to_sky
    sky_result = (
        np.array([0]),  # x
        np.array([0]),  # y
        np.array([0.0]),  # ra
        np.array([0.0]),  # dec
        np.array([1.0]),  # wave
        np.array([1]),  # slice_no
        np.array([0.1]),  # dwave
        corner_coords,  # corner_coord_all indexed along axis 0
    )

    # Override method using monkeypatch
    monkeypatch.setattr(cube, "map_miri_pixel_to_sky", lambda *args, **kwargs: sky_result)

    # 2. Run map_detector_to_outputframe for drizzle
    res = cube.map_detector_to_outputframe("MEDIUM", input_model)

    # Verify SATURATED DQ flag was preserved through map_detector_to_outputframe
    dq_out = res[7]
    assert (dq_out[0] & sat_flag) == sat_flag

    # 3. Simulate C-extension / C-drizzle output updating spaxel_dq
    drizzle_spaxel_dq = np.array([sat_flag], dtype=np.uint32)
    cube.spaxel_dq = np.bitwise_or(cube.spaxel_dq, drizzle_spaxel_dq)

    # Drizzle weight accumulation
    cube.spaxel_weight[0] = 0.85
    cube.spaxel_flux[0] = 425.0
    cube.spaxel_iflux[0] = 1.0

    # 4. Finalize spaxel flux division
    cube.find_spaxel_flux()
    cube.set_final_dq_flags()

    # Verify output flux and preserved DQ flag
    assert cube.spaxel_flux[0] == pytest.approx(500.0)
    assert (cube.spaxel_dq[0] & sat_flag) == sat_flag


def test_drizzle_saturated_nan_flux_sets_do_not_use(drizzle_cube_instance, monkeypatch):
    """
    Test that if drizzle interpolation results in NaN flux for a saturated spaxel,
    the DO_NOT_USE DQ flag is correctly appended during model setup.
    """
    cube = drizzle_cube_instance
    sat_flag = dqflags.pixel["SATURATED"]
    do_not_use_flag = dqflags.pixel["DO_NOT_USE"]

    # Assign NaN flux under drizzle combination with SATURATED DQ
    cube.spaxel_flux[0] = np.nan
    cube.spaxel_weight[0] = 0.5
    cube.spaxel_iflux[0] = 1.0
    cube.spaxel_dq[0] = sat_flag

    ref_model = FakeInputModel(
        data=np.array([[0.0]]), err=np.array([[0.0]]), dq=np.array([[0]], dtype=np.uint32)
    )

    # Intercept IFUCubeModel creation and create_fitswcs via monkeypatch
    last_created_cube = {}

    def fake_ifucube_constructor(**kwargs):
        inst = FakeIFUCubeModel(**kwargs)
        last_created_cube["instance"] = inst
        return inst

    import jwst.cube_build.ifu_cube as ifu_module

    monkeypatch.setattr(ifu_module.datamodels, "IFUCubeModel", fake_ifucube_constructor)
    # Return an AutoNamespace so setting bounding_box succeeds
    monkeypatch.setattr(
        ifu_module.pointing, "create_fitswcs", lambda *args, **kwargs: AutoNamespace()
    )

    # Execute model creation setup
    cube.setup_final_ifucube_model(ref_model)

    # Inspect final generated DQ array
    final_dq = last_created_cube["instance"].dq

    # Verify DO_NOT_USE was set alongside SATURATED for the NaN flux spaxel
    assert (final_dq[0, 0, 0] & sat_flag) == sat_flag
    assert (final_dq[0, 0, 0] & do_not_use_flag) == do_not_use_flag

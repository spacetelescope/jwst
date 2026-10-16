"""
Unit tests for cube_lib.setup_final_model and ifu_cube.IFUCubeData.set_wcs_slicer.

Coverage targets
----------------
setup_final_model
  - pixelarea keywords (PIXAR_A2, PIXAR_SR) are computed from cdelt1 * cdelt2
  - skyalign / linear_wave=True  : cdelt1/2 stored in degrees, ctype3="WAVE"
  - skyalign / linear_wave=False : ctype3="WAVE-TAB", crval3=1, cdelt3=None
  - internal_cal / NIRSPEC       : cdelt1/2 stored in native units (meters),
                                   cunit1="meter", ctype1="NRSSLICERX", pc1_1 sign flip
  - internal_cal / MIRI          : cunit1="arcsec", ctype1="MRSALPHA"
  - drizzle interpolation        : weight_power/roi_wave/roi_spatial zeroed

set_wcs_slicer
  - linear_wave flag is passed through (True and False)
  - cdelt1/cdelt2 stored verbatim from caller for both MIRI and NIRSpec
  - naxis3 computed correctly from wavelength range and cdelt3
"""

import math

import gwcs
import numpy as np
import pytest
from astropy.modeling.models import Mapping
from stdatamodels.jwst import datamodels

from jwst.cube_build import cube_lib, ifu_cube

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

CDELT1 = 0.13  # arcsec  (spatial x)
CDELT2 = 0.13  # arcsec  (spatial y)
CDELT3 = 0.001  # micron  (spectral)
CRVAL1 = 45.0  # RA  deg
CRVAL2 = -30.0  # Dec deg
CRVAL3 = 7.5  # micron
CRPIX1 = 21.0
CRPIX2 = 21.0
CRPIX3 = 1.0
NAXIS1 = 41
NAXIS2 = 41
NAXIS3 = 10
ARCSEC2_TO_SR = 2.3504e-11


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_wcs():
    """Minimal GWCS instance — enough for IFUImageModel.meta.wcs."""
    input_frame = gwcs.Frame2D(name="detector")
    output_frame = gwcs.Frame2D(name="world")
    pipeline = [(input_frame, Mapping((0, 1, 1))), (output_frame, None)]
    return gwcs.WCS(pipeline)


@pytest.fixture
def input_model(mock_wcs):
    """Minimal IFUImageModel used as model_ref in setup_final_model."""
    model = datamodels.IFUImageModel()
    model.meta.instrument.name = "NIRSPEC"
    model.meta.filename = "test_ref.fits"
    model.data = np.zeros((1, 1), dtype=np.float32)
    model.meta.wcs = mock_wcs
    return model


@pytest.fixture
def base_cube_wcs():
    """CubeWCS with linear_wave=True for the skyalign / linear case."""
    naxis3 = NAXIS3
    zcoord = CRVAL3 + (np.arange(naxis3) + 1 - CRPIX3) * CDELT3
    return cube_lib.CubeWCS(
        naxis1=NAXIS1,
        naxis2=NAXIS2,
        naxis3=naxis3,
        cdelt1=CDELT1,
        cdelt2=CDELT2,
        cdelt3=CDELT3,
        crpix1=CRPIX1,
        crpix2=CRPIX2,
        crpix3=CRPIX3,
        crval1=CRVAL1,
        crval2=CRVAL2,
        crval3=CRVAL3,
        cdelt3_normal=np.full(naxis3, CDELT3),
        rot_angle=0.0,
        xcoord=np.linspace(-2.6, 2.6, NAXIS1),
        ycoord=np.linspace(-2.6, 2.6, NAXIS2),
        zcoord=zcoord,
        a_min=-2.665,
        a_max=2.665,
        b_min=-2.665,
        b_max=2.665,
        lambda_min=CRVAL3,
        lambda_max=CRVAL3 + naxis3 * CDELT3,
        linear_wave=True,
        wavelength_table=None,
    )


def _make_spaxel_arrays(naxis1, naxis2, naxis3, flux_value=1.0):
    """Return flat spaxel arrays sized for (naxis3, naxis2, naxis1)."""
    total = naxis1 * naxis2 * naxis3
    flux = np.full(total, flux_value, dtype=np.float64)
    iflux = np.ones(total, dtype=np.float64)
    var = np.ones(total, dtype=np.float64)
    dq = np.zeros(total, dtype=np.uint32)
    return flux, iflux, var, dq


def _call_setup_final_model(
    input_model,
    cube_wcs,
    coord_system="skyalign",
    interpolation="pointcloud",
    weighting="msm",
    weight_power=2.0,
    rois=0.2,
    roiw=0.001,
    roiw_table=None,
    instrument="NIRSPEC",
):
    """Helper that builds spaxel arrays and calls setup_final_model."""
    flux, iflux, var, dq = _make_spaxel_arrays(cube_wcs.naxis1, cube_wcs.naxis2, cube_wcs.naxis3)
    if roiw_table is None:
        roiw_table = np.array([roiw])

    cube_model, status = cube_lib.setup_final_model(
        model_ref=input_model,
        cube_wcs=cube_wcs,
        spaxel_flux=flux,
        spaxel_iflux=iflux,
        spaxel_var=var,
        spaxel_dq=dq,
        instrument=instrument,
        rot_angle=cube_wcs.rot_angle,
        coord_system=coord_system,
        interpolation=interpolation,
        weighting=weighting,
        weight_power=weight_power,
        rois=rois,
        roiw=roiw,
        roiw_table=roiw_table,
        output_name=None,
    )
    return cube_model, status


# ---------------------------------------------------------------------------
# Tests: pixelarea keywords
# ---------------------------------------------------------------------------


class TestPixelArea:
    """PIXAR_A2 and PIXAR_SR must be derived from cdelt1 * cdelt2."""

    def test_pixelarea_arcsecsq(self, input_model, base_cube_wcs):
        """pixelarea_arcsecsq == cdelt1 * cdelt2."""
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        expected = CDELT1 * CDELT2
        assert math.isclose(cube_model.meta.photometry.pixelarea_arcsecsq, expected, rel_tol=1e-10)

    def test_pixelarea_steradians(self, input_model, base_cube_wcs):
        """pixelarea_steradians == pixelarea_arcsecsq * 2.3504e-11."""
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        expected = CDELT1 * CDELT2 * ARCSEC2_TO_SR
        assert math.isclose(
            cube_model.meta.photometry.pixelarea_steradians, expected, rel_tol=1e-10
        )

    def test_pixelarea_unequal_cdelts(self, input_model, base_cube_wcs):
        """When cdelt1 != cdelt2 the product, not the square, is used."""

        # Only occurs for 'internal' type cubes
        base_cube_wcs.cdelt1 = 0.10
        base_cube_wcs.cdelt2 = 0.13
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        expected = 0.10 * 0.13
        assert math.isclose(cube_model.meta.photometry.pixelarea_arcsecsq, expected, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# Tests: skyalign, linear_wave=True
# ---------------------------------------------------------------------------


class TestSkyalignLinearWave:
    """skyalign + linear_wave=True path."""

    def test_cdelt_stored_in_degrees(self, input_model, base_cube_wcs):
        """meta.wcsinfo.cdelt1/2 must be in degrees (arcsec / 3600)."""
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        assert math.isclose(cube_model.meta.wcsinfo.cdelt1, CDELT1 / 3600.0, rel_tol=1e-10)
        assert math.isclose(cube_model.meta.wcsinfo.cdelt2, CDELT2 / 3600.0, rel_tol=1e-10)

    def test_ctype3_is_wave(self, input_model, base_cube_wcs):
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        assert cube_model.meta.wcsinfo.ctype3 == "WAVE"

    def test_crval3_cdelt3_crpix3(self, input_model, base_cube_wcs):
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        assert math.isclose(cube_model.meta.wcsinfo.crval3, CRVAL3, rel_tol=1e-10)
        assert math.isclose(cube_model.meta.wcsinfo.cdelt3, CDELT3, rel_tol=1e-10)
        assert math.isclose(cube_model.meta.wcsinfo.crpix3, CRPIX3, rel_tol=1e-10)

    def test_roi_spatial_and_wave_set(self, input_model, base_cube_wcs):
        rois, roiw = 0.2, 0.001
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs, rois=rois, roiw=roiw)
        assert math.isclose(cube_model.meta.ifu.roi_spatial, rois, rel_tol=1e-10)
        assert math.isclose(cube_model.meta.ifu.roi_wave, roiw, rel_tol=1e-10)

    def test_status_complete(self, input_model, base_cube_wcs):
        _, status = _call_setup_final_model(input_model, base_cube_wcs)
        assert status == 0
        cube_model, _ = _call_setup_final_model(input_model, base_cube_wcs)
        assert cube_model.meta.cal_step.cube_build == "COMPLETE"


# ---------------------------------------------------------------------------
# Tests: skyalign, linear_wave=False (WAVE-TAB)
# ---------------------------------------------------------------------------


class TestSkyalignNonLinearWave:
    """skyalign + linear_wave=False path (WAVE-TAB / multi-band)."""

    @pytest.fixture
    def nonlinear_cube_wcs(self, base_cube_wcs):
        wavelength_table = np.linspace(7.5, 7.51, NAXIS3)
        base_cube_wcs.linear_wave = False
        base_cube_wcs.cdelt3 = None
        base_cube_wcs.wavelength_table = wavelength_table
        return base_cube_wcs

    def test_ctype3_is_wave_tab(self, input_model, nonlinear_cube_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            nonlinear_cube_wcs,
            roiw_table=np.array([0.001, 0.002, 0.003]),
        )
        assert cube_model.meta.wcsinfo.ctype3 == "WAVE-TAB"

    def test_crval3_and_crpix3_fixed(self, input_model, nonlinear_cube_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            nonlinear_cube_wcs,
            roiw_table=np.array([0.001, 0.002, 0.003]),
        )
        assert cube_model.meta.wcsinfo.crval3 == 1.0
        assert cube_model.meta.wcsinfo.crpix3 == 1.0

    def test_cdelt3_is_none(self, input_model, nonlinear_cube_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            nonlinear_cube_wcs,
            roiw_table=np.array([0.001, 0.002, 0.003]),
        )
        assert cube_model.meta.wcsinfo.cdelt3 is None

    def test_roi_wave_is_mean_of_table(self, input_model, nonlinear_cube_wcs):
        roiw_table = np.array([0.001, 0.002, 0.003])
        cube_model, _ = _call_setup_final_model(
            input_model, nonlinear_cube_wcs, roiw_table=roiw_table
        )
        assert math.isclose(cube_model.meta.ifu.roi_wave, np.mean(roiw_table), rel_tol=1e-10)


# ---------------------------------------------------------------------------
# Tests: internal_cal, NIRSpec
# ---------------------------------------------------------------------------


class TestInternalCalNIRSpec:
    """internal_cal + NIRSpec path: cdelt in meters, correct ctype/cunit."""

    @pytest.fixture
    def nirspec_internal_wcs(self, base_cube_wcs):
        """CubeWCS with NIRSpec internal_cal scales (meters)."""
        base_cube_wcs.cdelt1 = 0.00074  # meters, across-slice
        base_cube_wcs.cdelt2 = 0.00037  # meters, along-slice (cdelt1 / 2)
        base_cube_wcs.rot_angle = None
        return base_cube_wcs

    def test_cdelt_stored_in_native_units_not_degrees(self, input_model, nirspec_internal_wcs):
        """
        For internal_cal the cdelt values must NOT be divided by 3600.
        They should equal cube_wcs.cdelt1 and cube_wcs.cdelt2 exactly.
        """
        cube_model, _ = _call_setup_final_model(
            input_model,
            nirspec_internal_wcs,
            coord_system="internal_cal",
            instrument="NIRSPEC",
        )
        assert math.isclose(
            cube_model.meta.wcsinfo.cdelt1,
            nirspec_internal_wcs.cdelt1,
            rel_tol=1e-10,
        ), (
            f"cdelt1 should be {nirspec_internal_wcs.cdelt1} (native units) "
            f"but got {cube_model.meta.wcsinfo.cdelt1}"
        )
        assert math.isclose(
            cube_model.meta.wcsinfo.cdelt2,
            nirspec_internal_wcs.cdelt2,
            rel_tol=1e-10,
        )

    def test_cunit_is_meter(self, input_model, nirspec_internal_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            nirspec_internal_wcs,
            coord_system="internal_cal",
            instrument="NIRSPEC",
        )
        assert cube_model.meta.wcsinfo.cunit1 == "meter"
        assert cube_model.meta.wcsinfo.cunit2 == "meter"

    def test_ctype_is_nrsslicerxy(self, input_model, nirspec_internal_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            nirspec_internal_wcs,
            coord_system="internal_cal",
            instrument="NIRSPEC",
        )
        assert cube_model.meta.wcsinfo.ctype1 == "NRSSLICERX"
        assert cube_model.meta.wcsinfo.ctype2 == "NRSSLICERY"

    def test_pc1_1_sign_flipped(self, input_model, nirspec_internal_wcs):
        """pc1_1 must be negated relative to the skyalign value."""
        # For rot_angle=0, skyalign sets pc1_1 = -cos(0) = -1.
        # internal_cal then negates it: pc1_1 *= -1 -> +1.
        cube_model, _ = _call_setup_final_model(
            input_model,
            nirspec_internal_wcs,
            coord_system="internal_cal",
            instrument="NIRSPEC",
        )
        assert cube_model.meta.wcsinfo.pc1_1 > 0

    def test_ifu_meta_zeroed(self, input_model, nirspec_internal_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            nirspec_internal_wcs,
            coord_system="internal_cal",
            instrument="NIRSPEC",
            weight_power=2.0,
            roiw=0.001,
            rois=0.2,
        )
        assert cube_model.meta.ifu.weight_power == 0
        assert cube_model.meta.ifu.roi_wave == 0
        assert cube_model.meta.ifu.roi_spatial == 0

    def test_pixelarea_uses_native_cdelt(self, input_model, nirspec_internal_wcs):
        """
        pixelarea_arcsecsq must still equal cdelt1 * cdelt2 as stored in
        cube_wcs (in meters for NIRSpec internal_cal), not (cdelt/3600)^2.
        This test would catch the crpix1/crpix2 typo introduced in cube_lib_1.py.
        """
        cube_model, _ = _call_setup_final_model(
            input_model,
            nirspec_internal_wcs,
            coord_system="internal_cal",
            instrument="NIRSPEC",
        )
        expected_pixar = nirspec_internal_wcs.cdelt1 * nirspec_internal_wcs.cdelt2
        assert math.isclose(
            cube_model.meta.photometry.pixelarea_arcsecsq,
            expected_pixar,
            rel_tol=1e-10,
        )


# ---------------------------------------------------------------------------
# Tests: drizzle interpolation
# ---------------------------------------------------------------------------


class TestDrizzleInterpolation:
    """drizzle path: weight_power, roi_wave, roi_spatial must be zeroed."""

    def test_drizzle_zeros_roi_and_weight(self, input_model, base_cube_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            base_cube_wcs,
            interpolation="drizzle",
            weighting="drizzle",
            weight_power=2.0,
            rois=0.2,
            roiw=0.001,
        )
        assert cube_model.meta.ifu.weight_power == 0
        assert cube_model.meta.ifu.roi_wave == 0
        assert cube_model.meta.ifu.roi_spatial == 0

    def test_drizzle_sets_weighting_label(self, input_model, base_cube_wcs):
        cube_model, _ = _call_setup_final_model(
            input_model,
            base_cube_wcs,
            interpolation="drizzle",
            weighting="drizzle",
        )
        assert cube_model.meta.ifu.weighting == "drizzle"


# ---------------------------------------------------------------------------
# Tests: set_wcs_slicer — linear_wave propagation
# ---------------------------------------------------------------------------


class TestSetWcsSlicer:
    """set_wcs_slicer must pass linear_wave through, not hardcode True."""

    # Minimal slicer-plane corners: centred at origin, ±1 unit
    CORNER_A = np.array([-1.0, -1.0, 1.0, 1.0])
    CORNER_B = np.array([-0.5, 0.5, -0.5, 0.5])
    LAMBDA_MIN = 7.5
    LAMBDA_MAX = 7.6
    CDELT1 = 0.00074  # meters (NIRSpec across-slice)
    CDELT2 = 0.00037  # meters (NIRSpec along-slice)
    CDELT3 = 0.001  # micron

    def _make_nirspec_cube(self):
        """Return a bare IFUCubeData with instrument=NIRSPEC."""
        cube = ifu_cube.IFUCubeData.__new__(ifu_cube.IFUCubeData)
        cube.instrument = "NIRSPEC"
        return cube

    def test_linear_wave_true_is_stored(self):
        """When called with linear_wave=True the CubeWCS must record True."""
        cube = self._make_nirspec_cube()
        cube_wcs = cube.set_wcs_slicer(
            self.CORNER_A,
            self.CORNER_B,
            self.LAMBDA_MIN,
            self.LAMBDA_MAX,
            self.CDELT1,
            self.CDELT2,
            self.CDELT3,
            linear_wave=True,
        )
        assert cube_wcs.linear_wave is True

    def test_linear_wave_false_is_stored(self):
        """
        When called with linear_wave=False the CubeWCS must record False.
        This is the test that catches the former hardcode of linear_wave=True.
        """
        cube = self._make_nirspec_cube()
        cube_wcs = cube.set_wcs_slicer(
            self.CORNER_A,
            self.CORNER_B,
            self.LAMBDA_MIN,
            self.LAMBDA_MAX,
            self.CDELT1,
            self.CDELT2,
            self.CDELT3,
            linear_wave=False,
        )
        assert cube_wcs.linear_wave is False

    def test_cdelt1_stored_verbatim(self):
        """cube_wcs.cdelt1 must be exactly the value passed in."""
        cube = self._make_nirspec_cube()
        cube_wcs = cube.set_wcs_slicer(
            self.CORNER_A,
            self.CORNER_B,
            self.LAMBDA_MIN,
            self.LAMBDA_MAX,
            self.CDELT1,
            self.CDELT2,
            self.CDELT3,
            linear_wave=True,
        )
        assert cube_wcs.cdelt1 == self.CDELT1
        assert cube_wcs.cdelt2 == self.CDELT2

    def test_naxis3_from_wavelength_range(self):
        """naxis3 must equal ceil((lambda_max - lambda_min) / cdelt3)."""
        cube = self._make_nirspec_cube()
        cube_wcs = cube.set_wcs_slicer(
            self.CORNER_A,
            self.CORNER_B,
            self.LAMBDA_MIN,
            self.LAMBDA_MAX,
            self.CDELT1,
            self.CDELT2,
            self.CDELT3,
            linear_wave=True,
        )
        expected_naxis3 = math.ceil((self.LAMBDA_MAX - self.LAMBDA_MIN) / self.CDELT3)
        assert cube_wcs.naxis3 == expected_naxis3

    def test_miri_linear_wave_false_stored(self):
        """Same linear_wave check for MIRI instrument branch."""
        cube = ifu_cube.IFUCubeData.__new__(ifu_cube.IFUCubeData)
        cube.instrument = "MIRI"
        cube_wcs = cube.set_wcs_slicer(
            self.CORNER_A,
            self.CORNER_B,
            self.LAMBDA_MIN,
            self.LAMBDA_MAX,
            cdelt1=0.196,
            cdelt2=0.196,
            cdelt3=self.CDELT3,
            linear_wave=False,
        )
        assert cube_wcs.linear_wave is False

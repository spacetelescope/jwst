"""Unit tests for cube_lib functions."""

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import dqflags

from jwst.cube_build import cube_lib

DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
NON_SCIENCE = dqflags.pixel["NON_SCIENCE"]
SATURATED = dqflags.pixel["SATURATED"]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_cube_wcs_kwargs(
    naxis1=5,
    naxis2=5,
    naxis3=10,
    cdelt1=0.13,
    cdelt2=0.13,
    cdelt3=0.001,
    crval1=45.0,
    crval2=-30.0,
    lambda_min=7.5,
    lambda_max=None,
    linear_wave=True,
    wavelength_table=None,
):
    """Return keyword arguments suitable for cube_lib.set_wcs."""
    if lambda_max is None:
        lambda_max = lambda_min + naxis3 * cdelt3
    return dict(
        corner_ra=[crval1 - 0.01, crval1 + 0.01, crval1 - 0.01, crval1 + 0.01],
        corner_dec=[crval2 - 0.01, crval2 - 0.01, crval2 + 0.01, crval2 + 0.01],
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        ra_center=None,
        dec_center=None,
        rot_angle=None,
        nspax_x=None,
        nspax_y=None,
        cdelt1=cdelt1,
        cdelt2=cdelt2,
        cdelt3=cdelt3,
        linear_wave=linear_wave,
        wavelength_table=wavelength_table,
    )


# ---------------------------------------------------------------------------
# check_cube_size
# ---------------------------------------------------------------------------


class TestCheckCubeSize:
    """
    check_cube_size raises only when total greater than limit.

    Notes:
    tests for check_cube_size
        passes silently when total is within limit
        raises ValueError when total exceeds limit
        exact limit boundary (total == limit) passes
    """

    def test_within_limit_no_exception(self):
        """Total less than  limit must not raise."""
        cube_lib.check_cube_size(10, 10, 10, limit=2000)

    def test_at_exact_limit_no_exception(self):
        """Total equal limit must not raise (boundary passes)."""
        cube_lib.check_cube_size(10, 10, 10, limit=1000)

    def test_exceeds_limit_raises(self):
        """Total greater than  limit must raise ValueError with informative message."""
        with pytest.raises(ValueError, match="exceeds limit"):
            cube_lib.check_cube_size(10, 10, 10, limit=999)


# ---------------------------------------------------------------------------
# radec2std and rot_angle branch
# ---------------------------------------------------------------------------
# Test for radec2std and std2radec:
#  rot_angle != None rotates xi/eta by the expected angle
#  rot_angle=0 is equivalent to rot_angle=None
#  round-trip ra,dec -> xi,eta -> ra,dec is exact (no rot_angle)


class TestRadec2StdRotAngle:
    """Verify the rotation in radec2std works correctly."""

    def test_rot_angle_zero_equals_no_rot_angle(self):
        """rot_angle equal to 0 should produce the same output as rot_angle=None."""
        crval1, crval2 = 45.0, 30.0
        ra, dec = 45.1, 30.1

        xi_none, eta_none = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=None)
        xi_zero, eta_zero = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=0.0)

        np.testing.assert_allclose(xi_none, xi_zero, atol=1e-10)
        np.testing.assert_allclose(eta_none, eta_zero, atol=1e-10)

    def test_rot_angle_90_swaps_xi_eta(self):
        """
        A 90-degree rotation should swap xi and eta (up to sign).
        So xi_rot = eta_unrot, eta_rot = -xi_unrot.
        """
        crval1, crval2 = 0.0, 0.0
        ra, dec = 0.01, 0.0  # small offset purely in RA

        xi0, eta0 = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=None)
        xi90, eta90 = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=90.0)

        # After 90 degrees rotation: xi_rot ≈ eta0, eta_rot ≈ -xi0
        np.testing.assert_allclose(xi90, eta0, atol=1e-6)
        np.testing.assert_allclose(eta90, -xi0, atol=1e-6)

    def test_rot_angle_180_negates_both(self):
        """A 180-degree rotation should negate both xi and eta."""
        crval1, crval2 = 30.0, 20.0
        ra, dec = 30.05, 20.05

        xi0, eta0 = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=None)
        xi180, eta180 = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=180.0)

        np.testing.assert_allclose(xi180, -xi0, atol=1e-6)
        np.testing.assert_allclose(eta180, -eta0, atol=1e-6)

    def test_rot_angle_preserves_total_magnitude(self):
        """Rotation must not change the Euclidean distance from the origin."""
        crval1, crval2 = 60.0, -10.0
        ra, dec = 60.05, -9.95

        xi0, eta0 = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=None)
        magnitude0 = np.sqrt(xi0**2 + eta0**2)

        for angle in [30.0, 45.0, 90.0, 135.0, 270.0]:
            xi_r, eta_r = cube_lib.radec2std(crval1, crval2, ra, dec, rot_angle=angle)
            magnitude_r = np.sqrt(xi_r**2 + eta_r**2)
            np.testing.assert_allclose(
                magnitude_r,
                magnitude0,
                atol=1e-6,
                err_msg=f"Magnitude changed at rot_angle={angle}",
            )


# ---------------------------------------------------------------------------
# set_wcs: ra_center and dec_center override
# ---------------------------------------------------------------------------


class TestSetWcsCenterOverride:
    """User-supplied ra_center and dec_center must replace the computed mean."""

    def test_both_centers_override_computed_values(self):
        kw = _make_cube_wcs_kwargs()
        kw["ra_center"] = 12.34
        kw["dec_center"] = -56.78
        wcs = cube_lib.set_wcs(**kw)
        assert np.isclose(wcs.crval1, 12.34)
        assert np.isclose(wcs.crval2, -56.78)


# ---------------------------------------------------------------------------
# set_wcs: nspax_x and nspax_y override
# ---------------------------------------------------------------------------


class TestSetWcsNspax:
    """User-supplied nspax_x and nspax_y must control naxis1, naxis2 and crpix."""

    def test_wcs_nspax(self):
        kw = _make_cube_wcs_kwargs()
        kw["nspax_x"] = 11  # odd → half = 5, naxis1 = 5*2+1 = 11
        kw["nspax_y"] = 7
        wcs = cube_lib.set_wcs(**kw)
        assert wcs.naxis1 == 11
        assert wcs.naxis2 == 7
        assert np.isclose(wcs.crpix1, (wcs.naxis1 + 1) / 2)
        assert np.isclose(wcs.crpix2, (wcs.naxis2 + 1) / 2)


# ---------------------------------------------------------------------------
# set_wcs and non-linear wavelength path
# ---------------------------------------------------------------------------


class TestSetWcsNonLinearWave:
    """Non-linear wavelength path derives naxis3 from the table, not cdelt3."""

    @pytest.fixture
    def wave_table(self):
        return np.array([7.5, 7.6, 7.75, 8.0, 8.3, 8.7])

    def test_wcs(self, wave_table):
        kw = _make_cube_wcs_kwargs(linear_wave=False, wavelength_table=wave_table)
        wcs = cube_lib.set_wcs(**kw)
        assert wcs.naxis3 == len(wave_table)
        assert np.isclose(wcs.crval3, wave_table[0])
        assert wcs.crpix3 == 1.0
        assert wcs.linear_wave is False
        np.testing.assert_allclose(wcs.zcoord, wave_table)


# ---------------------------------------------------------------------------
# set_wcs and coordinate array lengths and cdelt3_normal
# ---------------------------------------------------------------------------


class TestSetWcsArrayLengths:
    """xcoord, ycoord, zcoord must have exactly naxis1, naxis2, naxis3 elements."""

    def test_xycoord_length_equals_naxis(self):
        kw = _make_cube_wcs_kwargs()
        wcs = cube_lib.set_wcs(**kw)
        assert len(wcs.xcoord) == wcs.naxis1
        assert len(wcs.ycoord) == wcs.naxis2

    def test_zcoord_length_equals_naxis3_linear(self):
        cdelt3 = 0.0001
        kw = _make_cube_wcs_kwargs(lambda_min=7.5, lambda_max=7.51, cdelt3=cdelt3)
        wcs = cube_lib.set_wcs(**kw)
        assert len(wcs.zcoord) == wcs.naxis3
        assert np.isclose(wcs.cdelt3_normal[-1], wcs.cdelt3_normal[-2])
        np.testing.assert_allclose(wcs.cdelt3_normal[:-1], cdelt3)


# ---------------------------------------------------------------------------
# find_spaxel_flux
# ---------------------------------------------------------------------------
# test for find_spaxel_flux:
#  spaxels with iflux == 0 are not modified (division-by-zero guard)
#  spaxels with iflux > 0 are normalized by weight
#  variance is normalized by weight**2
#  mixed array: zero-iflux spaxel unchanged, non-zero normalized


class TestFindSpaxelFlux:
    """Normailization logic in find_spaxel_flux."""

    def test_flux_normalized_by_weight(self):
        """flux[good] = flux[good] / weight[good]."""
        iflux = np.array([1.0])
        flux = np.array([500.0])
        weight = np.array([2.0])
        var = np.array([0.0])

        flux_out, _ = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        assert flux_out[0] == 250.0

    def test_variance_normalized_by_weight_squared(self):
        """var[good] = var[good] / weight[good]**2."""
        iflux = np.array([1.0])
        flux = np.array([0.0])
        weight = np.array([4.0])
        var = np.array([48.0])

        _, var_out = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        assert var_out[0] == 48.0 / 16.0

    def test_mixed_array_only_good_spaxels_modified(self):
        """Zero-iflux spaxels stay unchanged; non-zero-iflux spaxels are normalized."""
        iflux = np.array([0.0, 1.0, 1.0])
        flux = np.array([777.0, 100.0, 200.0])
        weight = np.array([0.0, 5.0, 10.0])
        var = np.array([555.0, 25.0, 400.0])

        flux_out, var_out = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        np.testing.assert_allclose(flux_out[1], 100.0 / 5.0)
        np.testing.assert_allclose(var_out[1], 25.0 / 5.0**2)
        np.testing.assert_allclose(flux_out[2], 200.0 / 10.0)
        np.testing.assert_allclose(var_out[2], 400.0 / 10.0**2)

    def test_unit_weight_leaves_flux_unchanged(self):
        """When weight == 1 and iflux > 0, flux and var are unscaled."""
        iflux = np.array([1.0, 1.0])
        flux = np.array([3.14, 2.72])
        weight = np.ones(2)
        var = np.array([1.0, 4.0])

        flux_out, var_out = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        np.testing.assert_allclose(flux_out, flux)
        np.testing.assert_allclose(var_out, var)


# ---------------------------------------------------------------------------
# set_final_dq_flags
# ---------------------------------------------------------------------------
# Test for set_final_dq_flags:
#  spaxels with weight > 0 keep their flux and are not flagged DO_NOT_USE
#  spaxels with weight == 0 get DO_NOT_USE | NON_SCIENCE added to dq
#  spaxels with weight == 0 have flux and var set to NaN
#  pre-existing DQ bits on zero-weight spaxels are preserved (OR semantics)
#  mixed array: only zero-weight spaxels are flagged


class TestSetFinalDqFlags:
    """DQ flag, NaN, and preservation behavior."""

    def test_positive_weight_spaxel_flux_unchanged(self):
        """Spaxels with weight > 0 must keep their flux value."""
        dq = np.array([0], dtype=np.uint32)
        flux = np.array([42.0])
        var = np.array([1.0])
        weight = np.array([1.5])
        flux_out, var_out, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert flux_out[0] == 42.0
        assert not np.isnan(var_out[0])

    def test_positive_weight_spaxel_not_flagged(self):
        """Spaxels with weight greater than 0 must NOT receive DO_NOT_USE or NON_SCIENCE."""
        dq = np.array([0], dtype=np.uint32)
        flux = np.array([1.0])
        var = np.array([1.0])
        weight = np.array([2.0])

        _, _, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert not (dq_out[0] & DO_NOT_USE)
        assert not (dq_out[0] & NON_SCIENCE)

    def test_zero_weight_gets_do_not_use(self):
        """Spaxels with weight equal 0 must receive DO_NOT_USE and NON_SCIENCE flag."""
        dq = np.array([0], dtype=np.uint32)
        flux = np.array([0.0])
        var = np.array([0.0])
        weight = np.array([0.0])

        flux_out, var_out, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert dq_out[0] == (DO_NOT_USE | NON_SCIENCE)
        assert np.isnan(flux_out[0])
        assert np.isnan(var_out[0])

    def test_pre_existing_dq_preserved_on_zero_weight(self):
        """
        Pre-existing DQ bits (e.g. SATURATED) on a zero-weight spaxel
        must still be present in the output.
        """
        dq = np.array([SATURATED], dtype=np.uint32)
        flux = np.array([0.0])
        var = np.array([0.0])
        weight = np.array([0.0])

        _, _, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert dq_out[0] & SATURATED
        assert dq_out[0] & DO_NOT_USE
        assert dq_out[0] & NON_SCIENCE

    def test_mixed_array_only_zero_weight_flagged(self):
        """
        In a mixed array, only the zero-weight spaxel gets DO_NOT_USE and NON_SCIENCE
        and NaN flux; the good-weight spaxel is untouched.
        """
        dq = np.array([0, 0], dtype=np.uint32)
        flux = np.array([5.0, 0.0])
        var = np.array([1.0, 0.0])
        weight = np.array([2.0, 0.0])

        flux_out, var_out, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        # good spaxel (index 0)
        assert flux_out[0] == 5.0
        assert not np.isnan(var_out[0])
        assert not (dq_out[0] & DO_NOT_USE)

        # bad spaxel (index 1)
        assert np.isnan(flux_out[1])
        assert np.isnan(var_out[1])
        assert dq_out[1] & DO_NOT_USE
        assert dq_out[1] & NON_SCIENCE

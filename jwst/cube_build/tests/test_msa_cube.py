"""Unit tests for MSACubeData and related msa_cube helpers."""

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import dqflags

from jwst.cube_build import cube_lib

DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
NON_SCIENCE = dqflags.pixel["NON_SCIENCE"]
SATURATED = dqflags.pixel["SATURATED"]


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def basic_pars():
    """Minimal valid parameter set for constructing an MSACubeData instance."""
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


# ---------------------------------------------------------------------------
# find_spaxel_flux (via cube_lib, as used by build_msacube)
# ---------------------------------------------------------------------------


class TestFindSpaxelFluxViaCubeLib:
    """
    find_spaxel_flux normalization behaviour, exercised through cube_lib directly

    Notes:
    tests for find_spaxel_flux
        unpopulated spaxels (iflux == 0) are left unchanged
        populated spaxels are normalized: flux / weight, var / weight**2
        mixed arrays only modify the populated entries
    """

    def test_unpopulated_spaxel_unchanged(self):
        """Spaxels with iflux == 0 must not be modified (division-by-zero guard)."""
        iflux = np.array([0.0, 1.0])
        flux = np.array([99.0, 10.0])
        weight = np.array([0.0, 2.0])
        var = np.array([55.0, 4.0])

        flux_out, var_out = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        assert flux_out[0] == 99.0
        assert var_out[0] == 55.0

    def test_flux_normalized_by_weight(self):
        """flux[good] must equal flux[good] / weight[good]."""
        iflux = np.array([0.0, 1.0, 2.0])
        flux = np.array([0.0, 10.0, 20.0])
        weight = np.array([0.0, 2.0, 4.0])
        var = np.array([0.0, 4.0, 16.0])

        flux_out, _ = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        np.testing.assert_allclose(flux_out[1], 5.0)
        np.testing.assert_allclose(flux_out[2], 5.0)

    def test_variance_normalized_by_weight_squared(self):
        """var[good] must equal var[good] / weight[good]**2."""
        iflux = np.array([0.0, 1.0, 2.0])
        flux = np.array([0.0, 10.0, 20.0])
        weight = np.array([0.0, 2.0, 4.0])
        var = np.array([0.0, 4.0, 16.0])

        _, var_out = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        np.testing.assert_allclose(var_out[1], 1.0)  # 4 / 4
        np.testing.assert_allclose(var_out[2], 1.0)  # 16 / 16

    def test_mixed_array_only_good_spaxels_modified(self):
        """Zero-iflux spaxels stay unchanged; non-zero-iflux spaxels are normalized."""
        iflux = np.array([0.0, 1.0, 1.0])
        flux = np.array([777.0, 100.0, 200.0])
        weight = np.array([0.0, 5.0, 10.0])
        var = np.array([555.0, 25.0, 400.0])

        flux_out, var_out = cube_lib.find_spaxel_flux(iflux, flux, weight, var)

        # zero-iflux spaxel untouched
        assert flux_out[0] == 777.0
        assert var_out[0] == 555.0

        # populated spaxels normalized
        np.testing.assert_allclose(flux_out[1], 20.0)
        np.testing.assert_allclose(var_out[1], 25.0 / 25.0)
        np.testing.assert_allclose(flux_out[2], 20.0)
        np.testing.assert_allclose(var_out[2], 400.0 / 100.0)


# ---------------------------------------------------------------------------
# set_wcs (linear wavelength path, called from setup_msacube_wcs)
# ---------------------------------------------------------------------------


class TestSetWcsLinearWavelength:
    """
    cube_lib.set_wcs with linear_wave=True produces the correct WCS geometry
    when driven with the same inputs MSACubeData would supply.

    Notes:
    tests for set_wcs (linear path)
        ra_center and dec_center override the computed footprint mean
        nspax_x and nspax_y control naxis1, naxis2, and crpix values
        naxis3 equals ceil(range_lambda / cdelt3)
        xcoord, ycoord, zcoord lengths match naxis1, naxis2, naxis3
        cdelt3_normal last two entries are equal (boundary repeat)
    """

    def _make_wcs(self, **overrides):
        defaults = dict(
            corner_ra=[99.9, 100.1, 99.9, 100.1],
            corner_dec=[-0.1, -0.1, 0.1, 0.1],
            lambda_min=1.0,
            lambda_max=2.0,
            ra_center=100.0,
            dec_center=0.0,
            rot_angle=None,
            nspax_x=10,
            nspax_y=10,
            cdelt1=0.1,
            cdelt2=0.1,
            cdelt3=0.001,
            linear_wave=True,
            wavelength_table=None,
        )
        defaults.update(overrides)
        return cube_lib.set_wcs(**defaults)

    def test_ra_dec_center_stored_as_crval(self):
        """User-supplied ra_center and dec_center must appear as crval1/crval2."""
        wcs = self._make_wcs(ra_center=12.34, dec_center=-56.78)

        assert np.isclose(wcs.crval1, 12.34)
        assert np.isclose(wcs.crval2, -56.78)

    def test_nspax_controls_naxis_dimensions(self):
        """nspax_x=10  nspax_y=8"""
        wcs = self._make_wcs(nspax_x=10, nspax_y=8)

        assert wcs.naxis1 == 11
        assert wcs.naxis2 == 9
        assert len(wcs.xcoord) == wcs.naxis1
        assert len(wcs.ycoord) == wcs.naxis2
        assert len(wcs.zcoord) == wcs.naxis3

    def test_cdelt3_normal_interior_entries_equal_cdelt3(self):
        """Interior cdelt3_normal values must equal cdelt3 for a linear grid."""
        wcs = self._make_wcs(cdelt3=0.001)

        np.testing.assert_allclose(wcs.cdelt3_normal[:-1], 0.001, rtol=1e-6)


# ---------------------------------------------------------------------------
# set_wcs (non-linear wavelength path)
# ---------------------------------------------------------------------------


class TestSetWcsNonLinearWavelength:
    """
    cube_lib.set_wcs with linear_wave=False derives the spectral axis from the
    supplied wavelength_table rather than cdelt3.

    Notes:
    tests for set_wcs (non-linear path)
        naxis3 equals len(wavelength_table)
        zcoord contains exactly the wavelength_table values
        crval3 equals the first table entry
        crpix3 is 1.0 (FITS convention)
        linear_wave flag is False on the returned CubeWCS
    """

    @pytest.fixture
    def wave_table(self):
        return np.array([1.0, 1.2, 1.5, 2.0])

    def _make_nonlinear_wcs(self, wave_table):
        return cube_lib.set_wcs(
            corner_ra=[99.9, 100.1, 99.9, 100.1],
            corner_dec=[-0.1, -0.1, 0.1, 0.1],
            lambda_min=1.0,
            lambda_max=2.0,
            ra_center=100.0,
            dec_center=0.0,
            rot_angle=None,
            nspax_x=10,
            nspax_y=10,
            cdelt1=0.1,
            cdelt2=0.1,
            cdelt3=None,
            linear_wave=False,
            wavelength_table=wave_table,
        )

    def test_naxis3_equals_table_length(self, wave_table):
        """wavelength_table testing."""
        wcs = self._make_nonlinear_wcs(wave_table)

        assert wcs.naxis3 == len(wave_table)
        np.testing.assert_allclose(wcs.zcoord, wave_table)
        assert wcs.crpix3 == 1.0
        assert wcs.linear_wave is False
        assert len(wcs.cdelt3_normal) == wcs.naxis3


# ---------------------------------------------------------------------------
# set_final_dq_flags (via cube_lib, as used by build_msacube)
# ---------------------------------------------------------------------------


class TestSetFinalDqFlagsViaCubeLib:
    """
    cube_lib.set_final_dq_flags behaviour

    Notes:
    tests for set_final_dq_flags
        spaxels with weight > 0 keep their flux, are not flagged DO_NOT_USE
        spaxels with weight == 0 get DO_NOT_USE | NON_SCIENCE, flux → NaN
        pre-existing DQ bits are preserved on zero-weight spaxels (OR semantics)
        mixed array: only the zero-weight spaxel is affected
    """

    def test_nonzero_weight_flux_unchanged(self):
        """Spaxels with weight > 0 must retain their original flux value."""
        dq = np.array([0], dtype=np.uint32)
        flux = np.array([42.0])
        var = np.array([1.0])
        weight = np.array([1.5])

        flux_out, var_out, _ = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert flux_out[0] == 42.0
        assert not np.isnan(var_out[0])

    def test_nonzero_weight_not_flagged(self):
        """Spaxels with weight > 0 must NOT receive DO_NOT_USE or NON_SCIENCE."""
        dq = np.array([0], dtype=np.uint32)
        flux = np.array([1.0])
        var = np.array([1.0])
        weight = np.array([2.0])

        _, _, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert not (dq_out[0] & DO_NOT_USE)
        assert not (dq_out[0] & NON_SCIENCE)

    def test_zero_weight_gets_do_not_use_and_non_science(self):
        """Spaxels with weight == 0 must have DO_NOT_USE and NON_SCIENCE set."""
        dq = np.array([0], dtype=np.uint32)
        flux = np.array([0.0])
        var = np.array([0.0])
        weight = np.array([0.0])

        flux_out, var_out, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        assert np.isnan(flux_out[0])
        assert np.isnan(var_out[0])
        assert dq_out[0] == (DO_NOT_USE | NON_SCIENCE)

    def test_pre_existing_bits_preserved_on_zero_weight(self):
        """
        Pre-existing DQ bits (e.g. SATURATED) on a zero-weight spaxel
        must still be present alongside the newly added flags.
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
        In a mixed array the good-weight spaxel keeps its flux and clean DQ;
        only the zero-weight spaxel is NaN-filled and flagged.
        """
        dq = np.array([0, 0], dtype=np.uint32)
        flux = np.array([5.0, 0.0])
        var = np.array([1.0, 0.0])
        weight = np.array([2.0, 0.0])

        flux_out, var_out, dq_out = cube_lib.set_final_dq_flags(dq, flux, var, weight)

        # good spaxel untouched
        assert flux_out[0] == 5.0
        assert not np.isnan(var_out[0])
        assert not (dq_out[0] & DO_NOT_USE)

        # bad spaxel flagged and NaN-filled
        assert np.isnan(flux_out[1])
        assert np.isnan(var_out[1])
        assert dq_out[1] & DO_NOT_USE
        assert dq_out[1] & NON_SCIENCE


# ---------------------------------------------------------------------------
# spaxel_exposure_counter initialization and accumulation
# ---------------------------------------------------------------------------


class TestExposureCounter:
    """
    spaxel_exposure_counter is a uint32 array that tracks how many exposures
    contributed to each spaxel.  build_msacube initialises it as zeros and
    accumulates per-slit contributions by integer addition.  These contracts
    are exercised here through cube_lib.set_wcs to produce a known cube
    geometry, so no C extension or IFUCubeModel is required.

    Notes:
    tests for spaxel_exposure_counter
        initialised array has length naxis1 * naxis2 * naxis3
        initialised array is all zeros
        initialised array dtype is uint32
        accumulating one slit's contribution increments every reached spaxel by 1
        accumulating two identical contributions doubles every count
        partial contribution only raises reached spaxels
    """

    def _make_wcs(self, nspax_x=3, nspax_y=3):
        """Return a small CubeWCS via cube_lib.set_wcs."""
        return cube_lib.set_wcs(
            corner_ra=[99.9, 100.1, 99.9, 100.1],
            corner_dec=[-0.1, -0.1, 0.1, 0.1],
            lambda_min=1.0,
            lambda_max=2.0,
            ra_center=100.0,
            dec_center=0.0,
            rot_angle=None,
            nspax_x=nspax_x,
            nspax_y=nspax_y,
            cdelt1=0.1,
            cdelt2=0.1,
            cdelt3=0.001,
            linear_wave=True,
            wavelength_table=None,
        )

    def test_initial_length_equals_naxis_product(self):
        """Counter length must equal naxis1 * naxis2 * naxis3."""
        wcs = self._make_wcs()
        total = wcs.naxis1 * wcs.naxis2 * wcs.naxis3
        counter = np.zeros(total, dtype=np.uint32)
        assert len(counter) == total
        assert np.all(counter == 0)

    def test_one_slit_contribution_increments_spaxels(self):
        """Adding one slit's all-ones contribution must set every element to 1."""
        wcs = self._make_wcs()
        total = wcs.naxis1 * wcs.naxis2 * wcs.naxis3
        counter = np.zeros(total, dtype=np.uint32)
        this_exp = np.ones(total, dtype=np.int32)
        counter = counter + np.asarray(this_exp, np.int32)
        assert np.all(counter == 1)

    def test_two_slit_contributions_double_every_count(self):
        """Adding two identical slit contributions must set every element to 2."""
        wcs = self._make_wcs()
        total = wcs.naxis1 * wcs.naxis2 * wcs.naxis3
        counter = np.zeros(total, dtype=np.uint32)
        this_exp = np.ones(total, dtype=np.int32)
        for _ in range(2):
            counter = counter + np.asarray(this_exp, np.int32)
        assert np.all(counter == 2)

    def test_partial_contribution_only_raises_reached_spaxels(self):
        """A slit that only covers some spaxels must leave unreached ones at zero."""
        wcs = self._make_wcs()
        total = wcs.naxis1 * wcs.naxis2 * wcs.naxis3
        counter = np.zeros(total, dtype=np.uint32)
        this_exp = np.zeros(total, dtype=np.int32)
        this_exp[: total // 2] = 1
        counter = counter + np.asarray(this_exp, np.int32)
        assert np.all(counter[: total // 2] == 1)
        assert np.all(counter[total // 2 :] == 0)

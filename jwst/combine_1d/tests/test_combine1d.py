"""Test operations in Combine1dStep."""

import logging

import numpy as np
import pytest

from jwst import datamodels
from jwst.combine_1d import Combine1dStep
from jwst.combine_1d.combine1d import InputSpectrumModel, check_exptime, check_monotonic
from jwst.datamodels.utils.tests.wfss_helpers import N_SOURCES, wfss_multi
from jwst.tests.helpers import LogWatcher


@pytest.fixture
def two_spectra():
    spec1 = create_spec_model(flux=1e-9)
    spec2 = create_spec_model(flux=1e-9)
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1.0
    ms.meta.exposure.integration_time = 2.0
    ms.spec.append(spec1)
    ms.spec.append(spec2)
    yield ms
    ms.close()


@pytest.fixture
def three_spectra():
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1
    for _ in range(3):
        spec = create_spec_model(flux=1e-9)
        spec.source_id = 0
        spec.spectral_order = 1
        ms.spec.append(spec)
    yield ms
    ms.close()


def test_dq(two_spectra):
    """Test that DQ exclusion works."""
    ms = two_spectra
    spec1 = ms.spec[0]
    spec2 = ms.spec[1]

    # Make one pixel bad by being large.
    # Result should be large
    bad_pix = 5
    spec1.spec_table["FLUX"][bad_pix] = 1.0
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table["FLUX"][bad_pix], 0.5)

    # Now mark that pixel bad.
    # Result should not just contain the second spectrum value.
    spec1.spec_table["DQ"][bad_pix] = datamodels.dqflags.pixel["DO_NOT_USE"]
    result_dq = Combine1dStep.call(ms)
    assert np.isclose(
        result_dq.spec[0].spec_table["FLUX"][bad_pix], spec2.spec_table["FLUX"][bad_pix]
    )


def test_err():
    """Test error propagation."""
    spec1 = create_spec_model(flux=1.0, error=0.1)
    spec2 = create_spec_model(flux=1.0, error=0.2)
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1
    ms.spec.append(spec1)
    ms.spec.append(spec2)

    result = Combine1dStep.call(ms)

    expected_error = np.sqrt(0.1**2 + 0.2**2) / 2
    assert np.allclose(result.spec[0].spec_table["FLUX"], 1.0)
    assert np.allclose(result.spec[0].spec_table["ERROR"], expected_error)
    assert np.allclose(result.spec[0].spec_table["SURF_BRIGHT"], 1.0)
    assert np.allclose(result.spec[0].spec_table["SB_ERROR"], expected_error)


def test_nan(two_spectra):
    """Test that nan exclusion works."""
    ms = two_spectra
    spec1 = ms.spec[0]
    spec2 = ms.spec[1]

    # Make one pixel bad by being nan, without being marked as DO_NOT_USE.
    # Result should just contain the second spectrum value
    bad_pix = 5
    spec1.spec_table["FLUX"][bad_pix] = np.nan
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table["FLUX"][bad_pix], spec2.spec_table["FLUX"][bad_pix])


def test_sigmaclip(three_spectra):
    """Test that sigma clipping works."""
    ms = three_spectra
    spec1 = ms.spec[0]
    spec2 = ms.spec[1]

    # Make one pixel bad by being large.
    # Result should be large
    bad_pix = 5
    spec1.spec_table["FLUX"][bad_pix] = 1.0
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table["FLUX"][bad_pix], 1 / 3)
    result.close()

    # Now process with sigma_clipping turned on
    # Result should not just contain the second and third spectra values.
    result_sc = Combine1dStep.call(ms, sigma_clip=3)
    assert np.isclose(
        result_sc.spec[0].spec_table["FLUX"][bad_pix], spec2.spec_table["FLUX"][bad_pix]
    )
    result_sc.close()


@pytest.mark.parametrize("casing", ["upper", "lower"])
@pytest.mark.parametrize(
    "exptime",
    ["exposure_time", "integration_time", "unit_weight", "effexptm", "effinttm", "unit weight"],
)
def test_exptime_keys(exptime, casing):
    if casing == "upper":
        exptime = exptime.upper()

    spec1 = create_spec_model(flux=1.0, error=0.1)
    spec2 = create_spec_model(flux=2.0, error=0.2)

    ms1 = datamodels.MultiSpecModel()
    ms1.meta.exposure.exposure_time = 1.0
    ms1.meta.exposure.integration_time = 2.0
    ms1.spec.append(spec1)

    ms2 = datamodels.MultiSpecModel()
    ms2.meta.exposure.exposure_time = 2.0
    ms2.meta.exposure.integration_time = 1.0
    ms2.spec.append(spec2)

    result = Combine1dStep.call(datamodels.ModelContainer([ms1, ms2]), exptime_key=exptime)
    if "exp" in exptime.lower():
        # closer to 2
        assert np.allclose(result.spec[0].spec_table["FLUX"], 1 + 2 / 3)
    elif "int" in exptime.lower():
        # closer to 1
        assert np.allclose(result.spec[0].spec_table["FLUX"], 1 + 1 / 3)
    else:
        # equidistant
        assert np.allclose(result.spec[0].spec_table["FLUX"], 1.5)


def test_bad_exptime(two_spectra):
    # Runtime error if bad key is passed to input model
    with pytest.raises(RuntimeError):
        InputSpectrumModel(two_spectra, two_spectra.spec[0], "bad")

    # Bad key is translated to unit_weight if checked
    assert check_exptime("bad") == "unit_weight"


def create_spec_model(npoints=10, flux=1e-9, error=1e-10, wave_range=(11, 13)):
    """Create a SpecModel."""
    wavelength = np.arange(*wave_range, step=(wave_range[1] - wave_range[0]) / npoints)
    flux = np.full(npoints, flux)
    error = np.full(npoints, error)

    surf_bright = np.full(npoints, flux)
    sb_error = np.full(npoints, error)

    var = np.zeros(npoints)
    dq = np.zeros(npoints)
    background = np.zeros(npoints)
    berror = np.zeros(npoints)
    npixels = np.zeros(npoints)

    # This data type is used for creating an output table.
    spec_dtype = datamodels.SpecModel().get_dtype("spec_table")

    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                var,
                var,
                var,
                surf_bright,
                sb_error,
                var,
                var,
                var,
                dq,
                background,
                berror,
                var,
                var,
                var,
                npixels,
                strict=False,
            ),
        ),
        dtype=spec_dtype,
    )

    spec_model = datamodels.SpecModel(spec_table=otab)

    return spec_model


def create_spec_model_nonmonotonic(npoints=10, flux=1e-9, error=1e-10, wave_range=(5, 12)):
    """Create a SpecModel."""
    wavelength = np.arange(*wave_range, step=(wave_range[1] - wave_range[0]) / npoints)

    # To break monotonicity, swap the last two elements.
    if npoints > 1:
        wavelength[-1], wavelength[-2] = wavelength[-2], wavelength[-1]

    flux = np.full(npoints, flux)
    error = np.full(npoints, error)

    surf_bright = np.full(npoints, flux)
    sb_error = np.full(npoints, error)

    var = np.zeros(npoints)
    dq = np.zeros(npoints)
    background = np.zeros(npoints)
    berror = np.zeros(npoints)
    npixels = np.zeros(npoints)

    # This data type is used for creating an output table.
    spec_dtype = datamodels.SpecModel().get_dtype("spec_table")

    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                var,
                var,
                var,
                surf_bright,
                sb_error,
                var,
                var,
                var,
                dq,
                background,
                berror,
                var,
                var,
                var,
                npixels,
                strict=False,
            ),
        ),
        dtype=spec_dtype,
    )

    spec_model = datamodels.SpecModel(spec_table=otab)

    return spec_model


@pytest.fixture
def wfss_multiexposure():
    return wfss_multi()


def test_wfss_multi_input(wfss_multiexposure):
    """Smoke test to ensure combine_1d works with WFSSMultiSpecModel"""
    result = Combine1dStep.call(wfss_multiexposure)
    assert isinstance(result, datamodels.WFSSMultiCombinedSpecModel)
    assert result is not wfss_multiexposure

    tab = result.spec[0].spec_table
    assert tab.shape == (N_SOURCES,)
    assert result.meta.cal_step.combine_1d == "COMPLETE"
    assert np.allclose(tab["FLUX"], 1.0)

    # check that metadata was passed through correctly
    assert np.allclose(tab["SOURCE_RA"], 0.0)
    assert np.allclose(tab["SOURCE_DEC"], 0.0)
    assert np.all(tab["SOURCE_TYPE"] == "POINT")

    assert result.spec[0].dispersion_direction == 3


def test_allnan_skip(wfss_multiexposure, monkeypatch):
    """Test that all-nan spectra are skipped."""
    # Set all flux values to NaN
    for spec in wfss_multiexposure.spec:
        spec.spec_table["FLUX"][:] = np.nan
        spec.spectral_order = 1

    # message when a single spectrum has no valid flux values
    watcher0 = LogWatcher(
        "Input spectrum 5 order 1 from group_id 1 has no valid flux values; skipping."
    )
    monkeypatch.setattr(logging.getLogger("jwst.combine_1d.combine1d"), "warning", watcher0)
    result = Combine1dStep.call(wfss_multiexposure)
    assert result.meta.cal_step.combine_1d == "SKIPPED"
    watcher0.assert_seen()

    # check for the other log messages
    # these must be done one at a time because the watcher can only be monkeypatched once
    # message when no valid input spectra are found for the source
    watcher1 = LogWatcher("No valid input spectra found for source. Skipping.")
    monkeypatch.setattr(logging.getLogger("jwst.combine_1d.combine1d"), "error", watcher1)
    result = Combine1dStep.call(wfss_multiexposure)
    assert result.meta.cal_step.combine_1d == "SKIPPED"
    watcher1.assert_seen()

    # message when no valid input spectra at all are found
    watcher2 = LogWatcher("No valid input spectra found in WFSSMultiSpecModel")
    monkeypatch.setattr(logging.getLogger("jwst.combine_1d.combine_1d_step"), "error", watcher2)
    result = Combine1dStep.call(wfss_multiexposure)
    assert result.meta.cal_step.combine_1d == "SKIPPED"
    watcher2.assert_seen()


def test_allnan_skip_non_wfss(caplog, three_spectra):
    """Test that all-nan spectra are skipped."""
    # Set all flux values to NaN
    for spec in three_spectra.spec:
        spec.spec_table["FLUX"][:] = np.nan

    result = Combine1dStep.call(three_spectra)
    assert result.meta.cal_step.combine_1d == "SKIPPED"

    # message when a single spectrum has no valid flux values
    msg1 = "Input spectrum 0 order 1 has no valid flux values; skipping."
    assert msg1 in caplog.text

    # message when no valid input spectra are found for the source
    msg2 = "No valid input spectra found for source. Skipping."
    assert msg2 in caplog.text

    # message when no valid input spectra at all are found
    msg3 = "No valid input spectra found in WFSSMultiSpecModel"
    assert msg2 in caplog.text


def test_output_is_not_input(two_spectra):
    """Test that input is not modified."""
    result = Combine1dStep.call(two_spectra)

    # Successful completion
    assert result.meta.cal_step.combine_1d == "COMPLETE"

    # Input is not modified
    assert result is not two_spectra
    assert two_spectra.meta.cal_step.combine_1d is None


def test_combination_error(caplog):
    """Test an error in combination."""
    bad_model = datamodels.ImageModel()
    result = Combine1dStep.call(bad_model)

    # Step is skipped
    assert result.meta.cal_step.combine_1d == "SKIPPED"
    assert result is not bad_model
    assert "Invalid input model" in caplog.text

    # Input is not modified
    assert result is not bad_model
    assert bad_model.meta.cal_step.combine_1d is None


def test_monotonic():
    # test_strictly_increasing
    assert check_monotonic([400, 500, 600, 700]) is True

    # test_strictly_decreasing
    assert check_monotonic([700, 600, 500, 400]) is True

    # test_not_monotonic
    assert check_monotonic([400, 600, 500, 700]) is False

    # test_duplicates_return_false
    # Documentation says duplicates should return False
    assert check_monotonic([400, 500, 500, 600]) is False

    # test_single_element
    assert check_monotonic([500]) is True

    # test_empty_list():
    assert check_monotonic([]) is True

    # test_numpy_array_input
    # Ensuring it handles actual numpy arrays, not just lists
    arr = np.array([1.1, 1.2, 1.3])
    assert check_monotonic(arr) is True


def test_combine1d_wavelength_merging_test():
    """Test that Combine1dStep ignores a spectrum with non-monotonic wavelengths."""

    # 1. Setup Input Spectra: monotonic
    spec1 = create_spec_model(flux=1e-9, wave_range=(5, 12))
    # 2. Set up non-monotonuc wavelength
    spec2 = create_spec_model_nonmonotonic(flux=1e-9, wave_range=(5, 12))

    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1.0
    ms.meta.exposure.integration_time = 2.0
    ms.spec.append(spec1)
    ms.spec.append(spec2)

    # 2. Initialize and Run the Step
    step = Combine1dStep()
    result = step.run(ms)

    # 3. get wavelengths
    combined_wave = result.spec[0].spec_table["wavelength"]
    spec1_wave = spec1.spec_table["wavelength"]

    # 4. Assert they are the same within a small tolerance
    assert np.allclose(combined_wave, spec1_wave, atol=1e-12)

    # 5. Assert monotonic
    assert np.all(np.diff(combined_wave) > 0)

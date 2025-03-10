import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.residual_fringe.residual_fringe import utils, ResidualFringeCorrection
from jwst.residual_fringe.residual_fringe_step import ResidualFringeStep
from jwst.residual_fringe.utils import fit_residual_fringes_1d as rf1d


@pytest.fixture()
def wave():
    x1, x2 = 10.0, 11.75
    wave_array = np.linspace(x1, x2, 1024)
    return wave_array


@pytest.fixture()
def linear_spectrum(wave):
    """Mock a spectrum with a linear continuum."""
    # linear flux signal between min and max wavelengths
    x1, x2 = wave.min(), wave.max()
    y1, y2 = 0.5, 0.4
    flux = (y2 - y1) / (x2 - x1) * (wave - x1) + y1

    return wave, flux


@pytest.fixture()
def fringed_spectrum(linear_spectrum):
    """Mock a spectrum with periodic fringe signature."""
    wave, flux = linear_spectrum

    # add a sinusoid signal on top of the linear flux
    amp = 0.01
    period = 0.04
    fringe = amp * np.sin(2 * np.pi * wave / period)
    return wave, flux + fringe


@pytest.fixture()
def miri_mrs_model_linear(monkeypatch, linear_spectrum):
    shape = (1024, 10)
    model = datamodels.IFUImageModel(shape)
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIFUSHORT"
    model.meta.instrument.channel = "12"
    model.meta.instrument.band = "SHORT"
    model.meta.exposure.type = "MIR_MRS"
    model.meta.observation.date = "2022-05-01"
    model.meta.observation.time = "01:01:01"
    model.meta.cal_step.fringe = "COMPLETE"

    wave, flux = linear_spectrum
    model.data[:, :] = flux[:, None]
    model.wavelength[:, :] = wave[:, None]
    model.err = model.data * 0.01

    return model


@pytest.fixture()
def miri_mrs_model_with_fringe(miri_mrs_model_linear, fringed_spectrum):
    wave, flux = fringed_spectrum
    model_copy = miri_mrs_model_linear.copy()
    model_copy.data[:, :] = flux[:, None]
    return model_copy


@pytest.fixture()
def mock_wavemap(monkeypatch, wave):
    # mock the wavelength function to avoid making a full test WCS
    def return_wavelength(*args):
        wavemap = np.zeros((1024, 10))
        wavemap[:, :] = wave[:, None]
        return wavemap

    monkeypatch.setattr(ResidualFringeCorrection, "_get_wave_map", return_wavelength)


@pytest.fixture()
def mock_wavemap_with_nans(monkeypatch, wave):
    # mock the wavelength function to avoid making a full test WCS
    def return_wavelength(*args):
        wavemap = np.zeros((1024, 10))
        wavemap[:, :] = wave[:, None]

        # add some scattered NaN values
        wavemap[::20, ::2] = np.nan
        return wavemap

    monkeypatch.setattr(ResidualFringeCorrection, "_get_wave_map", return_wavelength)


@pytest.fixture()
def mock_slice_info_short(monkeypatch):
    # mock a single slice to fit matching the test data, for testing speed
    def one_slice(*args):
        slices_in_channel = [101]
        xrange_channel = np.array([[0, 10]])
        slice_x_ranges = np.array([[101, 0, 10]])
        all_slice_masks = np.ones((1, 1024, 10))
        return slices_in_channel, xrange_channel, slice_x_ranges, all_slice_masks

    monkeypatch.setattr(utils, "slice_info", one_slice)


@pytest.fixture()
def mock_slice_info_long(monkeypatch):
    # mock a single slice to fit matching the test data, for testing speed
    def one_slice(*args):
        slices_in_channel = [301]
        xrange_channel = np.array([[0, 10]])
        slice_x_ranges = np.array([[301, 0, 10]])
        all_slice_masks = np.ones((1, 1024, 10))
        return slices_in_channel, xrange_channel, slice_x_ranges, all_slice_masks

    monkeypatch.setattr(utils, "slice_info", one_slice)


def test_rf1d(linear_spectrum, fringed_spectrum):
    """
    Test the performance of the 1d residual defringe routine.

    Input synthetic data mimics a Ch2C spectrum taken from observations
    of bright star 16 CygB.
    """
    wave, flux = fringed_spectrum
    outflux = rf1d(flux, wave, channel=2)

    # defringing won't remove the pure sinusoidal fringe completely, but
    # it should be reasonably close to linear and significantly better
    # than no correction
    expected_flux = linear_spectrum[1]

    # corrected output has small diffs from linear on average
    # (edge effects might be larger)
    relative_diff_output = np.abs(outflux - expected_flux) / expected_flux
    assert np.nanmean(relative_diff_output) < 0.005

    # input diffs from linear are much bigger
    relative_diff_input = np.abs(flux - expected_flux) / expected_flux
    assert np.nanmean(relative_diff_input) > 0.01


def test_get_wavemap():
    """
    Test the _get_wavemap function directly.

    A separate test is needed, since calls to the higher level correction method
    mock this function for synthetic data simplicity.
    """
    model = datamodels.IFUImageModel()

    # Mock a WCS that returns 1 for wavelengths
    def return_ones(x, y):
        return None, None, np.ones(x.shape)

    model.meta.wcs = return_ones

    rf = ResidualFringeCorrection(model, "N/A", "N/A", None)
    wavemap = rf._get_wave_map()
    assert wavemap.shape == model.data.shape
    assert np.all(wavemap == 1.0)


@pytest.mark.parametrize("band", ["SHORT", "MEDIUM", "LONG"])
def test_rf_step_short(
    miri_mrs_model_linear, miri_mrs_model_with_fringe, mock_slice_info_short, mock_wavemap, band
):
    model = miri_mrs_model_with_fringe
    model.meta.instrument.band = band
    result = ResidualFringeStep.call(model, skip=False)

    assert result.meta.cal_step.residual_fringe == "COMPLETE"

    # output should be closer to a linear spectrum than input,
    # correction will not be precise
    expected = miri_mrs_model_linear.data
    relative_diff_input = np.abs(model.data - expected) / expected
    relative_diff_output = np.abs(result.data - expected) / expected
    assert np.nanmean(relative_diff_output) < np.nanmean(relative_diff_input)


@pytest.mark.parametrize("band", ["SHORT", "MEDIUM", "LONG"])
def test_rf_step_long(
    miri_mrs_model_with_fringe, mock_slice_info_long, mock_wavemap, band, log_watcher
):
    model = miri_mrs_model_with_fringe
    model.meta.instrument.detector = "MIRIFULONG"
    model.meta.instrument.channel = "34"
    model.meta.instrument.band = band

    # Synthetic input data is reasonable for MIRIFUSHORT, but is expected
    # to fail with a warning when treated as MIRIFULONG.
    watcher = log_watcher("jwst.residual_fringe.residual_fringe", message="Skipping col")
    result = ResidualFringeStep.call(model, skip=False)
    watcher.assert_seen()

    # Output data should be identical to input, although step is complete
    assert result.meta.cal_step.residual_fringe == "COMPLETE"
    assert np.allclose(model.data, result.data)


def test_rf_step_nans_in_wavelength(
    miri_mrs_model_linear, miri_mrs_model_with_fringe, mock_slice_info_short, mock_wavemap_with_nans
):
    model = miri_mrs_model_with_fringe

    # wavelength array has some scattered NaNs:
    # they should be interpolated over and correction should succeed
    result = ResidualFringeStep.call(model, skip=False)

    # output should be closer to a linear spectrum than input,
    # correction will not be precise
    assert result.meta.cal_step.residual_fringe == "COMPLETE"
    expected = miri_mrs_model_linear.data
    relative_diff_input = np.abs(model.data - expected) / expected
    relative_diff_output = np.abs(result.data - expected) / expected
    assert np.nanmean(relative_diff_output) < np.nanmean(relative_diff_input)


def test_rf_step_save_intermediate(
    tmp_path, miri_mrs_model_with_fringe, mock_slice_info_short, mock_wavemap
):
    model = miri_mrs_model_with_fringe
    model.meta.filename = "test.fits"
    ResidualFringeStep.call(
        model,
        skip=False,
        output_dir=str(tmp_path),
        save_results=True,
        save_intermediate_results=True,
    )

    output_files = [
        "test_residual_fringe.fits",
        "test_stat_table.ecsv",
        "test_out_table.ecsv",
        "test_fit_results.fits",
    ]
    for output_file in output_files:
        assert (tmp_path / output_file).exists()


def test_rf_step_ignore_regions(miri_mrs_model_with_fringe, mock_slice_info_short, mock_wavemap):
    model = miri_mrs_model_with_fringe

    # ignore all the data
    ignore_region_min = [model.wavelength.min()]
    ignore_region_max = [model.wavelength.max()]
    result = ResidualFringeStep.call(
        model, skip=False, ignore_region_min=ignore_region_min, ignore_region_max=ignore_region_max
    )

    # output should be the same as input
    assert np.allclose(model.data, result.data)


def test_rf_step_low_snr(
    miri_mrs_model_with_fringe, mock_slice_info_short, mock_wavemap, log_watcher
):
    model = miri_mrs_model_with_fringe

    # set all the data to a very small value so SNR is too low to fit
    model.data[:] = 1e-6

    watcher = log_watcher("jwst.residual_fringe.residual_fringe", message="SNR too low")
    result = ResidualFringeStep.call(model, skip=False)
    watcher.assert_seen()

    # output should be the same as input
    assert np.allclose(model.data, result.data)


def test_rf_step_weights_gap(miri_mrs_model_with_fringe, mock_slice_info_short, mock_wavemap):
    model = miri_mrs_model_with_fringe

    # set some big emission lines in rows 1 and -2:
    # this triggers a weighting edge case that ignores the first
    # and last rows
    model.data[1, :] = 10
    model.data[-2, :] = 10

    result = ResidualFringeStep.call(model, skip=False)

    # Fit should complete
    assert not np.allclose(result.data, model.data)

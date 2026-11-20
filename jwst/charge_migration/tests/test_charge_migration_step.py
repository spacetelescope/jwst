import numpy as np
import pytest
from astropy.table import Table

from jwst.charge_migration.charge_migration_step import ChargeMigrationStep
from jwst.charge_migration.tests import helpers


def test_output_is_not_input():
    """Test that input is not modified."""
    ramp_model = helpers.niriss_ramp_model()

    result = ChargeMigrationStep.call(ramp_model, signal_threshold=0)

    # Successful completion
    assert result.meta.cal_step.charge_migration == "COMPLETE"

    # Output is not input
    assert result is not ramp_model
    assert ramp_model.meta.cal_step.charge_migration is None


def test_too_few_groups(caplog):
    """
    Test that processing for datasets having too few (<3) groups per integration
    are skipped.
    """
    shape = (1, 2, 1, 1)
    ramp_model = helpers.niriss_ramp_model(shape)
    result = ChargeMigrationStep.call(ramp_model, signal_threshold=0)

    assert result.meta.cal_step.charge_migration == "SKIPPED"
    assert "Too few groups" in caplog.text

    assert result is not ramp_model
    assert ramp_model.meta.cal_step.charge_migration is None


@pytest.mark.parametrize("signal_threshold", [-1, 0, 1])
@pytest.mark.parametrize("correction_threshold", [-1, 0, 1])
def test_thresholds(caplog, signal_threshold, correction_threshold):
    ramp_model = helpers.niriss_ramp_model()

    result = ChargeMigrationStep.call(
        ramp_model, signal_threshold=signal_threshold, correction_threshold=correction_threshold
    )

    if signal_threshold >= 0:
        # ran flagging
        assert "Flagging pixels" in caplog.text

    if correction_threshold >= 0:
        # ran correction
        assert "Fitting and correcting" in caplog.text

    if signal_threshold < 0 and correction_threshold < 0:
        # did not run either
        assert "All thresholds negative" in caplog.text
        assert result.meta.cal_step.charge_migration == "SKIPPED"
    else:
        # ran at least one
        assert result.meta.cal_step.charge_migration == "COMPLETE"


def test_save_fit_info(tmp_path):
    # Make a ramp with one pixel to correct
    threshold = 100
    xx, yy = 3, 2
    shape = (1, 3, 5, 5)
    ramp_model = helpers.niriss_ramp_model_with_signal(
        threshold, (yy, xx), shape=shape, start_group=1
    )

    result = ChargeMigrationStep.call(
        ramp_model,
        signal_threshold=-1,
        correction_threshold=threshold,
        mask_pad=1,
        save_fit_info=True,
        output_dir=str(tmp_path),
        output_file="test",
        suffix="chg",
        save_results=True,
    )

    # output file is saved
    assert (tmp_path / "test_chg.fits").exists()

    # additional fit info table is saved
    fit_file = tmp_path / "test_bfe_fit.asdf"
    assert fit_file.exists()

    # found one pixel to fit + 8 neighbors
    fit_info = Table.read(fit_file)
    assert len(fit_info) == 9
    np.testing.assert_array_equal(fit_info["x"], [2, 3, 4] * 3)
    np.testing.assert_array_equal(fit_info["y"], [1] * 3 + [2] * 3 + [3] * 3)

    # only one good fit, at the expected pixel
    assert np.sum(fit_info["fit_success"]) == 1
    assert fit_info["fit_success"][4]

    # corrected data at that pixel matches fit_info model
    model = fit_info["model"][4]
    slope = model[1] - model[0]
    corrected = (0, 2, yy, xx)
    np.testing.assert_allclose(result.data[corrected], ramp_model.data[corrected] - slope)

    # all other data not affected
    flat_index = np.ravel_multi_index(corrected, shape)
    np.testing.assert_allclose(
        np.delete(result.data, flat_index), np.delete(ramp_model.data, flat_index)
    )

    result.close()


def test_save_fit_info_no_fit(tmp_path, caplog):
    # Make a ramp that will attempt corrections but find no pixels to fit
    threshold = 100
    xx, yy = 3, 2
    shape = (1, 10, 5, 5)  # ngroup > 5 => high pixel is sigma clipped out
    ramp_model = helpers.niriss_ramp_model_with_signal(threshold, (yy, xx), shape=shape)

    result = ChargeMigrationStep.call(
        ramp_model,
        signal_threshold=-1,
        correction_threshold=threshold,
        mask_pad=1,
        save_fit_info=True,
        output_dir=str(tmp_path),
        output_file="test",
        suffix="chg",
        save_results=True,
    )

    # output file is saved
    assert (tmp_path / "test_chg.fits").exists()

    # additional fit info table is not saved
    fit_file = tmp_path / "test_bfe_fit.asdf"
    assert not fit_file.exists()

    assert "No pixels to correct" in caplog.text
    assert "No fit info to save" in caplog.text

    result.close()

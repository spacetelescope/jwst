import numpy as np
from stdatamodels.jwst import datamodels

from jwst.persistence import persistence
from jwst.persistence.persistence_step import PersistenceStep


def test_step_no_trapsfilled(create_sci_model, tmp_path):
    model = create_sci_model()
    result = PersistenceStep.call(model, output_dir=str(tmp_path), save_trapsfilled=True)

    # Step completed
    assert result.meta.cal_step.persistence == "COMPLETE"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.persistence is None

    # Without a trapsfilled image provided, it should create one,
    # output DQ is unmodified
    assert (tmp_path / "test_trapsfilled.fits").exists()
    np.testing.assert_array_equal(result.pixeldq, 0)


def test_step_with_trapsfilled(
    tmp_path, create_sci_model, create_traps_filled_model, create_trap_density_model
):
    sci = create_sci_model()
    traps = create_traps_filled_model()
    density = create_trap_density_model()

    trap_file = str(tmp_path / "trapsfilled.fits")
    traps.save(trap_file)

    result = PersistenceStep.call(
        sci,
        input_trapsfilled=trap_file,
        override_trapdensity=density,
        save_trapsfilled=True,
        save_persistence=True,
        output_dir=str(tmp_path),
    )

    # Step completed
    assert result.meta.cal_step.persistence == "COMPLETE"

    # Input is not modified
    assert result is not sci
    assert sci.meta.cal_step.persistence is None

    # Output trapsfilled and persistence files are saved, DQ is updated
    assert (tmp_path / "test_trapsfilled.fits").exists()
    assert (tmp_path / "test_output_pers.fits").exists()
    assert np.any(result.pixeldq & datamodels.dqflags.pixel["PERSISTENCE"] > 0)


def test_step_missing_one_ref(caplog, create_sci_model):
    sci = create_sci_model()
    result = PersistenceStep.call(sci, override_trapdensity="N/A")

    # Missing reference file is logged
    assert "Missing reference file type:  TRAPDENSITY" in caplog.text

    # Step is skipped
    assert result.meta.cal_step.persistence == "SKIPPED"

    # Input is not modified
    assert result is not sci
    assert sci.meta.cal_step.persistence is None


def test_step_missing_all_ref(caplog, create_sci_model):
    sci = create_sci_model()
    result = PersistenceStep.call(
        sci, override_trapdensity="N/A", override_trappars="N/A", override_persat="N/A"
    )

    # Missing reference files are logged
    assert "Missing reference file types:  PERSAT TRAPDENSITY TRAPPARS" in caplog.text

    # Step is skipped
    assert result.meta.cal_step.persistence == "SKIPPED"

    # Input is not modified
    assert result is not sci
    assert sci.meta.cal_step.persistence is None


def test_step_persistence_fails(monkeypatch, create_sci_model):
    # mock a known error condition in the persistence call
    monkeypatch.setattr(
        persistence.DataSet, "do_all", lambda *args: (datamodels.RampModel(), None, None, True)
    )

    sci = create_sci_model()
    result = PersistenceStep.call(sci)

    # Step is skipped
    assert result.meta.cal_step.persistence == "SKIPPED"

    # Input is not modified
    assert result is not sci
    assert sci.meta.cal_step.persistence is None

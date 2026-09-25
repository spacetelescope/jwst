import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file
from jwst.photom.photom_step import PhotomStep
from jwst.photom.tests.test_photom import create_input


@pytest.fixture()
def input_model():
    """
    Make an example model that can run successfully through the photom step.

    Yields
    ------
    MultiSlitModel
        The example model.
    """
    test_model = create_input(
        "NIRSPEC", "NRS1", "NRS_FIXEDSLIT", filter_used="F170LP", grating="G235M"
    )
    yield test_model
    test_model.close()


def test_step_complete(input_model):
    result = PhotomStep.call(input_model)

    # Step is complete
    assert result.meta.cal_step.photom == "COMPLETE"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.photom is None


def test_skip_unexpected_model(caplog, input_model):
    # make an unexpected datamodel type, with appropriate headers
    bad_model = datamodels.RampModel()
    bad_model.update(input_model)

    result = PhotomStep.call(bad_model)

    # Step is skipped
    assert "Input is not one of the supported model types" in caplog.text
    assert result.meta.cal_step.photom == "SKIPPED"

    # Input is not modified
    assert result is not bad_model
    assert bad_model.meta.cal_step.photom is None


def test_skip_no_reffile(caplog, input_model):
    result = PhotomStep.call(input_model, override_photom="N/A")

    # Step is skipped
    assert "No PHOTOM reference file found" in caplog.text
    assert result.meta.cal_step.photom == "SKIPPED"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.photom is None


def test_photom_fail(caplog):
    # Make a NIRISS cube model: it should fail photom, as an unexpected type
    nis_model = create_input("NIRISS", "NIS", "NIS_SOSS", filter_used="CLEAR", pupil="GR700XD")
    cube = datamodels.CubeModel()
    cube.update(nis_model)

    result = PhotomStep.call(cube)

    # Step is skipped
    assert "Unexpected data model type" in caplog.text
    assert result.meta.cal_step.photom == "SKIPPED"

    # Input is not modified
    assert result is not cube
    assert cube.meta.cal_step.photom is None

    nis_model.close()
    cube.close()


def test_photom_inverse(input_model):
    step = PhotomStep()
    result = step.run(input_model)

    # output data is not the same
    nnan = ~np.isnan(input_model.slits[0].data) & ~np.isnan(result.slits[0].data)
    assert not np.allclose(result.slits[0].data[nnan], input_model.slits[0].data[nnan])

    # run again but invert
    new_step = PhotomStep()
    new_step.inverse = True
    inverse_result = new_step.run(result)

    # output data is the same as input, except that NaNs do not invert
    assert np.allclose(inverse_result.slits[0].data[nnan], input_model.slits[0].data[nnan])
    assert np.allclose(inverse_result.slits[0].err[nnan], input_model.slits[0].err[nnan])
    assert np.allclose(
        inverse_result.slits[0].var_rnoise[nnan], input_model.slits[0].var_rnoise[nnan]
    )
    assert np.allclose(
        inverse_result.slits[0].var_poisson[nnan], input_model.slits[0].var_poisson[nnan]
    )
    assert np.allclose(inverse_result.slits[0].var_flat[nnan], input_model.slits[0].var_flat[nnan])

    result.close()
    inverse_result.close()


def test_photom_nrs_ifu():
    """
    Smoke test photom for NIRSpec IFU.

    Requires a realistic WCS, so testing this mode at the step
    interface level instead of with the other mocked modes in
    ``test_photom``.
    """
    # Make a model with reasonable metadata
    shape = (2048, 2048)
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    model = datamodels.IFUImageModel(hdul)
    hdul.close()

    # Assign some flat data and a WCS
    model.data = np.ones(shape, dtype=np.float32)
    model.err = np.full(shape, 0.1, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = model.get_default("var_poisson")
    model.var_rnoise = model.get_default("var_rnoise")
    model.var_flat = model.get_default("var_flat")
    model = AssignWcsStep.call(model)

    # Make one data pixel NaN, not matched with NaNs in the err/dq/var
    bad_idx = (1414, 690)
    model.data[bad_idx] = np.nan

    result = PhotomStep.call(model)
    assert result.meta.cal_step.photom == "COMPLETE"

    # result is not the same as input
    assert result is not model
    assert model.meta.cal_step.photom is None
    assert not np.allclose(result.data, model.data)
    assert not np.allclose(result.err, model.err)

    # NaNs and flags are matched for the bad pixel, input not modified
    assert not np.isnan(model.err[bad_idx])
    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        assert np.isnan(getattr(result, ext)[bad_idx])
    assert (result.dq[bad_idx] & datamodels.dqflags.pixel["DO_NOT_USE"]) > 0

    # inverting recovers the input, except that the new NaN is not invertible
    inverse_result = PhotomStep.call(result, inverse=True)
    np.testing.assert_allclose(inverse_result.data, model.data)
    mask = np.full(model.err.shape, True)
    mask[bad_idx] = False
    np.testing.assert_allclose(inverse_result.err[mask], model.err[mask])
    assert np.isnan(inverse_result.err[bad_idx])

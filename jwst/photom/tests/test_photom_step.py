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


def test_photom_correction_pars(input_model):
    step = PhotomStep()
    result = step.run(input_model)

    # output data is not the same
    nnan = ~np.isnan(input_model.slits[0].data) & ~np.isnan(result.slits[0].data)
    assert not np.allclose(result.slits[0].data[nnan], input_model.slits[0].data[nnan])

    # use the computed correction and invert
    new_step = PhotomStep()
    new_step.use_correction_pars = True
    new_step.correction_pars = step.correction_pars
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
    model = AssignWcsStep.call(model)

    result = PhotomStep.call(model)
    assert result.meta.cal_step.photom == "COMPLETE"

    # result is not the same as input
    assert result is not model
    assert model.meta.cal_step.photom is None
    assert not np.allclose(result.data, model.data)
    assert not np.allclose(result.err, model.err)

    # inverting recovers the input
    inverse_result = PhotomStep.call(result, inverse=True)
    np.testing.assert_allclose(inverse_result.data, model.data)
    np.testing.assert_allclose(inverse_result.err, model.err)

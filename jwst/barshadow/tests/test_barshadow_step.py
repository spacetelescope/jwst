import os
from copy import deepcopy

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

import jwst
from jwst.assign_wcs import AssignWcsStep
from jwst.barshadow import BarShadowStep
from jwst.extract_2d import Extract2dStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_mos_file


def create_nirspec_mos_model():
    hdul = create_nirspec_mos_file()
    msa_meta = os.path.join(
        jwst.__path__[0], *["assign_wcs", "tests", "data", "msa_configuration.fits"]
    )
    hdul[0].header["MSAMETFL"] = msa_meta
    hdul[0].header["MSAMETID"] = 12
    im = datamodels.ImageModel(hdul)
    hdul.close()

    im.data = np.full((2048, 2048), 1.0)
    im_wcs = AssignWcsStep.call(im)
    im_ex2d = Extract2dStep.call(im_wcs)

    # add error/variance arrays
    im_ex2d.slits[0].name = "0"
    im_ex2d.slits[0].err = im_ex2d.slits[0].data * 0.1
    im_ex2d.slits[0].var_rnoise = im_ex2d.slits[0].data * 0.01
    im_ex2d.slits[0].var_poisson = im_ex2d.slits[0].data * 0.01
    im_ex2d.slits[0].var_flat = im_ex2d.slits[0].data * 0.01

    # add a couple more slits
    for i in [1, 2]:
        slit_copy = deepcopy(im_ex2d.slits[0])
        slit_copy.name = str(i)
        im_ex2d.slits.append(slit_copy)

    return im_ex2d


@pytest.fixture(scope="module")
def nirspec_mos_model():
    return create_nirspec_mos_model()


def test_barshadow_step(nirspec_mos_model):
    model = nirspec_mos_model.copy()
    result = BarShadowStep.call(model)
    assert result.meta.cal_step.barshadow == "COMPLETE"

    # check all slits for appropriate correction
    for slit in result.slits:
        assert slit.barshadow_corrected is True

        # 5 shutter slitlet, correction should not be uniform
        shadow = slit.barshadow
        assert not np.all(shadow == 1)

        # maximum correction value should be 1.0, minimum should be above zero
        assert np.nanmax(shadow) <= 1.0
        assert np.nanmin(shadow) > 0.0

        # data should have been divided by barshadow
        nnan = ~np.isnan(slit.data)
        assert np.allclose(slit.data[nnan] * shadow[nnan], model.slits[0].data[nnan])
        assert np.allclose(slit.err[nnan] * shadow[nnan], model.slits[0].err[nnan])
        assert np.allclose(
            slit.var_rnoise[nnan] * shadow[nnan] ** 2,
            model.slits[0].var_rnoise[nnan],
        )
        assert np.allclose(
            slit.var_poisson[nnan] * shadow[nnan] ** 2,
            model.slits[0].var_poisson[nnan],
        )
        assert np.allclose(
            slit.var_flat[nnan] * shadow[nnan] ** 2,
            model.slits[0].var_flat[nnan],
        )

    result.close()


def test_barshadow_step_zero_length(nirspec_mos_model, log_watcher):
    model = nirspec_mos_model.copy()
    model.slits[0].shutter_state = ""

    watcher = log_watcher(
        "jwst.barshadow.bar_shadow", message="has zero length, correction skipped", level="info"
    )
    result = BarShadowStep.call(model)
    watcher.assert_seen()

    # correction step ran, but is all 1s
    assert result.meta.cal_step.barshadow == "COMPLETE"
    assert np.all(result.slits[0].barshadow == 1)

    # correction status is False
    assert result.slits[0].barshadow_corrected is False

    result.close()


def test_barshadow_step_not_uniform(nirspec_mos_model, log_watcher):
    model = nirspec_mos_model.copy()
    for slit in model.slits:
        slit.source_type = "POINT"

    watcher = log_watcher("jwst.barshadow.bar_shadow", message="source not uniform")
    result = BarShadowStep.call(model)
    watcher.assert_seen()

    # correction ran, but is all 1s
    assert result.meta.cal_step.barshadow == "COMPLETE"
    for slit in result.slits:
        assert np.all(slit.barshadow == 1)
        assert slit.barshadow_corrected is False
    result.close()


def test_barshadow_no_reffile(monkeypatch, nirspec_mos_model):
    model = nirspec_mos_model.copy()
    monkeypatch.setattr(BarShadowStep, "get_reference_file", lambda *args, **kwargs: "N/A")

    result = BarShadowStep.call(model)

    # correction did not run
    assert result.meta.cal_step.barshadow == "SKIPPED"
    assert result.slits[0].barshadow.size == 0
    assert result.slits[0].barshadow_corrected is None
    result.close()


def test_barshadow_wrong_exptype():
    model = datamodels.MultiSlitModel()
    model.meta.exposure.type = "ANY"
    model.slits.append(datamodels.SlitModel())
    result = BarShadowStep.call(model)

    # correction did not run
    assert result.meta.cal_step.barshadow == "SKIPPED"
    assert result.slits[0].barshadow.size == 0
    assert result.slits[0].barshadow_corrected is None

    result.close()


def test_barshadow_correction_pars(nirspec_mos_model):
    model = nirspec_mos_model.copy()
    step = BarShadowStep()
    result = step.run(model)

    # output data is not the same
    nnan = ~np.isnan(model.slits[0].data) & ~np.isnan(result.slits[0].data)
    assert not np.allclose(result.slits[0].data[nnan], model.slits[0].data[nnan])

    # use the computed correction and invert
    new_step = BarShadowStep()
    new_step.use_correction_pars = True
    new_step.correction_pars = step.correction_pars
    new_step.inverse = True
    inverse_result = new_step.run(result)

    # output data is the same as input, except that NaNs do not invert
    assert np.allclose(inverse_result.slits[0].data[nnan], model.slits[0].data[nnan])
    assert np.allclose(inverse_result.slits[0].err[nnan], model.slits[0].err[nnan])
    assert np.allclose(inverse_result.slits[0].var_rnoise[nnan], model.slits[0].var_rnoise[nnan])
    assert np.allclose(inverse_result.slits[0].var_poisson[nnan], model.slits[0].var_poisson[nnan])
    assert np.allclose(inverse_result.slits[0].var_flat[nnan], model.slits[0].var_flat[nnan])

    result.close()
    inverse_result.close()


def test_barshadow_step_missing_scale(nirspec_mos_model, log_watcher):
    model = nirspec_mos_model.copy()
    model.slits[0].slit_yscale = None

    watcher = log_watcher("jwst.barshadow.bar_shadow", message="Using default value")
    result = BarShadowStep.call(model)
    watcher.assert_seen()

    # correction ran and has an appropriate correction - the
    # default value is close enough for most purposes.
    assert result.meta.cal_step.barshadow == "COMPLETE"

    # maximum correction value should be 1.0, minimum should be above zero
    shadow = result.slits[0].barshadow
    assert np.nanmax(shadow) <= 1.0
    assert np.nanmin(shadow) > 0.0

    # data should have been divided by barshadow
    nnan = ~np.isnan(result.slits[0].data)
    assert np.allclose(result.slits[0].data[nnan] * shadow[nnan], model.slits[0].data[nnan])
    result.close()

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

import jwst
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import (
    create_nirspec_mos_file,
    create_nirspec_fs_file,
    create_nirspec_ifu_file,
)
from jwst.extract_2d import Extract2dStep
from jwst.pathloss import PathLossStep


def create_nirspec_fs_model(source_type="POINT"):
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    im = datamodels.ImageModel(hdul)
    hdul.close()

    im.data = np.full((2048, 2048), 1.0)
    im_wcs = AssignWcsStep.call(im)
    im_ex2d = Extract2dStep.call(im_wcs)

    # add error/variance arrays
    im_ex2d.slits[0].err = im_ex2d.slits[0].data * 0.1
    im_ex2d.slits[0].var_rnoise = im_ex2d.slits[0].data * 0.01
    im_ex2d.slits[0].var_poisson = im_ex2d.slits[0].data * 0.01
    im_ex2d.slits[0].var_flat = im_ex2d.slits[0].data * 0.01

    # set source type and position off center
    im_ex2d.slits[0].source_type = source_type
    im_ex2d.slits[0].source_xpos = 0.25
    im_ex2d.slits[0].source_ypos = 0.25

    return im_ex2d


def create_nirspec_mos_model(source_type="POINT"):
    hdul = create_nirspec_mos_file()
    msa_meta = Path(jwst.__path__[0]) / "assign_wcs" / "tests" / "data" / "msa_configuration.fits"
    hdul[0].header["MSAMETFL"] = str(msa_meta)
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

    # set source type and position off center
    im_ex2d.slits[0].source_type = source_type
    im_ex2d.slits[0].source_xpos = 0.25
    im_ex2d.slits[0].source_ypos = 0.25

    # add a couple more slits
    for i in [1, 2]:
        slit_copy = deepcopy(im_ex2d.slits[0])
        slit_copy.name = str(i)
        im_ex2d.slits.append(slit_copy)

    return im_ex2d


@pytest.fixture(scope="module")
def nirspec_mos_model_point():
    return create_nirspec_mos_model(source_type="POINT")


@pytest.fixture(scope="module")
def nirspec_mos_model_extended():
    return create_nirspec_mos_model(source_type="EXTENDED")


@pytest.fixture(scope="module")
def nirspec_fs_model_point():
    return create_nirspec_fs_model(source_type="POINT")


@pytest.fixture(scope="module")
def nirspec_fs_model_extended():
    return create_nirspec_fs_model(source_type="EXTENDED")


def test_pathloss_step_mos_point(nirspec_mos_model_point):
    model = nirspec_mos_model_point.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # check all slits for appropriate correction
    for slit in result.slits:
        # correction type should be "POINT"
        assert slit.pathloss_correction_type == "POINT"

        # uniform also present, but point is used
        pathloss = slit.pathloss_point
        pathloss_un = slit.pathloss_uniform
        assert not np.all(pathloss == 1)
        assert not np.all(pathloss == pathloss_un)

        # maximum correction value should be 1.0, minimum should be above zero
        assert np.nanmax(pathloss) <= 1.0
        assert np.nanmin(pathloss) > 0.0

        # data should have been divided by pathloss
        nnan = ~np.isnan(slit.data)
        assert np.allclose(slit.data[nnan] * pathloss[nnan], model.slits[0].data[nnan])
        assert np.allclose(slit.err[nnan] * pathloss[nnan], model.slits[0].err[nnan])
        assert np.allclose(
            slit.var_rnoise[nnan] * pathloss[nnan] ** 2,
            model.slits[0].var_rnoise[nnan],
        )
        assert np.allclose(
            slit.var_poisson[nnan] * pathloss[nnan] ** 2,
            model.slits[0].var_poisson[nnan],
        )
        assert np.allclose(
            slit.var_flat[nnan] * pathloss[nnan] ** 2,
            model.slits[0].var_flat[nnan],
        )

    result.close()


def test_pathloss_step_mos_uniform(nirspec_mos_model_extended):
    model = nirspec_mos_model_extended.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # check all slits for appropriate correction
    for slit in result.slits:
        # correction type should be "UNIFORM"
        assert slit.pathloss_correction_type == "UNIFORM"

        # point also present, but uniform is used
        pathloss = slit.pathloss_uniform
        pathloss_pt = slit.pathloss_point
        assert not np.all(pathloss == 1)
        assert not np.all(pathloss == pathloss_pt)

        # maximum correction value can be > 1.0 for uniform correction,
        # minimum should be above zero
        assert np.nanmin(pathloss) > 0.0

        # data should have been divided by pathloss
        nnan = ~np.isnan(slit.data)
        assert np.allclose(slit.data[nnan] * pathloss[nnan], model.slits[0].data[nnan])
        assert np.allclose(slit.err[nnan] * pathloss[nnan], model.slits[0].err[nnan])
        assert np.allclose(
            slit.var_rnoise[nnan] * pathloss[nnan] ** 2,
            model.slits[0].var_rnoise[nnan],
        )
        assert np.allclose(
            slit.var_poisson[nnan] * pathloss[nnan] ** 2,
            model.slits[0].var_poisson[nnan],
        )
        assert np.allclose(
            slit.var_flat[nnan] * pathloss[nnan] ** 2,
            model.slits[0].var_flat[nnan],
        )

    result.close()


def test_pathloss_step_fs_point(nirspec_fs_model_point):
    model = nirspec_fs_model_point.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # correction type should be "POINT" for the first slit,
    # with source type set; uniform otherwise
    for i, slit in enumerate(result.slits):
        if i == 0:
            assert slit.pathloss_correction_type == "POINT"
        else:
            assert slit.pathloss_correction_type == "UNIFORM"


def test_pathloss_step_fs_uniform(nirspec_fs_model_extended):
    model = nirspec_fs_model_extended.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # correction type should be "UNIFORM" for all slits
    for slit in result.slits:
        assert slit.pathloss_correction_type == "UNIFORM"


@pytest.mark.parametrize("mode", ["mos", "fs"])
def test_pathloss_correction_pars(mode, nirspec_mos_model_point, nirspec_fs_model_point):
    if mode == "mos":
        model = nirspec_mos_model_point.copy()
    else:
        model = nirspec_fs_model_point.copy()
    step = PathLossStep()
    result = step.run(model)

    # output data is not the same
    nnan = ~np.isnan(model.slits[0].data) & ~np.isnan(result.slits[0].data)
    assert not np.allclose(result.slits[0].data[nnan], model.slits[0].data[nnan])

    # use the computed correction and invert
    new_step = PathLossStep()
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

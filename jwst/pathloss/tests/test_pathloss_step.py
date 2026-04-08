from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

import jwst
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import (
    create_nirspec_fs_file,
    create_nirspec_ifu_file,
    create_nirspec_mos_file,
)
from jwst.extract_1d.tests.helpers import mock_miri_lrs_fs_func, mock_niriss_soss_96_func
from jwst.extract_2d import Extract2dStep
from jwst.pathloss import PathLossStep


def create_nirspec_fs_model(source_type="POINT"):
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    im = datamodels.ImageModel(hdul)
    hdul.close()

    im.data = np.full((2048, 2048), 1.0)
    im.dq = np.zeros((2048, 2048), dtype=np.uint32)
    im.err = im.data * 0.1
    im.var_rnoise = im.data * 0.01
    im.var_poisson = im.data * 0.01
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
    im.dq = np.zeros((2048, 2048), dtype=np.uint32)
    im.err = im.data * 0.1
    im.var_rnoise = im.data * 0.01
    im.var_poisson = im.data * 0.01
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


def create_nirspec_ifu_model(source_type="POINT"):
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    im = datamodels.IFUImageModel(hdul)
    hdul.close()

    im.meta.target.source_type = source_type
    im.data = np.full((2048, 2048), 1.0)
    im.dq = np.zeros((2048, 2048), dtype=np.uint32)
    im.err = im.data * 0.1
    im.var_rnoise = im.data * 0.01
    im.var_poisson = im.data * 0.01
    im.var_flat = im.data * 0.01
    im_wcs = AssignWcsStep.call(im)

    im.close()
    return im_wcs


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


@pytest.fixture(scope="module")
def nirspec_ifu_model_point():
    return create_nirspec_ifu_model(source_type="POINT")


@pytest.fixture(scope="module")
def nirspec_ifu_model_extended():
    return create_nirspec_ifu_model(source_type="EXTENDED")


@pytest.fixture(scope="module")
def miri_lrs_model_point():
    model = mock_miri_lrs_fs_func()

    # Mock some WCS data needed for pathloss correction
    model.meta.wcs.pipeline[0].transform.offset_1 = -10
    model.meta.wcs.pipeline[0].transform.offset_2 = -10
    model.meta.target.ra = 45.0
    model.meta.target.dec = 45.1
    model.meta.target.source_type = "POINT"
    return model


@pytest.fixture(scope="module")
def miri_lrs_model_extended():
    model = mock_miri_lrs_fs_func()
    return model


@pytest.fixture(scope="module")
def niriss_soss_model():
    cube_model = mock_niriss_soss_96_func()

    # Make a non-TSO rate image from the cube model
    model = datamodels.ImageModel()
    for attr in ["data", "dq", "err", "var_rnoise", "var_poisson", "var_flat"]:
        setattr(model, attr, getattr(cube_model, attr)[0, :, :])
    model.update(cube_model)
    cube_model.close()

    return model


def test_pathloss_step_mos_point(nirspec_mos_model_point):
    model = nirspec_mos_model_point.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

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

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

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

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

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

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

    # correction type should be "UNIFORM" for all slits
    for slit in result.slits:
        assert slit.pathloss_correction_type == "UNIFORM"


@pytest.mark.parametrize("dataset", ["nirspec_fs_model_point", "nirspec_mos_model_point"])
@pytest.mark.parametrize("all_slits", [True, False])
def test_pathloss_step_source_outside_slit(caplog, request, dataset, all_slits):
    # Datamodel with all source positions outside slit
    model = request.getfixturevalue(dataset).copy()
    if all_slits:
        for slit in model.slits:
            slit.source_xpos = 100
            slit.source_type = "POINT"
    else:
        # only the point source slit is outside the slit
        model.slits[0].source_xpos = 100

    result = PathLossStep.call(model)
    assert "outside slit" in caplog.text

    if all_slits:
        # No slits were corrected
        assert result.meta.cal_step.pathloss == "SKIPPED"
    else:
        # Some slits were corrected
        assert result.meta.cal_step.pathloss == "COMPLETE"


def test_pathloss_step_ifu_point(caplog, nirspec_ifu_model_point):
    model = nirspec_ifu_model_point.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

    # correction type should be "POINT"
    assert "Correction type used: POINT" in caplog.text


def test_pathloss_step_ifu_uniform(caplog, nirspec_ifu_model_extended):
    model = nirspec_ifu_model_extended.copy()
    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

    # correction type should be "UNIFORM"
    assert "Correction type used: UNIFORM" in caplog.text


def test_pathloss_step_miri_lrs_point(miri_lrs_model_point):
    model = miri_lrs_model_point.copy()

    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"
    assert isinstance(result, datamodels.ImageModel)
    assert not np.allclose(result.data, model.data)

    # LRS correction is multiplicative: it will be > 1
    assert np.mean(result.pathloss_point) > 1.0

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None
    assert not model.hasattr("pathloss_point")


def test_pathloss_step_miri_lrs_extended(miri_lrs_model_extended):
    model = miri_lrs_model_extended.copy()

    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "SKIPPED"
    assert isinstance(result, datamodels.ImageModel)
    assert np.allclose(result.data, model.data)
    assert not result.hasattr("pathloss_point")

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None
    assert not model.hasattr("pathloss_point")


def test_pathloss_step_miri_lrs_user_slit_loc(caplog, miri_lrs_model_point):
    model = miri_lrs_model_point.copy()

    result = PathLossStep.call(model, user_slit_loc=0.1)
    assert result.meta.cal_step.pathloss == "COMPLETE"
    assert not np.allclose(result.data, model.data)
    assert np.mean(result.pathloss_point) > 1.0

    # Check for the expected target location
    # Mock spatial scale is 0.1 arcsec, offset is along dispersion direction (y),
    # so expected location is (0, 1)
    assert "target center offset: 0.1 arcsec" in caplog.text
    assert "New target location = (0.000, 1.000)" in caplog.text


def test_pathloss_step_miri_lrs_outside_slit(caplog, miri_lrs_model_point):
    model = miri_lrs_model_point.copy()

    # Call the step once with default values: source is at the center
    result_center = PathLossStep.call(model)

    # Call again but place the target outside the slit
    result_outside = PathLossStep.call(model, user_slit_loc=5)

    # Results are the same: the correction falls back to using the center of the slit
    np.testing.assert_allclose(result_center.data, result_outside.data)
    assert "Source is outside slit. Correction defaulting to center of the slit" in caplog.text


def test_pathloss_step_niriss_soss_non_tso(niriss_soss_model):
    model = niriss_soss_model.copy()

    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"
    assert isinstance(result, datamodels.ImageModel)
    assert not np.allclose(result.data, model.data)

    # SOSS correction divides out: it will be < 1
    assert np.mean(result.pathloss_point) < 1.0

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None
    assert not model.hasattr("pathloss_point")


def test_pathloss_step_niriss_soss_pupil_out_of_range(niriss_soss_model):
    model = niriss_soss_model.copy()
    model.meta.instrument.pupil_position = 1

    result = PathLossStep.call(model)
    assert result.meta.cal_step.pathloss == "COMPLETE"
    assert isinstance(result, datamodels.ImageModel)

    # SOSS correction is set to 1.0
    assert np.allclose(result.data, model.data)
    assert np.allclose(result.pathloss_point, 1.0)

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None
    assert not model.hasattr("pathloss_point")


@pytest.mark.parametrize(
    "dataset", ["nirspec_fs_model_point", "nirspec_mos_model_point", "nirspec_ifu_model_point"]
)
def test_pathloss_inverse(request, dataset):
    model = request.getfixturevalue(dataset).copy()
    step = PathLossStep()
    result = step.run(model)

    # run again but invert
    new_step = PathLossStep()
    new_step.inverse = True
    inverse_result = new_step.run(result)

    if "ifu" in dataset:
        compare_image = model
        result_image = result
        inverse_image = inverse_result
    else:
        compare_image = model.slits[0]
        result_image = result.slits[0]
        inverse_image = inverse_result.slits[0]

    # output data is not the same, except for IFU,
    # for which the correction currently has no effect
    nnan = ~np.isnan(compare_image.data) & ~np.isnan(result_image.data)
    if "ifu" not in dataset:
        assert not np.allclose(result_image.data[nnan], compare_image.data[nnan])

    # output data is the same as input, except that NaNs do not invert
    assert np.allclose(inverse_image.data[nnan], compare_image.data[nnan])
    assert np.allclose(inverse_image.err[nnan], compare_image.err[nnan])
    assert np.allclose(inverse_image.var_rnoise[nnan], compare_image.var_rnoise[nnan])
    assert np.allclose(inverse_image.var_poisson[nnan], compare_image.var_poisson[nnan])
    assert np.allclose(inverse_image.var_flat[nnan], compare_image.var_flat[nnan])

    result.close()
    inverse_result.close()


def test_pathloss_no_reffile(caplog, nirspec_fs_model_point):
    model = nirspec_fs_model_point.copy()
    result = PathLossStep.call(model, override_pathloss="N/A")

    # step is skipped
    assert result.meta.cal_step.pathloss == "SKIPPED"
    assert "No PATHLOSS reference file" in caplog.text
    assert np.all(result.slits[0].data == model.slits[0].data)
    assert not result.slits[0].hasattr("pathloss_point")

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.pathloss is None

import os
import numpy as np
import pytest

from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer
from stdatamodels.jwst.datamodels.dqflags import pixel as flags

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file
from jwst.pixel_replace.pixel_replace_step import PixelReplaceStep
from glob import glob


def cal_data(shape, bad_idx, dispaxis=1, model="slit"):
    if model == "image":
        model = datamodels.ImageModel(shape)
    elif model == "ifu":
        model = datamodels.IFUImageModel(shape)
    else:
        model = datamodels.SlitModel(shape)
    model.meta.wcsinfo.dispersion_direction = dispaxis

    # Set the data and error arrays to all 1s except one bad pixel
    # to correct at the middle of the array
    model.data[:] = 1.0
    model.err[:] = 1.0
    model.var_poisson[:] = 1.0
    model.var_rnoise[:] = 1.0
    model.var_flat[:] = 1.0

    bad_flag = flags["DO_NOT_USE"] + flags["OTHER_BAD_PIXEL"]
    model.data[bad_idx] = np.nan
    model.err[bad_idx] = np.nan
    model.var_poisson[bad_idx] = np.nan
    model.var_rnoise[bad_idx] = np.nan
    model.var_flat[bad_idx] = np.nan
    model.dq[bad_idx] = bad_flag

    # Also add a non-science region in one row and one column
    non_science = flags["DO_NOT_USE"] + flags["NON_SCIENCE"]
    model.data[..., 1] = np.nan
    model.err[..., 1] = np.nan
    model.var_poisson[..., 1] = np.nan
    model.var_rnoise[..., 1] = np.nan
    model.var_flat[..., 1] = np.nan
    model.dq[..., 1] = non_science

    model.data[..., 1, :] = np.nan
    model.err[..., 1, :] = np.nan
    model.var_poisson[..., 1, :] = np.nan
    model.var_rnoise[..., 1, :] = np.nan
    model.var_flat[..., 1, :] = np.nan
    model.dq[..., 1, :] = non_science

    return model


def nirspec_tso():
    bad_idx = (1, 10, 10)
    model = cal_data(shape=(3, 20, 20), bad_idx=bad_idx, dispaxis=1)
    model.meta.instrument.name = "NIRSPEC"
    model.meta.exposure.type = "NRS_BRIGHTOBJ"
    return model, bad_idx


def nirspec_fs_slitmodel():
    bad_idx = (10, 10)
    model = cal_data(shape=(20, 20), bad_idx=bad_idx, dispaxis=1)
    model.meta.instrument.name = "NIRSPEC"
    model.meta.exposure.type = "NRS_FIXEDSLIT"
    return model, bad_idx


def nirspec_msa_multislit():
    bad_idx = (10, 10)
    slit_model = cal_data(shape=(20, 20), bad_idx=bad_idx, dispaxis=1)
    model = datamodels.MultiSlitModel()
    model.slits.append(slit_model)
    model.meta.instrument.name = "NIRSPEC"
    model.meta.exposure.type = "NRS_MSASPEC"
    return model, bad_idx


def nirspec_ifu():
    shape = (2048, 2048)
    bad_idx = (1414, 690)

    # IFU mode requires WCS information, so make a more realistic model
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    hdul["SCI"].data = np.ones(shape, dtype=float)

    model = datamodels.IFUImageModel(hdul)
    model = AssignWcsStep.call(model)

    test_data = cal_data(shape=shape, bad_idx=bad_idx, dispaxis=1, model="ifu")
    model.data = test_data.data
    model.dq = test_data.dq
    model.err = test_data.err
    model.var_poisson = test_data.var_poisson
    model.var_rnoise = test_data.var_rnoise
    model.var_flat = test_data.var_flat

    test_data.close()

    return model, bad_idx


def miri_lrs():
    bad_idx = (10, 10)
    model = cal_data(shape=(20, 20), bad_idx=bad_idx, dispaxis=2, model="image")
    model.meta.instrument.name = "MIRI"
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    return model, bad_idx


def miri_mrs():
    shape = (20, 20)

    def mock_transform(*args):
        return None, np.full(shape, 1), None

    bad_idx = (10, 10)
    model = cal_data(shape=shape, bad_idx=bad_idx, dispaxis=2, model="ifu")
    model.meta.instrument.name = "MIRI"
    model.meta.exposure.type = "MIR_MRS"
    model.meta.wcs = {"transform": mock_transform}
    return model, bad_idx


@pytest.mark.parametrize(
    "input_model_function", [nirspec_tso, nirspec_fs_slitmodel, miri_lrs, miri_mrs]
)
@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_pixel_replace_no_container(input_model_function, algorithm):
    """
    Test pixel replace for modes with no container.

    This includes ImageModel, SlitModel, and IFUImageModel.
    """
    input_model, bad_idx = input_model_function()

    # for this simple case, the results from either algorithm should
    # be the same
    result = PixelReplaceStep.call(input_model, skip=False, algorithm=algorithm)

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result, ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result, ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model, ext)[bad_idx])
        assert getattr(result, ext)[bad_idx] == 1.0

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.dq[bad_idx] == (
        input_model.dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.dq[..., :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.dq[..., 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    result.close()
    input_model.close()


@pytest.mark.parametrize("input_model_function", [nirspec_msa_multislit])
@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_pixel_replace_multislit(input_model_function, algorithm):
    """Test pixel replace for multislit modes."""
    input_model, bad_idx = input_model_function()

    # for this simple case, the results from either algorithm should
    # be the same
    result = PixelReplaceStep.call(input_model, skip=False, algorithm=algorithm)

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result.slits[0], ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result.slits[0], ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model.slits[0], ext)[bad_idx])
        assert getattr(result.slits[0], ext)[bad_idx] == 1.0

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.slits[0].dq[bad_idx] == (
        input_model.slits[0].dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.slits[0].dq[..., :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.slits[0].dq[..., 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    result.close()
    input_model.close()


@pytest.mark.slow
@pytest.mark.parametrize("input_model_function", [nirspec_ifu])
@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_pixel_replace_nirspec_ifu(tmp_cwd, input_model_function, algorithm):
    """
    Test pixel replacement for NIRSpec IFU.

    Larger data and more WCS operations required for testing make
    this test take more than a minute, so marking this test 'slow'.

    The test is otherwise the same as for other modes.
    """
    input_model, bad_idx = input_model_function()
    input_model.meta.filename = "jwst_nirspec_cal.fits"

    # for this simple case, the results from either algorithm should
    # be the same
    result = PixelReplaceStep.call(input_model, skip=False, algorithm=algorithm, save_results=True)

    assert result.meta.filename == "jwst_nirspec_pixelreplacestep.fits"
    assert result.meta.cal_step.pixel_replace == "COMPLETE"
    assert os.path.isfile(result.meta.filename)

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result, ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result, ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model, ext)[bad_idx])
        assert getattr(result, ext)[bad_idx] == 1.0

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.dq[bad_idx] == (
        input_model.dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.dq[..., :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.dq[..., 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    result.close()
    input_model.close()


@pytest.mark.parametrize("input_model_function", [nirspec_fs_slitmodel])
def test_pixel_replace_container_names(tmp_cwd, input_model_function):
    """Test pixel replace output names for input container."""
    input_model, _ = input_model_function()
    input_model.meta.filename = "jwst_nirspec_1_cal.fits"
    input_model2, _ = input_model_function()
    input_model2.meta.filename = "jwst_nirspec_2_cal.fits"
    cfiles = [input_model, input_model2]
    container = ModelContainer(cfiles)

    expected_name = ["jwst_nirspec_1_pixelreplacestep.fits", "jwst_nirspec_2_pixelreplacestep.fits"]

    result = PixelReplaceStep.call(container, skip=False, save_results=True)
    for i, model in enumerate(result):
        assert model.meta.filename == expected_name[i]
        assert model.meta.cal_step.pixel_replace == "COMPLETE"

    result_files = glob(os.path.join(tmp_cwd, "*pixelreplacestep.fits"))
    for i, file in enumerate(sorted(result_files)):
        basename = os.path.basename(file)
        assert expected_name[i] == basename
        with datamodels.open(file) as model:
            assert model.meta.cal_step.pixel_replace == "COMPLETE"
            assert model.meta.filename == expected_name[i]

    result.close()
    input_model.close()

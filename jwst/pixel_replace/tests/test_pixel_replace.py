import os
from glob import glob

import gwcs
import numpy as np
import pytest
from astropy.modeling.models import Const1D, Mapping
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.dqflags import pixel as flags

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file
from jwst.datamodels import ModelContainer
from jwst.pixel_replace.pixel_replace_step import PixelReplaceStep


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
    ones = np.ones(shape, dtype=float)
    model.data = ones.copy()
    model.err = ones.copy()
    model.var_poisson = ones.copy()
    model.var_rnoise = ones.copy()
    model.var_flat = ones.copy()

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
    bad_idx = (10, 10)
    model = cal_data(shape=shape, bad_idx=bad_idx, dispaxis=2, model="ifu")
    model.meta.instrument.name = "MIRI"
    model.meta.exposure.type = "MIR_MRS"

    # Mock a wcs that just returns 1 for alpha, beta, lam
    transform = Mapping((0, 1, 1), n_inputs=2) | Const1D(1) & Const1D(1) & Const1D(1)
    model.meta.wcs = gwcs.WCS([("detector", transform), ("alpha_beta", None)])

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

    # Step is recorded as complete
    assert result.meta.cal_step.pixel_replace == "COMPLETE"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.pixel_replace is None

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

    # Step is recorded as complete
    assert result.meta.cal_step.pixel_replace == "COMPLETE"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.pixel_replace is None

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

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.pixel_replace is None

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

    # Input is not modified
    for model in container:
        assert model.meta.cal_step.pixel_replace is None

    result.close()
    input_model.close()


def test_pixel_replace_no_valid_data(caplog):
    """Test pixel replace for no valid data."""
    input_model, bad_idx = nirspec_tso()

    # Set a middle region to NaN to test invalid data handling
    input_model.data[2, :, 5:15] = np.nan
    input_model.dq[2, 1:-1, 5:15] = 1

    # Set one pixel valid in the middle with no valid data next to
    # it to test missing adjacent data
    input_model.data[2, 10, 10] = 1.0
    input_model.dq[2, 10, 10] = 0
    input_model.dq[2, :, 9] = 1
    input_model.dq[2, :, 11] = 1

    result = PixelReplaceStep.call(input_model, algorithm="fit_profile", n_adjacent_cols=1)

    assert caplog.text.count("has no valid values - skipping") == 7
    assert caplog.text.count("has no valid adjacent values - skipping") == 1

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result, ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result, ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model, ext)[bad_idx])
        assert getattr(result, ext)[bad_idx] == 1.0

        # invalid data is left alone
        assert np.all(np.isnan(result.data[2, 1:-1, 5:10]))
        assert np.all(np.isnan(result.data[2, 1:-1, 11:15]))

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.dq[bad_idx] == (
        input_model.dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.dq[:2, :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.dq[:2, 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    # Invalid region is still marked DNU
    assert np.all(result.dq[2, 1:-1, 5:10] == flags["DO_NOT_USE"])
    assert np.all(result.dq[2, 1:-1, 11:15] == flags["DO_NOT_USE"])
    assert np.all(result.dq[2, 10, 10] == 0)
    assert np.all(result.dq[2, :, 9] == flags["DO_NOT_USE"])
    assert np.all(result.dq[2, :, 11] == flags["DO_NOT_USE"])

    result.close()
    input_model.close()


def test_skip_unexpected_type():
    bad_model = datamodels.RampModel()
    result = PixelReplaceStep.call(bad_model)

    # Step is skipped
    assert result.meta.cal_step.pixel_replace == "SKIPPED"

    # Input is not modified
    assert result is not bad_model
    assert bad_model.meta.cal_step.pixel_replace is None


def test_skip_unexpected_type_in_container():
    bad_model = datamodels.RampModel()
    container = ModelContainer([bad_model])
    result = PixelReplaceStep.call(container)

    # Step is skipped
    assert result[0].meta.cal_step.pixel_replace == "SKIPPED"

    # Input is not modified
    assert result[0] is not bad_model
    assert bad_model.meta.cal_step.pixel_replace is None

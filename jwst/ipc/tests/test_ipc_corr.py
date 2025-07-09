import numpy as np
import pytest
from collections import namedtuple

from stdatamodels.jwst import datamodels
from jwst.ipc import ipc_corr, x_irs2
from jwst.ipc.ipc_step import IPCStep

# Define the nref namedtuple
from jwst.ipc.ipc_corr import NumRefPixels as Nref


@pytest.fixture
def science_data():
    data = np.zeros((2048, 2048), dtype=np.float32)
    for i in range(2048):
        for j in range(2048):
            data[i, j] = 1000.0 + 2.0 * i + 3.0 * j
    return data


@pytest.fixture
def simple_kernel_2d():
    # 3x3 kernel that sums to 1.0
    return np.array([[0.0, -0.1, 0.0], [0.1, 1.0, -0.1], [0.0, 0.1, 0.0]], dtype=np.float32)


@pytest.fixture
def simple_kernel_4d(simple_kernel_2d):
    # 3x3 kernel repeated over 2048x2048 spatial grid
    kernel = np.zeros((3, 3, 2048, 2048), dtype=np.float32)
    kernel[:, :, :, :] = simple_kernel_2d[:, :, None, None]
    return kernel


@pytest.fixture
def miri_full_science_datamodel(science_data):
    miri_model = datamodels.RampModel()
    xsize = 1032
    ysize = 1024
    nints = 1
    ngroups = 5
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    data[:, :, :, :] = science_data[None, None, :ysize, :xsize]
    miri_model.data = data
    miri_model.meta.subarray.xstart = 1
    miri_model.meta.subarray.ystart = 1
    miri_model.meta.subarray.xsize = xsize
    miri_model.meta.subarray.ysize = ysize
    miri_model.meta.instrument.name = "MIRI"
    miri_model.meta.exposure.nints = nints
    miri_model.meta.exposure.ngroups = ngroups
    miri_model.meta.exposure.groupgap = 0
    miri_model.meta.exposure.nframes = 1
    miri_model.meta.instrument.detector = "MIRIMAGE"
    miri_model.meta.observation.date = "2022-06-01"
    miri_model.meta.observation.time = "00:00:00"
    return miri_model


@pytest.fixture
def miri_subarray_science_datamodel(science_data):
    miri_model = datamodels.RampModel()
    xsize = 100
    ysize = 100
    nints = 1
    ngroups = 6
    ystart = 924
    xstart = 932
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    data[:, :, :, :] = science_data[None, None, ystart : ystart + ysize, xstart : xstart + xsize]
    miri_model.data = data
    miri_model.meta.subarray.xstart = xstart + 1
    miri_model.meta.subarray.ystart = ystart + 1
    miri_model.meta.subarray.xsize = xsize
    miri_model.meta.subarray.ysize = ysize
    miri_model.meta.instrument.name = "MIRI"
    miri_model.meta.exposure.nints = nints
    miri_model.meta.exposure.ngroups = ngroups
    miri_model.meta.exposure.groupgap = 0
    miri_model.meta.exposure.nframes = 1
    miri_model.meta.instrument.detector = "MIRIMAGE"
    miri_model.meta.observation.date = "2022-03-01"
    miri_model.meta.observation.time = "00:00:00"
    return miri_model


@pytest.fixture
def nir_full_science_datamodel(science_data):
    nir_model = datamodels.RampModel()
    xsize = 2048
    ysize = 2048
    nints = 3
    ngroups = 4
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    data[:, :, :, :] = science_data[None, None, :ysize, :xsize]
    nir_model.data = data
    nir_model.meta.subarray.xstart = 1
    nir_model.meta.subarray.ystart = 1
    nir_model.meta.subarray.xsize = xsize
    nir_model.meta.subarray.ysize = ysize
    nir_model.meta.instrument.name = "NIRCAM"
    nir_model.meta.exposure.nints = nints
    nir_model.meta.exposure.ngroups = ngroups
    nir_model.meta.exposure.groupgap = 0
    nir_model.meta.exposure.nframes = 1
    nir_model.meta.instrument.detector = "NRCA3"
    nir_model.meta.observation.date = "2022-06-01"
    nir_model.meta.observation.time = "00:00:00"
    return nir_model


@pytest.fixture
def nir_subarray_science_datamodel(science_data):
    nir_model = datamodels.RampModel()
    xsize = 200
    ysize = 200
    nints = 2
    ngroups = 7
    ystart = 1848
    xstart = 1848
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    data[:, :, :, :] = science_data[None, None, ystart : ystart + ysize, xstart : xstart + xsize]
    nir_model.data = data
    nir_model.meta.subarray.xstart = xstart + 1
    nir_model.meta.subarray.ystart = ystart + 1
    nir_model.meta.subarray.xsize = xsize
    nir_model.meta.subarray.ysize = ysize
    nir_model.meta.instrument.name = "NIRCAM"
    nir_model.meta.exposure.nints = nints
    nir_model.meta.exposure.ngroups = ngroups
    nir_model.meta.exposure.groupgap = 0
    nir_model.meta.exposure.nframes = 1
    nir_model.meta.instrument.detector = "NRCA3"
    nir_model.meta.observation.date = "2022-06-01"
    nir_model.meta.observation.time = "00:00:00"
    return nir_model


@pytest.fixture
def irs2_science_datamodel():
    irs2_model = datamodels.RampModel()
    xsize = 2048
    ysize = 3200
    nints = 1
    ngroups = 5
    detector = "NRS1"
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32) + 200.0
    mask = x_irs2.make_mask(data)
    science_pixels = np.zeros((nints, ngroups, 2048, 2048), dtype=np.float32)
    for i in range(2048):
        for j in range(2048):
            science_pixels[:, :, i, j] = 1000.0 + 2.0 * i + 3.0 * j
    for integration in range(nints):
        for group in range(ngroups):
            x_irs2.to_irs2(
                data[integration, group], science_pixels[integration, group], mask, detector
            )
    irs2_model.data = data
    irs2_model.meta.subarray.xstart = 1
    irs2_model.meta.subarray.ystart = 1
    irs2_model.meta.subarray.xsize = 2048
    irs2_model.meta.subarray.ysize = 3200
    irs2_model.meta.instrument.name = "NIRSPEC"
    irs2_model.meta.exposure.nints = nints
    irs2_model.meta.exposure.ngroups = ngroups
    irs2_model.meta.instrument.detector = detector
    irs2_model.meta.exposure.nrs_normal = 16
    irs2_model.meta.exposure.nrs_reference = 4
    irs2_model.meta.exposure.groupgap = 0
    irs2_model.meta.exposure.nframes = 1
    return irs2_model


@pytest.fixture
def ipc_model_2d(simple_kernel_2d):
    ipcmodel = datamodels.IPCModel()
    ipcmodel.data = simple_kernel_2d
    return ipcmodel


@pytest.fixture
def ipc_miri_model_4d(ipc_model_2d):
    ipcmodel = datamodels.IPCModel()
    xsize = 1032
    ysize = 1024
    kernel = np.zeros((3, 3, ysize, xsize), dtype=np.float32)
    kernel[:, :, :, :] = ipc_model_2d.data[:, :, None, None]
    ipcmodel.data = kernel
    return ipcmodel


@pytest.fixture
def ipc_nir_model_4d(ipc_model_2d):
    ipcmodel = datamodels.IPCModel()
    xsize = 2048
    ysize = 2048
    kernel = np.zeros((3, 3, ysize, xsize), dtype=np.float32)
    kernel[:, :, :, :] = ipc_model_2d.data[:, :, None, None]
    ipcmodel.data = kernel
    return ipcmodel


def test_do_correction(miri_full_science_datamodel, ipc_model_2d):
    input_model_argument = miri_full_science_datamodel.copy()
    result = ipc_corr.do_correction(input_model_argument, ipc_model_2d)
    diff = result.data[:, :, 5:95, 10:90] - miri_full_science_datamodel.data[:, :, 5:95, 10:90]
    assert np.allclose(diff, 0.2, rtol=1.0e-3)


@pytest.mark.parametrize(
    "input_model, ipc_model",
    [
        ("miri_full_science_datamodel", "ipc_model_2d"),
        ("miri_full_science_datamodel", "ipc_miri_model_4d"),
        ("miri_subarray_science_datamodel", "ipc_miri_model_4d"),
        ("nir_full_science_datamodel", "ipc_model_2d"),
        ("irs2_science_datamodel", "ipc_model_2d"),
    ],
)
def test_ipc_correction(input_model, ipc_model, request):
    used_input_model = request.getfixturevalue(input_model)
    input_model_argument = used_input_model.copy()
    used_ipc_model = request.getfixturevalue(ipc_model)
    result = ipc_corr.ipc_correction(input_model_argument, used_ipc_model)
    analysis_slice = (slice(None), slice(None), slice(5, 95), slice(5, 95))
    diff = result.data[analysis_slice] - used_input_model.data[analysis_slice]
    if input_model != "irs2_science_datamodel":
        assert np.allclose(diff, 0.2, rtol=1.0e-3)
    else:
        detector = used_input_model.meta.instrument.detector
        mask = x_irs2.make_mask(used_input_model)
        input_science_pixels = x_irs2.from_irs2(used_input_model.data, mask, detector)
        result_science_pixels = x_irs2.from_irs2(result.data, mask, detector)
        diff = result_science_pixels[:, :, 5:95, 10:90] - input_science_pixels[:, :, 5:95, 10:90]
        assert np.allclose(diff, 0.2, rtol=1.0e-3)


@pytest.mark.parametrize(
    "input_model, expected",
    [
        ("miri_full_science_datamodel", Nref(0, 0, 4, 4)),
        ("nir_full_science_datamodel", Nref(4, 4, 4, 4)),
        ("miri_subarray_science_datamodel", Nref(0, 0, 0, 4)),
        ("nir_subarray_science_datamodel", Nref(0, 4, 0, 4)),
    ],
)
def test_get_num_ref_pixels(input_model, expected, request):
    used_input_model = request.getfixturevalue(input_model)
    nref = ipc_corr.get_num_ref_pixels(used_input_model)
    assert nref == expected


@pytest.mark.parametrize(
    "input_model",
    [
        "miri_full_science_datamodel",
        "nir_full_science_datamodel",
        "miri_subarray_science_datamodel",
        "nir_subarray_science_datamodel",
    ],
)
def test_get_ipc_slice_2d(input_model, ipc_model_2d, request):
    used_input_model = request.getfixturevalue(input_model)
    returned_kernel = ipc_corr.get_ipc_slice(input_model, ipc_model_2d.copy())
    # If the ipc model is 2d, return the 2d ipc model
    assert np.allclose(returned_kernel.data, ipc_model_2d.data, rtol=1.0e-7)


@pytest.mark.parametrize(
    "input_model, ipc_model",
    [
        ("miri_full_science_datamodel", "ipc_miri_model_4d"),
        ("nir_full_science_datamodel", "ipc_nir_model_4d"),
        ("miri_subarray_science_datamodel", "ipc_miri_model_4d"),
        ("nir_subarray_science_datamodel", "ipc_nir_model_4d"),
    ],
)
def test_get_ipc_slice_4d(input_model, ipc_model, request):
    used_input_model = request.getfixturevalue(input_model)
    used_ipc_model = request.getfixturevalue(ipc_model)
    returned_kernel = ipc_corr.get_ipc_slice(used_input_model, used_ipc_model)
    assert returned_kernel.shape[-2:] == used_input_model.shape[-2:]


def test_ipc_convolve_2d_with_reference(science_data, simple_kernel_2d):
    input_data = science_data.copy()
    ref_pixel_config = Nref(4, 4, 4, 4)
    ipc_corr.ipc_convolve(input_data, simple_kernel_2d, ref_pixel_config)

    # Check reference pixels remain the same
    br, tr, lc, rc = ref_pixel_config
    assert np.array_equal(input_data[:tr, :], science_data[:tr, :])
    assert np.array_equal(input_data[-br:, :], science_data[-br:, :])
    assert np.array_equal(input_data[:, :lc], science_data[:, :lc])
    assert np.array_equal(input_data[:, -rc:], science_data[:, -rc:])

    diffs = input_data - science_data
    # Check science data changed as expected
    assert np.allclose(diffs[5:95, 5:95], 0.2, rtol=3.0e-4)


def test_ipc_convolve_4d_with_reference(science_data, simple_kernel_4d):
    input_data = science_data.copy()
    ref_pixel_config = Nref(4, 4, 4, 4)
    ipc_corr.ipc_convolve(input_data, simple_kernel_4d, ref_pixel_config)

    # Check reference pixels unchanged
    br, tr, lc, rc = ref_pixel_config
    assert np.array_equal(input_data[:tr, :], science_data[:tr, :])
    assert np.array_equal(input_data[-br:, :], science_data[-br:, :])
    assert np.array_equal(input_data[:, :lc], science_data[:, :lc])
    assert np.array_equal(input_data[:, -rc:], science_data[:, -rc:])

    diffs = input_data - science_data
    assert np.allclose(diffs[5:2040, 5:2040], 0.2, rtol=3.0e-3)


def test_ipc_convolve_no_reference(science_data, simple_kernel_2d):
    input_data = science_data.copy()
    nref = Nref(0, 0, 0, 0)
    ipc_corr.ipc_convolve(input_data, simple_kernel_2d, nref)

    # All values should be 0.2 apart from the 1 pixel border
    diffs = input_data - science_data
    assert np.allclose(diffs[1:-1, 1:-1], 0.2, rtol=3.0e-3)


@pytest.mark.parametrize(
    "input_model",
    ["miri_full_science_datamodel", "nir_full_science_datamodel", "nir_subarray_science_datamodel"],
)
def test_ipc_step_run(input_model, request):
    step_instance = IPCStep()
    used_input_model = request.getfixturevalue(input_model)
    result = IPCStep.call(used_input_model)
    assert result.meta.cal_step.ipc == "COMPLETE"


def mock_crds_getref(*args, **kwargs):
    return "N/A"


def test_ipc_step_skip(miri_subarray_science_datamodel, monkeypatch):
    # This should fail the CRDS lookup, so would throw a CrdsLookupError and exit
    # Monkeypatch the CRDS lookup so we test the return for CRDS returning "N/A"
    monkeypatch.setattr("stpipe.crds_client.get_reference_file", mock_crds_getref)
    result = IPCStep.call(miri_subarray_science_datamodel)
    assert result.meta.cal_step.ipc == "SKIPPED"

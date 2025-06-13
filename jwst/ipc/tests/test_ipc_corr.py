import numpy as np
import pytest
from collections import namedtuple

from stdatamodels.jwst import datamodels
from jwst.ipc import ipc_corr

# Define the nref namedtuple
Nref = namedtuple("Nref", ["bottom_rows", "top_rows", "left_columns", "right_columns"])

@pytest.fixture
def science_data():
    data = np.zeros((100, 100), dtype=np.float32)
    for i in range(100):
        for j in range(100):
            data[i, j] = 1000.0 + 2.0 * i + 3.0 * j
    return data

@pytest.fixture
def simple_kernel_2d():
    # 3x3 kernel that sums to 1.0
    return np.array([
        [ 0.0,  -0.1,  0.0],
        [ 0.1,   1.0, -0.1],
        [ 0.0,   0.1,  0.0]
    ], dtype=np.float32)

@pytest.fixture
def simple_kernel_4d():
    # 3x3 kernel repeated over 100x100 spatial grid
    base_kernel = np.array([
        [ 0.0,  -0.1,  0.0],
        [ 0.1,   1.0, -0.1],
        [ 0.0,   0.1,  0.0]
    ], dtype=np.float32)
    kernel = np.zeros((3, 3, 100, 100), dtype=np.float32)
    for j in range(3):
        for i in range(3):
            kernel[j, i, :, :] = base_kernel[j, i]
    return kernel

@pytest.fixture
def miri_full_science_datamodel():
    miri_model = datamodels.RampModel()
    xsize = 1032
    ysize = 1024
    nints = 1
    ngroups = 5
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    for i in range(ysize):
        for j in range(xsize):
            data[:, :, i, j] = 1000.0 + 2.0 * i + 3.0 * j
    miri_model.data = data
    miri_model.meta.subarray.xstart = 1
    miri_model.meta.subarray.ystart = 1
    miri_model.meta.subarray.xsize = xsize
    miri_model.meta.subarray.ysize = ysize
    miri_model.meta.instrument.name = 'MIRI'
    miri_model.meta.exposure.nints = nints
    miri_model.meta.exposure.ngroups = ngroups
    miri_model.meta.exposure.groupgap = 0
    miri_model.meta.exposure.nframes = 1
    return miri_model

@pytest.fixture
def miri_subarray_science_datamodel():
    miri_model = datamodels.RampModel()
    xsize = 100
    ysize = 100
    nints = 1
    ngroups = 6
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    for i in range(ysize):
        for j in range(xsize):
            data[:, :, i, j] = 1000.0 + 2.0 * i + 3.0 * j
    miri_model.data = data
    miri_model.meta.subarray.xstart = 933
    miri_model.meta.subarray.ystart = 925
    miri_model.meta.subarray.xsize = xsize
    miri_model.meta.subarray.ysize = ysize
    miri_model.meta.instrument.name = 'MIRI'
    miri_model.meta.exposure.nints = nints
    miri_model.meta.exposure.ngroups = ngroups
    miri_model.meta.exposure.groupgap = 0
    miri_model.meta.exposure.nframes = 1
    return miri_model

@pytest.fixture
def nir_full_science_datamodel():
    nir_model = datamodels.RampModel()
    xsize = 2048
    ysize = 2048
    nints = 3
    ngroups = 4
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    for i in range(ysize):
        for j in range(xsize):
            data[:, :, i, j] = 1000.0 + 2.0 * i + 3.0 * j
    nir_model.data = data
    nir_model.meta.subarray.xstart = 1
    nir_model.meta.subarray.ystart = 1
    nir_model.meta.subarray.xsize = xsize
    nir_model.meta.subarray.ysize = ysize
    nir_model.meta.instrument.name = 'NIRCAM'
    nir_model.meta.exposure.nints = nints
    nir_model.meta.exposure.ngroups = ngroups
    nir_model.meta.exposure.groupgap = 0
    nir_model.meta.exposure.nframes = 1
    return nir_model

@pytest.fixture
def nir_subarray_science_datamodel():
    nir_model = datamodels.RampModel()
    xsize = 200
    ysize = 200
    nints = 2
    ngroups = 7
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    for i in range(ysize):
        for j in range(xsize):
            data[:, :, i, j] = 1000.0 + 2.0 * i + 3.0 * j
    nir_model.data = data
    nir_model.meta.subarray.xstart = 1849
    nir_model.meta.subarray.ystart = 1849
    nir_model.meta.subarray.xsize = xsize
    nir_model.meta.subarray.ysize = ysize
    nir_model.meta.instrument.name = 'NIRCAM'
    nir_model.meta.exposure.nints = nints
    nir_model.meta.exposure.ngroups = ngroups
    nir_model.meta.exposure.groupgap = 0
    nir_model.meta.exposure.nframes = 1
    return nir_model

@pytest.fixture
def irs2_science_datamodel():
    irs2_model = datamodels.RampModel()
    xsize = 2048
    ysize = 3200
    nints = 1
    ngroups = 5
    data = np.zeros((nints, ngroups, ysize, xsize), dtype=np.float32)
    for i in range(ysize):
        for j in range(xsize):
            data[:, :, i, j] = 1000.0 + 2.0 * i + 3.0 * j
    irs2_model.data = data
    irs2_model.meta.subarray.xstart = 1
    irs2_model.meta.subarray.ystart = 1
    irs2_model.meta.subarray.xsize = 2048
    irs2_model.meta.subarray.ysize = 3200
    irs2_model.meta.instrument.name = 'NIRSPEC'
    irs2_model.meta.exposure.nints = nints
    irs2_model.meta.exposure.ngroups = ngroups
    irs2_model.meta.instrument.detector = 'NRS1'
    irs2_model.meta.exposure.nrs_normal = 16
    irs2_model.meta.exposure.nrs_reference = 4
    irs2_model.meta.exposure.groupgap = 0
    irs2_model.meta.exposure.nframes = 1
    return irs2_model

@pytest.fixture
def ipc_model_2d():
    ipcmodel = datamodels.IPCModel()
    ipcmodel.data = np.array([
        [ 0.0,  -0.1,  0.0],
        [ 0.1,   1.0, -0.1],
        [ 0.0,   0.1,  0.0]
    ], dtype=np.float32)
    return ipcmodel

@pytest.fixture
def ipc_miri_model_4d():
    ipcmodel = datamodels.IPCModel()
    base_kernel = np.array([
        [ 0.0,  -0.1,  0.0],
        [ 0.1,   1.0, -0.1],
        [ 0.0,   0.1,  0.0]
    ], dtype=np.float32)
    xsize = 1032
    ysize = 1024
    kernel = np.zeros((3, 3, ysize, xsize), dtype=np.float32)
    for j in range(3):
        for i in range(3):
            kernel[j, i, :, :] = base_kernel[j, i]
    ipcmodel.data = kernel
    return ipcmodel

@pytest.fixture
def ipc_nir_model_4d():
    ipcmodel = datamodels.IPCModel()
    base_kernel = np.array([
        [ 0.0,  -0.1,  0.0],
        [ 0.1,   1.0, -0.1],
        [ 0.0,   0.1,  0.0]
    ], dtype=np.float32)
    xsize = 2048
    ysize = 2048
    kernel = np.zeros((3, 3, ysize, xsize), dtype=np.float32)
    for j in range(3):
        for i in range(3):
            kernel[j, i, :, :] = base_kernel[j, i]
    ipcmodel.data = kernel
    return ipcmodel

@pytest.mark.parametrize('input_model, ipc_model',
                         [('miri_full_science_datamodel', 'ipc_model_2d'),
                          ('miri_full_science_datamodel', 'ipc_miri_model_4d'),
                          ('miri_subarray_science_datamodel', 'ipc_miri_model_4d'),
                          ('nir_full_science_datamodel', 'ipc_model_2d'),
                          ('irs2_science_datamodel', 'ipc_model_2d')
                         ]
                        )
def test_do_correction(input_model, ipc_model, request):
    used_input_model = request.getfixturevalue(input_model)
    input_model_argument = used_input_model.copy()
    used_ipc_model = request.getfixturevalue(ipc_model)
    result = ipc_corr.do_correction(input_model_argument, used_ipc_model)
    diff = result.data[:, :, 5:95, 10:90] - used_input_model.data[:, :, 5:95, 10:90]
    assert result is input_model_argument
    if input_model != 'irs2_science_datamodel':
        assert np.allclose(diff, 0.2, rtol=1.0e-3)

@pytest.mark.parametrize('input_model, ipc_model',
                         [('miri_full_science_datamodel', 'ipc_model_2d'),
                          ('miri_full_science_datamodel', 'ipc_miri_model_4d'),
                          ('miri_subarray_science_datamodel', 'ipc_miri_model_4d'),
                          ('nir_full_science_datamodel', 'ipc_model_2d'),
                          ('irs2_science_datamodel', 'ipc_model_2d')
                         ]
                        )
def test_ipc_correction(input_model, ipc_model, request):
    used_input_model = request.getfixturevalue(input_model)
    input_model_argument = used_input_model.copy()
    used_ipc_model = request.getfixturevalue(ipc_model)
    result = ipc_corr.ipc_correction(input_model_argument, used_ipc_model)
    analysis_slice = (slice(None), slice(None), slice(5, 95), slice(5, 95))
    diff = result.data[analysis_slice] - used_input_model.data[analysis_slice]
    assert result is input_model_argument
    if input_model != 'irs2_science_datamodel':
        assert np.allclose(diff, 0.2, rtol=1.0e-3)

@pytest.mark.parametrize('input_model, expected',
                         [('miri_full_science_datamodel', Nref(0, 0, 4, 4)),
                          ('nir_full_science_datamodel', Nref(4, 4, 4, 4)),
                          ('miri_subarray_science_datamodel', Nref(0, 0, 0, 4)),
                          ('nir_subarray_science_datamodel', Nref(0, 4, 0, 4))
                         ]
                        )
def test_get_num_ref_pixels(input_model, expected, request):
    used_input_model = request.getfixturevalue(input_model)
    nref = ipc_corr.get_num_ref_pixels(used_input_model)
    assert nref == expected

@pytest.mark.parametrize('input_model', ['miri_full_science_datamodel',
                                         'nir_full_science_datamodel',
                                         'miri_subarray_science_datamodel',
                                         'nir_subarray_science_datamodel'])
def test_get_ipc_slice_2d(input_model, ipc_model_2d, request):
    used_input_model = request.getfixturevalue(input_model)
    returned_kernel = ipc_corr.get_ipc_slice(input_model, ipc_model_2d)
    assert returned_kernel is ipc_model_2d.data

@pytest.mark.parametrize('input_model, ipc_model',
                         [("miri_full_science_datamodel", "ipc_miri_model_4d"),
                          ("nir_full_science_datamodel", "ipc_nir_model_4d"),
                          ("miri_subarray_science_datamodel", "ipc_miri_model_4d"),
                          ("nir_subarray_science_datamodel", "ipc_nir_model_4d")])
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
    assert np.allclose(diffs[5:95, 5:95], 0.2, rtol=3.0e-4)

def test_ipc_convolve_no_reference(science_data, simple_kernel_2d):
    input_data = science_data.copy()
    nref = Nref(0, 0, 0, 0)
    ipc_corr.ipc_convolve(input_data, simple_kernel_2d, nref)

    # All values should be 0.2 apart from the 1 pixel border
    diffs = input_data - science_data
    assert np.allclose(diffs[1:-1, 1:-1], 0.2, rtol=3.0e-4)

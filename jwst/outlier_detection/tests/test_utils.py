import pytest
from jwst.outlier_detection.utils import DiskAppendableArray, OnDiskTimewiseOperation
from pathlib import Path
import numpy as np
import os


def test_disk_appendable_array(tmp_cwd):

    slice_shape = (8,7)
    dtype = 'float32'
    tempdir = tmp_cwd / Path("tmptest")
    os.mkdir(tempdir)

    arr = DiskAppendableArray(slice_shape, dtype=dtype, tempdir=tempdir)

    # check temporary file setup
    assert arr._filename.split("/")[-1] in os.listdir(tempdir)
    assert len(os.listdir(tempdir)) == 1
    assert arr.shape == (0,) + slice_shape

    # check cwd contains no files
    assert all([not os.path.isfile(f) for f in os.listdir(tmp_cwd)])

    # check expected append failures
    with pytest.raises(ValueError):
        candidate = np.zeros((7,7), dtype=dtype)
        arr.append(candidate)
    with pytest.raises(ValueError):
        candidate = np.zeros((8,7), dtype='float64')
        arr.append(candidate)

    # check append and read
    candidate0 = np.zeros(slice_shape, dtype=dtype)
    candidate1 = np.full(slice_shape, 2, dtype=dtype)
    candidate2 = np.full(slice_shape, np.nan, dtype=dtype)
    for candidate in [candidate0, candidate1, candidate2]:
        arr.append(candidate)
    
    arr_in_memory = arr.read()

    assert arr_in_memory.shape == (3,) + slice_shape
    assert np.all(arr_in_memory[0] == candidate0)
    assert np.all(arr_in_memory[1] == candidate1)
    assert np.allclose(arr_in_memory[2], candidate2, equal_nan=True)

    # check cleanup occurs after reassigning arr
    arr = None
    assert len(os.listdir(tempdir)) == 0



def test_on_disk_median(tmp_cwd):

    library_length = 3
    frame_shape = (21,20)
    dtype = 'float32'
    tempdir = tmp_cwd / Path("tmptest")
    os.mkdir(tempdir)
    shape = (library_length,) + frame_shape

    median_computer = OnDiskTimewiseOperation(shape, dtype=dtype, tempdir=tempdir)

    # test compute buffer indices
    # buffer size equals size of single input model by default
    # which means we expect same number of sections as library length
    # in reality there is often one more section than that because
    # of necessity of integer number of rows per section, but math is exact in this case
    expected_buffer_size = frame_shape[0] * frame_shape[1] * np.dtype(dtype).itemsize
    expected_section_nrows = frame_shape[0] // library_length
    assert median_computer.nsections == library_length
    assert median_computer.section_nrows == expected_section_nrows
    assert median_computer.buffer_size == expected_buffer_size

    # test temp file setup
    assert len(os.listdir(tempdir)) == 1
    assert str(median_computer._temp_path).startswith(str(tempdir))
    assert len(os.listdir(median_computer._temp_path)) == library_length
    # check cwd and parent tempdir contain no files
    assert all([not os.path.isfile(f) for f in os.listdir(tmp_cwd)])
    assert all([not os.path.isfile(f) for f in os.listdir(tempdir)])
    
    # test validate data
    with pytest.raises(ValueError):
        candidate = np.zeros((20,20), dtype=dtype)
        median_computer.add_image(candidate)
    with pytest.raises(ValueError):
        candidate = np.zeros((21,20), dtype='float64')
        median_computer.add_image(candidate)

    # test add and compute
    candidate0 = np.full(frame_shape, 3, dtype=dtype)
    candidate1 = np.full(frame_shape, 2, dtype=dtype)
    candidate2 = np.full(frame_shape, np.nan, dtype=dtype)
    for candidate in [candidate0, candidate1, candidate2]:
        median_computer.add_image(candidate)
    median = median_computer.compute_median()
    assert median.shape == frame_shape
    assert np.all(median == 2.5)

    # test expected error trying to add too many frames
    # for loop to ensure always happens, not just the first time
    candidate3 = np.zeros_like(candidate0)
    for _ in range(2):
        with pytest.raises(IndexError):
            median_computer.add_image(candidate3)

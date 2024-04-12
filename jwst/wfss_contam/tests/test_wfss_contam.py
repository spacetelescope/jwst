import pytest
from jwst.wfss_contam.wfss_contam import _determine_multiprocessing_ncores, _cut_frame_to_match_slit, build_common_slit
from jwst.datamodels import SlitModel
import numpy as np


@pytest.mark.parametrize("max_cores, num_cores, expected", 
                         [("none", 4, 1), 
                          ("quarter", 4, 1), 
                          ("half", 4, 2), 
                          ("all", 4, 4), 
                          ("none", 1, 1),])
def test_determine_multiprocessing_ncores(max_cores, num_cores, expected):
    assert _determine_multiprocessing_ncores(max_cores, num_cores) == expected  


@pytest.fixture(scope="module")
def contam():
    return np.ones((10, 10))*0.1

@pytest.fixture(scope="module")
def slit0():
    slit = SlitModel(data=np.ones((5, 3)))
    slit.xstart = 2
    slit.ystart = 3
    slit.xsize = 3
    slit.ysize = 5
    return slit


@pytest.fixture(scope="module")
def slit1():
    slit = SlitModel(data=np.ones((4, 4))*0.5)
    slit.xstart = 3
    slit.ystart = 2
    slit.xsize = 4
    slit.ysize = 4
    return slit


def test_cut_frame_to_match_slit(slit0, contam):
    cut_contam = _cut_frame_to_match_slit(contam, slit0)
    assert cut_contam.shape == (5, 3)
    assert np.all(cut_contam == 0.1)


def test_build_common_slit(slit0, slit1):
    slit0, slit1 = build_common_slit(slit0, slit1)

    # check indexing in metadata
    assert slit0.xstart == slit1.xstart
    assert slit0.ystart == slit1.ystart
    assert slit0.xsize == slit1.xsize
    assert slit0.ysize == slit1.ysize
    assert slit0.data.shape == slit1.data.shape

    # check data overlap
    assert np.count_nonzero(slit0.data) == 15
    assert np.count_nonzero(slit1.data) == 16
    assert np.count_nonzero(slit0.data * slit1.data) == 6

    # check data values
    assert np.all(slit0.data[1:6, 0:3] == 1)
    assert np.all(slit1.data[0:4, 1:5] == 0.5)

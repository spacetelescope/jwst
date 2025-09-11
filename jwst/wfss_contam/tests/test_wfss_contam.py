import numpy as np
import pytest
from stdatamodels.jwst.datamodels import SlitModel

from jwst.wfss_contam.wfss_contam import (
    SlitOverlapError,
    UnmatchedSlitIDError,
    _cut_frame_to_match_slit,
    _find_matching_simul_slit,
    determine_multiprocessing_ncores,
    match_backplane_encompass_both,
    match_backplane_prefer_first,
)


@pytest.mark.parametrize(
    "max_cores, num_cores, expected",
    [
        ("none", 4, 1),
        ("quarter", 4, 1),
        ("half", 4, 2),
        ("all", 4, 4),
        ("none", 1, 1),
        (
            None,
            1,
            1,
        ),
        (3, 5, 3),
        (100, 5, 5),
    ],
)
def test_determine_multiprocessing_ncores(max_cores, num_cores, expected):
    assert determine_multiprocessing_ncores(max_cores, num_cores) == expected


@pytest.fixture(scope="module")
def contam():
    return np.ones((10, 10)) * 0.1


@pytest.fixture(scope="module")
def slit0():
    slit = SlitModel(data=np.ones((5, 3)))
    slit.xstart = 2
    slit.ystart = 3
    slit.xsize = 3
    slit.ysize = 5
    slit.meta.wcsinfo.spectral_order = 1
    slit.source_id = 1
    return slit


@pytest.fixture(scope="module")
def slit1():
    slit = SlitModel(data=np.ones((4, 4)) * 0.5)
    slit.xstart = 3
    slit.ystart = 2
    slit.xsize = 4
    slit.ysize = 4
    return slit


@pytest.fixture(scope="module")
def slit2():
    slit = SlitModel(data=np.ones((3, 5)) * 0.1)
    slit.xstart = 300
    slit.ystart = 200
    slit.xsize = 5
    slit.ysize = 3
    return slit


def test_find_matching_simul_slit(slit0):
    sids = [0, 1, 1]
    orders = [1, 1, 2]
    idx = _find_matching_simul_slit(slit0, sids, orders)
    assert idx == 1


def test_find_matching_simul_slit_no_match(slit0):
    sids = [0, 1, 1]
    orders = [1, 2, 2]
    with pytest.raises(UnmatchedSlitIDError):
        _find_matching_simul_slit(slit0, sids, orders)


def test_cut_frame_to_match_slit(slit0, contam):
    cut_contam = _cut_frame_to_match_slit(contam, slit0)
    assert cut_contam.shape == (5, 3)
    assert np.all(cut_contam == 0.1)


def test_common_slit_encompass(slit0, slit1):
    slit0_final, slit1_final = match_backplane_encompass_both(slit0.copy(), slit1.copy())

    # check indexing in metadata
    assert slit0_final.xstart == slit1_final.xstart
    assert slit0_final.ystart == slit1_final.ystart
    assert slit0_final.xsize == slit1_final.xsize
    assert slit0_final.ysize == slit1_final.ysize
    assert slit0_final.data.shape == slit1_final.data.shape

    # check data overlap
    assert np.count_nonzero(slit0_final.data) == 15
    assert np.count_nonzero(slit1_final.data) == 16
    assert np.count_nonzero(slit0_final.data * slit1_final.data) == 6

    # check data values
    assert np.all(slit0_final.data[1:6, 0:3] == 1)
    assert np.all(slit1_final.data[0:4, 1:5] == 0.5)


def test_common_slit_prefer(slit0, slit1):
    slit0_final, slit1_final = match_backplane_prefer_first(slit0.copy(), slit1.copy())
    assert slit0_final.xstart == slit0.xstart
    assert slit0_final.ystart == slit0.ystart
    assert slit0_final.xsize == slit0.xsize
    assert slit0_final.ysize == slit0.ysize
    assert slit0_final.data.shape == slit0.data.shape
    assert np.all(slit0_final.data == slit0.data)

    assert slit1_final.xstart == slit0.xstart
    assert slit1_final.ystart == slit0.ystart
    assert slit1_final.xsize == slit0.xsize
    assert slit1_final.ysize == slit0.ysize
    assert slit1_final.data.shape == slit0.data.shape
    assert np.count_nonzero(slit1_final.data) == 6


def test_common_slit_prefer_expected_raise(slit0, slit2):
    with pytest.raises(SlitOverlapError):
        match_backplane_prefer_first(slit0.copy(), slit2.copy())

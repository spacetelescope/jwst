import numpy as np
import pytest
from stdatamodels.jwst.datamodels import ImageModel, RampModel, ReadnoiseModel, WfssBkgModel

from jwst.lib.reffile_utils import (
    detector_science_frame_transform,
    find_row,
    get_subarray_model,
    science_detector_frame_transform,
)


def generate_test_refmodel_metadata(refmodel):
    refmodel.meta.description = "filler"
    refmodel.meta.reftype = "filler"
    refmodel.meta.author = "Py Test"
    refmodel.meta.pedigree = "Pytest"
    refmodel.meta.useafter = "2015-01-01T01:00:00"


def test_find_row():
    filters = [
        {"column_offset": 1.0, "filter": "F277W", "pupil": "FLAT", "row_offset": 2.0},
        {"column_offset": 0.0, "filter": "F356W", "pupil": "FLAT", "row_offset": 0.0},
    ]
    match_keys = {"filter": "F277W", "pupil": "FLAT"}
    missing_key = {"filter": "F277H", "pupil": "FLAT"}

    result = find_row(filters, match_keys)
    assert result == {"column_offset": 1.0, "filter": "F277W", "pupil": "FLAT", "row_offset": 2.0}

    result = find_row(filters, missing_key)
    assert result is None


def test_multistripe_subarray_model():
    mock_rn = ReadnoiseModel(data=(np.ones((2048, 2048), dtype=int) * np.arange(2048)).T)
    mock_rn.meta.instrument.name = "NIRCAM"
    generate_test_refmodel_metadata(mock_rn)
    mock_sci = RampModel(data=np.ones((5, 5, 164, 2048)))
    mock_sci.meta.subarray = {
        "fastaxis": 1,
        "name": "SUB164STRIPE4_DHS",
        "slowaxis": -2,
        "xsize": 2048,
        "xstart": 1,
        "ysize": 164,
        "ystart": 1885,
        "multistripe_reads1": 1,
        "multistripe_skips1": 1549,
        "multistripe_reads2": 40,
        "multistripe_skips2": 82,
        "repeat_stripe": 1,
        "interleave_reads1": 1,
        "superstripe_step": 0,
        "num_superstripe": 0,
    }
    mock_rn_cutout = get_subarray_model(mock_sci, mock_rn)
    assert mock_rn_cutout.data.shape == mock_sci.shape[-2:]


@pytest.fixture
def mock_sci():
    """
    Create a mock science model in a subarray.

    Returns
    -------
    sci : ImageModel
        The mock science datamodel
    """
    sci = ImageModel(data=np.ones((64, 2048)), dq=np.zeros((64, 2048), dtype=int))
    sci.meta.subarray.xstart = 1
    sci.meta.subarray.ystart = 1985
    sci.meta.subarray.xsize = 2048
    sci.meta.subarray.ysize = 64
    generate_test_refmodel_metadata(sci)
    sci.meta.instrument.name = "NIRISS"
    return sci


@pytest.fixture
def mock_bkg():
    """
    Create a mock background model that is full detector size.

    Returns
    -------
    bkg : WfssBkgModel
        The mock background datamodel
    """
    bkg = WfssBkgModel(data=np.ones((2048, 2048)))
    bkg.meta.subarray.xstart = 1
    bkg.meta.subarray.ystart = 1
    bkg.meta.subarray.xsize = 2048
    bkg.meta.subarray.ysize = 2048
    generate_test_refmodel_metadata(bkg)
    bkg.meta.instrument.name = "NIRISS"
    return bkg


def test_get_subarray_model(mock_sci, mock_bkg):
    """Test that get_subarray_model returns a subarray model with correct shape and metadata."""
    sub_model = get_subarray_model(mock_sci, mock_bkg)
    assert sub_model.data.shape == (64, 2048)
    assert sub_model.dq.shape == (64, 2048)
    assert sub_model.meta.instrument.name == "NIRISS"


def test_get_subarray_model_typeerror(mock_sci, mock_bkg):
    """Test that TypeError is raised when non-model types are passed."""
    with pytest.raises(TypeError):
        get_subarray_model(mock_sci, "not_a_model")

    with pytest.raises(TypeError):
        get_subarray_model("not_a_model", mock_bkg)


@pytest.mark.parametrize(
    "fastaxis, slowaxis, result",
    [
        (1, 2, np.array([[1, 2], [3, 4]])),
        (1, -2, np.array([[3, 4], [1, 2]])),
        (-1, 2, np.array([[2, 1], [4, 3]])),
        (-1, -2, np.array([[4, 3], [2, 1]])),
        (2, 1, np.array([[1, 3], [2, 4]])),
        (2, -1, np.array([[2, 4], [1, 3]])),
        (-2, 1, np.array([[3, 1], [4, 2]])),
        (-2, -1, np.array([[4, 2], [3, 1]])),
    ],
)
def test_science_detector_frame_transform(fastaxis, slowaxis, result):
    detector_array = np.array([[1, 2], [3, 4]])
    returned = science_detector_frame_transform(detector_array.copy(), fastaxis, slowaxis)
    assert np.allclose(returned, result)


@pytest.mark.parametrize(
    "fastaxis, slowaxis, input",
    [
        (1, 2, np.array([[1, 2], [3, 4]])),
        (1, -2, np.array([[3, 4], [1, 2]])),
        (-1, 2, np.array([[2, 1], [4, 3]])),
        (-1, -2, np.array([[4, 3], [2, 1]])),
        (2, 1, np.array([[1, 3], [2, 4]])),
        (2, -1, np.array([[2, 4], [1, 3]])),
        (-2, 1, np.array([[3, 1], [4, 2]])),
        (-2, -1, np.array([[4, 2], [3, 1]])),
    ],
)
def test_detector_science_frame_transform(fastaxis, slowaxis, input):
    detector_array = np.array([[1, 2], [3, 4]])
    returned = detector_science_frame_transform(input.copy(), fastaxis, slowaxis)
    assert np.allclose(returned, detector_array)


@pytest.mark.parametrize(
    "fastaxis, slowaxis",
    [
        (1, 2),
        (1, -2),
        (-1, 2),
        (-1, -2),
        (2, 1),
        (2, -1),
        (-2, 1),
        (-2, -1),
    ],
)
def test_roundtrip(fastaxis, slowaxis):
    test_array = np.arange(100 * 100, dtype=np.float32).reshape((100, 100))
    forward = science_detector_frame_transform(test_array.copy(), fastaxis, slowaxis)
    reverse = detector_science_frame_transform(forward.copy(), fastaxis, slowaxis)
    assert np.allclose(reverse, test_array)

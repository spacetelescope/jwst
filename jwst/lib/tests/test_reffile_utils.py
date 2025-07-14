import pytest
import numpy as np
from jwst.lib.reffile_utils import (
    find_row,
    generate_stripe_array,
    get_subarray_model,
    science_detector_frame_transform,
)
from stdatamodels.jwst.datamodels import ReadnoiseModel, RampModel, ImageModel, WfssBkgModel


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


def test_generate_stripe():
    # Generate test array with pixel values
    # equal to row number in detector frame.
    test_array = (np.ones((2048, 2048), dtype=int) * np.arange(2048)).T

    # Use two NIRCam subarray cases to test.
    # Stripe params: xsize_sci, ysize_sci, nreads1, nreads2, nskips1,
    #                nskips2, repeat_stripe, interleave_reads1, fastaxis, slowaxis

    # SUB41STRIPE1_DHS nrca1 case
    stripe_params = (2048, 41, 1, 40, 1901, 0, 1, 1, -1, 2)

    # Function presumes input in science frame, so move test array to science frame
    # before supplying to function.
    stripe1_array = generate_stripe_array(
        science_detector_frame_transform(test_array, *stripe_params[-2:]), *stripe_params
    )
    assert stripe1_array.shape == (41, 2048)
    assert stripe1_array[0, 1024] == 0
    assert stripe1_array[1, 1024] == 1902  # nreads1 + nskips1

    # Test swapped axes
    stripe_params = (2048, 41, 1, 40, 1901, 0, 1, 1, 2, 1)
    stripe1swap_array = generate_stripe_array(
        science_detector_frame_transform(test_array, *stripe_params[-2:]), *stripe_params
    )
    assert stripe1swap_array.shape == (2048, 41)
    assert stripe1swap_array[1024, 1] == 1902  # nreads1 + nskips1

    # SUB82STRIPE2_DHS nrca2 case
    stripe_params = (2048, 82, 1, 40, 1662, 82, 1, 1, 1, -2)
    stripe2_array = generate_stripe_array(
        science_detector_frame_transform(test_array, *stripe_params[-2:]), *stripe_params
    )
    assert stripe2_array.shape == (82, 2048)
    # nrca2 has flipped row direction, so in science frame the row indices are flipped.
    assert stripe2_array[-1, 1024] == 0
    assert stripe2_array[-2, 1024] == 1663  # nreads1 + nskips1
    assert stripe2_array[-42, 1024] == 0
    assert stripe2_array[-43, 1024] == 1785  # nreads1 + nskips1 + nreads2 + nskips2

    # SUB164STRIPE4_DHS nrcalong case
    stripe_params = (2048, 164, 1, 40, 971, 0, 1, 0, -1, 2)
    stripe4_array = generate_stripe_array(
        science_detector_frame_transform(test_array, *stripe_params[-2:]), *stripe_params
    )
    assert stripe4_array.shape == (164, 2048)
    assert stripe4_array[0, 1024] == 0
    assert stripe4_array[1, 1024] == 972  # nreads1 + nskips1
    assert stripe4_array[42, 1024] == 972  # nreads1 + nskips1, stripe 2
    assert stripe4_array[83, 1024] == 972  # nreads1 + nskips1, stripe 3


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
    """Fixture to create a mock science model in a subarray."""
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
    """Fixture to create a mock background model that is full detector size."""
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

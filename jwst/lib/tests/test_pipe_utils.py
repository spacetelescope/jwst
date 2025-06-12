"""Test utilities"""

import pytest

import numpy as np

from stdatamodels.jwst import datamodels
from stdatamodels.properties import ObjectNode

from jwst.lib import pipe_utils
from jwst.associations.lib import dms_base
from jwst.lib.pipe_utils import is_tso

all_exp_types = (
    dms_base.IMAGE2_NONSCIENCE_EXP_TYPES
    + dms_base.IMAGE2_SCIENCE_EXP_TYPES
    + dms_base.SPEC2_SCIENCE_EXP_TYPES
)

exp_types = [
    (exp_type, False) for exp_type in all_exp_types if exp_type not in dms_base.TSO_EXP_TYPES
]
exp_types.extend([(exp_type, True) for exp_type in dms_base.TSO_EXP_TYPES])


@pytest.mark.parametrize("exp_type, expected", exp_types)
def test_is_tso_from_exptype(exp_type, expected):
    """Test is_tso integrity based on exp_type"""
    model = datamodels.JwstDataModel()
    model.meta.exposure.type = exp_type.upper()
    assert pipe_utils.is_tso(model) is expected


@pytest.mark.parametrize(
    "tsovisit, expected",
    [
        (True, True),
        (False, False),
    ],
)
def test_is_tso_from_tsoflag(tsovisit, expected):
    """Test is_tso integrity based on the TSO flag"""
    model = datamodels.JwstDataModel()
    model.meta.visit.tsovisit = tsovisit
    assert pipe_utils.is_tso(model) is expected


def test_is_tso_nrcgrism_nints1():
    """Test is_tso with NRC_TSGRISM and NINTS=1"""
    model = datamodels.JwstDataModel()
    model.meta.exposure.type = "NRC_TSGRISM"
    model.meta.visit.tsovisit = False
    model.meta.exposure.nints = 10

    # with NINTS>1, should be True
    assert pipe_utils.is_tso(model)

    # with NINTS=1, should be False
    model.meta.exposure.nints = 1
    assert not pipe_utils.is_tso(model)

    # with hardwired TSO EXP_TYPE's, should always be True
    assert is_tso(model) or model.meta.exposure.type.lower() in ["nrc_tsimage", "nrc_tsgrism"]


def test_is_irs2_1():
    """Test is_irs2 using a numpy array"""
    x = np.ones((2048, 2), dtype=np.float32)
    assert not pipe_utils.is_irs2(x)
    x = np.ones((3200, 2), dtype=np.float32)
    assert pipe_utils.is_irs2(x)


def test_is_irs2_2():
    """Test is_irs2 using a jwst data model"""
    x = np.ones((2048, 2), dtype=np.float32)
    model = datamodels.ImageModel(data=x)
    assert not pipe_utils.is_irs2(model)
    x = np.ones((3200, 2), dtype=np.float32)
    model = datamodels.ImageModel(data=x)
    assert pipe_utils.is_irs2(model)


@pytest.mark.parametrize("shape", [(10,), (10, 10), (10, 10, 10)])
@pytest.mark.parametrize(
    "extensions",
    [
        ("data",),
        ("data", "dq"),
        ("data", "err", "var_rnoise"),
        ("data", "err", "var_rnoise", "var_poisson", "var_flat", "dq"),
    ],
)
def test_match_nans_and_flags(shape, extensions):
    # Set up a model with matching data shapes, variable extensions
    dnu = datamodels.dqflags.pixel["DO_NOT_USE"]
    model = datamodels.SlitModel()
    invalid = np.full(shape, False)
    for i, extname in enumerate(extensions):
        if extname == "dq":
            data = np.zeros(shape, dtype=np.uint32)
            data.flat[i] = dnu
        else:
            data = np.ones(shape)
            data.flat[i] = np.nan
        setattr(model, extname, data)
        invalid.flat[i] = True

    # Match flags and NaNs across all extensions
    pipe_utils.match_nans_and_flags(model)

    # Check that all extensions have all invalid flags
    for extname in extensions:
        data = getattr(model, extname)
        if extname == "dq":
            assert np.all(data[invalid] == dnu)
            assert np.all(data[~invalid] == 0)
        else:
            assert np.all(np.isnan(data[invalid]))
            assert np.all(data[~invalid] == 1)

    model.close()


def test_match_nans_and_flags_dq_only():
    # Set up a model with only a DQ array, to test edge conditions
    dnu = datamodels.dqflags.pixel["DO_NOT_USE"]
    model = datamodels.MaskModel()
    shape = (10, 10)
    invalid = np.full(shape, False)

    model.dq = np.zeros(shape, dtype=np.uint32)
    model.dq.flat[5] = dnu
    invalid.flat[5] = True

    # Match flags and NaNs across all extensions
    dq_copy = model.dq.copy()
    pipe_utils.match_nans_and_flags(model)

    # Check that DQ extension is unchanged
    assert np.all(model.dq == dq_copy)

    model.close()


def test_match_nans_and_flags_no_data(caplog):
    # Set up a model with no data, to test edge conditions
    model = datamodels.JwstDataModel()
    model_copy = model.copy()
    pipe_utils.match_nans_and_flags(model)

    # nothing happens
    assert model == model_copy
    assert "WARNING" not in caplog.text
    assert "ERROR" not in caplog.text

    model.close()
    model_copy.close()


def test_match_nans_and_flags_not_a_model():
    # Make sure the function throws an error if input is not a model
    model = "not_a_model"
    with pytest.raises(TypeError, match="not a datamodel"):
        pipe_utils.match_nans_and_flags(model)

    # empty JwstDataModel is okay
    model = datamodels.JwstDataModel()
    pipe_utils.match_nans_and_flags(model)

    # a datamodel object node is also okay
    multislit = datamodels.MultiSlitModel()
    multislit.slits.append(datamodels.SlitModel())
    assert not isinstance(multislit.slits[0], datamodels.JwstDataModel)
    assert isinstance(multislit.slits[0], ObjectNode)
    pipe_utils.match_nans_and_flags(multislit.slits[0])


def test_match_nans_and_flags_shape_mismatch():
    # Set up a model with mismatched data shapes
    dnu = datamodels.dqflags.pixel["DO_NOT_USE"]
    model = datamodels.SlitModel()

    shape = (10, 10)
    bad_shape = (10, 5)
    invalid = np.full(shape, False)
    extensions = ["data", "err", "dq", "var_poisson"]
    for i, extname in enumerate(extensions):
        if extname == "dq":
            data = np.zeros(bad_shape, dtype=np.uint32)
            data.flat[i] = dnu
        elif extname == "var_poisson":
            data = np.ones(bad_shape, dtype=np.uint32)
            data.flat[i] = dnu
        else:
            data = np.ones(shape)
            data.flat[i] = np.nan
            invalid.flat[i] = True
        setattr(model, extname, data)

    # Match flags and NaNs across all extensions:
    # will warn for data mismatch
    model_copy = model.copy()
    pipe_utils.match_nans_and_flags(model)

    # Only data and err are updated
    for extname in extensions:
        data = getattr(model, extname)
        if extname == "data" or extname == "err":
            assert np.all(np.isnan(data[invalid]))
            assert np.all(data[~invalid] == 1)
        else:
            assert np.allclose(data, getattr(model_copy, extname))

    model.close()
    model_copy.close()

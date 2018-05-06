"""Test utilities"""
import inspect
import pytest

from .. import pipe_utils
from ... import datamodels
from ...associations.lib import dms_base

all_datamodels = [
    model
    for name, model in inspect.getmembers(datamodels, inspect.isclass)
    if issubclass(model, datamodels.DataModel)
    and model not in pipe_utils.TSO_MODEL_TYPES
]

model_list = [
    (model, True)
    for model in pipe_utils.TSO_MODEL_TYPES
]
model_list.extend([
    (model, False)
    for model in all_datamodels
])

all_exp_types = dms_base.IMAGE2_NONSCIENCE_EXP_TYPES + \
            dms_base.IMAGE2_SCIENCE_EXP_TYPES + \
            dms_base.SPEC2_SCIENCE_EXP_TYPES

exp_types = [
    (exp_type, False)
    for exp_type in all_exp_types
    if exp_type not in dms_base.TSO_EXP_TYPES
]
exp_types.extend([
    (exp_type, True)
    for exp_type in dms_base.TSO_EXP_TYPES
])


@pytest.mark.parametrize(
    'model_class, expected',
    model_list
)
def test_is_tso_from_datamodel(model_class, expected):
    """Test integrity of is_tso based on datamodels"""
    try:
        model = model_class()
    except Exception:
        # Can't generate a model. Mark as skipped
        pytest.skip(
            'Unable to generate a model from class'
            '{}'.format(model_class)
        )
    else:
        assert pipe_utils.is_tso(model) is expected


@pytest.mark.parametrize(
    'exp_type, expected',
    exp_types
)
def test_is_tso_from_exptype(exp_type, expected):
    """Test is_tso integrity based on exp_type"""
    model = datamodels.DataModel()
    model.meta.exposure.type = exp_type.upper()
    assert pipe_utils.is_tso(model) is expected


@pytest.mark.parametrize(
    'tsovisit, expected',
    [
        (True, True),
        (False, False),
    ]
)
def test_is_tso_from_tsoflag(tsovisit, expected):
    """Test is_tso integrity based on the TSO flag"""
    model = datamodels.DataModel()
    model.meta.observation.tsovisit = tsovisit
    assert pipe_utils.is_tso(model) is expected

"""Test utilities"""
import pytest

import numpy as np

from jwst.lib import pipe_utils
from jwst import datamodels
from jwst.associations.lib import dms_base
from jwst.lib.pipe_utils import is_tso

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
    'exp_type, expected',
    exp_types
)
def test_is_tso_from_exptype(exp_type, expected):
    """Test is_tso integrity based on exp_type"""
    model = datamodels.JwstDataModel()
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
    assert (is_tso(model) or model.meta.exposure.type.lower() in
            ['nrc_tsimage', 'nrc_tsgrism'])


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

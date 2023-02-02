"""Test that the stdatamodels.jwst.datamodels API is replicated in jwst.datamodels"""

import pytest

from stdatamodels.jwst import datamodels as stdm
from jwst import datamodels as jwstdm


@pytest.mark.parametrize("model", stdm.__all__)
def test_stdatamodels_api(model):
    assert hasattr(jwstdm, model)
    assert getattr(stdm, model) is getattr(jwstdm, model)


@pytest.mark.parametrize("model", set(jwstdm.__all__) - set(jwstdm._jwst_models))
def test_jwst_datamodels_api(model):
    assert hasattr(stdm, model)
    assert getattr(stdm, model) is getattr(jwstdm, model)


@pytest.mark.parametrize("model", jwstdm._jwst_models)
def test_jwst_datamodels(model):
    assert not hasattr(stdm, model)

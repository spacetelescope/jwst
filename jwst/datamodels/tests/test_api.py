"""Test that the stdatamodels.jwst.datamodels API is replicated in jwst.datamodels"""

import pytest
import pkgutil
import importlib

from stdatamodels.jwst import datamodels as stdm
from jwst import datamodels as jwstdm


STDM_MODULES = [mdl.name for mdl in pkgutil.iter_modules(stdm.__path__)
                if not mdl.ispkg and mdl.name not in stdm._private_modules]
JWST_MODULES = [mdl.name for mdl in pkgutil.iter_modules(jwstdm.__path__) 
                if not mdl.ispkg and mdl.name not in jwstdm._jwst_modules]


def assert_has_same_import(module_a, module_b, import_):
    assert hasattr(module_b, import_)
    assert getattr(module_a, import_) is getattr(module_b, import_)


@pytest.mark.parametrize("model", stdm.__all__)
def test_stdatamodels_api(model):
    assert_has_same_import(stdm, jwstdm, model)


@pytest.mark.parametrize("module", STDM_MODULES)
def test_stdatamodels_modules(module):
    stdm_module = importlib.import_module(f"stdatamodels.jwst.datamodels.{module}")
    jwst_module = importlib.import_module(f"jwst.datamodels.{module}")

    for import_ in stdm_module.__all__:
        assert_has_same_import(stdm_module, jwst_module, import_)

    for import_ in jwst_module.__all__:
        assert_has_same_import(jwst_module, stdm_module, import_)


@pytest.mark.parametrize("module", JWST_MODULES)
def test_jwst_datamodels_modules(module):
    stdm_module = importlib.import_module(f"stdatamodels.jwst.datamodels.{module}")
    jwst_module = importlib.import_module(f"jwst.datamodels.{module}")

    for import_ in stdm_module.__all__:
        assert_has_same_import(stdm_module, jwst_module, import_)

    for import_ in jwst_module.__all__:
        assert_has_same_import(jwst_module, stdm_module, import_)


@pytest.mark.parametrize("model", set(jwstdm.__all__) - set(jwstdm._jwst_models))
def test_jwst_datamodels_api(model):
    assert_has_same_import(jwstdm, stdm, model)


@pytest.mark.parametrize("model", jwstdm._jwst_models)
def test_jwst_datamodels(model):
    assert not hasattr(stdm, model)
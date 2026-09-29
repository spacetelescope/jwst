import pytest
from stdatamodels.jwst.datamodels import JwstDataModel, SlitModel

from jwst.lib.exposure_types import is_moving_target, is_point_source


@pytest.mark.parametrize("target_type", ["FIXED", "MOVING", None])
@pytest.mark.parametrize("model_type", ["meta", "datamodel"])
def test_is_moving_target(model_type, target_type):
    if model_type == "meta":
        model = {"meta.target.type": target_type}
    elif model_type == "datamodel":
        model = JwstDataModel()
        model.meta.target.type = target_type

    if target_type == "MOVING":
        assert is_moving_target(model)
    else:
        assert not is_moving_target(model)


def test_moving_target_typeerror():
    with pytest.raises(TypeError):
        is_moving_target("not a datamodel or dict")


@pytest.mark.parametrize("override", ["POINT", "EXTENDED", "UNKNOWN", None])
@pytest.mark.parametrize("level", ["top", "meta", None])
@pytest.mark.parametrize("value", ["POINT", "EXTENDED", "UNKNOWN", None])
def test_is_point_source(override, level, value):
    model = SlitModel()
    if level == "top":
        model.source_type = value
        model.meta.target.source_type = None
    elif level == "meta":
        model.source_type = None
        model.meta.target.source_type = value
    else:
        model.source_type = None
        model.meta.target.source_type = None

    check = is_point_source(model, override_srctype=override)
    if override is not None:
        if override == "POINT":
            assert check is True
        else:
            assert check is False
    elif level is not None and value == "POINT":
        assert check is True
    else:
        assert check is False

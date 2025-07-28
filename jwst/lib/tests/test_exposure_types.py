import pytest
from stdatamodels.jwst.datamodels import JwstDataModel

from jwst.lib.exposure_types import is_moving_target


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

import pytest

from jwst.ramp_fitting.ramp_fit_step import _validate_ramp_fit_outputs


def test_ramp_fit_output_validation_accepts_products():
    _validate_ramp_fit_outputs({"slope": object()}, {"slope": object()})


@pytest.mark.parametrize(
    "image_info, integ_info",
    [(None, None), ({"slope": object()}, None), (None, {"slope": object()})],
)
def test_ramp_fit_output_validation_rejects_missing_products(image_info, integ_info):
    with pytest.raises(RuntimeError, match="did not produce slope products"):
        _validate_ramp_fit_outputs(image_info, integ_info)

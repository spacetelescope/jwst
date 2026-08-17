import pytest

from jwst.pipeline.calwebb_spec2 import _validate_spec2_exposure_type


@pytest.mark.parametrize("exp_type", ["MIR_IMAGE", "NRC_IMAGE", "NIS_IMAGE", "FGS_IMAGE"])
def test_spec2_rejects_imaging_exposure_types(exp_type):
    with pytest.raises(ValueError, match="imaging data"):
        _validate_spec2_exposure_type(exp_type)


def test_spec2_accepts_spectroscopic_exposure_type():
    _validate_spec2_exposure_type("NRS_FIXEDSLIT")

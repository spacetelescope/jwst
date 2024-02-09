from synphot.spectrum import SpectralElement

import pytest

import jwst.ami.utils as utils


@pytest.mark.parametrize("filt", ["F277W", "F380M", "F430M", "F480M"])
def test_get_filt_spec(filt):
    spec = utils.get_filt_spec(filt)
    assert isinstance(spec, SpectralElement)


def test_get_filt_spec_fail():
    with pytest.raises(ValueError, match="not a known NIRISS AMI filter"):
        utils.get_filt_spec("unknown")

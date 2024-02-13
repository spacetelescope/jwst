import re

import numpy as np
import pytest
from synphot.spectrum import SourceSpectrum, SpectralElement

import jwst.ami.utils as utils


@pytest.mark.parametrize("filt, peak", [
    ("F277W", 30533.8),
    ("F380M", 37533.2),
    ("F430M", 42512.1),
    ("F480M", 48177.0),
])
def test_get_filt_spec(filt, peak):
    spec = utils.get_filt_spec(filt)
    assert isinstance(spec, SpectralElement)
    wpeak = spec.wpeak()
    assert wpeak.unit == "Angstrom"
    np.isclose(spec.wpeak().value, peak, rtol=1)


def test_get_filt_spec_fail():
    with pytest.raises(ValueError, match="not a known NIRISS AMI filter"):
        utils.get_filt_spec("unknown")


# only testing a few of the catalog entries
@pytest.mark.parametrize("spectype, teff", [
    ("O3V", 45000),
    ("M2I", 3500),
])
def test_get_src_spec(spectype, teff):
    spec = utils.get_src_spec(spectype)
    assert isinstance(spec, SourceSpectrum)
    # check model expression for matching teff
    re.match(r".*T_eff=(?P<teff>[0-9]+)", spec.meta["expr"]).group("teff") == str(teff)


def test_get_src_spec_default():
    with pytest.warns(UserWarning, match="Input spectral type missing did not match"):
        spec = utils.get_src_spec("missing")
    assert isinstance(spec, SourceSpectrum)
    # check model expression for matching teff
    re.match(r".*T_eff=(?P<teff>[0-9]+)", spec.meta["expr"]).group("teff") == "9500"
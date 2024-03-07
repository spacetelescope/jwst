import re

import numpy as np
import pytest
from synphot.spectrum import SourceSpectrum

import jwst.ami.utils as utils


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


@pytest.mark.parametrize("shape, center", [
    ((10, 10), (4.5, 4.5)),
    ((11, 11), (5, 5)),
])
def test_centerpoint(shape, center):
    assert utils.centerpoint(shape) == center


def test_find_centroid():
    arr = np.zeros((30, 30), dtype='f4')
    arr[15, 15] = 1
    thresh = 0.02
    assert np.allclose(utils.find_centroid(arr, thresh), (0.5, 0.5))


@pytest.mark.parametrize("mas, rad", [
    (206264.8062471, 0.001),
    (103132403.12355, 0.5),
])
def test_mas2rad(mas, rad):
    assert np.isclose(utils.mas2rad(mas), rad)

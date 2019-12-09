from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def patch_bigdata(tmpdir, monkeypatch):
    monkeypatch.setenv("TEST_BIGDATA", str(tmpdir))


def test_regtestdata_get_data(rtdata, tmpdir, _jail):
    rtdata._env = ''
    rtdata._inputs_root = ''
    rtdata.output = "foo.fits"
    Path(tmpdir.join("foo_uncal.fits")).touch()
    rtdata.get_data("foo_uncal.fits", docopy=False)

    assert rtdata.input == str(tmpdir.join("foo_uncal.fits"))


def test_regtestdata_get_truth(rtdata, tmpdir, _jail):
    rtdata._env = ''
    rtdata._inputs_root = ''
    Path("truth").mkdir()
    pt = Path("truth/foo_cal.fits").resolve()
    pt.touch()
    rtdata.get_truth("truth/foo_cal.fits", docopy=False)
    rtdata.output = "foo.fits"

    assert rtdata.truth == str(tmpdir.join("truth/foo_cal.fits"))
    assert rtdata.bigdata_root is not None


def test_fitsdiff_defaults(fitsdiff_default_kwargs):
    assert 'ASDF' in fitsdiff_default_kwargs['ignore_hdus']

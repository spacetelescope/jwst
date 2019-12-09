import os
from pathlib import Path

import pytest


def test_regtestdata_get_data(rtdata, tmpdir, _jail, monkeypatch):
    monkeypatch.setenv("TEST_BIGDATA", str(tmpdir))
    rtdata._env = ''
    rtdata._inputs_root = ''
    rtdata.output = "foo.fits"
    Path(tmpdir.join("foo_uncal.fits")).touch()
    rtdata.get_data("foo_uncal.fits", docopy=False)

    assert rtdata.input == str(tmpdir.join("foo_uncal.fits"))


def test_regtestdata_get_truth(rtdata, tmpdir, _jail, monkeypatch):
    monkeypatch.setenv("TEST_BIGDATA", str(tmpdir))
    rtdata._env = ''
    rtdata._inputs_root = ''
    Path("truth").mkdir()
    pt = Path("truth/foo_cal.fits").resolve()
    pt.touch()
    rtdata.get_truth("truth/foo_cal.fits", docopy=False)
    rtdata.output = "foo.fits"

    assert rtdata.truth == str(tmpdir.join("truth/foo_cal.fits"))
    assert rtdata.bigdata_root is not None


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    rtdata = rtdata_module
    rtdata.input = "foo_uncal.fits"

    # Pretend to run pipeline and return rtdata
    return rtdata


JAIL = None

@pytest.mark.parametrize("foo", ["foo1.fits", "foo2.fits"], ids=["foo1", "foo2"])
def test_jail_module(foo, request, run_pipeline):
    rtdata = run_pipeline
    rtdata.output = foo
    rtdata.truth = "truth/" + foo
    rtdata.truth_remote = "jwst-pipeline/truth/" + foo

    # Assert that both tests run in the same tmpdir provided by run_pipeline
    global JAIL
    if JAIL is None:
        JAIL = os.path.dirname(foo)
    else:
        assert JAIL == os.path.dirname(foo)

    # from pprint import pprint
    # pprint(request.node.__dict__)
    # assert 0


def test_fitsdiff_defaults(fitsdiff_default_kwargs):
    assert 'ASDF' in fitsdiff_default_kwargs['ignore_hdus']

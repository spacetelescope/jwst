import logging

import pytest
from stdatamodels.jwst import datamodels

import jwst.stpipe._cal_logs
from jwst.lib.basic_utils import LoggingContext
from jwst.pipeline import Detector1Pipeline, Image2Pipeline, Image3Pipeline
from jwst.pipeline.tests.helpers import make_miri_ramp_model
from jwst.stpipe._cal_logs import _scrub
from jwst.stpipe.tests.steps import CalLogsPipeline, CalLogsStep

_FAKE_HOSTNAME = "my_hostname"
_FAKE_USER = "tim"


@pytest.fixture(autouse=True)
def dont_want_no_scrubs(monkeypatch):
    """Fake hostname and user for consistent _scrub behavior."""
    monkeypatch.setattr(jwst.stpipe._cal_logs, "_HOSTNAME", _FAKE_HOSTNAME)
    monkeypatch.setattr(jwst.stpipe._cal_logs, "_USER", _FAKE_USER)
    yield


def test_cal_logs_step():
    # Set the log level to INFO, since it is not directly configured in `run`
    with LoggingContext(logging.getLogger("jwst"), level=logging.INFO):
        m = CalLogsStep().run("foo")
    assert any(("foo" in l for l in m.cal_logs.cal_logs_step))


def test_cal_logs_pipeline():
    # Set the log level to INFO, since it is not directly configured in `run`
    with LoggingContext(logging.getLogger("jwst"), level=logging.INFO):
        m = CalLogsPipeline().run("foo")
    assert not hasattr(m.cal_logs, "cal_logs_step")
    assert any(("foo" in l for l in m.cal_logs.cal_logs_pipeline))


@pytest.mark.parametrize(
    "msg, is_empty",
    [
        ("2025-02-21T19:16:07.219", False),  # our timestamp
        (_FAKE_HOSTNAME, True),
        (_FAKE_USER, True),
        (f" something from {_FAKE_USER}", True),
        (f":{_FAKE_USER},", True),  # scrub if surrounded by punctuation
        (f"run{_FAKE_USER}e", True),  # scrub whole line if part of word
        ("123.42.26.1", True),
        ("leading 123.42.26.1 trailing", True),
        ("123.42.26", False),
        ("2001:db8::ff00:42:8329", True),
        ("leading 2001:db8:4006:812::200e trailing", True),
    ],
)
def test_scrub(msg, is_empty):
    target = "" if is_empty else msg
    assert _scrub(msg) == target


@pytest.mark.parametrize(
    "msg, expected",
    [
        (f"/some/path/{_FAKE_USER}/subdir/file.txt", "file.txt"),
        ("before /path/file.txt after", "before file.txt after"),
        ("/filename/without/extension", "extension"),
        ("/hyphenated-path/hyphenated-file.tar.gz", "hyphenated-file.tar.gz"),
        (
            f"/some/path/file.txt and '/another/{_FAKE_USER}/file2.txt'",
            "file.txt and 'file2.txt'",
        ),
        ("before relative/file.txt after", "before relative/file.txt after"),
        ("relative/nested/file.txt", "relative/nested/file.txt"),
        ("../file.txt", "../file.txt"),
        ("file I/O error", "file I/O error"),
        (f"/{_FAKE_USER}/file3.txt, /some/path/file4.txt", "file3.txt, file4.txt"),
        (f"{_FAKE_HOSTNAME} /path/file5.txt", ""),
        (f"{_FAKE_USER} /{_FAKE_USER}/file6.txt", ""),
        ("123.42.26.1 /some/path/file7.txt", ""),
        (f"2001:db8::ff00:42:8329 /{_FAKE_USER}/file8.txt", ""),
    ],
)
def test_path_scrub(msg, expected):
    scrubbed = _scrub(msg)
    assert scrubbed == expected


def test_miri_pipeline_cal_logs(tmp_cwd):
    input_model = make_miri_ramp_model()
    input_model.meta.exposure.type = "MIR_IMAGE"
    Detector1Pipeline.call(input_model, save_results=True)

    par2_dict = {
        "assign_wcs": {"skip": True},
        "flat_field": {"skip": True},
        "photom": {"skip": True},
        "resample": {"skip": True},
    }
    Image2Pipeline.call(
        "test_miri_rate.fits", steps=par2_dict, output_file="test_miri", save_results=True
    )

    par3_dict = {
        "tweakreg": {"skip": True},
        "skymatch": {"skip": True},
        "outlier_detection": {"skip": True},
        "resample": {"skip": True},
        "source_catalog": {"skip": True},
    }
    Image3Pipeline.call(
        "test_miri_cal.fits", steps=par3_dict, output_file="test_miri", save_results=True
    )

    # https://github.com/spacetelescope/jwst/issues/9862
    n_image3_steps_found = 0
    expected_steps = ["outlier_detection", "resample", "skymatch", "tweakreg"]
    with datamodels.open("test_miri_i2d.fits") as im:
        assert sorted(im.cal_logs.instance.keys()) == ["calwebb_image2", "calwebb_image3"]
        image3_logs = im.cal_logs.calwebb_image3
        for step_name in expected_steps:
            for line in image3_logs:
                if f"INFO - Step {step_name} running" in line:
                    n_image3_steps_found += 1
                    break
        assert n_image3_steps_found == len(expected_steps)

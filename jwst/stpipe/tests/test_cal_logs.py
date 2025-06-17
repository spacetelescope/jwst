import pytest

from jwst.stpipe.tests.steps import CalLogsStep, CalLogsPipeline
import jwst.stpipe._cal_logs
from jwst.stpipe._cal_logs import _scrub


_FAKE_HOSTNAME = "my_hostname"
_FAKE_USER = "my_user"


@pytest.fixture(autouse=True)
def dont_want_no_scrubs(monkeypatch):
    """Fake hostname and user for consistent _scrub behavior."""
    monkeypatch.setattr(jwst.stpipe._cal_logs, "_HOSTNAME", _FAKE_HOSTNAME)
    monkeypatch.setattr(jwst.stpipe._cal_logs, "_USER", _FAKE_USER)
    yield


def test_cal_logs_step():
    m = CalLogsStep().run("foo")
    assert any(("foo" in l for l in m.cal_logs.cal_logs_step))


def test_cal_logs_pipeline():
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
        ("123.42.26.1", True),
        ("123.42.26", False),
        ("2001:db8::ff00:42:8329", True),
        ("2001:db8:4006:812::200e", True),
    ],
)
def test_scrub(msg, is_empty):
    target = "" if is_empty else msg
    assert _scrub(msg) == target

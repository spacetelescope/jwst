import logging

import stcal.tweakreg.tweakreg as twk

from jwst.tweakreg.tweakreg_step import _log_alignment_failure


def test_alignment_failure_logs_context_and_traceback(caplog):
    caplog.set_level(logging.WARNING, logger="jwst.tweakreg.tweakreg_step")
    try:
        raise twk.TweakregError("not enough matches")
    except twk.TweakregError as error:
        _log_alignment_failure("Relative", error)

    assert "Relative alignment failed: not enough matches" in caplog.text
    assert "TweakregError" in caplog.text

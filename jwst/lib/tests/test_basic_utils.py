import logging

from jwst.lib.basic_utils import LoggingContext, disable_logging

log = logging.getLogger(__name__)


def test_logging_context(caplog):
    with LoggingContext(logger=log, level=logging.CRITICAL):
        log.info("test 1")
    assert caplog.text == ""

    with LoggingContext(logger=log, level=logging.INFO):
        log.info("test 2")
    assert "test 2" in caplog.text


def test_disable_logging(caplog):
    with LoggingContext(logger=log, level=logging.INFO):
        # log is enabled
        log.info("test 1")
        assert "test 1" in caplog.text
        assert logging.root.manager.disable == logging.NOTSET

        with disable_logging():
            log.info("test 2")
            assert logging.root.manager.disable == logging.CRITICAL
        assert "test 2" not in caplog.text

        # log is re-enabled after exit
        log.info("test 3")
        assert "test 3" in caplog.text
        assert logging.root.manager.disable == logging.NOTSET

        # if the level is higher than the message, it will appear
        with disable_logging(level=logging.DEBUG):
            log.info("test 4")
            assert logging.root.manager.disable == logging.DEBUG
        assert "test 4" in caplog.text
        assert logging.root.manager.disable == logging.NOTSET


def test_disable_logging_already_disabled(caplog):
    with LoggingContext(logger=log, level=logging.INFO):
        with disable_logging():
            # log is disabled
            log.info("test 1")
            assert caplog.text == ""
            assert logging.root.manager.disable == logging.CRITICAL

            # log is still disabled inside the context manager
            with disable_logging(level=logging.DEBUG):
                log.info("test 2")
                assert caplog.text == ""
                assert logging.root.manager.disable == logging.CRITICAL

            # log is still disabled after exiting the context manager
            log.info("test 3")
            assert caplog.text == ""
            assert logging.root.manager.disable == logging.CRITICAL

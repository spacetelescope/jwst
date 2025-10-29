"""General utility objects."""

__all__ = ["LoggingContext", "disable_logging"]

import logging
from contextlib import contextmanager


class LoggingContext:
    """
    Logging context manager.

    Keep logging configuration within a context.
    Based on the Python 3 Logging Cookbook example.

    Parameters
    ----------
    logger : logging.Logger
        The logger to modify.

    level : int
        The log level to set.

    handler : logging.Handler
        The handler to use.

    close : bool
        Close the handler when done.
    """

    def __init__(self, logger, level=None, handler=None, close=True):
        self.logger = logger
        self.level = level
        self.handler = handler
        self.close = close

        self.old_level = None

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)
        if self.handler:
            self.logger.addHandler(self.handler)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        if self.handler:
            self.logger.removeHandler(self.handler)
        if self.handler and self.close:
            self.handler.close()
        # implicit return of None => don't swallow exceptions


@contextmanager
def disable_logging(level=logging.CRITICAL):
    """
    Disable logging within a context.

    Parameters
    ----------
    level : int, optional
        Logging level.  At this level and below, all logging is disabled.
        Defaults to `logging.CRITICAL`, which disables logging at all levels.

    Examples
    --------
    The context manager is used as::

        with disable_logging(level=logging.ERROR):
            # code containing logging to ignore, other than CRITICAL messages
            ...
    """
    current_level = logging.root.manager.disable
    if current_level < level:
        logging.disable(level)
        try:
            yield
        finally:
            logging.disable(current_level)
    else:
        yield

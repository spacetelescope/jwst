"""General utility objects"""

from stdatamodels.jwst.datamodels import dqflags
import numpy as np


def set_nans_to_donotuse(data, dq):
    """Set all NaN values in the data that have an even value to
    DO_NOT_USE.

    Parameters
    ----------
    data : numpy array
        The science data array to find NaN values and
        check of these have a DQ flag=DO_NOT_USE, or
        set it if not.

    dq : numpy array
        The DQ array to be checked.

    Returns
    -------
    dq : numpy array
        The updated DQ array.
    """
    dq[np.isnan(data)] |= dqflags.pixel['DO_NOT_USE']
    return dq


class LoggingContext:
    """Logging context manager

    Keep logging configuration within a context

    Based on the Python 3 Logging Cookbook example

    Parameters
    ==========
    logger: logging.Logger
        The logger to modify.

    level: int
        The log level to set.

    handler: logging.Handler
        The handler to use.

    close: bool
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

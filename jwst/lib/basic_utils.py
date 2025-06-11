"""General utility objects."""


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
    """

    def __init__(self, logger, level=None):
        self.logger = logger
        self.level = level
        self.old_level = None

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)

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

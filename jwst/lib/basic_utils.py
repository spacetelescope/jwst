"""General utility objects"""


def deprecate_class(
    new_class,
    message='"{old_class}" is deprecated and will be removed. Use {new_class}',
):
    """Deprecate a class in favor of another class

    Parameters
    ----------
    new_class : class
        The class/object that should be used instead

    message : str
        The deprecation warning message
    """

    # Implement the inner decorator
    def _decorator(old_class):
        def init(self, *args, **kwargs):
            import warnings

            warnings.simplefilter("default")
            warnings.warn(
                message.format(
                    old_class=old_class.__name__, new_class=new_class.__name__
                ),
                category=DeprecationWarning,
            )
            new_class.__init__(self, *args, **kwargs)

        # Create the class
        deprecated = type(
            old_class.__name__,
            (new_class,),
            {
                "__doc__": old_class.__doc__,
                "__init__": init,
                "__module__": old_class.__module__,
                "__name__": old_class.__name__,
            },
        )

        return deprecated

    return _decorator


def bytes2human(n):
    """Convert bytes to human-readable format

    Taken from the `psutil` library which references
    http://code.activestate.com/recipes/578019

    Parameters
    ----------
    n : int
        Number to convert

    Returns
    -------
    readable : str
        A string with units attached.

    Examples
    --------
    >>> bytes2human(10000)
        '9.8K'

    >>> bytes2human(100001221)
        '95.4M'
    """
    symbols = ("K", "M", "G", "T", "P", "E", "Z", "Y")
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return "%.1f%s" % (value, s)
    return "%sB" % n


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

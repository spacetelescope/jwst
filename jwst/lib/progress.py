"""Provide a visual progress bar facility

The actual functionality is provided by the external library `progress`

https://pypi.org/project/progress/

If the module is not available, then stub it out.
"""
import logging

__all__ = ['Bar']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class _BarStub:
    """Stub the Bar functionality"""
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    def next(self):
        pass


# Stub the Bar functionality if the actual module does not exist.
try:
    from progress.bar import Bar as _Bar
except ModuleNotFoundError:
    _Bar = _BarStub


def Bar(*args, log_level=logging.INFO, log_cutoff=logging.INFO, **kwargs):
    """Actually use Bar only if logging level is appropriate

    Parameters
    ----------
    *args : *args
        Positional arguments for the `progress.Bar` class

    log_level : int
        The current log level. Only produce a Bar if the `log_level`
        is less than or equal to `log_cutoff`

    log_cutoff : int
        Maximum logging level to allow Bar to be created.

    **kwargs : **kwargs
        Keyword arguments for `progress.Bar` class
    """
    if log_level <= log_cutoff:
        return _Bar(*args, **kwargs)
    return _BarStub(*args, **kwargs)

"""Provide a visual progress bar facility

The actual functionality is provided by the external library `progress`

https://pypi.org/project/progress/

If the module is not available, then stub it out.
"""
from contextlib import contextmanager

try:
    from progress.bar import Bar
except ModuleNotFoundError:


    class Bar:
        """Stub the Bar functionality"""
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args, **kwargs):
            pass

        def next(self):
            pass

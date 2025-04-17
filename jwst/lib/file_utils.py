"""File utility functions."""

import os
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def pushdir(directory):
    """
    Temporarily change to specified directory.

    Parameters
    ----------
    directory : File-like object
        Directory to change to.

    Returns
    -------
    new_directory : Path
        The directory changed to.
    """
    previous = Path.cwd()
    try:
        os.chdir(directory)
        yield Path.cwd()
    finally:
        os.chdir(previous)

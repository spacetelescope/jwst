"""File utility functions"""
from contextlib import contextmanager
import os
from pathlib import Path


@contextmanager
def pushdir(directory):
    """Temporarily change to specified directory

    Parameters
    ----------
    directory: File-like object
        Directory to change to

    Returns
    -------
    new_directory: Path
        The directory changed to.
    """
    previous = os.getcwd()
    try:
        os.chdir(directory)
        yield Path.cwd()
    finally:
        os.chdir(previous)

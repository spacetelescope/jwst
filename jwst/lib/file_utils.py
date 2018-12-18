"""File utility functions"""
from contextlib import contextmanager
import os


@contextmanager
def pushdir(directory):
    """Temporarily change to specified directory

    Parameters
    ----------
    directory: File-like object
        Directory to change to
    """
    previous = os.getcwd()
    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(previous)

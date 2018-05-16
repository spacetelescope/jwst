# Testing helpers
from contextlib import contextmanager
from glob import glob
from pathlib import Path
import os

from ..tests.helpers import abspath

# Python2 compatibility for TemporaryDirectory
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from .tempfile_py2 import TemporaryDirectory


INPUT_FILES_GLOB = 'data/jwst_nod?_cal.fits'


def t_path(partial_path):
    here_path = Path(os.path.dirname(__file__))
    path = Path(partial_path)
    """Construct the full path for test files"""
    if path.is_absolute():
        result = path
    else:
        result = here_path.joinpath(path)

    return str(result)


@contextmanager
def chdir(path):
    """A context manager which changes the working directory to the given
    path, and then changes it back to its previous value on exit.

    """
    prev_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)

# Calcuate some extra constants
INPUT_FILES_GLOB = t_path(INPUT_FILES_GLOB)
INPUT_FILES = glob(t_path(INPUT_FILES_GLOB))

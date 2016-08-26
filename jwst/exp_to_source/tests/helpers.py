# Testing helpers
from glob import glob
import os

INPUT_FILES_GLOB = 'data/jwst_nod?_cal.fits'


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)


# Python2 compatibility for TemporaryDirectory
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from .tempfile_py2 import TemporaryDirectory


# Calcuate some extra constants
INPUT_FILES_GLOB = t_path(INPUT_FILES_GLOB)
INPUT_FILES = glob(t_path(INPUT_FILES_GLOB))

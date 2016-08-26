# Provide Python2 compatibility
from contextlib import contextmanager
from glob import glob
import os
from shutil import rmtree
from tempfile import (mkdtemp, mkstemp)


@contextmanager
def TemporaryDirectory(suffix='', prefix='tmp', dir=None):
    """Provide Python3 functionality

    Notes
    -----
    See the python3 documentation
    """
    path = mkdtemp(suffix, prefix, dir)
    try:
        yield path
    finally:
        rmtree(path)


# #####
# Tests
# #####
def test_TemporaryDirectory():
    with TemporaryDirectory() as path:
        assert os.path.exists(path)
        os.chdir(path)
        files = glob('*')
        assert len(files) == 0

        (fd, file_path) = mkstemp(dir=path)

        assert os.path.exists(file_path)
        assert path in file_path
        files = glob('*')
        assert len(files) == 1

    assert not os.path.exists(path)
    assert not os.path.exists(file_path)

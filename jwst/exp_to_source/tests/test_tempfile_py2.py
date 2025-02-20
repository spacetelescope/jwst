# Provide Python2 compatibility
from contextlib import contextmanager
from pathlib import Path
from shutil import rmtree
from tempfile import (mkdtemp, mkstemp)


@contextmanager
def temporary_directory(suffix='', prefix='tmp', dir_nme=None):
    """Provide Python3 functionality.

    Notes
    -----
    See the python3 documentation
    """
    path = mkdtemp(suffix, prefix, dir_nme)
    try:
        yield path
    finally:
        rmtree(path)


# #####
# Tests
# #####
def test_temporary_directory():
    """Test the temporary directory function."""
    with temporary_directory() as path:
        assert Path.exists(path)
        files = path.glob('*')
        assert not list(files)

        (fd, file_path) = mkstemp(dir=path)

        assert Path.exists(file_path)
        assert path in file_path
        files = file_path.glob('*')
        assert not len(list(files)) == 1

    assert not Path.exists(path)
    assert not Path.exists(file_path)

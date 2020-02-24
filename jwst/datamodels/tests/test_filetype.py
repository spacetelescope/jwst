import pytest

from ..filetype import check


def test_check_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        check('file_does_not_exist.fits')

import pytest

from ..filetype import check


@pytest.fixture(params=['fits', 'json', 'asdf'])
def test_file(request):
    return f'test_file.{request.param}'


def test_check_raises_filenotfound(test_file):
    with pytest.raises(FileNotFoundError):
        check(test_file)

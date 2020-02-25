import pytest

from ..filetype import check

SUPPORTED_EXTS = (('fits', 'fits'), ('json', 'asn'), ('asdf', 'asdf'))  # (ext, expected filetype)


@pytest.fixture(params=SUPPORTED_EXTS)
def test_file(request):
    return f'test_file.{request.param[0]}', request.param[-1]


def test_check_on_str_init(test_file):
    filename, expected = test_file
    filetype = check(filename)

    assert filetype == expected


def test_check_fails_on_unsupported_ext():
    with pytest.raises(ValueError):
        check('test_file')

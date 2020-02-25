import pytest

from ..filetype import check

SUPPORTED_EXTS = (('fits', 'fits'), ('json', 'asn'), ('asdf', 'asdf'))  # (ext, expected filetype)


@pytest.fixture(params=SUPPORTED_EXTS)
def input_file(request):
    return f'test_file.{request.param[0]}', request.param[-1]


@pytest.fixture(params=['stpipe.MyPipeline.fits', 'stpipe.MyPipeline.fits.gz'])
def pipeline_file(request):
    return request.param


def test_check_on_str_init(input_file):
    filename, expected = input_file
    filetype = check(filename)

    assert filetype == expected


def test_check_fails_on_unsupported_ext():
    with pytest.raises(ValueError):
        check('test_file')


def test_check_works_for_zipped(input_file):
    filename, expected = input_file
    filename += '.gz'  # extra zip extension

    filetype = check(filename)

    assert filetype == expected


def test_check_works_for_pipeline_patters(pipeline_file):
    assert check(pipeline_file) == 'fits'

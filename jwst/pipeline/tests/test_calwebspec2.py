import pytest

from jwst.pipeline.calwebb_spec2 import Spec2Pipeline


def test_filenotfounderror_raised(capsys):
    # Verify the failure is in the traceback message
    with pytest.raises(RuntimeError, match="FileNotFoundError"):
        Spec2Pipeline().run('file_does_not_exist.fits')

    # Verify the failure is printed to stderr
    captured = capsys.readouterr()
    assert 'FileNotFoundError' in captured.err

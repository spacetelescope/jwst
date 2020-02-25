import pytest

from ..calwebb_spec2 import Spec2Pipeline


@pytest.fixture(scope='module')
def fake_pipeline():
    return Spec2Pipeline()


def test_filenotfounderror_raised(fake_pipeline, capsys):
    with pytest.raises(RuntimeError):
        fake_pipeline.run('file_does_not_extis.fits')

    captured = capsys.readouterr()
    assert 'FileNotFoundError' in captured.err

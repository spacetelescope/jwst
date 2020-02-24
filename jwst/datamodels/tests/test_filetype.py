import pytest

from ...pipeline import Spec2Pipeline


@pytest.fixture(scope='module')
def fake_pipeline():
    return Spec2Pipeline()


def test_pipeline_raises_filenotfound(fake_pipeline, capsys):
    with pytest.raises(RuntimeError):
        fake_pipeline.run('file_does_not_exist.fits')

    assert (
            "FileNotFoundError: [Errno 2] No such file or directory: 'file_does_not_exist.fits'"
            in capsys.readouterr()[-1]
    )

from glob import glob
import os

from jwst.exp_to_source.tests import helpers
from jwst.exp_to_source.main import Main


def test_help(capsys):
    Main(['-h'])
    out, err = capsys.readouterr()
    assert out.startswith('usage:')


def test_default_run(tmpdir, capsys):
    path = str(tmpdir)
    args = [
        '-o',
        path
    ]
    args.extend(helpers.INPUT_FILES)
    result = Main(args)
    assert len(result.sources) == 5
    files = glob(os.path.join(path, '*.fits'))
    assert len(files) == 5


def test_dry_run(_jail, capsys):
    no_files = glob('*.fits')
    assert len(no_files) == 0

    args = ['--dry-run']
    args.extend(helpers.INPUT_FILES)
    result = Main(args)
    assert len(result.sources) == 5
    no_files = glob('*.fits')
    assert len(no_files) == 0

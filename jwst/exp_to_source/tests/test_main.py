from glob import glob
import os
import pytest

from . import helpers
from ..main import Main


def test_help(capsys):
    result = Main(['-h'])
    out, err = capsys.readouterr()
    assert out.startswith('usage:')


def test_default_run(capsys):
    with helpers.TemporaryDirectory() as path:
        args = [
            '-o',
            path
        ]
        args.extend(helpers.INPUT_FILES)
        result = Main(args)
        assert len(result.sources) == 5
        files = glob(os.path.join(path, '*.fits'))
        assert len(files) == 5


def test_dry_run(capsys):
    with helpers.TemporaryDirectory() as path:
        with helpers.chdir(path):
            no_files = glob('*.fits')
            assert len(no_files) == 0

            args = ['--dry-run']
            args.extend(helpers.INPUT_FILES)
            result = Main(args)
            assert len(result.sources) == 5
            no_files = glob('*.fits')
            assert len(no_files) == 0

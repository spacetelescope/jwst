from pathlib import Path
import importlib
from pkgutil import iter_modules
import os

import pytest

from ci_watson.artifactory_helpers import get_bigdata_root

from jwst.regtest.regtestdata import (
    _data_glob_local,
    _data_glob_url
)
from jwst.tests.helpers import word_precision_check


def test_word_precision_check():
    """Test word_precision_check"""
    s1 = "a b c"
    s2 = "aa bb cc"
    s3 = "aa bb cc dd"
    s4 = "aazz bbzz cczz"

    assert word_precision_check(s1, s1)
    assert not word_precision_check(s1, s2)
    assert word_precision_check(s1, s2, length=1)
    assert not word_precision_check(s2, s3)
    assert word_precision_check(s2, s4, length=2)


@pytest.mark.parametrize(
    'glob_filter, nfiles',
    [
        ('*', 3),
        ('*.txt', 3),
        ('*.fits', 0)
    ], ids=['all', 'txt', 'fits']
)
def test_data_glob_local(glob_filter, nfiles, _jail):
    """Test working of local globbing

    Parameters
    ----------
    glob_filter: str
        The glob filter to use.

    nfiles: int
        The number of files expected to find.
    """
    path = Path('datadir')
    path.mkdir()
    for idx in range(3):
        with open(path / ('afile' + str(idx) + '.txt'), 'w') as fh:
            fh.write(f'I am file {idx}')

    files = _data_glob_local(path, glob_filter)
    assert len(files) == nfiles


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'glob_filter, nfiles',
    [
        ('*', 1),
        ('*.txt', 0),
        ('*.fits', 1)
    ]
)
def test_data_glob_url(glob_filter, nfiles, pytestconfig, request):
    """Test globbing over a URL

    Parameters
    ----------
    glob_filter: str
        The glob filter to use.

    nfiles: int
        The number of files expected to find.
    """
    inputs_root = pytestconfig.getini('inputs_root')[0]
    env = request.config.getoption('env')
    path = os.path.join(inputs_root, env, 'infrastructure/test_data_glob')

    files = _data_glob_url(path, glob_filter, root=get_bigdata_root())
    assert len(files) == nfiles


def test_submodules_can_be_imported():
    """Make sure all package submodules can be imported"""
    import jwst

    submodules = [mod for _, mod, ispkg in iter_modules(jwst.__path__) if ispkg]
    for module in submodules:
        importlib.import_module("jwst." + module)

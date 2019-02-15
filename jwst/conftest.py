"""Project default for pytest"""
import os
import tempfile
import pytest

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES

from jwst.associations import (AssociationRegistry, AssociationPool)
from jwst.associations.tests.helpers import t_path

try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['asdf'] = 'asdf'
    del PYTEST_HEADER_MODULES['h5py']
except (NameError, KeyError):
    pass


@pytest.fixture
def mk_tmp_dirs():
    """Create a set of temporary directorys and change to one of them."""
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()
    try:
        os.chdir(tmp_current_path)
        yield (tmp_current_path, tmp_data_path, tmp_config_path)
    finally:
        os.chdir(old_path)


@pytest.fixture(scope='session')
def full_pool_rules(request):
    """Setup to use the full example pool and registry"""
    pool_fname = t_path('data/mega_pool.csv')
    pool = AssociationPool.read(pool_fname)
    rules = AssociationRegistry()

    return (pool, rules, pool_fname)

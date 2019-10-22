"""Project default for pytest"""
import os
import tempfile
import pytest

from jwst.associations import (AssociationRegistry, AssociationPool)
from jwst.associations.tests.helpers import t_path
from jwst.lib.tests import helpers as lib_helpers
from jwst.lib import s3_utils


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


@pytest.fixture(autouse=True)
def monkey_patch_s3_client(monkeypatch):
    monkeypatch.setattr(s3_utils, "_CLIENT", lib_helpers.MockS3Client())

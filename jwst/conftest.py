"""Project default for pytest"""
import os
import tempfile
import pytest

from jwst.associations import (AssociationRegistry, AssociationPool)
from jwst.associations.tests.helpers import t_path
from jwst.lib.tests import helpers as lib_helpers
from jwst.lib import s3_utils


@pytest.fixture(scope='session')
def full_pool_rules(request):
    """Setup to use the full example pool and registry"""
    pool_fname = t_path('data/mega_pool.csv')
    pool = AssociationPool.read(pool_fname)
    rules = AssociationRegistry()

    return (pool, rules, pool_fname)


@pytest.fixture(scope="module")
def jail(request, tmpdir_factory):
    """Run test in a pristine temporary working directory, scoped to module.

    This fixture is the same as _jail in ci_watson, but scoped to module
    instead of function.  This allows a fixture using it to produce files in a
    temporary directory, and then have the tests access them.
    """
    old_dir = os.getcwd()
    path = request.module.__name__.split('.')[-1]
    if request._parent_request.fixturename is not None:
        path = path + "_" + request._parent_request.fixturename
    newpath = tmpdir_factory.mktemp(path)
    try:
        os.chdir(str(newpath))
        yield newpath
    finally:
        os.chdir(old_dir)


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


@pytest.fixture(autouse=True)
def monkey_patch_s3_client(monkeypatch):
    monkeypatch.setattr(s3_utils, "_CLIENT", lib_helpers.MockS3Client())


@pytest.fixture
def slow(request):
    """Setup slow fixture for tests to identify if --slow
    has been specified
    """
    return request.config.getoption('--slow')

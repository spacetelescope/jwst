"""Project default for pytest"""
import inspect
import os
import tempfile
import warnings

import pytest

from jwst.associations import (AssociationRegistry, AssociationPool)
from jwst.associations.tests.helpers import t_path


@pytest.fixture
def jail_environ():
    """Lock changes to the environment"""
    original = os.environ.copy()
    try:
        yield
    finally:
        os.environ = original


@pytest.fixture(scope='session')
def full_pool_rules(request):
    """Setup to use the full example pool and registry"""
    pool_fname = t_path('data/mega_pool.csv')
    pool = AssociationPool.read(pool_fname)
    rules = AssociationRegistry()

    return pool, rules, pool_fname


@pytest.fixture
def mk_tmp_dirs():
    """Create a set of temporary directories and change to one of them."""
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()
    try:
        os.chdir(tmp_current_path)
        yield (tmp_current_path, tmp_data_path, tmp_config_path)
    finally:
        os.chdir(old_path)


@pytest.fixture
def slow(request):
    """Setup slow fixture for tests to identify if --slow
    has been specified
    """
    return request.config.getoption('--slow')


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
    os.chdir(str(newpath))
    yield newpath
    os.chdir(old_dir)


@pytest.hookimpl(trylast=True)
def pytest_configure(config):
    terminal_reporter = config.pluginmanager.getplugin('terminalreporter')
    config.pluginmanager.register(TestDescriptionPlugin(terminal_reporter), 'testdescription')


class TestDescriptionPlugin:
    """Pytest plugin to print the test docstring when `pytest -vv` is used.

    This plug-in was added to support JWST instrument team testing and
    reporting for the JWST calibration pipeline.
    """

    def __init__(self, terminal_reporter):
        self.terminal_reporter = terminal_reporter
        self.desc = None

    def pytest_runtest_protocol(self, item):
        try:
            # Get the docstring for the test
            self.desc = inspect.getdoc(item.obj)
        except AttributeError:
            self.desc = None

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_logstart(self, nodeid, location):
        # When run as `pytest` or `pytest -v`, no change in behavior
        if self.terminal_reporter.verbosity <= 1:
            yield
        # When run as `pytest -vv`, `pytest -vvv`, etc, print the test docstring
        else:
            self.terminal_reporter.write('\n')
            yield
            if self.desc:
                self.terminal_reporter.write(f'\n{self.desc} ')


def pytest_collection_modifyitems(session, config, items):
    # don't order tests if not in xdist
    if not config.pluginmanager.hasplugin('xdist'):
        return
    with open('foo.txt', 'a+') as f:
        f.write(str(config.option.numprocesses) + '\n')
        f.write(str(config.option.tx) + '\n')
    # a few tests take a very long time (>1 hr)
    # schedule these first
    ordered_tests = [
        'test_spec3_multi[jw01249005001_03101_00001_nrs1_o005_crf.fits]',
        'test_residual_fringe_cal'
    ]
    # find indices of ordered tests
    ordered_test_indices = {k: None for k in ordered_tests}
    for (index, item) in enumerate(items):
        if item.name in ordered_test_indices:
            ordered_test_indices[item.name] = index

    for new_index, name in enumerate(ordered_tests):
        index = ordered_test_indices[name]
        if index is None:
            msg = f"Failed to find test {name} while ordering tests"
            warnings.warn(msg)
            continue
        if new_index == index:  # no need to move
            continue
        # swap this test (at index) with test at new_index
        items[new_index], items[index] = items[index], items[new_index]

        # update indices of ordered_tests for
        # the new index for the test we put in a specific order
        ordered_test_indices[name] = new_index
        # and the index for the test we moved (in case it's in the dict)
        if items[index].name in ordered_test_indices:
            ordered_test_indices[items[index].name] = index
    return items

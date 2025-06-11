"""Set project defaults and add fixtures for pytest."""

import inspect
import logging
import os
import tempfile
from pathlib import Path

import pytest
from astropy.utils.data import get_pkg_data_filename

from jwst.associations import AssociationRegistry, AssociationPool
from jwst.tests.helpers import LogWatcher


@pytest.fixture
def jail_environ():
    """Lock changes to the environment."""
    original = os.environ.copy()
    try:
        yield
    finally:
        os.environ = original


@pytest.fixture(scope="session")
def full_pool_rules(request):
    """
    Set up the full example pool and registry.

    Returns
    -------
    pool : AssociationPool
        The full example pool as read from data/mega_pool.csv.
    rules : AssociationRegistry
        The registry of available associations.
    pool_fname : str
        The full test path to mega_pool.csv.
    """
    pool_fname = get_pkg_data_filename("data/mega_pool.csv", package="jwst.associations.tests")
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
    """
    Set up slow fixture for tests to identify if --slow has been specified.

    Returns
    -------
    bool
        True if --slow has been specified, False otherwise.
    """
    return request.config.getoption("--slow")


@pytest.fixture(scope="module")
def tmp_cwd_module(request, tmp_path_factory):
    """
    Set up fixture to run test in a pristine temporary working directory, scoped to module.

    This allows a test using this fixture to produce files in a
    temporary directory, and then have the tests access them.

    Yields
    ------
    tmp_path
        The temporary directory path.
    """
    old_dir = os.getcwd()
    path = request.module.__name__.split(".")[-1]
    if request._parent_request.fixturename is not None:
        path = path + "_" + request._parent_request.fixturename
    newpath = tmp_path_factory.mktemp(path)
    os.chdir(str(newpath))
    yield newpath
    os.chdir(old_dir)


@pytest.fixture
def tmp_cwd(tmp_path):
    """Perform test in a pristine temporary working directory, scoped to function."""
    old_dir = Path.cwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(old_dir)


@pytest.hookimpl(trylast=True)
def pytest_configure(config):
    """Add the test description plugin to the pytest configuration."""
    terminal_reporter = config.pluginmanager.getplugin("terminalreporter")
    config.pluginmanager.register(TestDescriptionPlugin(terminal_reporter), "testdescription")


class TestDescriptionPlugin:
    """
    Pytest plugin to print the test docstring when `pytest -vv` is used.

    This plug-in was added to support JWST instrument team testing and
    reporting for the JWST calibration pipeline.
    """

    def __init__(self, terminal_reporter):
        self.terminal_reporter = terminal_reporter
        self.desc = None

    def pytest_runtest_protocol(self, item):
        """Get the docstring for the test."""
        try:
            self.desc = inspect.getdoc(item.obj)
        except AttributeError:
            self.desc = None

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_logstart(self, nodeid, location):
        """Print the test docstring when `pytest -vv` is used."""
        # When run as `pytest` or `pytest -v`, no change in behavior
        if self.terminal_reporter.verbosity <= 1:
            yield
        # When run as `pytest -vv`, `pytest -vvv`, etc, print the test docstring
        else:
            self.terminal_reporter.write("\n")
            yield
            if self.desc:
                self.terminal_reporter.write(f"\n{self.desc} ")


@pytest.fixture()
def log_watcher(monkeypatch):
    """
    Provide a fixture to watch for log messages.

    Parameters
    ----------
    monkeypatch : pytest.monkeypatch.MonkeyPatch
        Monkeypatch fixture.

    Returns
    -------
    _log_watcher : callable
        A function that when called, produces a LogWatcher object for
        a specified log name. Signature is:
        `_log_watcher(log_name, level=None, message="")`.
    """

    def _log_watcher(log_name, message="", level=None):
        """
        Set a log watcher to check for a log message in a specific module.

        To change the message to watch for, set the `message`
        attribute in the returned LogWatcher instance, prior to the
        call that is expected to trigger the message.

        Parameters
        ----------
        log_name : str
            Name of the log to watch.
        message : str, optional
            The message to watch for.  If not provided, the message
            is an empty string, which will match any log message.
        level : str or list of str None, optional
            The log level(s) to watch.  If not provided, the message may
            be raised at any log level.  Level options are:
            "debug", "info", "warning", "error", "critical".

        Returns
        -------
        LogWatcher
            The log watcher for the module.
        """
        # Set a log watcher to check for a log message at any level
        # in the barshadow module
        watcher = LogWatcher(message)
        logger = logging.getLogger(log_name)

        if level is None:
            level = ["debug", "info", "warning", "error", "critical"]
        elif isinstance(level, str):
            level = [level]

        for level_name in level:
            monkeypatch.setattr(logger, level_name, watcher)
        return watcher

    return _log_watcher

from datetime import datetime
import os
from urllib.parse import urlparse
from glob import glob

import pytest
import requests
from ci_watson.artifactory_helpers import (
    get_bigdata_root,
    get_bigdata,
    generate_upload_schema,
    BigdataError
)


TODAYS_DATE = datetime.now().strftime("%Y-%m-%d")


# def pytest_runtest_logreport(report):
#     test_name = report.nodeid.replace('[', '_').replace(']', '')
#     if report.when == 'call' and report.failed:
#         _func(test_name)

@pytest.fixture(scope="session")
def artifactory_repos(pytestconfig):
    inputs_root = pytestconfig.getini('inputs_root')[0]
    results_root = pytestconfig.getini('results_root')[0]
    return inputs_root, results_root


@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()

    # set a report attribute for each phase of a call, which can
    # be "setup", "call", "teardown"
    setattr(item, "rep_" + rep.when, rep)


@pytest.fixture(scope='function', autouse=True)
def generate_artifactory_json(request, artifactory_repos):
    inputs_root, results_root = artifactory_repos

    def _func(schema_pattern=[]):
        import getpass

        # Generate the Artifactory target path
        whoami = getpass.getuser() or 'nobody'
        user_tag = 'NOT_CI_{}'.format(whoami)
        build_tag = os.environ.get('BUILD_TAG', user_tag)
        build_matrix_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', '0')
        subdir = '{}_{}_{}'.format(TODAYS_DATE, build_tag, build_matrix_suffix)
        testname = request.node.originalname or request.node.name
        results_path = os.path.join(results_root, subdir, testname) + os.sep
        cwd = os.path.abspath(os.path.curdir)

        generate_upload_schema(schema_pattern, results_path, request.node.name)

    yield
    # Execute the following at test teardown
    if request.node.rep_setup.passed:
        if request.node.rep_call.failed:
            schema_pattern = []
            for prop in request.node.user_properties:
                name, pattern = prop
                if name == 'output':
                    schema_pattern.append(os.path.abspath(pattern))
            _func(schema_pattern)


@pytest.fixture(scope="module")
def jail(request, tmpdir_factory):
    """Run test in a pristine temporary working directory, scoped to module.

    This fixture is the same as _jail in ci_watson, but scoped to module
    instead of function.  This allows a fixture using it to produce files in a
    temporary directory, and then have the tests access them.
    """
    old_dir = os.getcwd()
    newpath = tmpdir_factory.mktemp(request.module.__name__.split('.')[-1])
    os.chdir(str(newpath))
    yield newpath
    os.chdir(old_dir)


class RegtestData:
    """Defines data paths on Artifactory and data retrieval methods"""

    def __init__(self, env="dev", inputs_root="jwst-pipeline",
        results_root="jwst-pipeline-results", docopy=True):
        self._env = env
        self._inputs_root = inputs_root
        self._results_root = results_root

        self.docopy = docopy

        self._input_remote = None
        self._truth_remote = None
        self._input = None
        self._truth = None
        self._output = None
        self._bigdata_root = get_bigdata_root()

    @property
    def input_remote(self):
        if self._input_remote is not None:
            return os.path.join(*self._input_remote)
        else:
            return None

    @input_remote.setter
    def input_remote(self, value):
        self._input_remote = value.split(os.sep)

    @property
    def truth_remote(self):
        if self._truth_remote is not None:
            return os.path.join(*self._truth_remote)
        else:
            return None

    @truth_remote.setter
    def truth_remote(self, value):
        try:
            dirname = os.path.dirname(os.path.join(*self._input_remote))
            self._truth_remote = os.path.join(dirname, 'truth', value)
        except TypeError as e:
            raise RuntimeError("Define self.input_remote first") from e

    @property
    def input(self):
        return self._input

    @property
    def truth(self):
        return self._truth

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, value):
        self._output = os.path.abspath(value)


    @property
    def bigdata_root(self):
        return self._bigdata_root

    @bigdata_root.setter
    def bigdata_root(self, value):
        return NotImplementedError("Set TEST_BIGDATA environment variable "
            "to change this value.")

    # The methods
    def get_data(self, path):
        """Copy data from Artifactory remote resource to the CWD

        Updates self.input on completion
        """
        local_path = get_bigdata(self._inputs_root, self._env,
            os.path.dirname(path), os.path.basename(path), docopy=self.docopy)
        self._input = local_path
        return local_path

    def get_truth(self, path=None):
        """Copy truth data from Artifactory remote resource to the CWD/truth

        Updates self.truth on completion
        """
        if path is None:
            pass
        else:
            os.makedirs('truth', exist_ok=True)
            os.chdir('truth')
            try:
                local_truth = get_bigdata(self._inputs_root, self._env,
                    os.path.dirname(path), os.path.basename(path), docopy=self.docopy)
            except BigdataError:
                os.chdir('..')
                raise
            os.chdir('..')
            self._truth = local_truth
            return local_truth

    def get_association(self, asn):
        return NotImplementedError

    def generate_artifactory_json(self):
        return NotImplementedError

    def okify_truth(self):
        return NotImplementedError

    def okify_input(self):
        return NotImplementedError


@pytest.fixture(scope='module')
def rtdata_module(artifactory_repos, envopt):
    """Provides the RemoteResource class"""
    inputs_root, results_root = artifactory_repos
    resource = RegtestData(env=envopt, inputs_root=inputs_root,
        results_root=results_root)

    return resource


@pytest.fixture(scope='function')
def rtdata(artifactory_repos, envopt):
    """Provides the RemoteResource class"""
    inputs_root, results_root = artifactory_repos
    resource = RegtestData(env=envopt, inputs_root=inputs_root,
        results_root=results_root)

    return resource


@pytest.fixture
def fitsdiff_defaults():
    return dict(
        ignore_hdus=['ASDF'],
        ignore_keywords=['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
        rtol=0.00001,
        atol=0.0000001,
    )

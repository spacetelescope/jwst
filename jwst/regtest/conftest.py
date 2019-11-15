from datetime import datetime
import os
import copy
import json

import getpass
import pytest
from ci_watson.artifactory_helpers import (
    get_bigdata_root,
    get_bigdata,
    BigdataError,
    UPLOAD_SCHEMA,
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
    """Create request.node.report_setup and request.node.report_call to be
    used in generate_artifactory_json() fixture.
    """
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()

    # set a report attribute for each phase of a call, which can
    # be "setup", "call", "teardown"
    setattr(item, "report_" + rep.when, rep)


def postmortem(request, want_property):
    if isinstance(want_property, str):
        want_property = [want_property]

    result = {}
    if request.node.report_setup.passed:
        if request.node.report_call.failed:
            for prop in request.node.user_properties:
                name, data = prop
                for wanted in want_property:
                    if name == wanted:
                        result.update([prop])
    return result


@pytest.fixture(scope='function', autouse=True)
def generate_artifactory_json(request, artifactory_repos):
    inputs_root, results_root = artifactory_repos

    def _func(schema_pattern):
        # Generate the Artifactory target path
        whoami = getpass.getuser() or 'nobody'
        user_tag = 'NOT_CI_{}'.format(whoami)
        build_tag = os.environ.get('BUILD_TAG', user_tag)
        build_matrix_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', '0')
        subdir = '{}_{}_{}'.format(TODAYS_DATE, build_tag, build_matrix_suffix)
        testname = request.node.originalname or request.node.name
        remote_results_path = os.path.join(results_root, subdir, testname) + os.sep

        # Generate an upload schema
        return generate_upload_schema(schema_pattern, remote_results_path)

    yield
    # Execute the following at test teardown
    schema_pattern = []

    props = postmortem(request, 'output')
    if props:
        path = os.path.abspath(props['output'])
        schema_pattern.append(path)
        cwd, _ = os.path.split(path)

        upload_schema = _func(schema_pattern)

        # Write the schema to JSON
        jsonfile = os.path.join(cwd, "{}_results.json".format(request.node.name))
        with open(jsonfile, 'w') as outfile:
            json.dump(upload_schema, outfile, indent=2)


# @pytest.fixture(scope='function', autouse=True)
def generate_artifactory_okify_json(request, artifactory_repos):
    inputs_root, results_root = artifactory_repos

    yield
    # Execute the following at test teardown
    schema_pattern = []
    props = postmortem(request, 'output')
    if props:
        schema_pattern.append(os.path.abspath(props['output']))
        generate_okify_schema(schema_pattern)


def generate_upload_schema(pattern, target, recursive=False):
    """
    Write out JSON file to upload Jenkins results from test to
    Artifactory storage area.

    This function relies on the JFROG JSON schema for uploading data into
    artifactory using the Jenkins plugin.  Docs can be found at
    https://www.jfrog.com/confluence/display/RTF/Using+File+Specs

    Parameters
    ----------
    pattern : str or list of strings
        Specifies the local file system path to test results which should be
        uploaded to Artifactory. You can specify multiple artifacts by using
        wildcards or a regular expression as designated by the regexp property.

    target : str
        Specifies the target path in Artifactory in the following format::

            [repository_name]/[repository_path]

    recursive : bool, optional
        Specify whether or not to identify files listed in sub-directories
        for uploading.  Default: `False`

    """
    recursive = repr(recursive).lower()

    if not isinstance(pattern, str):
        # Populate schema for this test's data
        upload_schema = {"files": []}

        for p in pattern:
            temp_schema = copy.deepcopy(UPLOAD_SCHEMA["files"][0])
            temp_schema.update({"pattern": p, "target": target,
                                "recursive": recursive})
            upload_schema["files"].append(temp_schema)
    else:
        # Populate schema for this test's data
        upload_schema = copy.deepcopy(UPLOAD_SCHEMA)
        upload_schema["files"][0].update({"pattern": pattern, "target": target,
                                          "recursive": recursive})
    return upload_schema


# The following function can eventually be moved to ci-watson perhaps
def generate_okify_schema(pattern, target, testname, recursive=False):
    """
    Build an upload schema for updating truth files on Artifactory.

    This function uses Artifactory JSON spec files for uploading data using
    the Artifactory Jenkins plugin.  Docs can be found at
    https://www.jfrog.com/confluence/display/RTF/Using+File+Specs

    Parameters
    ----------
    pattern : str or list of strings
        Specifies the local file system path to test results which should be
        uploaded to Artifactory. You can specify multiple artifacts by using
        wildcards or a regular expression as designated by the regexp property.
    target : str
        Specifies the target path in Artifactory in the following format::
            [repository_name]/[repository_path]
    testname : str
        Name of test that generate the results. This will be used to create the
        name of the JSON file to enable these results to be uploaded to
        Artifactory.
    recursive : bool, optional
        Specify whether or not to identify files listed in sub-directories
        for uploading.  Default: `False`
    """
    jsonfile = "{}_update_truth.json".format(testname)
    recursive = repr(recursive).lower()

    if not isinstance(pattern, str):
        # Populate schema for this test's data
        upload_schema = {"files": []}

        for p in pattern:
            temp_schema = copy.deepcopy(UPLOAD_SCHEMA["files"][0])
            temp_schema.update({"pattern": p, "target": target,
                                "recursive": recursive})
            upload_schema["files"].append(temp_schema)

    else:
        # Populate schema for this test's data
        upload_schema = copy.deepcopy(UPLOAD_SCHEMA)
        upload_schema["files"][0].update({"pattern": pattern, "target": target,
                                          "recursive": recursive})

    # Write out JSON file with description of test results
    with open(jsonfile, 'w') as outfile:
        json.dump(upload_schema, outfile, indent=2)


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
        self._truth_remote = value.split(os.sep)

    @property
    def input(self):
        return self._input

    @input.setter
    def input(self, value):
        self._input = os.path.abspath(value)

    @property
    def truth(self):
        return self._truth

    @truth.setter
    def truth(self, value):
        self._truth = os.path.abspath(value)

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
    def get_data(self, path=None):
        """Copy data from Artifactory remote resource to the CWD

        Updates self.input and self.input_remote upon completion
        """
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path
        self.input = get_bigdata(self._inputs_root, self._env,
            os.path.dirname(path), os.path.basename(path), docopy=self.docopy)

        return self.input

    def get_truth(self, path=None):
        """Copy truth data from Artifactory remote resource to the CWD/truth

        Updates self.truth and self.truth_remote on completion
        """
        if path is None:
            path = self.truth_remote
        else:
            self.truth_remote = path
        os.makedirs('truth', exist_ok=True)
        os.chdir('truth')
        try:
            self.truth = get_bigdata(self._inputs_root, self._env,
                os.path.dirname(path), os.path.basename(path), docopy=self.docopy)
        except BigdataError:
            os.chdir('..')
            raise
        os.chdir('..')

        return self.truth

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

    yield resource


@pytest.fixture(scope='function')
def rtdata(artifactory_repos, envopt):
    """Provides the RemoteResource class"""
    inputs_root, results_root = artifactory_repos
    resource = RegtestData(env=envopt, inputs_root=inputs_root,
        results_root=results_root)

    yield resource


@pytest.fixture
def fitsdiff_default_kwargs():
    return dict(
        ignore_hdus=['ASDF'],
        ignore_keywords=['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
        rtol=0.00001,
        atol=0.0000001,
    )

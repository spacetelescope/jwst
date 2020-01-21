from datetime import datetime
import os
import copy
import json

import getpass
import pytest
from ci_watson.artifactory_helpers import UPLOAD_SCHEMA
from astropy.table import Table
from numpy.testing import assert_allclose

from .regtestdata import RegtestData


TODAYS_DATE = datetime.now().strftime("%Y-%m-%d")


@pytest.fixture(scope="session")
def artifactory_repos(pytestconfig):
    """Provides Artifactory inputs_root and results_root"""
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


def postmortem(request, fixturename):
    """Retrieve a fixture object if a test failed
    """
    if request.node.report_setup.passed:
        if request.node.report_call.failed:
            return request.node.funcargs.get(fixturename, None)


@pytest.fixture(scope='function', autouse=True)
def generate_artifactory_json(request, artifactory_repos):
    """Pytest fixture that leaves behind JSON upload and okify specfiles
    if the rtdata or rtdata_module fixtures are used in either a test or a
    module-scoped fixture that runs a pipeline and provides the results to a
    series of test.
    """
    inputs_root, results_root = artifactory_repos

    def artifactory_result_path():
        # Generate the Artifactory result path
        whoami = getpass.getuser() or 'nobody'
        user_tag = 'NOT_CI_{}'.format(whoami)
        build_tag = os.environ.get('BUILD_TAG', user_tag)
        build_matrix_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', '0')
        subdir = '{}_{}_{}'.format(TODAYS_DATE, build_tag, build_matrix_suffix)
        testname = request.node.originalname or request.node.name

        return os.path.join(results_root, subdir, testname) + os.sep

    yield
    # Execute the following at test teardown
    upload_schema_pattern = []
    okify_schema_pattern = []

    rtdata = postmortem(request, 'rtdata') or postmortem(request, 'rtdata_module')
    if rtdata:
        try:
            # The _jail fixture from ci_watson sets tmp_path
            cwd = str(request.node.funcargs['tmp_path'])
        except KeyError:
            # The jail fixture (module-scoped) returns the path
            cwd = str(request.node.funcargs['jail'])
        rtdata.remote_results_path = artifactory_result_path()
        rtdata.test_name = request.node.name
        # Dump the failed test traceback into rtdata
        rtdata.traceback = str(request.node.report_call.longrepr)

        upload_schema_pattern.append(rtdata.input)
        upload_schema_pattern.append(rtdata.output)
        upload_schema = generate_upload_schema(upload_schema_pattern,
            rtdata.remote_results_path)

        # Write the upload schema to JSON file
        jsonfile = os.path.join(cwd, f"{request.node.name}_results.json")
        with open(jsonfile, 'w') as outfile:
            json.dump(upload_schema, outfile, indent=2)

        pattern = os.path.join(rtdata.remote_results_path, os.path.basename(rtdata.output))
        okify_schema_pattern.append(pattern)
        okify_schema = generate_upload_schema(okify_schema_pattern, rtdata.truth_remote)

        # Write the okify schema to JSON file
        jsonfile = os.path.join(cwd, f"{request.node.name}_okify.json")
        with open(jsonfile, 'w') as outfile:
            json.dump(okify_schema, outfile, indent=2)

        # Write the rtdata class out as an ASDF file
        path = os.path.join(cwd, f"{request.node.name}_rtdata.asdf")
        rtdata.to_asdf(path)


def generate_upload_schema(pattern, target, recursive=False):
    """
    Generate JSON schema for Artifactory upload specfile using JFROG.

    Docs can be found at
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

    Returns
    -------
    upload_schema : dict
        Dictionary specifying the upload schema
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


def _rtdata_fixture_implementation(artifactory_repos, envopt, request):
    """Provides the RemoteResource class"""
    inputs_root, results_root = artifactory_repos
    rtdata = RegtestData(env=envopt, inputs_root=inputs_root,
        results_root=results_root)

    yield rtdata


@pytest.fixture(scope='function')
def rtdata(artifactory_repos, envopt, request, _jail):
    yield from _rtdata_fixture_implementation(artifactory_repos, envopt, request)


@pytest.fixture(scope='module')
def rtdata_module(artifactory_repos, envopt, request, jail):
    yield from _rtdata_fixture_implementation(artifactory_repos, envopt, request)


@pytest.fixture
def fitsdiff_default_kwargs():
    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX',
        'NAXIS1', 'TFORM*']
    return dict(
        ignore_hdus=['ASDF'],
        ignore_keywords=ignore_keywords,
        ignore_fields=ignore_keywords,
        rtol=1e-5,
        atol=1e-7,
    )


@pytest.fixture
def diff_astropy_tables():
    """Compare astropy tables with tolerances for float columns."""

    def _diff_astropy_tables(result_path, truth_path, rtol=1e-5, atol=1e-7):
        result = Table.read(result_path)
        truth = Table.read(truth_path)

        diffs = []

        if result.colnames != truth.colnames:
            diffs.append("Column names (or order) do not match")

        if len(result) != len(truth):
            diffs.append("Row count does not match")

        # If either the columns or the row count is mismatched, then don't
        # bother checking the individual column values.
        if len(diffs) > 0:
            return diffs

        if result.meta != truth.meta:
            diffs.append("Metadata does not match")

        for col_name in truth.colnames:
            try:
                if result[col_name].dtype != truth[col_name].dtype:
                    diffs.append(f"Column '{col_name}' dtype does not match")
                    continue

                dtype = truth[col_name].dtype
                if dtype.kind == "f":
                    try:
                        assert_allclose(result[col_name], truth[col_name],
                            rtol=rtol, atol=atol)
                    except AssertionError as err:
                        diffs.append(
                            f"Column '{col_name}' values do not match (within tolerances) \n{str(err)}"
                        )
                else:
                    if not (result[col_name] == truth[col_name]).all():
                        diffs.append(f"Column '{col_name}' values do not match")
            except AttributeError:
                # Ignore case where a column does not have a dtype, as in the case
                # of SkyCoord objects
                pass

        return diffs

    return _diff_astropy_tables

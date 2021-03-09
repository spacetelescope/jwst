from datetime import datetime
import copy
import json
import os
from pathlib import Path

import getpass
import pytest
from ci_watson.artifactory_helpers import UPLOAD_SCHEMA
from astropy.table import Table
from numpy.testing import assert_allclose, assert_equal
from astropy.io.fits import conf

from jwst.regtest.regtestdata import RegtestData
from jwst.regtest.sdp_pools_source import SDPPoolsSource


TODAYS_DATE = datetime.now().strftime("%Y-%m-%d")

# Turn of FITS memmap for all regtests (affects FITSDiff)
conf.use_memmap = False


@pytest.fixture(scope="session")
def artifactory_repos(pytestconfig):
    """Provides Artifactory inputs_root and results_root"""
    try:
        inputs_root = pytestconfig.getini('inputs_root')[0]
        results_root = pytestconfig.getini('results_root')[0]
    except IndexError:
        inputs_root = "jwst-pipeline"
        results_root = "jwst-pipeline-results"
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
    """Retrieve a fixture object if a test failed, else return None
    """
    try:
        if request.node.report_setup.passed:
            try:
                if request.node.report_call.failed:
                    return request.node.funcargs.get(fixturename, None)
            # Handle case where `report_call` hasn't been generated yet
            # because the test hasn't finished running, as in the case of
            # a user interrupt
            except AttributeError:
                return None
    except AttributeError:
        return None


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

        # Upload and allow okify of truth by rtdata.output, if the test did not
        # fail before producing rtdata.output
        if rtdata.output and os.path.exists(rtdata.output):
            # Write the rtdata class out as an ASDF file
            path_asdf = os.path.join(cwd, f"{request.node.name}_rtdata.asdf")
            rtdata.to_asdf(path_asdf)

            # Generate an OKify JSON file
            pattern = os.path.join(rtdata.remote_results_path,
                os.path.basename(rtdata.output))
            okify_schema_pattern.append(pattern)
            okify_schema = generate_upload_schema(okify_schema_pattern,
                f"{os.path.dirname(rtdata.truth_remote)}/")

            jsonfile = os.path.join(cwd, f"{request.node.name}_okify.json")
            with open(jsonfile, 'w') as fd:
                json.dump(okify_schema, fd, indent=2)

            # Generate an upload JSON file, including the OKify, asdf file
            upload_schema_pattern.append(rtdata.output)
            upload_schema_pattern.append(os.path.abspath(jsonfile))
            upload_schema_pattern.append(path_asdf)
            upload_schema = generate_upload_schema(upload_schema_pattern,
                rtdata.remote_results_path)

            jsonfile = os.path.join(cwd, f"{request.node.name}_results.json")
            with open(jsonfile, 'w') as fd:
                json.dump(upload_schema, fd, indent=2)


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
        __tracebackhide__ = True
        result = Table.read(result_path)
        truth = Table.read(truth_path)

        diffs = []

        try:
            assert result.colnames == truth.colnames
        except AssertionError as err:
            diffs.append(f"Column names (or order) do not match\n{err}")

        try:
            assert len(result) == len(truth)
        except AssertionError as err:
            diffs.append(f"Row count does not match\n{err}")

        # If either the columns or the row count is mismatched, then don't
        # bother checking the individual column values.
        if len(diffs) > 0:
            raise AssertionError("\n".join(diffs))

        # Disable meta comparison for now, until we're able to specify
        # individual entries for comparison
        #if result.meta != truth.meta:
        #    diffs.append("Metadata does not match")

        for col_name in truth.colnames:
            try:
                try:
                    assert result[col_name].dtype == truth[col_name].dtype
                except AssertionError as err:
                    diffs.append(f"Column '{col_name}' dtype does not match\n{err}")
                    continue

                if truth[col_name].dtype.kind == "f":
                    try:
                        assert_allclose(result[col_name], truth[col_name],
                            rtol=rtol, atol=atol)
                    except AssertionError as err:
                        diffs.append("\n----------------------------------\n"
                            + f"Column '{col_name}' values do not "
                            + f"match (within tolerances) \n{err}"
                        )
                else:
                    try:
                        assert_equal(result[col_name], truth[col_name])
                    except AssertionError as err:
                        diffs.append(f"Column '{col_name}' values do not match\n{err}")
            except AttributeError:
                # Ignore case where a column does not have a dtype, as in the case
                # of SkyCoord objects
                pass

        if len(diffs) != 0:
            raise AssertionError("\n".join(diffs))

        # No differences
        return True

    return _diff_astropy_tables

# Add option to specify a single pool name
def pytest_addoption(parser):
    parser.addoption(
        '--sdp-pool', metavar='sdp_pool', default=None,
        help='SDP test pool to run. Specify the name only, not extension or path'
    )
    parser.addoption(
        '--standard-pool', metavar='standard_pool', default=None,
        help='Standard test pool to run. Specify the name only, not extension or path'
    )


@pytest.fixture
def sdp_pool(request):
    """Retrieve a specific SDP pool to test"""
    return request.config.getoption('--sdp-pool')


@pytest.fixture
def standard_pool(request):
    """Retrieve a specific standard pool to test"""
    return request.config.getoption('--standard-pool')


def pytest_generate_tests(metafunc):
    """Prefetch and parametrize a set of test pools"""
    if 'pool_path' in metafunc.fixturenames:
        try:
            SDPPoolsSource.inputs_root = metafunc.config.getini('inputs_root')[0]
            SDPPoolsSource.results_root = metafunc.config.getini('results_root')[0]
            SDPPoolsSource.env = metafunc.config.getoption('env')
        except IndexError:
            SDPPoolsSource.inputs_root = "jwst-pipeline"
            SDPPoolsSource.results_root = "jwst-pipeline-results"
            SDPPoolsSource.env = "dev"

        pools = SDPPoolsSource()

        try:
            pool_paths = pools.pool_paths
        except Exception:
            pool_paths = []

        ids = [
            Path(pool_path).stem
            for pool_path in pool_paths
        ]

        metafunc.parametrize('pool_path', pool_paths, ids=ids)

import copy
import getpass
import json
import os
from datetime import datetime

import pytest
from astropy.table import Table
from astropy.io.fits import conf
from ci_watson.artifactory_helpers import UPLOAD_SCHEMA
from numpy.testing import assert_allclose, assert_equal

from jwst.regtest.regtestdata import RegtestData

TODAYS_DATE = datetime.now().strftime("%Y-%m-%d")

# Turn of FITS memmap for all regtests (affects FITSDiff)
conf.use_memmap = False


@pytest.fixture(scope="session")
def artifactory_repos(pytestconfig):
    """
    Provide Artifactory inputs_root and results_root.

    Parameters
    ----------
    pytestconfig : pytest.config.Config
        Pytest configuration object.

    Returns
    -------
    inputs_root : str
        The Artifactory inputs root, set to "jwst-pipeline" unless found in the config.
    results_root : str
        The Artifactory results root, set to "jwst-pipeline-results" unless found in the config.
    """
    inputs_root = pytestconfig.getini("inputs_root")
    # in pytest 8 inputs_root will be None
    # in pytest <8 inputs_root will be []
    # see: https://github.com/pytest-dev/pytest/pull/11594
    # using "not inputs_root" will handle both cases
    if not inputs_root:
        inputs_root = "jwst-pipeline"
    else:
        inputs_root = inputs_root[0]

    results_root = pytestconfig.getini("results_root")
    # see not above about inputs_root
    if not results_root:
        results_root = "jwst-pipeline-results"
    else:
        results_root = results_root[0]
    return inputs_root, results_root


@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """
    Create request.node.report_setup and request.node.report_call.

    These are used in the generate_artifactory_json() fixture.

    Parameters
    ----------
    item : pytest.Item
        The pytest item object.
    call : pytest.CallInfo
        The pytest call object.
    """
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()

    # set a report attribute for each phase of a call, which can
    # be "setup", "call", "teardown"
    setattr(item, "report_" + rep.when, rep)


def postmortem(request, fixturename):
    """
    Retrieve a fixture object if a test failed, else return None.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The pytest fixture request object.
    fixturename : str
        The name of the fixture to retrieve.

    Returns
    -------
    fixture : object or None
        The fixture object if the test failed, else None.
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


def pytest_collection_modifyitems(config, items):
    """
    Allow artifactory_result_path to be unique for each test.

    Add the index of each item in the list of items.
    This is used below for artifactory_result_path
    to produce a unique result subdirectory for
    each test (even if that test shares a name with
    another test which is the case for parametrized tests).

    Parameters
    ----------
    config : _pytest.config.Config
        Pytest configuration object.
    items : list
        List of pytest items.
    """
    for i, item in enumerate(items):
        item.index = i


@pytest.fixture(scope="function", autouse=True)
def generate_artifactory_json(request, artifactory_repos):
    """
    Create JSON files for Artifactory upload and OKify.

    Pytest fixture that leaves behind JSON upload and okify specfiles
    if the rtdata or rtdata_module fixtures are used in either a test or a
    module-scoped fixture that runs a pipeline and provides the results to a
    series of test.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The pytest fixture request object.
    artifactory_repos : tuple(str, str)
        Tuple of the Artifactory inputs_root and results_root.

    Yields
    ------
    None
    """
    _, results_root = artifactory_repos

    def artifactory_result_path():
        # Generate the Artifactory result path
        whoami = getpass.getuser() or "nobody"
        user_tag = f"NOT_CI_{whoami}"
        build_tag = os.environ.get("BUILD_TAG", user_tag)
        build_matrix_suffix = os.environ.get("BUILD_MATRIX_SUFFIX", "0")
        subdir = f"{TODAYS_DATE}_{build_tag}_{build_matrix_suffix}"
        testname = request.node.originalname or request.node.name
        basename = f"{request.node.index}_{testname}"

        return os.path.join(results_root, subdir, basename) + os.sep

    yield
    # Execute the following at test teardown
    upload_schema_pattern = []
    okify_schema_pattern = []

    rtdata = postmortem(request, "rtdata") or postmortem(request, "rtdata_module")
    if rtdata:
        try:
            # The tmp_cwd fixture sets tmp_path
            cwd = str(request.node.funcargs["tmp_path"])
        except KeyError:
            # The tmp_cwd_module fixture (module-scoped) returns the path
            cwd = str(request.node.funcargs["tmp_cwd_module"])
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
            if rtdata.okify_op == "sdp_pool_copy":
                pattern = (
                    os.path.join(rtdata.remote_results_path, os.path.basename(rtdata.output))
                    + "/jw*.json"
                )
            else:
                pattern = os.path.join(rtdata.remote_results_path, os.path.basename(rtdata.output))
            okify_schema_pattern.append(pattern)
            okify_schema = generate_upload_schema(
                okify_schema_pattern, f"{os.path.dirname(rtdata.truth_remote)}/"
            )

            jsonfile = os.path.join(cwd, f"{request.node.name}_okify.json")
            with open(jsonfile, "w") as fd:
                json.dump(okify_schema, fd, indent=2)

            # Generate an upload JSON file, including the OKify, asdf file
            upload_schema_pattern.append(os.path.abspath(jsonfile))
            upload_schema_pattern.append(path_asdf)
            upload_schema = generate_upload_schema(
                upload_schema_pattern, rtdata.remote_results_path
            )

            if rtdata.okify_op == "file_copy":
                upload_schema = generate_upload_schema(
                    rtdata.output, rtdata.remote_results_path, schema=upload_schema
                )
            elif rtdata.okify_op in ("folder_copy", "sdp_pool_copy"):
                output = rtdata.output + "/"
                target = rtdata.remote_results_path + os.path.basename(rtdata.output) + "/"
                upload_schema = generate_upload_schema(output, target, schema=upload_schema)
            else:
                raise RuntimeError(f"Unknown artifactory operation: {rtdata.okify_op}")

            jsonfile = os.path.join(cwd, f"{request.node.name}_results.json")
            with open(jsonfile, "w") as fd:
                json.dump(upload_schema, fd, indent=2)


def generate_upload_schema(pattern, target, recursive=False, schema=None):
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
    schema : dict
        Existing schema to append to.

    Returns
    -------
    upload_schema : dict
        Dictionary specifying the upload schema
    """
    recursive = repr(recursive).lower()
    if schema is None:
        upload_schema = {"files": []}
    else:
        upload_schema = copy.deepcopy(schema)
    if isinstance(pattern, str):
        pattern = [pattern]

    for p in pattern:
        temp_schema = copy.deepcopy(UPLOAD_SCHEMA["files"][0])
        temp_schema.update({"pattern": p, "target": target, "recursive": recursive})
        upload_schema["files"].append(temp_schema)

    return upload_schema


def _rtdata_fixture_implementation(artifactory_repos, envopt, request):
    """
    Provide the RegtestData class.

    Parameters
    ----------
    artifactory_repos : tuple(str, str)
        Tuple of the Artifactory inputs_root and results_root.
    envopt : str
        The Artifactory environment, e.g, "dev".
    request : pytest.FixtureRequest
        The pytest fixture request object.

    Returns
    -------
    `jwst.regtest.regtestdata.RegtestData`
        RegtestData class instance.
    """
    inputs_root, results_root = artifactory_repos
    return RegtestData(env=envopt, inputs_root=inputs_root, results_root=results_root)


@pytest.fixture(scope="function")
def rtdata(artifactory_repos, envopt, request, tmp_cwd):
    """
    Provide access to regression test data via a function-scoped fixture.

    Parameters
    ----------
    artifactory_repos : tuple(str, str)
        Tuple of the Artifactory inputs_root and results_root.
    envopt : str
        The Artifactory environment, e.g, "dev".
    request : pytest.FixtureRequest
        The pytest fixture request object.
    tmp_cwd : pytest fixture
        Pytest fixture that creates a function-scoped temporary working directory.

    Returns
    -------
    `jwst.regtest.regtestdata.RegtestData`
        A RegtestData object providing access to regression test data.
    """
    return _rtdata_fixture_implementation(artifactory_repos, envopt, request)


@pytest.fixture(scope="module")
def rtdata_module(artifactory_repos, envopt, request, tmp_cwd_module):
    """
    Provide access to regression test data via a module-scoped fixture.

    Parameters
    ----------
    artifactory_repos : tuple(str, str)
        Tuple of the Artifactory inputs_root and results_root.
    envopt : str
        The Artifactory environment, e.g, "dev".
    request : pytest.FixtureRequest
        The pytest fixture request object.
    tmp_cwd_module : pytest fixture
        Pytest fixture that creates a module-scoped temporary working directory.

    Returns
    -------
    `jwst.regtest.regtestdata.RegtestData`
        A RegtestData object providing access to regression test data.
    """
    return _rtdata_fixture_implementation(artifactory_repos, envopt, request)


@pytest.fixture
def fitsdiff_default_kwargs():
    """
    Provide default tolerances and ignores for FITSDiff comparisons.

    Returns
    -------
    dict
        Keyword arguments to pass into FITSDiff.
    """
    ignore_keywords = ["DATE", "CAL_VER", "CAL_VCS", "CRDS_VER", "CRDS_CTX", "NAXIS1", "TFORM*"]
    return {
        "ignore_hdus": ["ASDF"],
        "ignore_keywords": ignore_keywords,
        "ignore_fields": ignore_keywords,
        "rtol": 1e-5,
        "atol": 1e-7,
    }


@pytest.fixture
def diff_astropy_tables():
    """
    Compare astropy tables with tolerances for float columns.

    Returns
    -------
    function
        Function to compare two astropy tables.
    """

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
        # if result.meta != truth.meta:
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
                        assert_allclose(result[col_name], truth[col_name], rtol=rtol, atol=atol)
                    except AssertionError as err:
                        diffs.append(
                            "\n----------------------------------\n"
                            f"Column '{col_name}' values do not "
                            f"match (within tolerances) \n{err}"
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

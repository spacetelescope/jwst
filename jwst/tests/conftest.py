from datetime import datetime
import os
from urllib.parse import urlparse
from glob import glob

import pytest
import requests
from ci_watson.artifactory_helpers import (
    generate_upload_schema,
    get_bigdata_root,
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
        results_path = os.path.join(results_root, subdir, request.node.originalname) + os.sep

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


@pytest.fixture(scope="session")
def update_local_bigdata_from_artifactory(artifactory_repos):
    """Update files in locally-defined TEST_BIGDATA from Artifactory source"""
    inputs_root, _ = artifactory_repos
    params = {'name': "*.*", 'repos': inputs_root}
    search_url = os.path.join(get_bigdata_root(), 'api/search/artifact')
    with requests.get(search_url, params=params) as r:
        artifacts = [a['uri'] for a in r.json()['results']]
    paths = [urlparse(a).path.replace('/artifactory/api/storage/' + inputs_root + '/', '') for a in artifacts]
    local_paths = [os.path.join(_bigdata, path) for path in paths]

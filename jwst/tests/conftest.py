import os
from urllib.parse import urlparse

import pytest
import requests


RT_SERVER = "https://bytesalad.stsci.edu/artifactory/"
REPO = "scsb-jwst-pipeline"


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
def update_local_bigdata_from_artifactory(_bigdata):
    "Update files in locally-defined TEST_BIGDATA from Artifactory source"
    params = {'name': "*.*", 'repos': REPO}
    search_url = os.path.join(RT_SERVER, 'api/search/artifact')
    with requests.get(search_url, params=params) as r:
        artifacts = [a['uri'] for a in r.json()['results']]
    paths = [urlparse(a).path.replace('/artifactory/api/storage/' + REPO + '/', '') for a in artifacts]
    local_paths = [os.path.join(_bigdata, path) for path in paths]


@pytest.fixture(scope="session")
def artifactory_repos(pytestconfig):
    inputs_root = pytestconfig.getini('inputs_root')[0]
    results_root = pytestconfig.getini('results_root')[0]

    return inputs_root, results_root


@pytest.fixture(scope="function")
def artifactory_upload_json(request):
    yield
    

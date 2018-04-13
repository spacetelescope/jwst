import os
import pytest
import re
import requests


RE_URL = re.compile('\w+://\S+')


def pytest_addoption(parser):
    parser.addoption('--bigdata', action='store_true',
                     help='Big local datasets')


def check_url(url):
    if RE_URL.match(url) is None:
        return False

    r = requests.head(url, allow_redirects=True)
    if r.status_code >= 400:
        return False
    return True


class BigdataError(Exception):
    pass


@pytest.fixture
def _bigdata():
    origins = [os.environ.get('TEST_BIGDATA', ''),
               '/data4/jwst_test_data',
               'https://bytesalad.stsci.edu/artifactory/jwst_pipeline']
    for path in origins:
        if os.path.exists(path) or check_url(path):
            return path

    raise BigdataError('Data files are not available.')


@pytest.fixture(scope='function')
def _jail(tmpdir):
    path = str(tmpdir)
    os.chdir(path)
    yield

"""Project default for pytest"""
import os
import tempfile
import pytest
import re
import requests

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions

# Uncomment the following line to treat all DeprecationWarnings as exceptions
enable_deprecations_as_exceptions()

try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['asdf'] = 'asdf'
    del PYTEST_HEADER_MODULES['h5py']
except (NameError, KeyError):
    pass

pytest_plugins = [
    'asdf.tests.schema_tester'
]

RE_URL = re.compile('\w+://\S+')

def check_url(url):
    """ Determine if `url` can be resolved without error
    """
    if RE_URL.match(url) is None:
        return False

    r = requests.head(url, allow_redirects=True)
    if r.status_code >= 400:
        return False
    return True


@pytest.fixture
def mk_tmp_dirs():
    """Create a set of temporary directorys and change to one of them."""
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()
    try:
        os.chdir(tmp_current_path)
        yield (tmp_current_path, tmp_data_path, tmp_config_path)
    finally:
        os.chdir(old_path)


class BigdataError(Exception):
    pass

import os

# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions

# Uncomment the following line to treat all DeprecationWarnings as
# exceptions
enable_deprecations_as_exceptions()

try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['asdf'] = 'asdf'
    del PYTEST_HEADER_MODULES['h5py']
except (NameError, KeyError):
    pass

# This is to figure out the affiliated package version, rather than
# using Astropy's
#from . import version

#try:
#    packagename = os.path.basename(os.path.dirname(__file__))
#    TESTED_VERSIONS[packagename] = version.version
#except NameError:   # Needed to support Astropy <= 1.0.0
#    pass

pytest_plugins = [
    'asdf.tests.schema_tester'
]

"""Project default for pytest"""
from astropy.tests.helper import enable_deprecations_as_exceptions

enable_deprecations_as_exceptions()

pytest_plugins = [
    'asdf.tests.schema_tester'
]

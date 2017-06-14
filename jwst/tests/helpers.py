"""Handy helpful pytest helpers helping pytest test"""

import os
import pytest


# Decorator to indicate slow tests
runslow = pytest.mark.skipif(
    not pytest.config.getoption("--runslow"),
    reason="need --runslow option to run"
)


# Decorator to indicate TEST_BIGDATA required
require_bigdata = pytest.mark.skipif(
    'TEST_BIGDATA' not in os.environ,
    reason='"TEST_BIGDATA" environmental not defined. Cannot access test data.'
)


# Decorator to skip test if running under a TravisCI
not_under_travis = pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason='Temporarily disable due to performance issues'
)

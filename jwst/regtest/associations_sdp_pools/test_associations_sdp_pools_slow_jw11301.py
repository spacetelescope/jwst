"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw11301_pool(_jail, rtdata, resource_tracker, request):
    """Check that subarray name does not include underscores or dashes."""
    pool_args = ("jw11301_20261004t142009_pool", [])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

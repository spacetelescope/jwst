"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw01192_o008_pool(_jail, rtdata, resource_tracker, request):
    """Check imprint behavior."""
    pool_args = ("jw01192_o008_pool", ["--DMS", "-i", "o008"])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

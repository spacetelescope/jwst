"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw03855_20241103t042455_pool(_jail, rtdata, resource_tracker, request):
    pool_args = ("jw03855_20241103t042455_pool", [])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

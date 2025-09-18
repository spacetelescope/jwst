"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw01288_c1005_mostilno12_pool(_jail, rtdata, resource_tracker, request):
    """See JP-3230."""
    pool_args = (
        "jw01288_2025_c1005_mostilno12_pool",
        ["-i", "o003", "c1001", "c1005"],
    )
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

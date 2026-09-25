"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw01480_20250319t173819_pool(_jail, rtdata, resource_tracker, request):
    """NRC_GRISM, NRC_TSIMAGE, NRC_TSGRISM."""
    pool_args = ("jw01480_20250319t173819_pool", [])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

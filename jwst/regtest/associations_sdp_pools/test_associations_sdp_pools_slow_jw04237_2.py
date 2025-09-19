"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw04237_20250321t192812_dms_pool(_jail, rtdata, resource_tracker, request):
    """Lvl 3 MIRI MRS BKG with DMS."""
    pool_args = ("jw04237_20250321t192812_dms_pool", ["--DMS"])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

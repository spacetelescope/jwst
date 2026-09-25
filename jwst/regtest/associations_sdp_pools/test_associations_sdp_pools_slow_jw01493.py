"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_jw01493_20230307t040130_pool(_jail, rtdata, resource_tracker, request):
    pool_args = ("jw01493_20230307t040130_pool", [])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

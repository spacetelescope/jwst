"""
Test using SDP-generated pools (slow).

.. note::
    These are inflight equivalent approximate replacements for
    test_associations_standards.py test module.
"""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


def test_imprint_single_obs_1(_jail, rtdata, resource_tracker, request):
    pool_args = ("jw01928_20250318t105523_pool", [])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

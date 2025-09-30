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


def test_pool_021_tso_NIRCAM_TSIMAGE_TSGRISM(_jail, rtdata, resource_tracker, request):
    pool_args = ("jw02084_20250320t081615_pool", [])
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

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


def test_pool_034_wfss_parallel_NIRCAM_3ROW_2DIRECT(_jail, rtdata, resource_tracker, request):
    """Three row grism images and two direct images."""
    pool_args = (
        "jw05398_20250312t000159_3r2d_pool",
        ["-i", "o030", "o031", "o032", "o033", "o035", "c1012", "c1013", "c1014"],
    )
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

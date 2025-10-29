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


def test_pool_034_wfss_parallel_NIRCAM_ROW_ONLY(_jail, rtdata, resource_tracker, request):
    """Only row grism and direct image, but no column grism image."""
    pool_args = (
        "jw05398_20250312t000159_row_pool",
        ["-i", "o036", "o037", "c1029"],
    )
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)

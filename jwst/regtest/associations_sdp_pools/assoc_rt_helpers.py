"""Helper functions for associations regression testing."""

import os
import re
from glob import glob

from jwst.associations.lib.diff import compare_asn_files
from jwst.associations.main import Main as asn_generate

__all__ = ["parfunc", "assoc_sdp_against_standard"]

# Decompose pool name to retrieve proposal and version id.
pool_regex = re.compile(r"(?P<proposal>jw.+?)_(?P<versionid>.+)_pool")


def parfunc(x):
    """
    For readable param ids.

    Parameters
    ----------
    x : tuple
        Tuple of ``(pool, args)``.

    Returns
    -------
    param_id : str
        Just the pool.
    """
    return x[0]


def assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args):
    """Worker test function for associations regression tests."""
    pool, args = pool_args
    proposal, version_id = pool_regex.match(pool).group("proposal", "versionid")

    input_csv = rtdata.get_data(f"associations/sdp/pools/{pool}.csv")
    rtdata.output = os.curdir  # This test is jailed and parametrized
    rtdata.okify_op = "sdp_pool_copy"  # failures as folder content replacements
    # Create the associations
    with resource_tracker.track(log=request):
        asn_generate.cli(args + ["-p", rtdata.output, "--version-id", version_id, input_csv])
    out_paths = sorted(glob("*.json"))

    # Compare to the truth associations.
    truth_pool_path = f"truth/test_associations_sdp_pools/{pool}/"
    rtdata.truth_remote = truth_pool_path
    truth_paths = sorted(
        [rtdata.get_truth(p) for p in rtdata.data_glob(truth_pool_path, glob="*.json")]
    )
    if truth_paths == []:  # truth dir does not exist
        rtdata.truth_remote = f"{rtdata._inputs_root}/{rtdata.env}/{truth_pool_path}"
    compare_asn_files(out_paths, truth_paths)

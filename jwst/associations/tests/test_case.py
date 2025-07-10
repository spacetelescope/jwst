"""Test case insensitivity"""

import pytest
from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests.helpers import (
    level3_rule_path,
    combine_pools,
)
from jwst.associations.main import Main


@pytest.mark.parametrize(
    "partial_args, n_asns",
    [
        # Invalid ACID
        (["-i", "nosuchid"], 0),
        # Basic observation ACIDs
        (["-i", "o001"], 2),
        (["-i", "o002"], 2),
        # Specifying multiple ACIDs
        (["-i", "o001", "o002"], 4),
        # Candidate ID's
        (["-i", "c1000"], 2),
        (["-i", "c1000", "c1001"], 4),
        # Whole program
        ([], 22),
        # Discovered only
        (["--discover"], 2),
    ],
)
def test_candidate_observation_caseagnostic(partial_args, n_asns, tmp_path):
    """Use the extensive candidate test as a test for case"""
    pool_path = str(tmp_path / "pool.csv")
    pool = combine_pools(
        get_pkg_data_filename(
            "data/pool_001_candidates_lower.csv", package="jwst.associations.tests"
        )
    )
    pool.write(pool_path, format="ascii", delimiter="|")
    cmd_args = [
        pool_path,
        "--dry-run",
        "-r",
        level3_rule_path(),
        "--ignore-default",
    ]
    cmd_args.extend(partial_args)
    generated = Main.cli(cmd_args)
    assert len(generated.associations) == n_asns

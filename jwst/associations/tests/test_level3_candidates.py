"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import
import pytest

from .helpers import (
    generate_params,
    generate_pool_paths,
    level3_rule_path,
    t_path,
)

from ..main import Main

pool_path = pytest.yield_fixture(
    scope='module',
    params=['data/pool_001_candidates.csv']
)(generate_pool_paths)


pool_params = pytest.fixture(
    scope='module',
    params=[
        # Invalid ACID
        (
            ['-i', 'nosuchid'],
            0
        ),
        # Basic observation ACIDs
        (
            ['-i', 'o001'],
            2
        ),
        (
            ['-i', 'o002'],
            2
        ),
        # Specifying multiple ACIDs
        (
            ['-i', 'o001', 'o002'],
            4
        ),
        # Candidate ID's
        (
            ['-i', 'c1000'],
            2
        ),
        (
            ['-i', 'c1000', 'c1001'],
            4
        ),
        # Whole program
        (
            [],
            22
        ),
    ]
)(generate_params)


def test_candidate_observation(pool_path, pool_params):
    partial_args, n_asns = pool_params
    cmd_args = [
        pool_path,
        '--dry-run',
        '-r', level3_rule_path(),
        '--ignore-default',
    ]
    cmd_args.extend(partial_args)
    generated = Main(cmd_args)
    assert len(generated.associations) == n_asns

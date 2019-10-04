"""Test Level2 candidate operation"""
import pytest

from .helpers import (
    level2_rule_path,
    mkstemp_pool_file,
    t_path,
)

from ..main import Main


@pytest.mark.parametrize(
    "partial_args, n_asns",
    [
        # Invalid ACID
        (
            ['-i', 'nosuchid'],
            0
        ),
        # Basic observation ACIDs
        (
            ['-i', 'o001'],
            10
        ),
        # Whole program
        (
            [],
            70
        ),
        # Discovered only
        (
            ['--discover'],
            0
        ),
        # Candidates only
        (
            ['--all-candidates'],
            70
        ),
    ]
)
def test_candidate_observation(partial_args, n_asns):
    with mkstemp_pool_file(t_path('data/pool_001_candidates.csv')) as pool_path:
        cmd_args = [
            pool_path,
            '--dry-run',
            '-r', level2_rule_path(),
            '--ignore-default',
        ]
        cmd_args.extend(partial_args)
        generated = Main(cmd_args)
        assert len(generated.associations) == n_asns

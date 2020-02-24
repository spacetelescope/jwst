"""test_associations: Test of general Association functionality."""
import pytest

from jwst.associations import AssociationPool
from jwst.associations.main import Main


@pytest.mark.parametrize(
    'args',
    [
        [
            '--dry-run',
            '--discover',
            '--all-candidates',
            '-i', 'o001',
        ],
        [
            '--dry-run',
            '--discover',
            '--all-candidates',
        ],
        [
            '--dry-run',
            '--discover',
            '-i', 'o001',
        ],
        [
            '--dry-run',
            '--all-candidates',
            '-i', 'o001',
        ]
    ]
)
def test_toomanyoptions(args):
    """Test argument parsing for failures"""
    pool = AssociationPool()

    with pytest.raises(SystemExit):
        Main(args, pool=pool)


@pytest.mark.xfail(
    reason='Need to further investigate',
    run=False
)
def test_discovered(full_pool_rules):
    pool, rules, pool_fname = full_pool_rules

    full = Main([pool_fname, '--dry-run'])
    candidates = Main([pool_fname, '--dry-run', '--all-candidates'])
    discovered = Main([pool_fname, '--dry-run', '--discover'])
    assert len(full.associations) == len(candidates.associations) + len(discovered.associations)

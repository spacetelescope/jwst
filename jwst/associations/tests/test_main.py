"""test_associations: Test of general Association functionality."""
import pytest

from ..main import Main


def test_toomanyoptions(full_pool_rules):
    """Test argument parsing"""
    pool, rules, pool_fname = full_pool_rules

    with pytest.raises(SystemExit):
        Main([
            pool_fname,
            '--dry-run',
            '--discover',
            '--all-candidates',
            '-i', 'o001',
        ])
    with pytest.raises(SystemExit):
        Main([
            pool_fname,
            '--dry-run',
            '--discover',
            '--all-candidates',
        ])
    with pytest.raises(SystemExit):
        Main([
            pool_fname,
            '--dry-run',
            '--discover',
            '-i', 'o001',
        ])
    with pytest.raises(SystemExit):
        Main([
            pool_fname,
            '--dry-run',
            '--all-candidates',
            '-i', 'o001',
        ])


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

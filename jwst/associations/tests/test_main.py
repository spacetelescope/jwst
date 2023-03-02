"""test_associations: Test of general Association functionality."""
import os
import re
import subprocess

import pytest

from jwst.associations.tests.helpers import combine_pools, t_path

from jwst.associations import AssociationPool
from jwst.associations.main import Main


# Basic pool
POOL_PATH = 'pool_018_all_exptypes.csv'


@pytest.mark.parametrize(
    'case', [
        (None, 62),  # Don't re-run, just compare to the generator fixture results
        (['-i', 'o001'], 2),
        (['-i', 'o001', 'o002'], 3),
        (['-i', 'c1001'], 1),
        (['-i', 'o001', 'c1001'], 3),
        (['-i', 'c1001', 'c1002'], 2),
    ]
)
def test_asn_candidates(pool, all_candidates, case):
    """Test candidate selection option"""
    args, n_expected = case

    if args:
        generated = Main.cli(['--dry-run'] + args, pool=pool)
        n_actual = len(generated.associations)
    else:
        n_actual = len(all_candidates.associations)

    assert n_actual == n_expected


@pytest.mark.parametrize(
    'args, expected', [
        ([t_path(os.path.join('data', POOL_PATH))], 0),
        (['nosuchpool.csv'], 1),
    ]
)
def test_cmdline_status(args, expected, _jail):
    """Ensure command line status are as expected."""
    full_args = ['asn_generate'] + args
    status = subprocess.run(full_args)
    assert status.returncode == expected


@pytest.mark.parametrize(
    'version_id, expected', [
        ('', r'\d{3}t\d{6}'),
        ('mytestid', r'mytestid'),
    ]
)
def test_generate_version_id(version_id, expected, pool):
    """Check that an association has been given the appropriate version id"""
    regex = re.compile(expected)
    args = ['--dry-run', '-i', 'o001', '--version-id']
    if version_id:
        args.append(version_id)
    generated = Main.cli(args, pool=pool)
    for asn in generated.associations:
        assert regex.search(asn.asn_name)


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
        Main.cli(args, pool=pool)


# ########
# Fixtures
# ########
@pytest.fixture(scope='module')
def pool():
    """Retrieve pool path"""
    pool_path = t_path(os.path.join('data', POOL_PATH))
    pool = combine_pools(pool_path)

    return pool


@pytest.fixture(scope='module')
def all_candidates(pool):
    """"Retrieve the all exposure pool"""
    all_candidates = Main.cli(['--dry-run', '--all-candidates'], pool=pool)
    return all_candidates

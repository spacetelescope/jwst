from glob import glob
import os

import pytest
from astropy.table import Table

from jwst.associations.tests.helpers import TemporaryDirectory, combine_pools, t_path

from jwst.associations import load_asn
from jwst.associations.main import Main


# Basic pool
POOL_PATH = 'pool_018_all_exptypes.csv'


@pytest.fixture(scope='module')
def pool():
    """Retrieve pool path"""
    pool_path = t_path(os.path.join('data', POOL_PATH))
    pool = combine_pools(pool_path)

    return pool


@pytest.fixture(
    scope='module',
    params=['yaml', 'json']
)
def make_asns(pool, request):
    asn_format = request.param
    with TemporaryDirectory() as path:
        generated = Main.cli([
            '-p', path,
            '-i', 'o001',
            '--save-orphans',
            '--format', asn_format
        ], pool=pool)
        yield generated, path, asn_format


def test_roundtrip(make_asns):
    generated, path, asn_format = make_asns
    asn_files = glob(os.path.join(path, '*.' + asn_format))
    assert len(asn_files) == len(generated.associations)

    for asn_file in asn_files:
        with open(asn_file, 'r') as asn_fp:
            load_asn(asn_fp)

    orphaned_files = glob(os.path.join(path, '*.csv'))
    assert len(orphaned_files) == 1
    orphaned = Table.read(
        orphaned_files[0],
        format='ascii',
        delimiter='|'
    )
    assert len(orphaned) == len(generated.orphaned)


def test_load_asn_all(make_asns):
    generated, path, asn_format = make_asns
    asn_files = glob(os.path.join(path, '*.' + asn_format))
    assert len(asn_files) == len(generated.associations)

    for asn_file in asn_files:
        with open(asn_file, 'r') as asn_fp:
            asns = load_asn(asn_fp, registry=generated.rules, first=False)
        assert len(asns) > 1

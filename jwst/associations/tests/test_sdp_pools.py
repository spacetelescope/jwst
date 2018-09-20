"""Test using SDP-generated pools

Notes
-----
Most of the standard associations which are compared
against are built in the jupyter notebook

./notebooks/make_tests.ipynb
"""
from glob import glob
from os import path
import pytest

from .helpers import (
    compare_asns,
    runslow,
    t_path,
)

from .. import (AssociationPool, load_asn)
from ..main import Main

# Main test args
TEST_ARGS = ['--dry-run']

pool_paths = glob(t_path(path.join(
    'data', 'sdp', 'pools', '*.csv'
)))


@runslow
@pytest.mark.parametrize(
    'pool_path',
    pool_paths
)
def test_against_standard(pool_path):
    """Compare a generated assocaition against a standard
    """
    generated, standards = generate_asns(pool_path)
    assert len(generated) == len(standards)
    for asn in generated:
        for idx, standard in enumerate(standards):
            try:
                compare_asns(asn, standard)
            except AssertionError as e:
                last_err = e
            else:
                del standards[idx]
                break
        else:
            assert False, '{}'.format(last_err)


def generate_asns(pool_path):
    pool_dir, pool_root = path.split(pool_path)
    pool_root, pool_ext = path.splitext(pool_root)
    pool = AssociationPool.read(pool_path)

    standards = []
    asn_paths = glob(t_path(path.join(
        'data', 'sdp', 'asns', pool_root + '*.json'
    )))
    for asn_path in asn_paths:
        with open(asn_path) as fp:
            asn = load_asn(fp)
        standards.append(asn)

    results = Main(TEST_ARGS, pool=pool)

    return results.associations, standards

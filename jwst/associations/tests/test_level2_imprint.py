"""Test IMPRINT inclusion"""
from __future__ import absolute_import
import pytest

from .helpers import (
    combine_pools,
    compare_membership,
    t_path,
)

from .. import load_asn
from ..main import Main


@pytest.mark.parametrize(
    'pool_path, asn_standard_path',
    [
        (
            t_path('data/pool_017_spec_nirspec_lv2imprint.csv'),
            t_path('data/pool_017_spec2_001_asn.json')
        ),
    ]
)
def test_imprint(pool_path, asn_standard_path):
    """Test to ensure backgrounds get added"""
    with open(asn_standard_path) as fp:
        asn_standard = load_asn(fp)
    pool = combine_pools([pool_path])
    results = Main(
        [
            '--dry-run',
            '-r',
            t_path('../lib/rules_level2b.py'),
            '--ignore-default'
        ],
        pool=pool
    )
    assert compare_membership(results.associations[0], asn_standard)

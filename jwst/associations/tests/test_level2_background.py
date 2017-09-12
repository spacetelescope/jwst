"""test_level3_dithers: Test of dither rules."""
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
            t_path('data/pool_009_spec_miri_lv2bkg.csv'),
            t_path('data/pool_009_spec2_001_asn.json')
        ),
        (
            t_path('data/pool_010_spec_nirspec_lv2bkg.csv'),
            t_path('data/pool_010_spec2_001_asn.json')
        ),
        (
            t_path('data/pool_011_spec_miri_lv2bkg_lrs.csv'),
            t_path('data/pool_011_spec2_001_asn.json')
        ),
        (
            t_path('data/pool_015_spec_nirspec_lv2bkg_reversed.csv'),
            t_path('data/pool_015_spec2_001_asn.json')
        ),
        (
            t_path('data/pool_016_spec_nirspec_lv2bkg_double.csv'),
            t_path('data/pool_016_spec2_001_asn.json')
        ),
    ]
)
def test_background(pool_path, asn_standard_path):
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


@pytest.mark.xfail(
    reason='Not determined yet',
    run=False
)
def test_level2_bkgnod():
    assert False

"""
class TestLevel2Bkgnod(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/jw95070_20150615T143412_pool.csv'),
            n_asns=2,
            n_orphaned=2,
            kwargs={'delimiter': ','}
        ),
    ]

    valid_rules = [
        'Asn_MIRI_LRS_BKGNOD',
    ]
"""

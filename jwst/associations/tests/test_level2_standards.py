"""Test against Level2 standard associations"""
from __future__ import absolute_import
import pytest

from .helpers import (
    combine_pools,
    compare_asns,
    t_path,
)

from .. import load_asn
from ..main import Main


@pytest.yield_fixture(
    scope='module',
    params=[
        (
            t_path('data/pool_009_spec_miri_lv2bkg.csv'),
            [
                t_path('data/pool_009_spec2_001_asn.json'),
                t_path('data/pool_009_image2_001_asn.json')
            ]
        ),
        (
            t_path('data/pool_010_spec_nirspec_lv2bkg.csv'),
            [
                t_path('data/pool_010_spec2_001_asn.json'),
                t_path('data/pool_010_image2_001_asn.json')
            ]
        ),
        (
            t_path('data/pool_011_spec_miri_lv2bkg_lrs.csv'),
            [
                t_path('data/pool_011_spec2_001_asn.json'),
            ]
        ),
        (
            t_path('data/pool_015_spec_nirspec_lv2bkg_reversed.csv'),
            [
                t_path('data/pool_015_spec2_001_asn.json'),
                t_path('data/pool_015_image2_001_asn.json')
            ]
        ),
        (
            t_path('data/pool_016_spec_nirspec_lv2bkg_double.csv'),
            [
                t_path('data/pool_016_spec2_001_asn.json'),
                t_path('data/pool_016_image2_001_asn.json')
            ]
        ),
        (
            t_path('data/pool_017_spec_nirspec_lv2imprint.csv'),
            [
                t_path('data/pool_017_spec2_001_asn.json'),
                t_path('data/pool_017_image2_001_asn.json')
            ]
        ),
        (
            t_path('data/pool_018_all_exptypes.csv'),
            [
                t_path('data/pool_018_spec2_001_asn.json'),
                t_path('data/pool_018_image2_001_asn.json')
            ]
        ),
    ]
)
def generate_asns(request):
    """Test exp_type inclusion based on standard associations"""
    pool_path, standards_paths = request.param

    standards = {}
    for standard_path in standards_paths:
        with open(standard_path) as fp:
            asn = load_asn(fp)
        standards[asn['asn_type']] = asn

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

    asns = results.associations
    assert len(asns) == len(standards)
    yield asns, standards


def test_against_standard(generate_asns):
    """Compare a generated assocaition against a standard
    """
    generated, standards = generate_asns
    for asn in generated:
        assert compare_asns(asn, standards[asn['asn_type']])

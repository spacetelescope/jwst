"""Test basic usage of Level2 associations"""
from __future__ import absolute_import

import pytest

from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from .. import (
    generate,
    load_asn
)

NONSSCIENCE = ['BACKGROUND']


def from_level2_schema():
    with open(t_path('data/asn_level2.json')) as asn_file:
        asn = load_asn(asn_file)
    return [asn]


def from_level2_image():
    """Test creation of a Level 2 image association"""
    rules = registry_level2_only()
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns, orphaned = generate(pool, rules)
    return asns


def from_level2_spec():
    """Test creation of a Level 2 spectral association"""
    rules = registry_level2_only()
    pool = combine_pools(t_path('data/pool_007_spec_miri.csv'))
    asns, orphaned = generate(pool, rules)
    return asns


@pytest.mark.parametrize(
    'asns, n_asns, n_products, asn_type, asn_rule',
    [
        (from_level2_schema(), 1, 6, 'image2', 'Asn_Lv2Image'),
        (from_level2_image(), 1, 8, 'image2', 'Asn_Lv2Image'),
        (from_level2_spec(), 1, 11, 'spec2', 'Asn_Lv2Spec'),
    ]
)
def test_level2_basics(asns, n_asns, n_products, asn_type, asn_rule):
    assert len(asns) == n_asns
    for asn in asns:
        assert asn['asn_type'] == asn_type
        assert asn['asn_rule'] == asn_rule
        products = asn['products']
        assert len(products) == n_products
        for product in products:
            members = product['members']
            assert len(members) > 0
            science = [
                member
                for member in members
                if member['exptype'] == 'SCIENCE'
            ]
            assert len(science) == 1
            if len(members) > 1:
                nonscience = [
                    member
                    for member in members
                    if member['exptype'] != 'SCIENCE'
                ]
                assert len(nonscience) > 0
                exptypes = set(
                    member['exptype']
                    for member in nonscience
                )
                assert exptypes.issubset(NONSSCIENCE)

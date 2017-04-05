"""Test basic usage of Level2 associations"""
from __future__ import absolute_import

from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from .. import (
    generate,
    load_asn
)


def test_level2_schema():
    with open(t_path('data/asn_level2.json')) as asn_file:
        asn = load_asn(asn_file)
    assert len(asn['members']) == 4
    member = asn['members'][0]
    assert member['expname'] == 'test_lrs1_rate.fits'
    assert member['exptype'] == 'SCIENCE'
    assert isinstance(member['bkgexps'], list)
    bkg = member['bkgexps'][0]
    assert bkg['expname'] == 'test_lrsbkg_rate.fits'


def test_level2_image():
    """Test creation of a Level 2 image association"""
    rules = registry_level2_only()
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns, orphaned = generate(pool, rules)
    assert len(asns) == 8
    len(orphaned) == 0
    asn = asns[0]
    assert asn['asn_rule'] == 'Asn_Lv2Image'
    assert asn['asn_type'] == 'image2'
    assert len(asn['members']) == 1
    member = asn['members'][0]
    base_keys = {'expname', 'exptype'}
    assert base_keys.issubset(member.keys())
    assert member['expname'] == 'jw_00001_rate.fits'
    assert member['exptype'] == 'SCIENCE'


def test_level2_spec():
    """Test creation of a Level 2 spectral association"""
    rules = registry_level2_only()
    pool = combine_pools(t_path('data/pool_007_spec_miri.csv'))
    asns, orphaned = generate(pool, rules)
    assert len(asns) == 5
    len(orphaned) == 0
    asn = asns[0]
    assert asn['asn_rule'] == 'Asn_Lv2Spec'
    assert asn['asn_type'] == 'spec2'
    assert len(asn['members']) == 1
    member = asn['members'][0]
    base_keys = {'expname', 'exptype'}
    assert base_keys.issubset(member.keys())
    assert member['expname'] == 'jw_00001_rate.fits'
    assert member['exptype'] == 'SCIENCE'

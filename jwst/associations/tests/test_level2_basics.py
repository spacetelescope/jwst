"""Test basic usage of Level2 associations"""
from __future__ import absolute_import

from .helpers import (combine_pools, t_path)
from .. import (
    AssociationRegistry,
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
    rules = AssociationRegistry(
        definition_files=[t_path('../lib/rules_level2b.py')],
        include_default=False
    )
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns, orphaned = generate(pool, rules)
    assert len(asns) > 0

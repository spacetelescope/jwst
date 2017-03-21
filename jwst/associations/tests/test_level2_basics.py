"""Test basic usage of Level2 associations"""
from __future__ import absolute_import

from .helpers import t_path
from .. import load_asn


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

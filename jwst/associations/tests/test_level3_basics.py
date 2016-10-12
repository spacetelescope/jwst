"""Test general Level 3 rules environment"""
from __future__ import absolute_import

from .helpers import (combine_pools, t_path)

from .. import (AssociationRegistry, generate)


def test_meta():
    rules = AssociationRegistry()
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns, orphaned = generate(pool, rules)
    assert len(asns) > 0
    asn = asns[0]
    data = asn.data
    assert data['program'] == '99009'
    assert data['target'] == '1'
    assert data['asn_type'] == 'image'
    assert data['asn_id'] == 'a3001'
    assert data['asn_pool'] == 'pool_002_image_miri'
    assert data['asn_rule'] == 'Asn_Image'
    assert data['degraded_status'] == 'No known degraded exposures in association.'
    assert data['version_id'] is None
    assert data['constraints'] is not None

"""test_level3_dithers: Test of WFS rules."""
from __future__ import absolute_import

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)


class TestLevel3WFS(helpers.BasePoolRule):

    pools = [
        helpers.PoolParams(
            path=helpers.t_path('data/jw82600_002_20151107t165901_pool.csv'),
            n_asns=1,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_WFSCMB',
    ]

    def test_wfs_product_name(self):
        rules = AssociationRegistry()
        pool = AssociationPool.read(
            helpers.t_path('data/jw82600_002_20151107t165901_pool.csv')
        )
        (asns, orphaned) = generate(pool, rules)
        name = asns[0].data['products'][0]['name']
        assert name == 'jw82600_t001_nircam_clear_wfscmb-01.fits'

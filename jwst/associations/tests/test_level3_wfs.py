"""test_level3_dithers: Test of WFS rules."""
from __future__ import absolute_import

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)


class TestLevel3WFS(helpers.BasePoolRule):

    pools = [
        helpers.PoolParams(
            path=helpers.t_path('data/pool_004_wfs.csv'),
            n_asns=1,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_WFSCMB',
    ]

    def test_wfs_product_name(self):
        rules = AssociationRegistry()
        pool = helpers.combine_pools(
            helpers.t_path('data/pool_004_wfs.csv')
        )
        (asns, orphaned) = generate(pool, rules)
        name = asns[0]['products'][0]['name']
        assert name == 'jw99009-a3001_t001_nircam_clear_{product_type}-01.fits'

"""test_level3_dithers: Test of WFS rules."""

from . import helpers

from .. import generate


class TestLevel3WFS(helpers.BasePoolRule):

    pools = [
        helpers.PoolParams(
            path=helpers.t_path('data/pool_004_wfs.csv'),
            n_asns=35,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_WFSCMB',
    ]


def test_wfs_product_name():
    rules = helpers.registry_level3_only()
    pool = helpers.combine_pools(
        helpers.t_path('data/pool_004_wfs.csv')
    )
    asns = generate(pool, rules)
    name = asns[0]['products'][0]['name']
    assert name == 'jw99009-c1000_t001_nircam_f150w2'
    assert asns[0]['asn_type'] == 'wfs'

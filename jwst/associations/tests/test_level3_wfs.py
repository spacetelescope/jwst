"""test_level3_dithers: Test of WFS rules."""

from jwst.associations.tests import helpers

from jwst.associations import generate
from jwst.associations.main import constrain_on_candidates

# Generate Level3 associations
all_candidates = constrain_on_candidates(None)
rules = helpers.registry_level3_only(global_constraints=all_candidates)
pool = helpers.combine_pools(
    helpers.t_path('data/pool_004_wfs.csv')
)
level3_asns = generate(pool, rules)


class TestLevel3WFS(helpers.BasePoolRule):

    pools = [
        helpers.PoolParams(
            path=helpers.t_path('data/pool_004_wfs.csv'),
            n_asns=42,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_Lv3WFSCMB',
    ]


def test_wfs_duplicate_product_names():
    """Test for duplicate product names"""
    global level3_asns

    name_list = [
        product['name']
        for asn in level3_asns
        for product in asn['products']
    ]
    assert len(name_list)
    name_set = set(name_list)
    assert len(name_set) == len(name_list)

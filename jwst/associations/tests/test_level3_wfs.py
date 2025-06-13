"""test_level3_dithers: Test of WFS rules."""

from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests import helpers
from jwst.associations import generate
from jwst.associations.lib.utilities import constrain_on_candidates


class TestLevel3WFS(helpers.BasePoolRule):
    pools = [
        helpers.PoolParams(
            path=get_pkg_data_filename("data/pool_004_wfs.csv", package="jwst.associations.tests"),
            n_asns=42,
            n_orphaned=0,
        ),
    ]

    valid_rules = [
        "Asn_Lv3WFSCMB",
    ]


def test_wfs_duplicate_product_names():
    """Test for duplicate product names"""

    # Generate Level3 associations
    rules = helpers.registry_level3_only(global_constraints=constrain_on_candidates(None))
    pool = helpers.combine_pools(
        get_pkg_data_filename("data/pool_004_wfs.csv", package="jwst.associations.tests")
    )
    level3_asns = generate(pool, rules)

    name_list = [product["name"] for asn in level3_asns for product in asn["products"]]
    assert len(name_list)
    name_set = set(name_list)
    assert len(name_set) == len(name_list)

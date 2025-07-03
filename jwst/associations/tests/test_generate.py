"""Test basic generate operations"""

from astropy.utils.data import get_pkg_data_filename

from jwst.associations import AssociationPool, AssociationRegistry, generate, load_asn


def test_simple():
    """Test generate on simple registry"""
    registry = AssociationRegistry(
        [get_pkg_data_filename("data/rules_basic.py", package="jwst.associations.tests")],
        include_default=False,
    )
    pool = AssociationPool()
    pool["value"] = ["row1", "row2"]

    asns = generate(pool, registry)
    assert len(asns) == 1
    assert len(asns[0]["members"]) == 2


def test_unserialize():
    """Test basic unserializing"""
    asn_file = get_pkg_data_filename("data/asn_mosaic.json", package="jwst.associations.tests")
    with open(asn_file, "r") as asn_fp:
        asn = load_asn(asn_fp)
    assert isinstance(asn, dict)

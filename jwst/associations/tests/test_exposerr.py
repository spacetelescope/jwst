"""Test degraded exposure info"""

from astropy.utils.data import get_pkg_data_filename

from jwst.associations.lib.dms_base import _DEGRADED_STATUS_NOTOK, _DEGRADED_STATUS_OK, _EMPTY
from jwst.associations.main import Main
from jwst.associations.tests.helpers import combine_pools


def test_exposerr():
    pool = combine_pools(
        get_pkg_data_filename("data/pool_008_exposerr.csv", package="jwst.associations.tests")
    )
    generated = Main.cli(
        [
            "--dry-run",
            "-i",
            "o001",
        ],
        pool=pool,
    )
    asns = generated.associations
    assert len(asns) == 1
    for asn in asns:
        any_degraded = False
        for product in asn["products"]:
            any_degraded = any_degraded or any(
                [member["exposerr"] not in _EMPTY for member in product["members"]]
            )
        if any_degraded:
            assert asn["degraded_status"] == _DEGRADED_STATUS_NOTOK
        else:
            assert asn["degraded_status"] == _DEGRADED_STATUS_OK

"""Test Image rules

Image rules produce associations meant to be
processed by CALIMAGE3.

Such associations are mosaics and dithers.
"""

from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests.helpers import BasePoolRule, PoolParams

_miri_data = get_pkg_data_filename(
    "data/pool_002_image_miri.csv", package="jwst.associations.tests"
)
_nircam_data = get_pkg_data_filename(
    "data/pool_003_image_nircam.csv", package="jwst.associations.tests"
)


class TestLevel3Image(BasePoolRule):
    pools = [
        PoolParams(path=_miri_data, n_asns=1, n_orphaned=0),
        PoolParams(path=_nircam_data, n_asns=2, n_orphaned=0),
        # Below tested cannot be run due to an obscure numpy.ma bug.
        # PoolParams(
        #    path=[
        #        _miri_data,
        #        _nircam_data,
        #    ],
        #    n_asns=3,
        #    n_orphaned=0
        # ),
    ]

    valid_rules = [
        "Asn_Lv3Image",
    ]

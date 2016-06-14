"""test_level3_dithers: Test of spectrographic rules."""
from __future__ import absolute_import

from .helpers import BasePoolRule, PoolParams, t_path


class TestLevel3Spectrographic(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/jw00034_20150804T190826_pool.csv'),
            n_asns=1,
            n_orphaned=0
        ),
        PoolParams(
            path=t_path('data/jw84500_20150722T215143_pool.csv'),
            n_asns=3,
            n_orphaned=0
        ),
        PoolParams(
            path=t_path('data/jw80500_20150722T215139_pool.csv'),
            n_asns=3,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_MIRI_LRS_FIXEDSLIT',
        'Asn_NIR_SO_SLITLESS',
        'Asn_NIR_FIXEDSLIT',
        'Asn_NIR_MSA',
        'Asn_MIRI_MRS',
        'Asn_NIR_IFU'
    ]

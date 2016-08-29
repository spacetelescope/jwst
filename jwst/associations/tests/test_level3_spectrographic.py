"""test_level3_dithers: Test of spectrographic rules."""
from __future__ import absolute_import
from glob import glob

from .helpers import BasePoolRule, PoolParams, t_path


class TestLevel3Spectrographic(BasePoolRule):

    pools = [
        PoolParams(
            path=glob(t_path('data/pool_*_spec_*.csv')),
            n_asns=7,
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

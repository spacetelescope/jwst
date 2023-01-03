"""Test degraded exposure info"""
from jwst.associations.tests.helpers import (
    combine_pools,
    t_path
)

from jwst.associations.lib.dms_base import (
    _DEGRADED_STATUS_OK,
    _DEGRADED_STATUS_NOTOK,
    _EMPTY
    )
from jwst.associations.main import Main


def test_exposerr():
    pool = combine_pools(t_path('data/pool_008_exposerr.csv'))
    generated = Main.cli(
        [
            '--dry-run', '-i', 'o001',
        ],
        pool=pool
    )
    asns = generated.associations
    assert len(asns) > 1
    for asn in asns:
        any_degraded = False
        for product in asn['products']:
            any_degraded = any_degraded or any([
                member['exposerr'] not in _EMPTY
                for member in product['members']
            ])
        if any_degraded:
            assert asn['degraded_status'] == _DEGRADED_STATUS_NOTOK
        else:
            assert asn['degraded_status'] == _DEGRADED_STATUS_OK

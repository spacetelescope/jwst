"""Test degraded exposure info"""
from __future__ import absolute_import

from .helpers import (
    combine_pools,
    level3_rule_path,
    t_path
)

from ..lib.rules_level3_base import (
    _DEGRADED_STATUS_OK,
    _DEGRADED_STATUS_NOTOK,
    _EMPTY
    )
from ..main import Main


def test_exposerr():
    pool = combine_pools(t_path('data/pool_008_exposerr.csv'))
    generated = Main(
        [
            '--dry-run', '-i', 'o001',
            '-r', level3_rule_path(),
            '--ignore-default',
        ],
        pool=pool
    )
    asns = generated.associations
    assert len(asns) == 1
    asn = asns[0]
    any_degraded = any([
        member['exposerr'] not in _EMPTY
        for member in asn['products'][0]['members']
    ])
    assert any_degraded
    assert asn['degraded_status'] == _DEGRADED_STATUS_NOTOK

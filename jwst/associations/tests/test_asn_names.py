from __future__ import absolute_import

import pytest
import re

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)

LEVEL3_ASN_ACID_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>(o|c)\d{3,4})'
    '_(?P<stamp>.+)'
    '_(?P<asn_type>.+)'
    '_(?P<sequence>\d{3})'
    '_asn'
)
LEVEL3_ASN_DISCOVERED_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>a\d{4})'
    '_(?P<stamp>.+)'
    '_(?P<asn_type>.+)'
    '_(?P<sequence>\d{3})'
    '_asn'
)


pool_params = pytest.fixture(
    scope='module',
    params=[
        'data/pool_002_image_miri.csv'
    ]
)(helpers.generate_params)


class TestASNtNames():

    def test_level3_asn_names(self, pool_params):
        pool_path = helpers.t_path(pool_params)
        pool = helpers.combine_pools(pool_path)
        rules = AssociationRegistry()
        asns, orphaned = generate(pool, rules)
        assert len(asns) > 0
        for asn in asns:
            name = asn.asn_name
            if any(
                    c.get('is_acid', False)
                    for _, c in asn.constraints.items()
            ):
                m = re.match(LEVEL3_ASN_ACID_NAME_REGEX, name)
            else:
                m = re.match(LEVEL3_ASN_DISCOVERED_NAME_REGEX, name)
            assert m is not None

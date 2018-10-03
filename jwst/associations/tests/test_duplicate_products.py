"""Check for duplicate product names"""


from collections import(
    Counter,
    defaultdict,
)
from glob import glob
from os import path as op
import pytest

from .helpers import (
    t_path,
)

from ..main import Main as asn_generate


# Main test args
TEST_ARGS = ['--dry-run', '--no-merge']

# Pools to examine
SDP_POOL_PATHS = glob(t_path(op.join(
    'data', 'sdp', 'pools', '*.csv'
)))


@pytest.mark.slow
@pytest.mark.parametrize(
    'pool_path',
    SDP_POOL_PATHS
)
def test_dup_product_names(pool_path):
    """Check for duplicate product names for a pool"""

    results = asn_generate([pool_path] + TEST_ARGS)
    asns = results.associations

    product_names = Counter(
        product['name']
        for asn in asns
        for product in asn['products']
    )

    multiples = [
        product_name
        for product_name, count in product_names.items()
        if count > 1
    ]

    assert not len(multiples), 'Multiple product names: {}'.format(multiples)

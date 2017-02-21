"""Test file_to_asn"""

import pytest

from .. import (Association, load_asn)
from ..asn_from_list import asn_from_list


@pytest.mark.xfail
def test_base_association():
    items = ['a', 'b', 'c']
    asn = asn_from_list(items, rule=Association)
    assert asn['asn_rule'] == 'Assocation'
    assert asn['asn_type'] == 'None'
    assert asn['members'] == items


def test_default_simple():
    product_name = 'test_product'
    items = ['a', 'b', 'c']
    asn = asn_from_list(items, product_name=product_name)
    assert asn['asn_rule'] == 'DMS_Level3_Base'
    assert asn['asn_type'] == 'None'
    assert len(asn['products']) == 1
    product = asn['products'][0]
    assert product['name'] == product_name
    assert len(product['members']) == len(items)
    for member in product['members']:
        assert member['expname'] in items
        assert member['exptype'] == 'SCIENCE'


def test_default_with_type():
    product_name = 'test_product'
    items = {
        'a': 'SCIENCE',
        'b': 'TARGET_ACQ',
        'c': 'SOMETHINGELSE'
    }
    asn = asn_from_list(
        [(item, type_) for item, type_ in items.items()],
        product_name=product_name,
        with_type=True
    )
    assert asn['asn_rule'] == 'DMS_Level3_Base'
    assert asn['asn_type'] == 'None'
    assert len(asn['products']) == 1
    product = asn['products'][0]
    assert product['name'] == product_name
    assert len(product['members']) == len(items)
    for member in product['members']:
        assert member['expname'] in items
        assert member['exptype'] == items[member['expname']]


def test_default_fail():
    """Test default DMS_Level3_Base fail

    A product name needs to be included, but is not.
    """
    items = ['a']
    with pytest.raises((KeyError, TypeError)):
        asn = asn_from_list(items)


def test_default_roundtrip():
    product_name = 'test_product'
    items = {
        'a': 'SCIENCE',
        'b': 'TARGET_ACQ',
        'c': 'SOMETHINGELSE'
    }
    asn = asn_from_list(
        [(item, type_) for item, type_ in items.items()],
        product_name=product_name,
        with_type=True
    )
    name, serialized = asn.dump()
    reloaded = load_asn(serialized)
    assert asn['asn_rule'] == reloaded['asn_rule']
    assert asn['asn_type'] == reloaded['asn_type']
    assert len(asn['products']) == len(reloaded['products'])

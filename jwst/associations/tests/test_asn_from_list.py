"""Test file_to_asn"""

import pytest

from .. import (Association, AssociationRegistry, load_asn)
from ..asn_from_list import (Main, asn_from_list)


def test_base_association():
    items = ['a', 'b', 'c']
    asn = asn_from_list(items, rule=Association)
    assert asn['asn_rule'] == 'Association'
    assert asn['asn_type'] == 'None'
    assert asn['members'] == items


def test_base_roundtrip():
    items = ['a', 'b', 'c']
    asn = asn_from_list(items, rule=Association)
    name, serialized = asn.dump()
    reloaded = load_asn(serialized, registry=None)
    assert asn['asn_rule'] == reloaded['asn_rule']
    assert asn['asn_type'] == reloaded['asn_type']
    assert asn['members'] == reloaded['members']


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


def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        Main([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        Main(['-o', 'test_asn.json'])


@pytest.mark.parametrize(
    "format",
    ['json', 'yaml']
)
def test_cmdline_success(format, tmpdir):
    path = tmpdir.join('test_asn.json')
    product_name = 'test_product'
    inlist = ['a', 'b', 'c']
    args = [
        '-o', path.strpath,
        '--product-name', product_name,
        '--format', format
    ]
    args = args + inlist
    Main(args)
    with open(path.strpath, 'r') as fp:
        asn = load_asn(fp, format=format)
    assert len(asn['products']) == 1
    assert asn['products'][0]['name'] == product_name
    members = asn['products'][0]['members']
    expnames = [
        member['expname']
        for member in members
    ]
    assert inlist == expnames


def test_cmdline_change_rules(tmpdir):
    rule = 'Association'
    path = tmpdir.join('test_asn.json')
    inlist = ['a', 'b', 'c']
    args = [
        '-o', path.strpath,
        '-r', rule,
    ]
    args = args + inlist
    Main(args)
    with open(path.strpath, 'r') as fp:
        asn = load_asn(fp, registry=AssociationRegistry(include_bases=True))
    assert inlist == asn['members']


def test_api_list():
    """Test api call with simple list"""
    product_name = 'test_product'
    inlist = ['a', 'b', 'c']

    asn = asn_from_list(inlist, product_name=product_name)
    assert len(asn['products']) == 1
    assert asn['products'][0]['name'] == product_name
    members = asn['products'][0]['members']
    expnames = [
        member['expname']
        for member in members
    ]
    assert inlist == expnames


def test_api_with_type():
    """Test api call with type tuple"""
    product_name = 'test_product'
    inlist = [
        ('a', 'science'),
        ('b', 'psf'),
        ('c', 'science')
    ]

    asn = asn_from_list(inlist, product_name=product_name, with_type=True)
    assert len(asn['products']) == 1
    assert asn['products'][0]['name'] == product_name
    members = asn['products'][0]['members']
    members_dict = {
        member['expname']: member['exptype']
        for member in members
    }
    for name, type_ in inlist:
        assert members_dict[name] == type_

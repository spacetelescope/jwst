"""Test utility update_path"""
from ..asn_from_list import asn_from_list
from ..lib.rules_level2_base import DMSLevel2bBase
from ..lib.update_path import update_path


def test_update_path_level2():
    members = ['a', 'b', 'c']
    new_path = 'new_path'
    asn = asn_from_list(members, rule=DMSLevel2bBase)
    update_path(asn, new_path)
    for product in asn['products']:
        for member in product['members']:
            assert member['expname'].startswith(new_path)


def test_update_path_level3():
    members = ['a', 'b', 'c']
    new_path = 'new_path'
    asn = asn_from_list(members, product_name='test')
    update_path(asn, new_path)
    for product in asn['products']:
        for member in product['members']:
            assert member['expname'].startswith(new_path)

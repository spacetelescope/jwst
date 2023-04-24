"""Test association diffing"""
import json
import pytest

import jwst.associations.lib.diff as asn_diff

# Test cases for individual associations
badcandidate_asn = {
    'asn_type': 'test',
    'asn_id': 'o001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_a', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_b',
            'members': [
                {'expname': 'member_b_a', 'exptype': 'science'},
                {'expname': 'member_b_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_d',
            'members': [
                {'expname': 'member_d_a', 'exptype': 'science'},
                {'expname': 'member_d_b', 'exptype': 'science'}
            ]
        },
    ]
}

badexptype_asn = {
    'asn_type': 'test',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_b', 'exptype': 'background'}
            ]
        },
        {
            'name': 'product_b',
            'members': [
                {'expname': 'member_b_a', 'exptype': 'science'},
                {'expname': 'member_b_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_c',
            'members': [
                {'expname': 'member_c_a', 'exptype': 'science'},
                {'expname': 'member_c_b', 'exptype': 'science'}
            ]
        },
    ]
}

badmember_asn = {
    'asn_type': 'test',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_a', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_b',
            'members': [
                {'expname': 'member_b_a', 'exptype': 'science'},
                {'expname': 'member_b_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_d',
            'members': [
                {'expname': 'member_d_a', 'exptype': 'science'},
                {'expname': 'member_d_b', 'exptype': 'science'}
            ]
        },
    ]
}

badtype_asn = {
    'asn_type': 'badtype',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_b',
            'members': [
                {'expname': 'member_b_a', 'exptype': 'science'},
                {'expname': 'member_b_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_d',
            'members': [
                {'expname': 'member_d_a', 'exptype': 'science'},
                {'expname': 'member_d_b', 'exptype': 'science'}
            ]
        },
    ]
}

disjoint_product_asn = {
    'asn_type': 'test',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_b',
            'members': [
                {'expname': 'member_b_a', 'exptype': 'science'},
                {'expname': 'member_b_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_d',
            'members': [
                {'expname': 'member_d_a', 'exptype': 'science'},
                {'expname': 'member_d_b', 'exptype': 'science'}
            ]
        },
    ]
}

dup_product_asn = {
    'asn_type': 'test',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_c',
            'members': [
                {'expname': 'member_c_a', 'exptype': 'science'},
                {'expname': 'member_c_b', 'exptype': 'science'}
            ]
        },
    ]
}

standard_asn = {
    'asn_type': 'test',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_a', 'exptype': 'science'},
                {'expname': 'member_a_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_b',
            'members': [
                {'expname': 'member_b_a', 'exptype': 'science'},
                {'expname': 'member_b_b', 'exptype': 'science'}
            ]
        },
        {
            'name': 'product_c',
            'members': [
                {'expname': 'member_c_a', 'exptype': 'science'},
                {'expname': 'member_c_b', 'exptype': 'science'}
            ]
        },
    ]
}

# Test cases for association comparison.
# The products in each association should be considered separate
# associations to be compared.
asns_subset = {
    'asn_type': 'test',
    'asn_id': 'c1001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_1', 'exptype': 'science'},
                {'expname': 'member_a_2', 'exptype': 'science'},
            ]
        },
        {
            'name': 'product_a',
            'members': [
                {'expname': 'member_a_1', 'exptype': 'science'},
                {'expname': 'member_a_2', 'exptype': 'science'},
                {'expname': 'member_a_3', 'exptype': 'science'},
            ]
        }
    ]
}


def test_duplicate_members():
    """Test duplicate members"""
    product = badcandidate_asn['products'][0]
    with pytest.raises(asn_diff.DuplicateMembersError):
        asn_diff.check_duplicate_members(product)


def test_equivalency():
    """Test that two associations equivalency

    Success is the fact that no errors happen.
    """
    asns = asn_diff.separate_products(standard_asn)
    asn_diff.compare_asn_lists(asns, asns)


@pytest.mark.parametrize(
    'mismatched',
    [
        dup_product_asn,
        disjoint_product_asn,
        badtype_asn,
        badmember_asn,
        badexptype_asn,
        badcandidate_asn,
    ]
)
def test_fails(mismatched, standard=standard_asn):
    """Test cases of failures

    Parameters
    ----------
    mismatched: Association
        The association that should not match to the standard.

    standard: Association
        The standard association.
    """
    with pytest.raises(AssertionError):
        left_asns = asn_diff.separate_products(mismatched)
        right_asns = asn_diff.separate_products(standard)
        asn_diff.compare_asn_lists(left_asns, right_asns)


@pytest.mark.usefixtures('_jail')
def test_fromfiles():
    """Test from files

    Success is the fact that no errors happen.
    """
    with open('test.json', 'w') as fh:
        json.dump(standard_asn, fh)

    asn_diff.compare_asn_files(['test.json'], ['test.json'])


def test_separate_products():
    new_asns = asn_diff.separate_products(standard_asn)

    assert len(new_asns) == 3
    for asn in new_asns:
        assert len(asn['products']) == 1
        idx = asn['products'][0]['name'][-1]
        for member in asn['products'][0]['members']:
            assert member['expname'][:-2] == 'member' + '_' + idx


def test_subset_product():
    """Test subset identification"""
    left, right = asn_diff.separate_products(asns_subset)
    try:
        asn_diff.compare_product_membership(left['products'][0], right['products'][0])
    except asn_diff.MultiDiffError as exception:
        if len(exception) > 1 or not isinstance(exception[0], asn_diff.SubsetError):
            raise

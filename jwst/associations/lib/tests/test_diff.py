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


# Test case when comparing groups of associations.
# There should be no difference because of the `rate/rateints`
# equivalence. Tests will try different ordering due to a bug
# which made the comparison order-dependent.
NODIFF_ASNS = [
    {
        "asn_type": "image2",
        "asn_id": "o006",
        "products": [
            {
                "name": "jw00016006001_04105_00001_nrca4",
                "members": [
                    {
                        "expname": "jw00016006001_04105_00001_nrca4_rate.fits",
                        "exptype": "science",
                        "exposerr": "null"
                    }
                ]
            }
        ]
    },
    {
        "asn_type": "image2",
        "asn_id": "o006",
        "products": [
            {
                "name": "jw00016006001_04105_00001_nrca4",
                "members": [
                    {
                        "expname": "jw00016006001_04105_00001_nrca4_rateints.fits",
                        "exptype": "science",
                        "exposerr": "null"
                    }
                ]
            }
        ]
    },
]


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


@pytest.mark.parametrize(
    'left_indexes, right_indexes',
    [
        ([0, 1], [0, 1]),
        ([0, 1], [1, 0]),
        ([1, 0], [0, 1]),
        ([1, 0], [1, 0]),
    ]
)
def test_nodiff_order(left_indexes, right_indexes):
    """Associations should have no difference, regardless of order"""

    # Order the associations as specified in the indexes.
    left_asns = [NODIFF_ASNS[idx] for idx in left_indexes]
    right_asns = [NODIFF_ASNS[idx] for idx in right_indexes]

    # Execute test. Success will not generate any exceptions.
    asn_diff.compare_asn_lists(left_asns, right_asns)


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

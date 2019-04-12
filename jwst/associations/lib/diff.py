""""Diff and compare associations"""

from collections import Counter
from copy import copy
import logging
import re

from ..load_asn import load_asn

logger = logging.getLogger(__name__)


def compare_asn_files(left_paths, right_paths):
    """Compare association files

    Parameters
    ----------
    left_paths: [Path or str[, ...]]
        Set of association files

    right_paths: [Path or str[, ...]]
        Set of association files to compare against

    Raises
    ------
    AssertionError
        If there are differences. The message will contain
        all the differences.
    """
    # Read in all the associations, separating out the products into separate associations.
    left_asns = []
    for path in left_paths:
        with open(path, 'r') as fh:
            asn = load_asn(fh)
        single_product_asns = separate_products(asn)

        # Collect the associations.
        left_asns += single_product_asns

    right_asns = []
    for path in right_paths:
        with open(path, 'r') as fh:
            asn = load_asn(fh)
        single_product_asns = separate_products(asn)

        # Collect the associations.
        right_asns += single_product_asns

    # Compare the associations.
    compare_asn_lists(left_asns, right_asns)


def compare_asn_lists(left_asns, right_asns):
    """Compare to lists of associations

    Both association lists must contain associations that only have
    single products. Use `separate_products` prior to calling this
    function.

    Parameters
    ----------
    left_asns: [`Association`[, ...]]
        Group of associations

    right_asns: [`Association`[, ...]]
        Group of associations to compare

    Raises
    ------
    AssertionError
        If there are differences. The message will contain
        all the differences.
    """

    # Ensure that product names are unique
    left_product_names = no_duplicate_products(left_asns)
    right_product_names = no_duplicate_products(right_asns)

    # Ensure that the product name lists are the same.
    name_diff = left_product_names ^ right_product_names
    assert not len(name_diff), (
        'Associations do not share a common set of products: {}'
        ''.format(name_diff)
    )

    # Compare like product associations
    left_asns_by_product = {
        asn['products'][0]['name']: asn
        for asn in left_asns
    }
    right_asns_by_product = {
        asn['products'][0]['name']: asn
        for asn in right_asns
    }
    for product_name in left_product_names:
        compare_asns(left_asns_by_product[product_name], right_asns_by_product[product_name])


def compare_asns(left, right):
    """Compare two associations

    This comparison will include metadata such as
    `asn_type` and membership

    Parameters
    ---------
    left, right : dict
        Two, individual, associations to compare

    Returns
    -------
    equality : boolean
        True if matching.

    Raises
    ------
    AssertionError
        If there is a difference.
    """

    try:
        _compare_asns(left, right)
    except AssertionError as excp:
        message = (
            'Associations do not match. Mismatch because:'
            '\n{reason}'
            '\nLeft association = {left}'
            '\nRight association = {right}'
            ''.format(left=left, right=right, reason=excp)
        )
        raise AssertionError(message) from excp

    return True


def _compare_asns(left, right):
    """Compare two associations

    This comparison will include metadata such as
    `asn_type` and membership

    Parameters
    ---------
    left, right : dict
        Two, individual, associations to compare

    Returns
    -------
    equality : boolean

    Raises
    ------
    AssertionError
    """

    # Metadata
    assert left['asn_type'] == right['asn_type'], \
        'Type mismatch {} != {}'.format(left['asn_type'], right['asn_type'])

    # Membership
    return compare_membership(left, right)


def compare_membership(left, right):
    """Compare two associations' membership

    Parameters
    ---------
    left, right : dict
        Two, individual, associations to compare

    Returns
    -------
    equality : boolean

    Raises
    ------
    AssertionError
    """
    products_left = left['products']
    products_right = copy(right['products'])

    assert len(products_left) == len(products_right), (
        '# products differ: {left_len} != {right_len}'
        ''.format(left_len=len(products_left), right_len=len(products_right))
    )

    for left_idx, left_product in enumerate(products_left):
        left_product_name = components(left_product['name'])
        for right_idx, right_product in enumerate(products_right):
            if components(right_product['name']) != left_product_name:
                continue

            assert len(right_product['members']) == len(left_product['members']), (
                'Product Member length differs:'
                ' Left Product #{left_idx} len {left_len} !=  '
                ' Right Product #{right_idx} len {right_len}'
                ''.format(left_idx=left_idx, left_len=len(left_product['members']),
                          right_idx=right_idx, right_len=len(right_product['members'])
                )
            )

            members_right = copy(right_product['members'])
            for left_member in left_product['members']:
                for right_member in members_right:
                    if left_member['expname'] != right_member['expname']:
                        continue

                    assert left_member['exptype'] == right_member['exptype'], (
                        'Left {left_expname}:{left_exptype}'
                        ' != Right {right_expname}:{right_exptype}'
                        ''.format(left_expname=left_member['expname'], left_exptype=left_member['exptype'],
                                  right_expname=right_member['expname'], right_exptype=right_member['exptype'])
                    )

                    members_right.remove(right_member)
                    break
                else:
                    raise AssertionError(
                        'Left {left_expname}:{left_exptype} has no counterpart in right'
                        ''.format(left_expname=left_member['expname'], left_exptype=left_member['exptype'])
                    )

            assert len(members_right) == 0, (
                'Right has {len_over} unaccounted for members starting with'
                ' {right_expname}:{right_exptype}'
                ''.format(len_over=len(members_right),
                          right_expname=members_right[0]['expname'],
                          rigth_exptype=members_right[0]['exptype']
                )
            )

            products_right.remove(right_product)
            break
        else:
            raise AssertionError(
                'Left has {n_products} left over'
                ''.format(n_products=len(products_left))
            )


    assert len(products_right) == 0, (
        'Right has {n_products} left over'
        ''.format(n_products=len(products_right))
    )

    return True


def components(s):
    """split string into its components"""
    return set(re.split('[_-]', s))


def separate_products(asn):
    """Seperate products into individual associations

    Parameters
    ----------
    asn: `Association`
        The association to split

    Returns
    -------
    separated: [`Association`[, ...]]
        The list of separated associations
    """
    separated = []
    for product in asn['products']:
        new_asn = copy(asn)
        new_asn['products'] = [product]
        separated.append(new_asn)
    return separated


def no_duplicate_products(asns):
    """Ensure there are no duplicate product names in a list of associations

    Parameters
    ----------
    asns: [`Association`[, ...]]

    Returns
    -------
    product_names: set(str[, ...])
        The set of product names.

    Raises
    ------
    AssertionError
        If there are differences. The message will contain
        all the differences.
    """
    product_names = [
        asn['products'][0]['name']
        for asn in asns
    ]

    dups = [
        name
        for name, count in Counter(product_names).items()
        if count > 1
    ]
    assert not len(dups), (
        'Associations have the following product name duplication'
        '{}'.format(dups)
    )

    return set(product_names)

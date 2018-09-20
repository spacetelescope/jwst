""""Diff and compare associations"""

from copy import copy
import logging
import re

logger = logging.getLogger(__name__)


def compare_asns(left, right):
    """Compare two associations

    This comparison will include metadata such as
    `asn_type` and membership

    Parameters
    ---------
    left, right: dict
        Two, individual, associations to compare

    Returns
    -------
    equality: boolean
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
    left, right: dict
        Two, individual, associations to compare

    Returns
    -------
    equality: boolean

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
    left, right: dict
        Two, individual, associations to compare

    Returns
    -------
    equality: boolean

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
                                  rigth_expname=right_member['expname'], right_exptype=right_member['exptype'])
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



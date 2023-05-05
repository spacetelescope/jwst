""""Diff and compare associations"""

from collections import Counter, UserList, defaultdict
from copy import copy
import logging
from pathlib import Path
import re

from .product_utils import sort_by_candidate
from ..load_asn import load_asn
from ...lib.suffix import remove_suffix

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = [
    'compare_asn_files',
    'compare_asn_lists',
    'compare_asns',
    'compare_membership',
    'compare_product_membership',
]

# #########################
# Define the types of diffs
# #########################
class DiffError(AssertionError):
    """Base Class for difference errors

    Attributes
    ----------
    args : tuple
        Standard exception arguments

    asns : [Association[,...]]
        List of associations that generated the exception
    """
    def __init__(self, *args, asns=None, **kwargs):
        super().__init__(*args, **kwargs)
        if asns is None:
            asns = list()
        self.asns = asns


class CandidateLevelError(DiffError):
    """Candidate level mismatch"""


class DuplicateMembersError(DiffError):
    """Duplicate members within a product"""


class DuplicateProductError(DiffError):
    """Duplicate products found"""


class DifferentProductSetsError(DiffError):
    """Different product sets between groups of associations"""


class MemberLengthDifference(DiffError):
    """Difference in number of members between associations"""


class MemberMismatchError(DiffError):
    """Membership does not match"""


class SubsetError(DiffError):
    """One product is a subset of another"""


class TypeMismatchError(DiffError):
    """Association type mismatch"""


class UnaccountedMembersError(DiffError):
    """Members not present in other"""


class MultiDiffError(UserList, DiffError):
    """List of diff errors"""
    def __init__(self, *args, **kwargs):
        super(MultiDiffError, self).__init__(*args, **kwargs)

    @property
    def asns(self):
        """Return the list of affected associations for all contained errors"""
        asns = list()
        for err in self:
            asns.extend(err.asns)
        return asns

    @property
    def err_types(self):
        """Return list of error types"""
        err_types = [type(err) for err in self]
        return err_types

    def __str__(self):
        message = ['Following diffs found:\n']
        for diff in self:
            message.extend([
                '\n****\n', str(diff), '\n'
            ])
        return ''.join(message)


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
    MultiDiffError
        If there are differences. The message will contain
        all the differences.
    """
    __tracebackhide__ = True
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
    MultiDiffError
        If there are differences. The message will contain
        all the differences.
    """
    __tracebackhide__ = True
    diffs = MultiDiffError()

    # Ensure that product names are unique
    left_product_names, left_duplicates = get_product_names(left_asns)
    right_product_names, right_duplicates = get_product_names(right_asns)
    if left_duplicates:
        try:
            check_duplicate_products(left_asns)
        except MultiDiffError as dup_errors:
            diffs.extend(dup_errors)
    if right_duplicates:
        try:
            check_duplicate_products(right_asns)
        except MultiDiffError as dup_errors:
            diffs.extend(dup_errors)

    # Ensure that the product name lists are the same.
    name_diff = left_product_names ^ right_product_names
    if name_diff:
        diffs.append(DifferentProductSetsError(
            'Left and right associations do not share a common set of products: {}'
            ''.format(name_diff)
        ))

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
        try:
            compare_asns(left_asns_by_product[product_name], right_asns_by_product[product_name])
        except MultiDiffError as compare_diffs:
            diffs.extend(compare_diffs)
        except KeyError:
            # Most likely due to a previous error. Ignore
            pass

    if diffs:
        raise diffs


def compare_asns(left, right):
    """Compare two associations

    This comparison will include metadata such as
    `asn_type` and membership

    Parameters
    ---------
    left, right : dict
        Two, individual, associations to compare

    Raises
    ------
    MultiDiffError
        If there is a difference.
    """
    _compare_asns(left, right)


def _compare_asns(left, right):
    """Compare two associations

    This comparison will include metadata such as
    `asn_type` and membership

    Parameters
    ---------
    left, right : dict
        Two, individual, associations to compare

    Raises
    ------
    MultiDiffError
        If there are differences. The message will contain
        all the differences.

    Note
    ----
    This comparison is dependent on the associations being JWST-like associations.
    The attributes that are compared are as follows:

        - key `asn_type`
        - key `products`. Specifically the following are compared:
            - Length of the list
            - key `name` for each product
            - key `members` for each product
        - For the member lists of each product, the following are compared:
            - Length of the list
            - key `expname` for each member
            - key 'exptype` for each member
    """
    diffs = MultiDiffError()

    # Assert that the same result type is the same.
    if left['asn_type'] != right['asn_type']:
        diffs.append(TypeMismatchError(
            'Type mismatch {} != {}'.format(left['asn_type'], right['asn_type'])
        ))

    # Assert that the level of association candidate is the same.
    # Cannot guarantee value, but that the 'a'/'c'/'o' levels are similar.
    if left['asn_id'][0] != right['asn_id'][0]:
        diffs.append(CandidateLevelError(
            f"Candidate level mismatch left '{left['asn_id'][0]}' != right '{right['asn_id'][0]}'"
        ))

    # Membership
    try:
        compare_membership(left, right)
    except MultiDiffError as compare_diffs:
        diffs.extend(compare_diffs)

    if diffs:
        raise diffs


def compare_membership(left, right):
    """Compare two associations' membership

    Parameters
    ---------
    left, right : dict
        Two, individual, associations to compare

    Raises
    ------
    MultiDiffError
        If there are differences. The message will contain
        all the differences.
    """
    diffs = MultiDiffError()
    products_left = left['products']
    products_right = copy(right['products'])

    if len(products_left) != len(products_right):
        diffs.append(DifferentProductSetsError(
            '# products differ: {left_len} != {right_len}'
            ''.format(left_len=len(products_left), right_len=len(products_right))
        ))

    for left_idx, left_product in enumerate(products_left):
        left_product_name = components(left_product['name'])
        for right_idx, right_product in enumerate(products_right):
            if components(right_product['name']) != left_product_name:
                continue
            try:
                compare_product_membership(left_product, right_product)
            except MultiDiffError as compare_diffs:
                diffs.extend(compare_diffs)
            products_right.remove(right_product)
            break
        else:
            diffs.append(DifferentProductSetsError(
                'Left has {n_products} left over'
                ''.format(n_products=len(products_left))
            ))


    if len(products_right) != 0:
        diffs.append(DifferentProductSetsError(
            'Right has {n_products} left over'
            ''.format(n_products=len(products_right))
        ))

    if diffs:
        raise diffs


def compare_nosuffix(left, right):
    """Check if the only difference is in rate vs rateints suffixes

    A valid situation is to have two associations be exactly the same except
    for the suffix used on all the science inputs. If one association uses
    "rate" and the other uses "rateints", this is OK.

    If no errors are raised, the two associations are similar except
    for the suffixes used in the `expname`s, which is valid.

    Otherwise, the following errors can be raised:

    - DifferentProductSetsError
      Products do not match.

    Parameters
    ----------
    left, right : Association, Association
        The associations to compare.

    Raises
    ------
    MultiDiffError
    """
    if len(left['products']) != len(right['products']):
        raise MultiDiffError([DifferentProductSetsError(
            f'Different products between {left.asn_name}({len(left["products"])}) and {right.asn_name}({len(right["products"])})',
            asns=[left, right]
        )])

    diffs = MultiDiffError()
    for left_product in left['products']:
        left_members = set(exposure_name(member['expname'])[0]
                         for member in left_product['members'])
        for right_product in right['products']:
            right_members = set(exposure_name(member['expname'])[0]
                              for member in right_product['members'])
            if left_members == right_members:
                break
        else:
            # No right product matches the left product.
            # This is a fail.
            diffs.append(DifferentProductSetsError(
                f'Left product {left_product["name"]} has not counterpart in right products',
                asns=[left, right]
            ))

    if diffs:
        raise diffs


def compare_product_membership(left, right):
    """Compare membership between products

    If the membership is exactly the same, or other special conditions,
    no error is raised.

    Otherwise, the following errors will be raised:

    - DuplicateMemberError
      A member exists multiple times in the member list.

    - MemberLengthDifference
      Number of members in the two products differ.

    - MemberMismatchError
      Members with the same `expname` do not share the other attributes.

    - SubsetError
      A member list is a subset of another member list.

    - UnaccountedMembersError
      Members exist in one association that do not exist in the other.

    Parameters
    ---------
    left, right : dict
        Two, individual, association products to compare. Note that
        these are not associations, just products from associations.

    Raises
    ------
    MultiDiffError
        If there are differences. The message will contain
        all the differences.
    """
    diffs = MultiDiffError()

    # Check for duplicate members.
    try:
        check_duplicate_members(left)
    except DuplicateMembersError as dup_member_error:
        diffs.append(dup_member_error)
    try:
        check_duplicate_members(right)
    except DuplicateMembersError as dup_member_error:
        diffs.append(dup_member_error)

    if len(right['members']) != len(left['members']):
        diffs.append(MemberLengthDifference(
            'Product Member length differs:'
            ' Left Product {left_product_name} len {left_len} !=  '
            ' Right Product {right_product_name} len {right_len}'
            ''.format(left_product_name=left['name'], left_len=len(left['members']),
                      right_product_name=right['name'], right_len=len(right['members']))
        ))

    members_right = copy(right['members'])
    left_unaccounted_members = []
    for left_member in left['members']:
        for right_member in members_right:
            if left_member['expname'] != right_member['expname']:
                continue

            if left_member['exptype'] != right_member['exptype']:
                diffs.append(MemberMismatchError(
                    'Left {left_expname}:{left_exptype}'
                    ' != Right {right_expname}:{right_exptype}'
                    ''.format(left_expname=left_member['expname'], left_exptype=left_member['exptype'],
                              right_expname=right_member['expname'], right_exptype=right_member['exptype'])
                ))

            members_right.remove(right_member)
            break
        else:
            left_unaccounted_members.append(left_member)

    if len(left_unaccounted_members):
        diffs.append(UnaccountedMembersError(
            f'Left has {len(left_unaccounted_members)} unaccounted members. Members are {left_unaccounted_members}'
        ))

    if len(members_right) != 0:
        diffs.append(UnaccountedMembersError(
            f'Right has {len(members_right)} unaccounted members. Members are {members_right}'
        ))

    # Check if one is a subset of the other.
    err_types = diffs.err_types
    is_subset = (len(diffs) == 2) and (MemberLengthDifference in err_types) and (UnaccountedMembersError in err_types)
    if is_subset:
        diffs = MultiDiffError([SubsetError(f'Products are subsets: {left["name"]} {right["name"]}')])

    if diffs:
        raise diffs


def check_duplicate_members(product):
    """Check for duplicate members in an association product

    The check is based solely on `expname`.

    Parameters
    ----------
    product : dict
        Association product to check.

    Raises
    ------
    MultiDiffError
        If the product has duplicate members.
    """
    seen = set()
    dups = []
    for expname in [member['expname'] for member in product['members']]:
        if expname in seen:
            dups.append(expname)
        else:
            seen.add(expname)

    if dups:
        raise DuplicateMembersError(
            f'Product {product["name"]} has duplicate members {dups}'
        )


def check_duplicate_products(asns):
    """Check for duplicate products in a list of associations

    Duplicate products are defined as any products that share the same name.
    The errors flagged are listed below. If no `MultiDiffError` is raised,
    There are no duplicate products.

    - DuplicateProductError
      The general error for two products that share name and otherwise do
      not fall into any other category

    - SubsetError
      When one or more products are full subsets of another.

    There is a special case where products may have the same name. For Level 2
    processing, exposures can be processed either as time-series, indicated by
    the "rateints" suffix, or as single exposures, indicated by the "rate" suffix.
    Usually any particular exposure is treated one way or the other. However, for
    coronagraphic data, exposures are processed in both ways. In this case, the product
    name and membership will be the same, except that the members will have the different
    suffixes.

    Parameters
    ----------
    asns : [Associations[,...]]
        The associations to compare. Each association should only have one
        product. Use `separate_products` prior to calling if necessary.

    Raises
    ------
    MultiDiffError
    """
    product_names, dup_names = get_product_names(asns)
    if not dup_names:
        return

    # Order associations by candidate. This will allow choosing which
    # association should be considered "primary" when duplicates are found.
    # Categorize by duplicated product names.
    ordered_asns = sort_by_candidate(asns)
    asn_by_product = defaultdict(list)
    for asn in ordered_asns:
        asn_by_product[asn['products'][0]['name']].append(asn)

    diffs = MultiDiffError()
    for product in dup_names:
        dup_asns = asn_by_product[product]
        current_asn = dup_asns.pop()
        for asn in dup_asns:
            try:
                compare_product_membership(current_asn['products'][0], asn['products'][0])
            except MultiDiffError as compare_diffs:
                # Check that the associations differ only in suffix.
                # If so, the associations are not duplicates
                try:
                    compare_nosuffix(current_asn, asn)
                except MultiDiffError:
                    # Not interested in the diffs from `compare_nosuffix`.
                    # The diffs from `compare_product_membership` will be more informative.
                    diffs.extend(compare_diffs)
            else:
                # Associations are exactly the same. Pure duplicate.
                diffs.append(DuplicateProductError(
                    f'Associations share product name {product}',
                    asns=[current_asn, asn]
                ))
    if diffs:
        raise diffs


# #########
# Utilities
# #########
def components(s):
    """split string into its components"""
    return set(re.split('[_-]', s))


def exposure_name(path):
    """Extract the exposure name from a Stage 2 file name

    Parameters
    ----------
    path : Path or str
        The file name or path.

    Returns
    -------
    exposure : str
        The exposure name
    """
    path = Path(path)
    exposure = remove_suffix(path.stem)
    return exposure


def get_product_names(asns):
    """Return product names from associations and flag duplicates

    Parameters
    ----------
    asns: [`Association`[, ...]]

    Returns
    -------
    product_names, duplicates: set(str[, ...]), [str[,...]]
        2-tuple consisting of the set of product names and the list of duplicates.
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
    if dups:
        logger.debug(
            'Duplicate product names: %s', dups
        )

    return set(product_names), dups


def separate_products(asn):
    """Separate products into individual associations

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

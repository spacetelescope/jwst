""" Utilities for product manipulation."""

from collections import defaultdict, Counter
import copy
import logging

from .diff import compare_product_membership

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def sort_by_candidate(asns):
    """Sort associations by candidate

    Parameters
    ----------
    asns: [Association[,...]]
        List of associations

    Returns
    -------
    sorted_by_candidate: [Associations[,...]]
        New list of the associations sorted.

    Notes
    -----
    The current definition of candidates allows strictly lexigraphical
    sorting:
    aXXXX > cXXXX > oXXX

    If this changes, a comparision function will need be implemented
    """
    return sorted(asns, key=lambda asn: asn['asn_id'])


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


def prune_duplicate_associations(asns):
    """Remove duplicate associations in favor of lower level versions

    For Level 3 associations, multiple associations with the same membership, but different
    levels, can be created. Remove duplicate associations of higher level.

    The assumption is that there is only one product per association, before merging.

    Parameters
    ----------
    asns: [Association[,...]]
        Associations to prune

    Returns
    pruned: [Association[,...]]
        Pruned list of associations
    """
    ordered_asns = sort_by_candidate(asns)
    pruned = list()
    while True:
        try:
            original = ordered_asns.pop()
        except IndexError:
            break
        pruned.append(original)
        to_prune = list()
        for asn in ordered_asns:
            try:
                compare_product_membership(original['products'][0], asn['products'][0])
            except AssertionError:
                continue
            to_prune.append(asn)
        for prune in to_prune:
            ordered_asns.remove(prune)

    return pruned


def prune_duplicate_products(asns):
    """Remove duplicate products in favor of higher level versions

    For Level 2 associations, since the products are always just the input
    exposures, different candidates can be created for each exposure. Remove
    those associations of lesser candidates.

    The assumption is that there is only one product per association, before merging

    Parameters
    ----------
    asns: [Association[,...]]
        Associations to prune

    Returns
    pruned: [Association[,...]]
        Pruned list of associations
    """
    product_names, dups = get_product_names(asns)
    if not dups:
        return asns

    pruned = copy.copy(asns)
    to_prune = defaultdict(list)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if product_name in dups:
            to_prune[product_name].append(asn)

    for product_name, asns_to_prune in to_prune.items():
        asns_to_prune = sort_by_candidate(asns_to_prune)
        for asn in asns_to_prune[1:]:
            pruned.remove(asn)

    return pruned

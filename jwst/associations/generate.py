import logging

from astropy.table import Table

from .association import (
    AssociationError,
    make_timestamp
)

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def generate(pool, rules):
    """Generate associations in the pool according to the rules.

    Parameters
    ----------
    pool: AssociationPool
        The pool to generate from.

    rules: Associations
        The associaton rule set.

    Returns
    -------
    ([association,], orphans)
        A 2-tuple consisting of:
           * List of associations
           * Table of members from the pool that
             do not belong to any association.
    """
    associations = []
    orphaned = Table(dtype=pool.dtype)
    timestamp = make_timestamp()
    logger.debug('Starting...')
    for member in pool:
        logger.debug('Working member="{}"'.format(member))
        member_associations = generate_from_member(
            member,
            timestamp,
            associations,
            rules
        )
        if len(member_associations) == 0:
            orphaned.add_row(member)

    return associations, orphaned


def generate_from_member(member, timestamp, associations, rules):
    """Either match or generate a new assocation

    Parameters
    ----------
    member: dict
        The member to match to existing associations
        or generate new associations from

    timestamp: str
        Timestamp to use with association creation.

    associations: [association, ...]
        List of already existing associations.
        If the member matches any of these, it will be added
        to them.

    rules: AssociationRegistry
        List of rules to create new associations

    Returns
    -------
    [association,...]
        List of associations member belongs to. Empty if none match
    """

    # Check membership in existing associations.
    member_associations = match_member(member, associations)

    # Now see if this member will create new associatons.
    # By default, a member will not be allowed to create
    # an association based on rules of existing associations.
    ignore_asns = set([type(asn) for asn in member_associations])
    try:
        new_asns = rules.match(
            member,
            timestamp=timestamp,
            ignore=ignore_asns,
        )
    except AssociationError as error:
        logger.debug('Did not match any rule.')
        logger.debug('error="{}"'.format(error))
    else:
        member_associations.extend(new_asns)
        logger.debug(
            'Member created new associations "{}"'.format(new_asns)
        )

    return member_associations


def match_member(member, associations):
    """Match member to a list of associations

    Parameters
    ----------
    member: dict
        The member to match to the associations.

    associations: [association, ...]
        List of already existing associations.
        If the member matches any of these, it will be added
        to them.

    Returns
    -------
    [association,...]
        List of associations member belongs to. Empty if none match
    """
    member_associations = []
    for asn in associations:
        if asn in member_associations:
            continue
        try:
            asn.add(member)
        except AssociationError as error:
            logger.debug(
                'Did not match association "{}"'.format(asn)
            )
            logger.debug('error="{}"'.format(error))
            continue
        else:
            logger.debug('Matched association "{}"'.format(asn))
            member_associations.append(asn)
    return member_associations

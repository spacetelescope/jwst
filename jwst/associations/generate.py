import logging

from astropy.table import Table

from .association import (
    make_timestamp
)
from .exceptions import (
    AssociationError,
    AssociationProcessMembers
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
    logger.debug('Starting...')

    associations = []
    orphaned = Table(dtype=pool.dtype)
    timestamp = make_timestamp()
    process_list = [
        AssociationProcessMembers(
            members=pool,
            allowed_rules=[rule for _, rule in rules.items()]
        )
    ]

    for process_event in process_list:
        for member in process_event.members:
            logger.debug('Working member="{}"'.format(member))
            existing_asns, new_asns, to_process = generate_from_member(
                member,
                timestamp,
                associations,
                rules,
                process_event.allowed_rules
            )
            associations.extend(new_asns)
            process_list.extend(to_process)

    logger.debug('Number of processes: "{}"'.format(len(process_list)))
    return associations, orphaned


def generate_from_member(
        member,
        timestamp,
        associations,
        rules,
        allowed_rules):
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

    allwed_rules: [rule, ...]
        Only compare existing associations and rules if the
        the rule is in this list. If none, all rules are valid.

    Returns
    -------
    (associations, process_list): 3-tuple where
        existing_asns: [association,...]
            List of existing associations member belongs to.
            Empty if none match
        new_asns: [association,...]
            List of new associations member creates. Empty if none match
        process_list: [AssociationProcess, ...]
            List of process events.
    """

    # Check membership in existing associations.
    associations = [
        asn
        for asn in associations
        if type(asn) in allowed_rules
    ]
    existing_asns, process_list = match_member(member, associations)

    # Now see if this member will create new associatons.
    # By default, a member will not be allowed to create
    # an association based on rules of existing associations.
    ignore_asns = set([type(asn) for asn in existing_asns])
    new_asns, to_process = rules.match(
        member,
        timestamp=timestamp,
        allow=allowed_rules,
        ignore=ignore_asns,
    )
    logger.debug(
        'Member created new associations "{}"'.format(new_asns)
    )

    process_list.extend(to_process)
    return existing_asns, new_asns, process_list


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
    (associations, process_list): 2-tuple where
        associations: [association,...]
            List of associations member belongs to. Empty if none match
        process_list: [AssociationProcess, ...]
            List of process events.
    """
    member_associations = []
    process_list = []
    for asn in associations:
        if asn in member_associations:
            continue
        try:
            asn.add(member)
        except AssociationError as error:
            logger.debug(
                'Did not match association "{}"'.format(asn)
            )
            logger.debug('Reason="{}"'.format(error))
        except AssociationProcessMembers as process_event:
            logger.debug('Process event "{}"'.format(process_event))
            process_list.append(process_event)
        else:
            logger.debug('Matched association "{}"'.format(asn))
            member_associations.append(asn)
    return member_associations, process_list

import logging

from astropy.table import Table

from jwst.associations.association import (AssociationError,
                                                 make_timestamp)

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
    ([association, ], [orphans,])
        A 2-tuple consisting of:
           * List of associations
           * List of members from the pool that
             do not belong to any association.
    """
    asns = []
    orphaned = Table(dtype=pool.dtype)
    timestamp = make_timestamp()
    logger.debug('Starting...')
    for member in pool:
        logger.debug('Working member="{}"'.format(member))
        member_asns = set()
        is_a_member = False
        for asn in asns:
            try:
                asn.add(member)
                member_asns.add(type(asn))
                is_a_member = True
                logger.debug('Matched association "{}"'.format(asn))
            except AssociationError as error:
                logger.debug(
                    'Did not match association "{}"'.format(asn)
                )
                logger.debug('error="{}"'.format(error))
                continue
        try:
            asns.extend(rules.match(member, timestamp, member_asns))
            is_a_member = True
            logger.debug(
                'Member created new association "{}"'.format(asns[-1])
            )
        except AssociationError as error:
            logger.debug('Did not match any rule.')
            logger.debug('error="{}"'.format(error))
        if not is_a_member:
            orphaned.add_row(member)
    return asns, orphaned

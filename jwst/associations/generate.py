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
    asns = []
    orphaned = Table(dtype=pool.dtype)
    timestamp = make_timestamp()
    logger.debug('Starting...')
    for member in pool:
        logger.debug('Working member="{}"'.format(member))
        member_asns = set()
        for asn in asns:
            try:
                asn.add(member)
            except AssociationError as error:
                logger.debug(
                    'Did not match association "{}"'.format(asn)
                )
                logger.debug('error="{}"'.format(error))
                continue
            else:
                member_asns.add(type(asn))
                logger.debug('Matched association "{}"'.format(asn))

        try:
            new_asns = rules.match(member, timestamp, member_asns)
        except AssociationError as error:
            logger.debug('Did not match any rule.')
            logger.debug('error="{}"'.format(error))
        else:
            asns.extend(new_asns)
            member_asns.update(
                type(asn)
                for asn in new_asns
            )
            logger.debug(
                'Member created new associations "{}"'.format(new_asns)
            )

        if len(member_asns) == 0:
            orphaned.add_row(member)

    return asns, orphaned

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
    def reprocess_member(member):
        """A member which should be reprocessed."""
        if member is None:
            return
        logger.debug('Reprocess request for member="{}"'.format(member))
        working_pool.add_row(member)

    asns = []
    orphaned = Table(dtype=pool.dtype)
    timestamp = make_timestamp()
    working_pool = pool.copy(copy_data=False)
    n_orignal = len(working_pool)
    logger.debug('Starting...')
    for member in working_pool:
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
            new_asns = rules.match(
                member,
                timestamp=timestamp,
                ignore=member_asns,
                reprocess_cb=reprocess_member
            )
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

        # If no associations were matched or made
        # AND
        # if the member is not reprocessed member
        # THEN
        # orphan it.
        if len(member_asns) == 0 and member.index < n_orignal:
            orphaned.add_row(member)

    return asns, orphaned

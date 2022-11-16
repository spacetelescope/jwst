import logging
from timeit import default_timer as timer

from .association import make_timestamp
from .lib.process_list import (
    ListCategory,
    ProcessList,
    ProcessQueueSorted,
    workover_filter
)
from .pool import PoolRow
from ..lib.progress import Bar

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['generate']


def generate(pool, rules, version_id=None):
    """Generate associations in the pool according to the rules.

    Parameters
    ----------
    pool : AssociationPool
        The pool to generate from.

    rules : AssociationRegistry
        The association rule set.

    version_id : None, True, or str
        The string to use to tag associations and products.
        If None, no tagging occurs.
        If True, use a timestamp
        If a string, the string.

    Returns
    -------
    associations : [Association[,...]]
        List of associations

    Notes
    -----
    Refer to the :ref:`Association Generator <design-generator>`
    documentation for a full description.
    """
    associations = []
    if type(version_id) is bool:
        version_id = make_timestamp()
    process_queue = ProcessQueueSorted([
        ProcessList(
            items=pool,
            rules=[rule for _, rule in rules.items()]
        )
    ])

    logger.debug('Initial process queue: %s', process_queue)
    for process_list in process_queue:
        logger.debug('** Working process list: %s', process_list)
        time_start = timer()
        total_mod_existing = 0
        total_new = 0
        total_reprocess = 0
        with Bar('Processing items', log_level=logger.getEffectiveLevel(),
                 max=len(process_list.items)) as bar:
            for item in process_list.items:
                item = PoolRow(item)

                existing_asns, new_asns, to_process = generate_from_item(
                    item,
                    version_id,
                    associations,
                    rules,
                    process_list
                )
                total_mod_existing += len(existing_asns)
                total_new += len(new_asns)
                associations.extend(new_asns)

                # If working on a process list EXISTING
                # remove any new `to_process` that is
                # also EXISTING. Prevent infinite loops.
                to_process_modified = []
                for next_list in to_process:
                    next_list = workover_filter(next_list, process_list.work_over)
                    if next_list:
                        to_process_modified.append(next_list)
                process_queue.extend(to_process_modified)
                total_reprocess += len(to_process_modified)
                bar.next()

        logger.debug('Existing associations modified: %d New associations created: %d', total_mod_existing, total_new)
        logger.debug('New process lists: %d', total_reprocess)
        logger.debug('Updated process queue: %s', process_queue)
        logger.debug('# associations: %d', len(associations))
        logger.debug('Seconds to process: %.2f\n', timer() - time_start)

    # Finalize found associations
    logger.debug('# associations before finalization: %d', len(associations))
    try:
        finalized_asns = rules.callback.reduce('finalize', associations)
    except KeyError:
        finalized_asns = associations

    return finalized_asns


def generate_from_item(
        item,
        version_id,
        associations,
        rules,
        process_list):
    """Either match or generate a new association

    Parameters
    ----------
    item : dict
        The item to match to existing associations
        or generate new associations from

    version_id : str or None
        Version id to use with association creation.
        If None, no versioning is used.

    associations : [association, ...]
        List of already existing associations.
        If the item matches any of these, it will be added
        to them.

    rules : AssociationRegistry or None
        List of rules to create new associations

    process_list : ProcessList
        The `ProcessList` from which the current item belongs to.

    Returns
    -------
    (associations, process_list): 3-tuple where
        existing_asns : [association,...]
            List of existing associations item belongs to.
            Empty if none match
        new_asns : [association,...]
            List of new associations item creates. Empty if none match
        process_list : [ProcessList, ...]
            List of process events.
    """

    # Setup the rules allowed to be examined.
    if process_list.rules is None or len(process_list.rules) == 0:
        allowed_rules = list(rules.values())
    else:
        allowed_rules = process_list.rules

    # Check membership in existing associations.
    existing_asns = []
    reprocess_list = []
    if process_list.work_over in (
            ListCategory.BOTH,
            ListCategory.EXISTING,
            ListCategory.NONSCIENCE,
    ):
        associations = [
            asn
            for asn in associations
            if type(asn) in allowed_rules
        ]
        existing_asns, reprocess_list = match_item(
            item, associations
        )

    # Now see if this item will create new associations.
    # By default, a item will not be allowed to create
    # an association based on rules of existing associations.
    reprocess = []
    new_asns = []
    if process_list.work_over in (
            ListCategory.BOTH,
            ListCategory.RULES,
    ) and rules is not None:
        ignore_asns = set([type(asn) for asn in existing_asns])
        new_asns, reprocess = rules.match(
            item,
            version_id=version_id,
            allow=allowed_rules,
            ignore=ignore_asns,
        )
    reprocess_list.extend(reprocess)

    return existing_asns, new_asns, reprocess_list


def match_item(item, associations):
    """Match item to a list of associations

    Parameters
    ----------
    item : dict
        The item to match to the associations.

    associations : [association, ...]
        List of already existing associations.
        If the item matches any of these, it will be added
        to them.

    Returns
    -------
    (associations, process_list): 2-tuple where
        associations : [association,...]
            List of associations item belongs to. Empty if none match
        process_list : [ProcessList, ...]
            List of process events.
    """
    item_associations = []
    process_list = []
    for asn in associations:
        if asn in item_associations:
            continue
        matches, reprocess = asn.add(item)
        process_list.extend(reprocess)
        if matches:
            item_associations.append(asn)
    return item_associations, process_list

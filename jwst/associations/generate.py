import logging

from .association import (
    make_timestamp
)
from .lib.process_list import (
    ProcessList,
    ProcessQueueSorted
)

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
        The associaton rule set.

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

    for process_list in process_queue:
        # logger.debug(
        #     'Working process list:'
        #     f'\n\t#items {len(process_list.items)}'
        #     f' working over {process_list.work_over}'
        #     f' matching on {process_list.only_on_match}'
        #     f'\n\trules {process_list.rules}'
        # )
        for item in process_list.items:
            # logger.debug(f'Processing item {item}')
            existing_asns, new_asns, to_process = generate_from_item(
                item,
                version_id,
                associations,
                rules,
                process_list
            )
            # logger.debug(f'Associations updated: {existing_asns}')
            # logger.debug(f'New associations: {new_asns}')
            associations.extend(new_asns)

            # If working on a process list EXISTING
            # remove any new `to_process` that is
            # also EXISTING. Prevent infinite loops.
            if process_list.work_over in (ProcessList.EXISTING, ProcessList.NONSCIENCE):
                to_process = [
                    to_process_list
                    for to_process_list in to_process
                    if to_process_list.work_over != process_list.work_over
                ]
            process_queue.extend(to_process)

    # Finalize found associations
    logger.debug('# associations before finalization: %s', len(associations))
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
    """Either match or generate a new assocation

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
            ProcessList.BOTH,
            ProcessList.EXISTING,
            ProcessList.NONSCIENCE,
    ):
        associations = [
            asn
            for asn in associations
            if type(asn) in allowed_rules
        ]
        # logger.debug(
        #     f'Checking against {len(associations)} existing associations'
        # )
        existing_asns, reprocess_list = match_item(
            item, associations
        )

    # Now see if this item will create new associatons.
    # By default, a item will not be allowed to create
    # an association based on rules of existing associations.
    reprocess = []
    new_asns = []
    if process_list.work_over in (
            ProcessList.BOTH,
            ProcessList.RULES,
    ) and rules is not None:
        ignore_asns = set([type(asn) for asn in existing_asns])
        # logger.debug(f'Ignore asns {ignore_asns}')
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

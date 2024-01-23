import collections
import logging

from .generate import generate
from .lib.utilities import evaluate
from .main import CANDIDATE_RULESET, DISCOVER_RULESET, constrain_on_candidates
from .registry import AssociationRegistry

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def generate_per_candidate(pool, user_rules, cids=None, version_id=None, finalize=True, ignore_default=False, discover=False):
    """Generate associations in the pool according to the rules.

    Parameters
    ----------
    pool : AssociationPool
        The pool to generate from.

    user_rules : File-like
        The rule definitions to use. None to use the defaults if `ignore_default` is False.

    cids : [str,[...]] or None
        List of candidates to produce for. If None, do all possible candidates

    version_id : None, True, or str
        The string to use to tag associations and products.
        If None, no tagging occurs.
        If True, use a timestamp
        If a string, the string.

    finalize : bool
        Run all rule methods marked as 'finalized'.

    ignore_default : bool
        Ignore the default rules. Use only the user-specified ones.

    discover : bool
        Find associations that are not candidate-based.

    Returns
    -------
    associations : [Association[,...]]
        List of associations

    Notes
    -----
    Refer to the :ref:`Association Generator <design-generator>`
    documentation for a full description.
    """
    # Get the candidates
    cids_by_type = ids_by_ctype(pool)
    if cids is None:
        cids_ctypes = [(cid, ctype)
                       for ctype in cids_by_type
                       for cid in cids_by_type[ctype]
                       ]
    else:
        cids_ctypes = []
        for cid in cids:
            for ctype in cids_by_type:
                if cid in cids_by_type[ctype]:
                    cids_ctypes.append((cid, ctype))
                    break
            else:
                logger.warning('cid %s not found in pool', cid)

    associations = []
    for cid_ctype in cids_ctypes:
        cid, ctype = cid_ctype
        logger.info(f'Working on {cid_ctype}')

        # Get the pool
        pool_cid = pool_from_candidate(pool, cid)
        pool_cid['asn_candidate'] = [f"[('{cid}', '{ctype}')]"] * len(pool_cid)
        logger.info(f'Len(pool_{cid}): {len(pool_cid)}')

        # Create the rules with the simplified asn_candidate constraint
        asn_constraint = constrain_on_candidates([cid])
        rules = AssociationRegistry(user_rules, include_default=not ignore_default, global_constraints=asn_constraint, name=CANDIDATE_RULESET)

        # Get the associations
        associations_cid = generate(pool_cid, rules, version_id=version_id, finalize=False)

        # Add to the list
        associations.extend(associations_cid)

    # The ruleset has been generated on a per-candidate case.
    # Here, need to do a final rebuild of the ruleset to get the finalization
    # functions. This ruleset does not need any of the candidate specifications.
    # This ruleset is also used if discovery is in play.
    rules = AssociationRegistry(user_rules, include_default=not ignore_default, name=DISCOVER_RULESET)
    if discover:
        logger.info('Discovering associations...')
        associations_discover = generate(pool, rules, version_id=version_id, finalize=False)
        logger.info('# discovered associations: %d', len(associations_discover))
        associations.extend(associations_discover)

    # Finalize found associations
    logger.debug('# associations before finalization: %d', len(associations))
    finalized_asns = associations
    if finalize and len(associations):
        logger.debug('Performing association finalization.')

        try:
            finalized_asns = rules.callback.reduce('finalize', associations)
        except KeyError as exception:
            logger.debug('Finalization failed for reason: %s', exception)

    logger.info('Associations generated: %s', len(finalized_asns))
    return finalized_asns


def ids_by_ctype(pool):
    """Groups candidate ids by the candidate type

    Parameters
    ----------
    pool : AssociationPool
        The association pool

    Returns
    -------
    ids_by_ctype : {ctype: counter}
        Dict with the key of the candidate type. Value is a
        `collections.Counter` object of the ids and their counts.
    """
    ids_by_ctype = collections.defaultdict(list)
    for exp_candidates in pool['asn_candidate']:
        candidates = evaluate(exp_candidates)
        if isinstance(candidates, int):
            ids_by_ctype['unknown'].append(str(candidates))
            continue
        try:
            for id, ctype in candidates:
                ids_by_ctype[ctype].append(id)
        except ValueError:
            logger.debug('Cannot parse asn_candidate field: %s', candidates)

    for ctype in ids_by_ctype:
        ids_by_ctype[ctype] = collections.Counter(ids_by_ctype[ctype])

    return ids_by_ctype


def pool_from_candidate(pool, candidate):
    """Create a pool containing only the candidate

    Parameters
    ----------
    pool : AssociationPool
        The pool to filter from.

    candidate : str
        The candidate id to filter.

    Returns
    -------
    candidate_pool : AssociationPool
        Pool containing only the candidate
    """
    candidate_pool = pool[[candidate in row['asn_candidate'] for row in pool]]
    return candidate_pool

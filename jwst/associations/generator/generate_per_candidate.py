import collections
import logging
import numpy as np
from timeit import default_timer as timer

from .generate import generate
from .generate_per_pool import CANDIDATE_RULESET, DISCOVER_RULESET, constrain_on_candidates
from jwst.associations.lib.utilities import evaluate, filter_discovered_only
from jwst.associations.registry import AssociationRegistry

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def generate_per_candidate(
    pool,
    rule_defs,
    candidate_ids=None,
    all_candidates=True,
    discover=False,
    version_id=None,
    finalize=True,
    merge=False,
    ignore_default=False,
    dms_enabled=False,
):
    """
    Generate associations in the pool according to the rules.

    Parameters
    ----------
    pool : AssociationPool
        The pool to generate from.

    rule_defs : [File-like[,...]] or None
        The rule definitions to use. None to use the defaults if `ignore_default` is False.

    candidate_ids : [str,[...]] or None
        List of candidates to produce for. If None, do all possible candidates

    all_candidates : bool
        Keep associations generated from candidates when in discovery mode.

    discover : bool
        Find associations that are not candidate-based.

    version_id : None, True, or str
        The string to use to tag associations and products.
        If None, no tagging occurs.
        If True, use a timestamp
        If a string, the string.

    finalize : bool
        Run all rule methods marked as 'finalized'.

    merge : bool
        Merge single associations into a common association with separate products.

    ignore_default : bool
        Ignore the default rules. Use only the user-specified ones.

    dms_enabled : bool
        Flag for DMS processing, true if command-line argument '--DMS' was used.

    Returns
    -------
    associations : [Association[,...]]
        List of associations

    Notes
    -----
    Refer to the :ref:`Association Generator <design-generator>`
    documentation for a full description.
    """
    logger.info("Generating based on the per-candidate algorithm.")

    # If DMS flag is present, add all c-type ids associated with input ids to generation
    if dms_enabled and candidate_ids is not None:
        # Generate a selection mask to find all cids associated with DMS-provided ids
        row_mask = np.zeros((len(pool),), dtype=bool)
        for cid in candidate_ids:
            row_mask |= [cid in x for x in pool["asn_candidate"]]
        input_candidate_ids = list(candidate_ids)
        bkg_cids = ids_by_ctype(pool[row_mask]).get("background", None)
        if bkg_cids is not None:
            for key in bkg_cids:
                if key not in candidate_ids:
                    candidate_ids.append(key)

    # Get the candidates
    cids_by_type = ids_by_ctype(pool)
    if candidate_ids is None:
        cids_ctypes = [(cid, ctype) for ctype in cids_by_type for cid in cids_by_type[ctype]]
    else:
        cids_ctypes = []
        for cid in candidate_ids:
            for ctype in cids_by_type:
                if cid in cids_by_type[ctype]:
                    cids_ctypes.append((cid, ctype))
                    break
            else:
                logger.warning("Candidate id %s not found in pool", cid)

    associations = []
    for cid_ctype in cids_ctypes:
        time_start = timer()
        # Generate the association for the given candidate
        associations_cid = generate_on_candidate(
            cid_ctype,
            pool,
            rule_defs,
            version_id=version_id,
            ignore_default=ignore_default,
        )

        # Add to the list
        associations.extend(associations_cid)

        logger.info("Time to process candidate %s: %.2f", cid_ctype[0], timer() - time_start)

    # The ruleset has been generated on a per-candidate case.
    # Here, need to do a final rebuild of the ruleset to get the finalization
    # functions. This ruleset does not need any of the candidate specifications.
    # This ruleset is also used if discovery is in play.
    rules = AssociationRegistry(
        rule_defs, include_default=not ignore_default, name=DISCOVER_RULESET
    )
    if discover:
        logger.info("Discovering associations...")
        associations_discover = generate(pool, rules, version_id=version_id, finalize=False)
        associations.extend(associations_discover)
        logger.info("# associations found before discover filter: %d", len(associations_discover))
        associations = filter_discovered_only(
            associations,
            DISCOVER_RULESET,
            CANDIDATE_RULESET,
            keep_candidates=all_candidates,
        )
        rules.Utility.resequence(associations)

    # Finalize found associations
    logger.debug("# associations before finalization: %d", len(associations))
    finalized_asns = associations
    if finalize and len(associations):
        logger.debug("Performing association finalization.")

        try:
            finalized_asns = rules.callback.reduce("finalize", associations)
        except KeyError as exception:
            logger.debug("Finalization failed for reason: %s", exception)

    # Do a grand merging. This is done particularly for
    # Level2 associations.
    if merge:
        try:
            finalized_asns = rules.Utility.merge_asns(finalized_asns)
        except AttributeError:
            pass

    if dms_enabled and candidate_ids is not None:
        # We now must remove the asns that were used to duplicate-check the input candidates.
        finalized_asns = [
            asn
            for asn in finalized_asns
            if any(cid in asn["asn_id"] for cid in input_candidate_ids)
        ]

    logger.info("Total associations generated: %s", len(finalized_asns))
    return finalized_asns


def generate_on_candidate(cid_ctype, pool, rule_defs, version_id=None, ignore_default=False):
    """
    Generate associations based on a candidate ID.

    Parameters
    ----------
    cid_ctype : (str, str)
        2-tuple of candidate ID and the candidate type

    pool : AssociationPool
        The pool to generate from.

    rule_defs : [File-like[,...]] or None
        The rule definitions to use. None to use the defaults if `ignore_default` is False.

    version_id : None, True, or str
        The string to use to tag associations and products.
        If None, no tagging occurs.
        If True, use a timestamp
        If a string, the string.

    ignore_default : bool
        Ignore the default rules. Use only the user-specified ones.

    Returns
    -------
    associations : [Association[,...]]
        List of associations
    """
    cid, ctype = cid_ctype
    logger.info(f"Generating associations on candidate {cid_ctype}")

    # Get the pool
    pool_cid = pool_from_candidate(pool, cid)

    pool_cid["asn_candidate"] = [f"[('{cid}', '{ctype}')]"] * len(pool_cid)
    logger.info(f"Length of pool for {cid}: {len(pool_cid)}")

    # Create the rules with the simplified asn_candidate constraint
    asn_constraint = constrain_on_candidates([cid])
    rules = AssociationRegistry(
        rule_defs,
        include_default=not ignore_default,
        global_constraints=asn_constraint,
        name=CANDIDATE_RULESET,
    )

    # Get the associations
    associations = generate(pool_cid, rules, version_id=version_id, finalize=False)

    return associations


def ids_by_ctype(pool):
    """
    Group candidate ids by the candidate type.

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
    for exp_candidates in pool["asn_candidate"]:
        candidates = evaluate(exp_candidates)
        if isinstance(candidates, int):
            ids_by_ctype["unknown"].append(str(candidates))
            continue
        try:
            for cand_id, ctype in candidates:
                ids_by_ctype[ctype].append(cand_id)
        except ValueError:
            logger.debug("Cannot parse asn_candidate field: %s", candidates)

    for ctype in ids_by_ctype:
        ids_by_ctype[ctype] = collections.Counter(ids_by_ctype[ctype])

    return ids_by_ctype


def pool_from_candidate(pool, candidate):
    """
    Create a pool containing only the candidate.

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
    candidate_pool = pool[[candidate in row["asn_candidate"] for row in pool]]
    return candidate_pool

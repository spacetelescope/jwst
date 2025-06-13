# Generate associations per-pool
import logging

from .generate import generate
from jwst.associations.lib.utilities import constrain_on_candidates, filter_discovered_only
from jwst.associations.registry import AssociationRegistry

__all__ = ["generate_per_pool"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Ruleset names
DISCOVER_RULESET = "discover"
CANDIDATE_RULESET = "candidate"


def generate_per_pool(
    pool,
    rule_defs=None,
    candidate_ids=None,
    all_candidates=True,
    discover=False,
    version_id=None,
    finalize=True,
    merge=False,
    ignore_default=False,
):
    """
    Generate associations on a specified pool.

    Association candidates are filtered based on global constraints added to the rules.

    Parameters
    ----------
    pool : AssociationPool
        The pool to generate from.

    rule_defs : [File-like[,...] or None
        List of rule definitions to use. None to use the defaults if `ignore_default` is False.

    candidate_ids : [str,[...]] or None
        List of candidates to produce for.

    all_candidates : bool
        Find associations for all possible candidates

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

    Returns
    -------
    associations : [Association[,...]]
        List of associations

    Notes
    -----
    Refer to the :ref:`Association Generator <design-generator>`
    documentation for a full description.
    """
    logger.info("Generating based on the per-pool algorithm")

    # Setup the rule registry
    global_constraints = None
    if discover or all_candidates:
        global_constraints = constrain_on_candidates(None)
    elif candidate_ids is not None:
        global_constraints = constrain_on_candidates(candidate_ids)

    rules = AssociationRegistry(
        rule_defs,
        include_default=not ignore_default,
        global_constraints=global_constraints,
        name=CANDIDATE_RULESET,
    )

    if discover:
        rules.update(
            AssociationRegistry(
                rule_defs, include_default=not ignore_default, name=DISCOVER_RULESET
            )
        )

    # Generate the associations
    associations = generate(pool, rules, version_id=version_id, finalize=finalize)

    # If doing discover, filter out all specified candidates
    if discover:
        logger.debug(f"# asns found before discover filtering={len(associations)}")
        associations = filter_discovered_only(
            associations,
            DISCOVER_RULESET,
            CANDIDATE_RULESET,
            keep_candidates=all_candidates,
        )
        rules.Utility.resequence(associations)

    # Do a grand merging. This is done particularly for
    # Level2 associations.
    if merge:
        try:
            associations = rules.Utility.merge_asns(associations)
        except AttributeError:
            pass

    # That's all folks
    return associations

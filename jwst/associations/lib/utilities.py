"""General Utilities."""

from ast import literal_eval
from functools import wraps
import logging

from numpy.ma import masked

from jwst.associations import config

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def constrain_on_candidates(candidates):
    """
    Create a constraint based on a list of candidates.

    Parameters
    ----------
    candidates : (str, ...) or None
        List of candidate id's.
        If None, then all candidates are matched.

    Returns
    -------
    DMSAttrConstraint
        The constraint built off the candidate list.
    """
    from .dms_base import DMSAttrConstraint

    if candidates is not None and len(candidates):
        c_list = "|".join(candidates)
        values = "".join([".+(", c_list, ").+"])
    else:
        values = None
    constraint = DMSAttrConstraint(
        name="asn_candidate",
        sources=["asn_candidate"],
        value=values,
        force_unique=True,
        is_acid=True,
        evaluate=True,
    )

    return constraint


def evaluate(value):
    """
    Evaluate a value.

    Parameters
    ----------
    value : str
        The string to evaluate.

    Returns
    -------
    type or str
        The evaluation. If the value cannot be
        evaluated, the value is simply returned
    """
    try:
        evaled = literal_eval(value)
    except (ValueError, SyntaxError):
        evaled = value
    return evaled


def filter_discovered_only(
    associations,
    discover_ruleset,
    candidate_ruleset,
    keep_candidates=True,
):
    """
    Return only those associations that have multiple candidates.

    Parameters
    ----------
    associations : iterable
        The list of associations to check. The list
        is that returned by the `generate` function.
    discover_ruleset : str
        The name of the ruleset that has the discover rules
    candidate_ruleset : str
        The name of the ruleset that finds just candidates
    keep_candidates : bool
        Keep explicit candidate associations in the list.

    Returns
    -------
    iterable
        The new list of just cross candidate associations.

    Notes
    -----
    This utility is only meant to run on associations that have
    been constructed. Associations that have been Association.dump
    and then Association.load will not return proper results.
    """
    from .prune import identify_dups

    # Split the associations along discovered/not discovered lines
    dups, valid = identify_dups(associations)
    asn_by_ruleset = {candidate_ruleset: [], discover_ruleset: []}
    for asn in valid:
        asn_by_ruleset[asn.registry.name].append(asn)
    candidate_list = asn_by_ruleset[candidate_ruleset]
    discover_list = asn_by_ruleset[discover_ruleset]

    # Filter out the non-unique discovered.
    for candidate in candidate_list:
        if len(discover_list) == 0:
            break
        unique_list = []
        for discover in discover_list:
            if discover != candidate:
                unique_list.append(discover)

        # Reset the discovered list to the new unique list
        # and try the next candidate.
        discover_list = unique_list

    if keep_candidates:
        discover_list.extend(candidate_list)

    if config.DEBUG:
        discover_list += dups
    return discover_list


def getattr_from_list(adict, attributes, invalid_values=None):
    """
    Retrieve value from dict using a list of attributes.

    Parameters
    ----------
    adict : dict
        Dictionary to retrieve from.
    attributes : list
        List of attributes.
    invalid_values : set
        A set of values that essentially mean the
        attribute does not exist.

    Returns
    -------
    (attribute, value)
        Returns the value and the attribute from
        which the value was taken.

    Raises
    ------
    KeyError
        None of the attributes are found in the dict.
    """
    if invalid_values is None:
        invalid_values = set()

    for attribute in attributes:
        try:
            result = adict[attribute]
        except KeyError:
            continue
        else:
            if result is masked:
                continue
            if result not in invalid_values:
                return attribute, result
            else:
                continue
    else:
        raise KeyError(f"Object has no attributes in {attributes}")


def return_on_exception(exceptions=(Exception,), default=None):
    """
    Force functions raising exceptions to return a value.

    This function returns a decorator to accomplish the value return.

    Parameters
    ----------
    exceptions : (Exception(,...))
        Tuple of exceptions to catch.

    default : obj
        The value to return when a specified exception occurs.

    Returns
    -------
    decorator
        The decorator to wrap functions that will return on certain exceptions.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except exceptions as err:
                logger.debug(
                    "Caught exception %s in function %s, forcing return value of %s",
                    err,
                    func,
                    default,
                )
                return default

        return wrapper

    return decorator


@return_on_exception(exceptions=(KeyError,), default=None)
def getattr_from_list_nofail(*args, **kwargs):
    """
    Call getattr_from_list without allowing exceptions.

    If the specified exceptions are caught, return `default`
    instead.

    Parameters
    ----------
    *args, **kwargs : dicts
        Arguments passed to getattr_from_list.

    Returns
    -------
    (attribute, value)
        Returns the value and the attribute from
        which the value was taken.

    Raises
    ------
    KeyError
        None of the attributes are found in the dict.
    """
    return getattr_from_list(*args, **kwargs)


def is_iterable(obj):
    """
    General iterator check.

    Parameters
    ----------
    obj : object
        The object to be checked for defined `__iter__` method.

    Returns
    -------
    bool
        True if iterable, false otherwise.
    """
    return not isinstance(obj, str) and not isinstance(obj, tuple) and hasattr(obj, "__iter__")

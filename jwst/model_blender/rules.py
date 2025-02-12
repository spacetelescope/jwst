import numpy as np


__all__ = ["RULE_FUNCTIONS", "AttributeBlender", "make_blender"]


def _multi(vals):
    """
    Return either the common value from a list of identical values or 'MULTIPLE'.

    Parameters
    ----------
    vals : list
        List of values to check.

    Returns
    -------
    value : str or None
        The common value from the list of identical values or 'MULTIPLE'.
    """
    uniq_vals = list(set(vals))
    num_vals = len(uniq_vals)
    if num_vals == 0:
        return None
    if num_vals == 1:
        return uniq_vals[0]
    if num_vals > 1:
        return "MULTIPLE"


RULE_FUNCTIONS = {
    "multi": _multi,
    "mean": np.mean,
    "sum": np.sum,
    "max": np.max,
    "min": np.min,
    # retained date/time names for backwards compatibility
    # as these all assume ISO8601 format the lexical and
    # chronological sorting match
    "mintime": min,
    "maxtime": max,
    "mindate": min,
    "maxdate": max,
    "mindatetime": min,
    "maxdatetime": max,
}
"""
Mapping of rule names to functions.

Used for `make_blender`.

The following rules are considered deprecated
and should not be used for new schemas.

  - mintime
  - maxtime
  - mindate
  - maxdate
  - mindatetime
  - maxdatetime
"""


class AttributeBlender:
    """Single attribute metadata blender."""

    def __init__(self, blend_function):
        """
        Create a new metadata attribute blender.

        Parameters
        ----------
        blend_function : callable
            Function to blend accumulated metadata values
        """
        self.blend_function = blend_function
        self.values = []

    def accumulate(self, value):
        """
        Add a metadata value for blending.

        Parameters
        ----------
        value : any
            Value for this metadata attribute to use when blending.
        """
        self.values.append(value)

    def finalize(self):
        """
        Blend the accumulated metadata values.

        Returns
        -------
        value :
            The blended result.
        """
        if not self.values:
            return None
        return self.blend_function(self.values)


def make_blender(rule):
    """
    Make an `AttributeBlender` instance using the provided rule.

    Parameters
    ----------
    rule : str
        Name of the blending rule. Must be in `RULE_FUNCTIONS`.

    Returns
    -------
    attr_blender : `AttrBlender`
        Blender instance using the provided rule.
    """
    return AttributeBlender(RULE_FUNCTIONS[rule])

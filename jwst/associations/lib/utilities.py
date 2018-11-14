"""General Utilities"""

from ast import literal_eval

from numpy.ma import masked


def evaluate(value):
    """Evaluate a value

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


def getattr_from_list(adict, attributes, invalid_values=None):
    """Retrieve value from dict using a list of attributes

    Parameters
    ----------
    adict : dict
        dict to retrieve from

    attributes : list
        List of attributes

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
        raise KeyError('Object has no attributes in {}'.format(attributes))


def is_iterable(obj):
    """General iterator check"""
    return not isinstance(obj, str) and \
        not isinstance(obj, tuple) and \
        hasattr(obj, '__iter__')

"""Association Candidate Identifier"""
from ast import literal_eval

__all__ = [
    'ACID'
]


class ACID(object):
    """Association Candidate Identifer

    Parameters
    ----------
    in: str or 2-tuple
        The string representation or 2-tuple containing the
        candidate ID and TYPE. The string should be itself
        the 2-tuple representation when evaluated. The 2-tuple
        will have the form:
            (id, type)

    Attributes
    ----------
    id: str
        The id number.

    type: str
        The type of candidate. Some, but not all
        possibilities include 'OBSERVATION',
        'MOSAIC', 'DISCOVERED'

    __str__: str
        The DMS-specified string representation of
        a candidate identifier. Possibilities include
        'oXXX' for OBSERVATION, 'c1XXX' for 'MOSAIC' or
        other PPS-defined candidates, and 'a3XXX' for
        'DISCOVERED' associations.
    """
    def __init__(self, input):
        try:
            self.id, self.type = literal_eval(input)
        except (ValueError, SyntaxError):
            self.id, self.type = input

    def __str__(self):
        return self.id

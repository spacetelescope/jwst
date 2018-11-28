"""Association Candidate Identifier"""
from ast import literal_eval
import re

from .counter import Counter

__all__ = [
    'ACID'
]

# Start of the discovered association ids.
_DISCOVERED_ID_START = 3001


class ACID():
    """Association Candidate Identifer

    Parameters
    ----------
    input : str or 2-tuple
        The string representation or 2-tuple containing the
        candidate ID and TYPE. The string should be itself
        the 2-tuple representation when evaluated. The 2-tuple
        will have the form:
            (id, type)

    Attributes
    ----------
    id : str
        The id number.

    type : str
        The type of candidate. Some, but not all
        possibilities include 'OBSERVATION',
        'MOSAIC', 'DISCOVERED'

    __str__ : str
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


class ACIDMixin():
    """Enable ACID for rules"""
    def __init__(self, *args, **kwargs):

        # Initialize discovered association ID
        self.discovered_id = Counter(_DISCOVERED_ID_START)

        super(ACIDMixin, self).__init__(*args, **kwargs)

    def acid_from_constraints(self):
        """Determine ACID from constraints"""
        for constraint in self.constraints:
            if getattr(constraint, 'is_acid', False):
                value = re.sub(
                    '\\\\',
                    '',
                    '-'.join(constraint.found_values)
                )
                try:
                    acid = ACID(value)
                except ValueError:
                    pass
                else:
                    break
        else:
            id = 'a{:0>3}'.format(self.discovered_id.value)
            acid = ACID((id, 'DISCOVERED'))

        return acid

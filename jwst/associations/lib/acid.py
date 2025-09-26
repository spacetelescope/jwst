"""Association Candidate Identifier."""

import re
from ast import literal_eval

from jwst.associations.lib.counter import Counter

__all__ = ["ACID", "ACIDMixin"]

# Start of the discovered association ids.
_DISCOVERED_ID_START = 3001


class ACID:
    """
    Association Candidate Identifier.

    Parameters
    ----------
    inits : str or tuple of str
        String representation or 2-tuple containing the candidate ID
        and TYPE. The string should be itself the 2-tuple representation
        when evaluated. The 2-tuple will have the form ``(id, type)``.

    Attributes
    ----------
    id : str
        The ID number.

    type : str
        The type of candidate. Some, but not all
        possibilities include 'OBSERVATION',
        'MOSAIC', 'DISCOVERED'

    Notes
    -----
    ``str(obj)`` returns the DMS-specified string representation of
    a candidate identifier. Possibilities include:

    * 'oXXX' for OBSERVATION
    * 'c1XXX' for 'MOSAIC'
    * other PPS-defined candidates
    * 'a3XXX' for 'DISCOVERED' associations
    """

    def __init__(self, inits):
        try:
            self.id, self.type = literal_eval(inits)
        except (ValueError, SyntaxError):
            self.id, self.type = inits

    def __str__(self):
        return self.id


class ACIDMixin:
    """Enable ACID for rules."""

    def __init__(self, *args, **kwargs):
        # Initialize discovered association ID
        self.discovered_id = Counter(_DISCOVERED_ID_START)

        super(ACIDMixin, self).__init__(*args, **kwargs)

    def acid_from_constraints(self):
        """
        Determine ACID from constraints.

        Returns
        -------
        acid : `~jwst.associations.lib.acid.ACID`
            Association candidate identifier from constraints.
        """
        for constraint in self.constraints:
            if getattr(constraint, "is_acid", False):
                value = re.sub("\\\\", "", "-".join(constraint.found_values))
                try:
                    acid = ACID(value)
                except ValueError:
                    pass
                else:
                    break
        else:
            cand_id = f"a{self.discovered_id.value:0>3}"
            acid = ACID((cand_id, "DISCOVERED"))

        return acid

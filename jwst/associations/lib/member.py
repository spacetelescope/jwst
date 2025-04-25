"""Association Member definition - store exposure metadata as dict-like objects."""

from collections import UserDict
from copy import copy


class Member(UserDict):
    """
    Member of an association.

    Attributes
    ----------
    item : object
        The original item that created this member.
    """

    def __init__(self, initialdata=None, item=None):
        """
        Initialize a Member.

        Parameters
        ----------
        initialdata : Dict-like or Member
            Initialization data. Any type of initialization that
            `collections.UserDict` allows or `Member` itself.

        item : object
            The item to initialize with. This will override
            any `Member.item` given in `initialdata`. Most common
            object type is a ~jwst.associations.pool.PoolRow instance.
        """
        self.item = None

        if isinstance(initialdata, Member):
            self.data = copy(initialdata.data)
            self.item = copy(initialdata.item)
        else:
            super(Member, self).__init__(initialdata)

        if item is not None:
            self.item = copy(item)

    def __eq__(self, other):
        """
        Compare members.

        If both Members have attributes `expname` and `exptype`,
        compare only those attributes. Otherwise, use the default
        comparison.

        Parameters
        ----------
        other : object
            The comparison object.

        Returns
        -------
        bool
            True if deemed equal/equivalent.
        """
        hasexpkeys = all(k in data for k in ("expname", "exptype") for data in (self, other))
        if hasexpkeys:
            return all(self[k] == other[k] for k in ("expname", "exptype"))
        else:
            return super().__eq__(other)

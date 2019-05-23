"""Association Member"""
from collections import UserDict
from copy import copy


class Member(UserDict):
    """Member of an association

    Parameters
    ----------
    initialdata: Dict-like or Member
        Initialization data. Any type of initialization that
        `collections.UserDict` allows or `Member` itself.

    item: obj
        The item to initialize with. This will override
        any `Member.item` given in `initialdata`.

    Attributes
    ----------
    item: obj
        The original item that created this member.
    """
    def __init__(self, initialdata=None, item=None):
        self.item = None

        if isinstance(initialdata, Member):
            self.data = copy(initialdata.data)
            self.item = copy(initialdata.item)
        else:
            super(Member, self).__init__(initialdata)

        if item is not None:
            self.item = copy(item)

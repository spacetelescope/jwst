"""Association Member"""
from collections import UserDict
from copy import deepcopy


class Member(UserDict):
    """Member of an association

    Parameters
    ----------
    initialdata: Dict-like or Member
        Initialization data. Any type of initialization that
        `collections.UserDict` allows or `Member` itself.
        Either way, copy of the initialization data is done.

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
            self.data = deepcopy(initialdata.data)
            self.item = deepcopy(initialdata.item)
        else:
            super(Member, self).__init__(initialdata)
            
        if item is not None:
            self.item = deepcopy(item)

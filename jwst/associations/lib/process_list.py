"""Reprocessing List"""
from collections import (defaultdict, deque)
from functools import reduce

__all__ = [
    'ProcessList',
    'ProcessItem',
    'ProcessQueue',
    'ProcessQueueSorted'
]


class ProcessItem:
    """Items to be processed

    Parameters
    ----------
    obj : object
        The object to make a `ProcessItem`.
        Objects must be equatable.
    """
    def __init__(self, obj):
        self.obj = obj

    @classmethod
    def to_process_items(cls, iterable):
        """Iterable to convert a list to ProcessItem's

        Parameters
        ----------
        iterable : iterable
            A source of objects to convert

        Returns
        -------
        An iterable where the object has been
        converted to a `ProcessItem`
        """
        for obj in iterable:
            yield cls(obj)

    def __hash__(self):
        try:
            hash_value = self.obj.__hash__()
        except (AttributeError, TypeError):
            hash_value = hash(repr(self))
        return hash_value

    def __eq__(self, other):
        try:
            equality = self.obj == other.obj
        except AttributeError:
            equality = self.__hash__() == other.__hash__()
        return equality


class ProcessList:
    """A Process list

    Parameters
    ----------
    items : [item[, ...]]
        The list of items to process

    rules : [Association[, ...]]
        List of rules to process the items against.

    work_over : int
        What the reprocessing should work on:
        - `ProcessList.EXISTING`:   Only existing associations
        - `ProcessList.RULES`:      Only on the rules to create new associations
        - `ProcessList.BOTH`:       Compare to both existing and rules
        - `ProcessList.NONSCIENCE`: Only on non-science items

    only_on_match : bool
        Only use this object if the overall condition
        is True.
    """

    (
        RULES,
        BOTH,
        EXISTING,
        NONSCIENCE,
    ) = range(0, 4)
    """
    Categories of different sets of associations or items
    to process
    """

    _str_attrs = ('rules', 'work_over', 'only_on_match')

    def __init__(self,
                 items=None,
                 rules=None,
                 work_over=BOTH,
                 only_on_match=False):
        self.items = items
        self.rules = rules
        self.work_over = work_over
        self.only_on_match = only_on_match

    def __str__(self):
        result = '{}(n_items: {}, {})'.format(
            self.__class__.__name__,
            len(self.items),
            {
                str_attr: getattr(self, str_attr)
                for str_attr in self._str_attrs
            }
        )
        return result


class ProcessQueue(deque):
    """Make a deque iterable and mutable"""
    def __iter__(self):
        while True:
            try:
                yield self.popleft()
            except:
                break


class ProcessQueueNoDups:
    """First-In-First-Out queue

    Checks on whether the objects are already in
    the queue

    Parameters
    ----------
    init : [obj[,...]] or None
        List of objects to put on the queue
    """
    def __init__(self, init=None):

        self._members = set()
        self._queue = deque()

        if init is not None:
            self.extend(init)

    def extend(self, iterable):
        """Add objects if not already in the queue"""
        for obj in iterable:
            self.append(obj)

    def append(self, obj):
        """Add object if not already in the queue"""
        if obj not in self._members:
            self._queue.append(obj)
            self._members.add(obj)

    def popleft(self):
        """Pop the first-in object"""
        obj = self._queue.popleft()
        self._members.remove(obj)
        return obj

    def __len__(self):
        return len(self._queue)

    def __iter__(self):
        while True:
            try:
                yield self.popleft()
            except:
                break


class ProcessQueueSorted:
    """Sort ProcessItem based on work_over

    `ProcessList`s are handled in order of `RULES`, `BOTH`,
    `EXISTING`, and `NONSCIENCE`.

    Parameters
    ----------
    init : [ProcessList[,...]]
        List of `ProcessList` to start the queue with.
    """
    def __init__(self, init=None):
        self.queues = [
            ProcessQueue(),
            ProcessQueue(),
            ProcessQueue(),
            ProcessQueue(),
        ]

        if init is not None:
            self.extend(init)

    def extend(self, process_lists):
        """Add the list of process items to their appropriate queues"""
        for process_list in process_lists:
            self.queues[process_list.work_over].append(process_list)

    def __len__(self):
        return reduce(lambda x, y: x + len(y), self.queues, 0)

    def __iter__(self):
        """Return the queues in order"""
        while len(self) > 0:
            for queue in self.queues:
                for process_list in queue:
                    yield process_list
                    break
                else:
                    continue
                break

"""Reprocessing List"""
from collections import deque
from enum import Enum
from functools import reduce

class ListCategory(Enum):
    RULES      = 0
    BOTH       = 1
    EXISTING   = 2
    NONSCIENCE = 3

__all__ = [
    'ListCategory',
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

    _str_attrs = ('rules', 'work_over', 'only_on_match')

    def __init__(self,
                 items=None,
                 rules=None,
                 work_over=ListCategory.BOTH,
                 only_on_match=False):
        self.items = items
        self.rules = rules
        self.work_over = work_over
        self.only_on_match = only_on_match

    @property
    def hash(self):
        """Create a unique hash"""
        return (tuple(self.rules), self.work_over, self.only_on_match)

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


class ProcessListQueue:
    """First-In-First-Out queue of ProcessLists

    Parameters
    ----------
    init : [ProcessList[,...]] or None
        List of ProcessLists to put on the queue.
    """
    def __init__(self, init=None):
        self._queue = dict()
        if init is not None:
            self.extend(init)

    def append(self, process_list):
        """Add object if not already in the queue"""
        plhash = process_list.hash
        if plhash not in self._queue:
            self._queue[plhash] = process_list
        else:
            self._queue[plhash].items += process_list.items

    def extend(self, iterable):
        """Add objects if not already in the queue"""
        for process_list in iterable:
            self.append(process_list)

    def items(self):
        """Return list generator of all items"""
        for plhash in self._queue:
            for item in self._queue[plhash].items:
                yield item

    def popleft(self):
        """Pop the first-in object"""
        plhash = next(iter(self._queue))
        process_list = self._queue[plhash]
        del self._queue[plhash]
        return process_list

    def __len__(self):
        return len(self._queue)

    def __iter__(self):
        while True:
            try:
                yield self.popleft()
            except:
                break

    def __str__(self):
        result = f'{self.__class__.__name__}: rulesets {len(self)} items {len(list(self.items()))}'
        return result


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
        self.queues = {
            list_category: ProcessListQueue()
            for list_category in ListCategory
        }

        if init is not None:
            self.extend(init)

    def extend(self, process_lists):
        """Add the list of process items to their appropriate queues"""
        for process_list in process_lists:
            self.queues[process_list.work_over].append(process_list)

    def __iter__(self):
        """Return the queues in order"""
        while len(self) > 0:
            for category in ListCategory:
                for process_list in self.queues[category]:
                    yield process_list
                    break
                else:
                    continue
                break

    def __len__(self):
        return reduce(lambda x, y: x + len(y), self.queues.values(), 0)

    def __str__(self):
        result = f'{self.__class__.__name__}:'
        for queue in self.queues:
            result += f'\n\tQueue {queue}: {self.queues[queue]}'
        return result

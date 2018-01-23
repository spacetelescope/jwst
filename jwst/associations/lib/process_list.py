"""Reprocessing List"""
from collections import (defaultdict, deque)
from functools import reduce

__all__ = [
    'ProcessList',
    'ProcessQueue',
    'ProcessQueueSorted'
]


class ProcessList:
    """A Process list

    Parameters
    ----------
    items: [item[, ...]]
        The list of items to process

    rules: [Association[, ...]]
        List of rules to process the items against.

    work_over: int
        What the reprocessing should work on:
        - `ProcessList.EXISTING`: Only existing associations
        - `ProcessList.RULES`: Only on the rules to create new associations
        - `ProcessList.BOTH`: Compare to both existing and rules

    only_on_match: bool
        Only use this object if the overall condition
        is True.
    """

    (
        RULES,
        BOTH,
        EXISTING,
    ) = range(0, 3)

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


class ProcessQueueSorted:
    """Sort ProcessItem based on work_over

    `ProcessItem`s are handled in order of `RULES`, `BOTH`, and
    `EXISTING`.

    Parameters
    ----------
    init: [ProcessList[,...]]
        List of `ProcessList` to start the queue with.
    """
    def __init__(self, init=None):
        self.queues = [
            ProcessQueue(),
            ProcessQueue(),
            ProcessQueue()
        ]

        if init is not None:
            self.extend(init)

    def extend(self, process_items):
        """Add the list of process items to their appropriate queues"""
        for item in process_items:
            self.queues[item.work_over].append(item)

    def __iter__(self):
        """Return the queues in order"""
        while reduce(lambda x, y: x + len(y), self.queues, 0) > 0:
            for queue in self.queues:
                for item in queue:
                    yield item
                    break
                else:
                    continue
                break

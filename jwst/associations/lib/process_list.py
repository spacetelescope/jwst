"""
Reprocessing Lists and Queues.

This modules defines what process lists are and queues of process lists.

A process list, `ProcessList`, is a list of (items, rules) and meta information
, most notably `work_over`. `work_over` is one of the values of `ListCategory`.
A `ListCategory` defines which stage of association processing the list is
relevant to. In other words, the order, or priority, of when a list should be processed
is defined by its `ListCategory`. The priority is the value of each `ListCategory`,
starting with zero.

ProcessLists are primarily put into queues for processing. There are two
queues for handling ProcessLists. `ProcessListQueue` is a basic
First-In-First-Out (FIFO) queue that can be used as a generator.

The second queue is `ProcessQueueSorted`, which returns ProcessLists according to
their priority as defined by each ProcessList's `work_over`. An important aspect of
ProcessQueueSorted is that it is mutable: New ProcessLists can be added to the queue
while iterating over the queue.
"""

from collections import deque
from enum import Enum
from functools import reduce


__all__ = ["ListCategory", "ProcessList", "ProcessItem", "ProcessQueue", "ProcessQueueSorted"]


class ListCategory(Enum):
    """The work_over categories for ProcessLists."""

    RULES = 0  # Operate over rules only
    BOTH = 1  # Operate over both rules and existing associations
    EXISTING = 2  # Operate over existing associations only
    NONSCIENCE = 3  # Items that are not science specific that should be applied to only
    # existing associations


class ProcessItem:
    """
    Items to be processed.

    Create hashable objects from a list of arbitrary objects.

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
        """
        Convert a list to ProcessItems.

        Parameters
        ----------
        iterable : iterable
            A source of objects to convert

        Returns
        -------
        iterable
            An iterable where the object has been converted to a `ProcessItem`.
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
    """A Process list."""

    _str_attrs = ("rules", "work_over", "only_on_match", "trigger_constraints", "trigger_rules")

    def __init__(
        self,
        items=None,
        rules=None,
        work_over=ListCategory.BOTH,
        only_on_match=False,
        trigger_constraints=None,
        trigger_rules=None,
    ):
        """
        Initialize a ProcessList.

        Parameters
        ----------
        items : [item[, ...]]
            The list of items to process

        rules : [Association[, ...]]
            List of rules to process the items against.

        work_over : int
            What the reprocessing should work on:
            - `ProcessList.RULES`:      Only on the rules to create new associations
            - `ProcessList.EXISTING`:   Only existing associations
            - `ProcessList.BOTH`:       Compare to both existing and rules
            - `ProcessList.NONSCIENCE`: Only on non-science items

        only_on_match : bool
            Only use this object if the overall condition
            is True.

        trigger_constraints : [Constraint[,...]]
            The constraints that created the ProcessList

        trigger_rules : [Association[,...]]
            The association rules that created the ProcessList
        """
        self.items = items
        self.rules = rules
        self.work_over = work_over
        self.only_on_match = only_on_match
        self.trigger_constraints = set(trigger_constraints) if trigger_constraints else set()
        self.trigger_rules = set(trigger_rules) if trigger_rules else set()

    @property
    def hash(self):
        """
        Create a unique hash.

        Returns
        -------
        Tuple(Rule, ...), int, bool
            Tuple of: tuple of rule objects, integer and bool. Used as a unique hash.
        """
        return (tuple(self.rules), self.work_over, self.only_on_match)

    def update(self, process_list, full=False):
        """
        Update with information from ProcessList.

        Attributes from `process_list` are added to self's attributes. If `not
        full`, the attributes `rules`, 'work_over`, and `only_on_match` are not
        taken.

        Note that if `full`, destructive action will occur with respect to
        `work_over` and `only_on_match`.

        Parameters
        ----------
        process_list : ProcessList
            The source process list to absorb.

        full : bool
            Include the hash attributes `rules`, `work_over`, and `only_on_match`.
        """
        self.items += process_list.items
        self.trigger_constraints.update(process_list.trigger_constraints)
        self.trigger_rules.update(process_list.trigger_rules)
        if full:
            self.rules += process_list.rules
            self.work_over = process_list.work_over
            self.only_on_match = process_list.only_on_match

    def __str__(self):
        result = (
            f"{self.__class__.__name__}(n_items: {len(self.items)}, "
            f"{ ({str_attr: getattr(self, str_attr) for str_attr in self._str_attrs}) })"
        )
        return result


class ProcessQueue(deque):
    """Make a deque iterable and mutable."""

    def __iter__(self):
        while True:
            try:
                yield self.popleft()
            except Exception:
                break


class ProcessListQueue:
    """
    First-In-First-Out queue of ProcessLists.

    ProcessLists can be added either individually using `append` method, or
    a list of ProcessLists can be added through object initialization or
    the `extend` method.

    There are two generators implement. The first is the ProcessListQueue
    object itself. When the object is used as a generator, the generator will
    return the earliest ProcessList added to the queue (FIFO), popping that
    ProcessList from the queue, hence draining the queue.

    The second generator is returned by the `items` method. This method will
    return all the items from all the ProcessLists in the queue,
    non-destructively. The ProcessLists are accessed in their order in the
    queue, and then each item is retrieved from their ProcessList in the list
    order of the ProcessList.

    A final feature of ProcessListQueue is that it is mutable: New items can
    be added to the queue while items are being popped from the queue.

    Notes
    -----
    The FIFO operations depends on the fact that, inherently,
    `dict` preserves order in which key/value pairs are added to the
    dictionary.
    """

    def __init__(self, init=None):
        """
        Initialize a ProcessListQueue.

        Parameters
        ----------
        init : [ProcessList[,...]] or None
            List of ProcessLists to put on the queue.
        """
        self._queue = {}
        if init is not None:
            self.extend(init)

    def append(self, process_list):
        """Add ProcessList to queue, if not already in the queue."""
        plhash = process_list.hash
        if plhash not in self._queue:
            self._queue[plhash] = process_list
        else:
            self._queue[plhash].update(process_list)

    def extend(self, iterable):
        """Add lists of ProcessLists if not already in the queue."""
        for process_list in iterable:
            self.append(process_list)

    def items(self):
        """Return list generator of all items."""
        for plhash in self._queue:
            yield from self._queue[plhash].items

    def popleft(self):
        """
        Pop the first-in object.

        Returns
        -------
        [ProcessList, ...]
            The queue of ProcessList objects with the first one popped.
        """
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
            except Exception:
                break

    def __str__(self):
        result = f"{self.__class__.__name__}: rulesets {len(self)} items {len(list(self.items()))}"
        return result


class ProcessQueueSorted:
    """
    Sort ProcessItem based on work_over.

    Create a generator that implements a First-In-First-Out (FIFO) queue, with the one
    modification that the queues are handled in order of their `work_over` priority.
    For example, even if a ProcessList with work_over of ListCategory.EXISTING had
    been added to the queue before a ProcessList with work_over of ListCategory.RULES,
    the second ProcessList will be returned before the first.

    ProcessQueueSorted is also mutable: ProcessLists can be added to the queue
    while the lists are being popped from the queue. When doing so, it is
    important to remember that the order of return, as described above, still
    pertains. For example, if the queue only has ProcessLists of work_over
    ListCategory.EXISTING, and a new ProcessList of work_over
    ListCategory.RULES is added during iteration, the next list returned will
    be the RULES one, because the RULES lists have priority over EXISTING
    lists, regardless of when the list was added.
    """

    def __init__(self, init=None):
        """
        Initialize a ProcessQueueSorted.

        Parameters
        ----------
        init : [ProcessList[,...]]
            List of `ProcessList` to start the queue with.
        """
        self.queues = {list_category: ProcessListQueue() for list_category in ListCategory}

        if init is not None:
            self.extend(init)

    def extend(self, process_lists):
        """Add the list of process items to their appropriate queues."""
        for process_list in process_lists:
            self.queues[process_list.work_over].append(process_list)

    def __iter__(self):
        """Return the queues in order."""
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
        result = f"{self.__class__.__name__}:"
        for queue in self.queues:
            result += f"\n\tQueue {queue}: {self.queues[queue]}"
        return result


def workover_filter(process_list, work_over):
    """
    Determine and modify workover of input process list.

    Parameters
    ----------
    process_list : ProcessList
        The process list under consideration

    work_over : ListCategory
        The ListCategory to compare against.

    Returns
    -------
    process_list : ProcessList or None
        The input process_list with work_over modified.
        None if the process list should not be continued.
    """
    result = process_list
    if process_list.work_over in [ListCategory.RULES, ListCategory.BOTH]:
        if work_over in [ListCategory.RULES, ListCategory.BOTH]:
            result.work_over = ListCategory.BOTH
        else:
            result = None
    elif work_over not in [ListCategory.RULES, ListCategory.BOTH]:
        result = None
    return result

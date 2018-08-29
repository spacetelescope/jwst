"""Test ProcessList, ProcessQueue, ProcessQueueSorted"""

from .helpers import (
    combine_pools,
    t_path
)

from ..lib.process_list import *


def test_item():
    pool = combine_pools(t_path('data/pool_013_coron_nircam.csv'))
    item1 = ProcessItem(pool[0])
    item2 = ProcessItem(pool[1])
    assert item1 == item1
    assert item1 != item2
    s = set([item1, item2])
    assert len(s) == 2


def test_item_iterable():
    pool = combine_pools(t_path('data/pool_013_coron_nircam.csv'))
    process_items = ProcessItem.to_process_items(pool)
    for process_item in process_items:
        assert isinstance(process_item, ProcessItem)


def test_process_queue():
    queue = ProcessQueue()
    items_to_add = [
        [1],
        [2, 3, 4],
        # [3, 4],  # Neither should get added but not implemented
        [5]
    ]
    standard = [1, 2, 3, 4, 5]

    idx = 0
    queue.extend(items_to_add[0])
    results = []
    for item in queue:
        results.append(item)
        idx += 1
        try:
            to_add = items_to_add[idx]
        except:
            pass
        else:
            queue.extend(to_add)

    assert results == standard


def test_process_queue_sorted():
    """Test the sorted process queue"""
    items_to_add = [
        [
            ProcessList([1], work_over=ProcessList.RULES),
            ProcessList([2], work_over=ProcessList.BOTH),
            ProcessList([3], work_over=ProcessList.EXISTING),
        ],
        [
            ProcessList([4], work_over=ProcessList.EXISTING),
        ],
        [
            ProcessList([5], work_over=ProcessList.BOTH),
        ],
        [
            ProcessList([6], work_over=ProcessList.RULES),
        ],
    ]

    # This is the order the items should be retrieved in:
    # RULES -> BOTH -> EXISTING
    standard = [1, 2, 5, 6, 3, 4]

    queue = ProcessQueueSorted()
    idx = 0
    queue.extend(items_to_add[idx])
    assert len(queue) == 3
    results = []
    for item in queue:
        results.append(item.items[0])
        idx += 1
        try:
            to_add = items_to_add[idx]
        except:
            pass
        else:
            queue.extend(to_add)

    assert results == standard

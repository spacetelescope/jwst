"""Test ProcessList, ProcessQueue, ProcessQueueSorted"""

from ..lib.process_list import *


def test_process_queue():
    queue = ProcessQueue()
    items_to_add = [
        [1],
        [2, 3, 4],
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

from asdf.tags.core import HistoryEntry

def _iterable(values):
    if isinstance(values, str) or not hasattr(values, '__iter__'):
        values = (values,)
    return values

class HistoryList:
    """
    A list that coerces a new value into a HistoryEntry.
    Only a subset of the list interface is implemented.
    """
    def __init__(self, asdf):
        self._context = asdf
        if len(self._context.get_history_entries()):
            self._entries = self._context.get_history_entries()
        else:
            self._context.add_history_entry("fake entry")
            self._entries = self._context.get_history_entries()
            self._entries.clear()

    def __len__(self):
        return len(self._entries)

    def __getitem__(self, key):
        return self._entries[key]

    def __setitem__(self, key, value):
        self.append(value)
        value = self._entries.pop()
        self._entries[key] = value

    def __delitem__(self, key):
        del self._entries[key]

    def __iter__(self):
        return iter(self._entries)

    def __repr__(self):
        return repr(self._entries)

    def __str__(self):
        return str(self._entries)

    def __eq__(self, other):
        if isinstance(other, HistoryList):
            other = other._entries
        else:
            other = _iterable(other)

        if len(self) != len(other):
            return False

        for self_entry, other_entry in zip(self._entries, other):
            if isinstance(other_entry, str):
                if self_entry.get('description') != other_entry:
                    return False
            elif isinstance(other_entry, dict):
                for key in other_entry.keys():
                    if self_entry.get(key) != other_entry.get(key):
                        return False
        return True

    def append(self, value):
        if isinstance(value, HistoryEntry):
            self._entries.append(value)
        else:
            self._context.add_history_entry(value)

    def clear(self):
        self._entries.clear()

    def extend(self, values):
        values = _iterable(values)
        for value in values:
            self.append(value)

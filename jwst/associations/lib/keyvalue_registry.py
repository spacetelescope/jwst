"""Key/Value Registry"""

import pytest

try:
    from collections import UserDict
except ImportError:
    from UserDict import IterableUserDict

    class UserDict(IterableUserDict, object):
        pass

__all__ = [
    'KeyValueRegistry',
    'KeyValueRegistryError',
    'KeyValueRegistryNoKeyFound',
    'KeyValueRegistryNotSingleItemError'
]


class KeyValueRegistry(UserDict):
    """Provide a dict-like registry

    Differences from just a `dict`:
        - Can be given single item or a 2-tuple.
          If an item, attempts to read the `__name__` attribute
          and use that as the key.

        - If None is given as a key, a default key can
          be specified.

        - Instances can be used as decorators.

    Parameters
    ----------
    items : object or (str, object) or dict
        Initializing items.

    default : str or object
        The default to use when key is `None`
    """

    def __init__(self, items=None, default=None):
        super_args = ()
        if items is not None:
            super_args = (make_dict(items), )
        super(KeyValueRegistry, self).__init__(*super_args)

        self.default = None
        if default is not None:
            default_dict = make_dict(default)
            if len(default_dict) > 1:
                raise KeyValueRegistryNotSingleItemError
            default_dict = make_dict(default)
            self.update(default_dict)
            self.default = next(iter(default_dict.keys()))
            self.update({None: default_dict[self.default]})

    def update(self, item):
        """Add item to registry"""
        item_dict = make_dict(item)
        super(KeyValueRegistry, self).update(item_dict)

    def __call__(self, item):
        """Add item by calling instance

        This allows an instance to be used as a decorator.
        """
        self.update(item)
        return item


# ******
# Errors
# ******
class KeyValueRegistryError(Exception):
    def __init__(self, *args):
        if len(args) == 0:
            args = (self.msg, )
        super(KeyValueRegistryError, self).__init__(*args)


class KeyValueRegistryNotSingleItemError(KeyValueRegistryError):
    msg = 'Item cannot be a list'


class KeyValueRegistryNoKeyFound(KeyValueRegistryError):
    msg = 'Cannot deduce key from given value'


# *********
# Utilities
# *********
def make_dict(item):
    """Create a dict from an item

    Parameters
    ----------
    item : object or (name, object) or dict
        If dict, just return dict.
        If 2-tuple, return dict with the key/value pair
        If just object, use `__name__` as key
    """
    try:
        item_dict = dict(item)
    except (TypeError, ValueError):
        try:
            key, value = item
        except (TypeError, ValueError):
            try:
                key = item.__name__
            except (AttributeError, SyntaxError):
                raise KeyValueRegistryNoKeyFound
            else:
                value = item

        # At they point we have a key/value pair
        item_dict = {key: value}

    return item_dict


# **********
# Test suite
# **********
def test_dict_like():
    """Test the basic to ensure similar to a dict"""
    data = {'a': 1, 'b': 2}
    kvr = KeyValueRegistry(data)
    assert data.keys() == kvr.keys()
    assert set(data.values()) == set(kvr.values())

    assert kvr.get('a') == 1
    assert kvr.get('c', 3) == 3

    keys, values = zip(*kvr.items())
    assert set(data.keys()) == set(keys)
    assert set(data.values()) == set(values)

    kvr_copy = kvr.copy()
    assert set(kvr_copy) == set(kvr)

    assert kvr.pop('a') == 1
    assert kvr.pop('a', 3) == 3
    assert kvr.popitem() == ('b', 2)
    with pytest.raises(KeyError):
        kvr.popitem()

    kvr = KeyValueRegistry()
    kvr.update(data)
    assert data.keys() == kvr.keys()
    assert set(data.values()) == set(kvr.values())

    kvr.clear()
    assert len(kvr) == 0


def test_default():
    kvr = KeyValueRegistry(default={'a': 1})
    assert kvr.default == 'a'
    assert kvr[None] == 1


def test_tuple():
    data = ('a', 1)
    kvr = KeyValueRegistry(data)
    assert kvr['a'] == 1

    kvr = KeyValueRegistry(default=data)
    assert kvr['a'] == 1
    assert kvr[None] == 1


def test_fn():
    def fn():
        return 1

    kvr = KeyValueRegistry(fn)
    assert kvr[fn.__name__] == fn

    kvr = KeyValueRegistry(default=fn)
    assert kvr[fn.__name__] == fn
    assert kvr[None] == fn


def test_decorator():
    kvr = KeyValueRegistry()

    @kvr
    def fn():
        return 1

    assert kvr[fn.__name__] == fn
    assert kvr[fn.__name__]() == fn()

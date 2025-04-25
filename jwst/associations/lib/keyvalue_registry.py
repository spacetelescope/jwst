"""Key/Value Registry."""

from collections import UserDict


__all__ = [
    "KeyValueRegistry",
    "KeyValueRegistryError",
    "KeyValueRegistryNoKeyFoundError",
    "KeyValueRegistryNotSingleItemError",
]


class KeyValueRegistry(UserDict):
    """
    Provide a dict-like registry.

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
            super_args = (make_dict(items),)
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
        """Add item to registry."""
        item_dict = make_dict(item)
        super(KeyValueRegistry, self).update(item_dict)

    def __call__(self, item):
        """
        Add item by calling instance.

        This allows an instance to be used as a decorator.

        Parameters
        ----------
        item : object or (str, object) or dict
            Item used for decoration.

        Returns
        -------
        item : object or (str, object) or dict
            The item used to update self.
        """
        self.update(item)
        return item


# Errors
class KeyValueRegistryError(Exception):
    """Exception class for key value in registry."""

    def __init__(self, *args):
        if len(args) == 0:
            args = (self.msg,)
        super(KeyValueRegistryError, self).__init__(*args)


class KeyValueRegistryNotSingleItemError(KeyValueRegistryError):
    """Exception class for when passed item is a list and not a single item."""

    msg = "Item cannot be a list"


class KeyValueRegistryNoKeyFoundError(KeyValueRegistryError):
    """Exception class for no key in registry."""

    msg = "Cannot deduce key from given value"


# Utilities
def make_dict(item):
    """
    Create a dict from an item.

    Items may be a dict, a 2-tuple or an object. Objects are most
    often a file format class from ~jwst.associations.association_io - `json` or `yaml`.

    Parameters
    ----------
    item : object or (name, object) or dict
        If dict, just return dict.
        If 2-tuple, return dict with the key/value pair
        If just object, use `__name__` as key

    Returns
    -------
    dict
        The dictionary created from the item.
    """
    try:
        item_dict = dict(item)
    except (TypeError, ValueError):
        try:
            key, value = item
        except (TypeError, ValueError):
            try:
                key = item.__name__
            except (AttributeError, SyntaxError) as err:
                raise KeyValueRegistryNoKeyFoundError from err
            else:
                value = item

        # At they point we have a key/value pair
        item_dict = {key: value}

    return item_dict

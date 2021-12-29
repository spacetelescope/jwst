"""Update path of members in an association"""
from os.path import (join, basename)


def update_path(asn, file_path, target='expname'):
    """Update path of members in an association

    Parameters
    ----------
    asn : Association
        An association. The association is modified in-place.

    file_path : str
        New path to prepend to each member

    target : str
        Key to replace
    """
    update_key_value(asn, target, (file_path, ), mod_func=_replace_path)


def update_key_value(obj, target, func_args, mod_func=None):
    """Update all instances of a key using a modifier

    Parameters
    ----------
    obj : object
        Hashable object to modify

    target : str
        The target key to modify

    func_args : (arg(, ...))
        Arguments to pass to the modification function

    mod_func : function
        Function to modify the target key with. If `None`,
        the key will be replaced by the arg list.

    Notes
    -----
    The first argument to `mode_func` will always be the value of the
    target key. Any other arguments given will then be passed to the
    the function.
    """
    if mod_func is None:
        mod_func = lambda value, args: args

    if hasattr(obj, 'items'):
        for key, value in obj.items():
            if key == target:
                obj[key] = mod_func(value, *func_args)
            if isinstance(value, dict):
                update_key_value(value, target, func_args, mod_func=mod_func)
            elif isinstance(value, list):
                for item in value:
                    update_key_value(item, target, func_args, mod_func=mod_func)


def _replace_path(old_path, new_path):
    """Replace the path prefix of the basename

    Parameters
    ----------
    old_path : str
        A string with a file path.

    new_path : str
        The new path to prepend to the basename.

    Returns
    -------
    new_full_path : str
        The basename of the path with the
    """
    file_name = basename(old_path)
    new_full_path = join(new_path, file_name)
    return new_full_path


# #####
# Tests
# #####
_test_obj = {
    'a': 'change',
    'b': 'nochange',
    'c': {
        'a': 'change',
        'b': 'nochange'
    },
    'd': [
        {
            'a': 'change',
            'b': 'nochange'
        },
        {
            'a': 'change',
            'b': 'nochange'
        },
    ]
}


def test_update_key_value_default():
    from copy import deepcopy
    obj = deepcopy(_test_obj)
    target = 'a'
    new_value = 'changed'
    update_key_value(obj, target, (new_value,))
    for value in _gen_dict_extract(target, obj):
        assert value == new_value


def test_update_key_value_mod_func():
    from copy import deepcopy
    obj = deepcopy(_test_obj)
    target = 'a'
    new_value = 'changed'
    mod_func = lambda v, suffix: v + suffix
    update_key_value(obj, target, (new_value,), mod_func=mod_func)
    assert_values = [
        mod_func(value, new_value)
        for value in _gen_dict_extract(target, _test_obj)
    ]
    for value in _gen_dict_extract(target, obj):
        assert assert_values.pop(0) == value


def test_replace_path():
    new_path = join('hello', 'there')
    base = 'base.me'
    assert_value = join(new_path, base)
    assert assert_value == _replace_path(base, new_path)
    assert assert_value == _replace_path(
        join('somepath', base),
        new_path
    )


def _gen_dict_extract(key, var):
    if hasattr(var, 'items'):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in _gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in _gen_dict_extract(key, d):
                        yield result

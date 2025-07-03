from copy import deepcopy
from jwst.associations.lib.update_path import update_key_value, replace_path

_test_obj = {
    "a": "change",
    "b": "nochange",
    "c": {"a": "change", "b": "nochange"},
    "d": [
        {"a": "change", "b": "nochange"},
        {"a": "change", "b": "nochange"},
    ],
}


def _gen_dict_extract(key, var):
    if hasattr(var, "items"):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                yield from _gen_dict_extract(key, v)
            elif isinstance(v, list):
                for d in v:
                    yield from _gen_dict_extract(key, d)


def test_update_key_value_default():
    from copy import deepcopy

    obj = deepcopy(_test_obj)
    target = "a"
    new_value = "changed"
    update_key_value(obj, target, (new_value,))
    for value in _gen_dict_extract(target, obj):
        assert value == new_value


def test_update_key_value_mod_func():
    obj = deepcopy(_test_obj)
    target = "a"
    new_value = "changed"

    def mod_func(v, suffix):
        return v + suffix

    update_key_value(obj, target, (new_value,), mod_func=mod_func)
    assert_values = [mod_func(value, new_value) for value in _gen_dict_extract(target, _test_obj)]
    for value in _gen_dict_extract(target, obj):
        assert assert_values.pop(0) == value


def test_replace_path():
    new_path = "hello/there"
    base = "base.me"
    assert_value = "hello/there/base.me"
    assert assert_value == replace_path(base, new_path)
    assert assert_value == replace_path(("somepath/" + base), new_path)

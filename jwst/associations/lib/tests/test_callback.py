"""Test callback_registry"""
from jwst.associations.lib.callback_registry import CallbackRegistry


def fn_single_1(arg):
    return arg + 1


def fn_single_2(arg):
    return arg + arg


def fn_double_1(arg1, arg2):
    return arg1 + 1, arg2 - 1


def fn_double_2(*args):
    return args[0] + args[0], args[1] + args[1]


def test_basics():
    cbr = CallbackRegistry()
    cbr.add('event1', fn_single_1)
    assert cbr.reduce('event1', 1) == 2


def test_reduce():
    cbr = CallbackRegistry()
    cbr.add('event1', fn_single_1)
    cbr.add('event1', fn_single_2)
    assert cbr.reduce('event1', 1) == 4


def test_double():
    cbr = CallbackRegistry()
    cbr.add('event1', fn_double_1)
    assert cbr.reduce('event1', 1, 2) == (2, 1)

    cbr.add('event2', fn_double_2)
    assert cbr.reduce('event2', 2, 1) == (4, 2)


def test_double_reduce():
    cbr = CallbackRegistry()
    cbr.add('event1', fn_double_1)
    cbr.add('event1', fn_double_2)
    result = cbr.reduce('event1', 1, 2)
    assert result in [(4, 2), (3, 3)]


def test_decorator():
    cbr = CallbackRegistry()

    @cbr('event1')
    def fn_single_1(arg):
        return arg + 1

    @cbr('event1')
    def fn_single_2(arg):
        return arg + arg

    assert cbr.reduce('event1', 1) in (3, 4)

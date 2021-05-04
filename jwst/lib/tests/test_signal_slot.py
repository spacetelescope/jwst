"""Test signal_slot.py"""

from jwst.lib import signal_slot as ss


REFLECT_CALLED = 'reflect: called.'


def reflect(*args):
    """Handler function that simply reflects the input args"""
    print(REFLECT_CALLED)
    return args


def list_append(inlist, new_item):
    """Append new_item to inlist"""
    inlist.append(new_item)
    return inlist, new_item


class AClass:

    def __init__(self):
        self.kwargs = None

    def set_and_reflect(self, *args, **kwargs):
        """Handler function that simply reflects the input args

        Keyword arguments become an attribute
        """
        self.kwargs = kwargs
        return args


def test_basic_structure():
    """Test initialization and basic structure"""
    signal = ss.Signal()

    assert len(signal._slots) == 0


def test_emit(capsys):
    """Test emitting a signal"""
    signal = ss.Signal()
    signal.connect(reflect)
    signal.emit()
    out, err = capsys.readouterr()
    assert REFLECT_CALLED in out


def test_implicit_emit(capsys):
    """Test implicating emission of signal"""
    signal = ss.Signal()
    signal.connect(reflect)
    signal()
    out, err = capsys.readouterr()
    assert REFLECT_CALLED in out


def test_call():
    """Test calling a signal and getting results"""
    def another_reflect(*args):
        return reflect(*args)

    signal = ss.Signal()
    signal.connect(reflect)
    signal.connect(another_reflect)
    result = list(signal.call('hello'))
    assert len(result) == 2
    assert result[0] == ('hello',)
    assert result[1] == ('hello',)


def test_reduce():
    """Test reducing results to a single result"""
    def another_list_append(inlist, new_item):
        return list_append(inlist, new_item)

    signal = ss.Signal(list_append, another_list_append)
    a_list = []
    result = signal.reduce(a_list, 'hello')
    assert len(result) == 2
    assert result[0] == ['hello', 'hello']
    assert result[1] == 'hello'


def test_disable():
    """Test disable/enable of a signal"""
    signal = ss.Signal()
    signal.connect(reflect)

    signal.enabled = False
    result = list(signal.call())
    assert len(result) == 0

    signal.enabled = True
    result = list(signal.call())
    assert len(result) == 1

    signal.set_enabled(False, push=True)
    result = list(signal.call())
    assert len(result) == 0
    signal.reset_enabled()
    result = list(signal.call())
    assert len(result) == 1


def test_single_shot():
    """Test single shot signals"""
    signal = ss.Signal()
    signal.connect(reflect, single_shot=True)
    result = list(signal.call())
    assert len(result) == 1
    result = list(signal.call())
    assert len(result) == 0


def test_disconnect():
    """Disconnect a slot"""
    signal = ss.Signal()
    signal.connect(reflect)
    result = list(signal.call())
    assert len(result) == 1
    signal.disconnect(reflect)
    result = list(signal.call())
    assert len(result) == 0


def test_clear():
    """Test clearing all slots"""
    signal = ss.Signal()
    signal.connect(reflect)
    result = list(signal.call())
    assert len(result) == 1
    signal.clear()
    result = list(signal.call())
    assert len(result) == 0


def test_call_method():
    """Test with method slots"""
    signal = ss.Signal()
    aclass = AClass()
    signal.connect(aclass.set_and_reflect)
    results = list(signal.call('hello', akeyword='look a keyword'))
    assert len(results) == 1
    assert results[0] == ('hello', )
    assert aclass.kwargs == {'akeyword': 'look a keyword'}

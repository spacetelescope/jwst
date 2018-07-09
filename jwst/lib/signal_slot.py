""" A signal/slot implementation

Original: See below
File:    signal.py
Author:  Thiago Marcos P. Santos
Author:  Christopher S. Case
Author:  David H. Bronke
Created: August 28, 2008
Updated: December 12, 2011
License: MIT

"""
from __future__ import print_function
from collections import namedtuple
import inspect
import logging


__all__ = ['Signal',
           'Signals',
           'SignalsNotAClass']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

Slot = namedtuple('Slot', ['func', 'single_shot'])


class Signal(object):
    def __init__(self, *funcs):
        """Setup a signal

        Parameters
        ----------
        logger: logging.Logger
            Logger to use. If None, one will be created.

        *funcs: (func[, ...))
            Remaining arguments will be functions to connect
            to this signal.
        """
        self._slots = set()
        self._methods = dict()
        self._enabled = True
        self._states = []

        for func in funcs:
            self.connect(func)

    def __call__(self, *args, **kwargs):
        # Call handler functions
        logger.debug(
            'Signal {}: Emitting with args:"{}", kwargs:"{}"'.format(
                self.__class__.__name__,
                args,
                kwargs
            )
        )

        if not self.enabled:
            logger.debug(
                'Signal {}: Disabled, exiting...'.format(
                    self.__class__.__name__
                )
            )
            return

        # No recursive signalling
        self.set_enabled(False, push=True)

        # Call the slots.
        try:
            to_be_removed = []
            for slot in self._slots.copy():
                try:
                    slot.func(*args, **kwargs)
                except RuntimeError:
                    logger.warning(
                        'Signal {}: Signals func->RuntimeError: '
                        'func "{}" will be removed.'.format(
                            self.__class__.__name_,
                            slot.func
                        )
                    )
                    to_be_removed.append(slot)
                finally:
                    if slot.single_shot:
                        to_be_removed.append(slot)

            for remove in to_be_removed:
                self._slots.discard(remove)

            # Call handler methods
            to_be_removed = []
            emitters = self._methods.copy()
            for obj, slots in emitters.items():
                for slot in slots.copy():
                    try:
                        slot.func(obj, *args, **kwargs)
                    except RuntimeError:
                        logger.warning(
                            'Signal {}: Signals methods->RuntimeError, '
                            'obj.func "{}.{}" will be removed'.format(
                                self.__class__.__new__,
                                obj,
                                slot.func
                            )
                        )
                        to_be_removed.append((obj, slot))
                    finally:
                        if slot.single_shot:
                            to_be_removed.append((obj, slot))

            for obj, slot in to_be_removed:
                self._methods[obj].discard(slot)
        finally:
            self.reset_enabled()

    @property
    def enabled(self):
        return self._enabled

    @enabled.setter
    def enabled(self, state):
        self.set_enabled(state, push=False)

    def set_enabled(self, state, push=False):
        """Set whether signal is active or not

        Parameters
        ----------
        state: boolean
            New state of signal

        push: boolean
            If True, current state is saved.
        """
        if push:
            self._states.append(self._enabled)
        self._enabled = state

    def reset_enabled(self):
            self._enabled = self._states.pop()

    def connect(self, func, single_shot=False):
        """Connect a function to the signal
        Parameters
        ----------
        func: function or method
            The function/method to call when the signal is activated

        single_shot: bool
            If True, the function/method is removed after being called.
        """
        logger.debug(
            'Signal {}: Connecting function:"{}"'.format(
                self.__class__.__name__,
                func
            )
        )
        if inspect.ismethod(func):
            if func.__self__ not in self._methods:
                self._methods[func.__self__] = set()

            slot = Slot(
                func=func.__func__,
                single_shot=single_shot
            )
            self._methods[func.__self__].add(slot)

        else:
            slot = Slot(
                func=func,
                single_shot=single_shot
            )
            self._slots.add(slot)

    def disconnect(self, func):
        logger.debug(
            'Signal {}: Disconnecting func:"{}"'.format(
                self.__class__.__name__,
                func
            )
        )
        if inspect.ismethod(func):
            logger.debug(
                'func is a method: "{}"'.format(func)
            )
            if func.__self__ in self._methods:
                logger.debug(
                    'class "{}" is in list'.format(func.__self__)
                )
                logger.debug(
                    'methods="{}"'.format(self._methods[func.__self__])
                )
                slots = [
                    slot
                    for slot in self._methods[func.__self__]
                    if slot.func == func.__func__
                ]
                logger.debug(
                    'slots="{}"'.format(slots)
                )
                try:
                    self._methods[func.__self__].remove(slots[0])
                except IndexError:
                    logger.debug('slot not found.')
                    pass
        else:
            slots = [
                slot
                for slot in self._slots
                if slot.func == func
            ]
            try:
                self._slots.remove(slots[0])
            except IndexError:
                pass

    def clear(self, single_shot=False):
        """Clear slots

        Parameters
        ----------
        single_shot: bool
            If True, only remove single shot
            slots.
        """
        logger.debug(
            'Signal {}: Clearing slots'.format(
                self.__class__.__name__
            )
        )
        if not single_shot:
            self._slots.clear()
            self._methods.clear()
        else:
            to_be_removed = []
            for slot in self._slots.copy():
                if slot.single_shot:
                    to_be_removed.append(slot)
            for remove in to_be_removed:
                self._slots.discard(remove)

            to_be_removed = []
            emitters = self._methods.copy()
            for obj, slots in emitters.items():
                for slot in slots.copy():
                    if slot.single_shot:
                        to_be_removed.append((obj, slot))
            for obj, slot in to_be_removed:
                self._methods[obj].discard(slot)


class SignalsErrorBase(Exception):
    '''Base Signals Error'''

    default_message = ''

    def __init__(self, *args):
        if len(args):
            super(SignalsErrorBase, self).__init__(*args)
        else:
            super(SignalsErrorBase, self).__init__(self.default_message)


class SignalsNotAClass(SignalsErrorBase):
    '''Must add a Signal Class'''
    default_message = 'Signal must be a class.'


class Signals(dict):
    '''Manage the signals.'''

    def __setitem__(self, key, value):
        if key not in self:
            super(Signals, self).__setitem__(key, value)
        else:
            logger.warning('Signals: signal "{}" already exists.'.format(key))

    def __getattr__(self, key):
        for signal in self:
            if signal.__name__ == key:
                return self[signal]
        raise KeyError('{}'.format(key))

    def add(self, signal_class, *args, **kwargs):
        if inspect.isclass(signal_class):
            self.__setitem__(signal_class, signal_class(*args, **kwargs))
        else:
            raise SignalsNotAClass

# ------
# Tests
# ------


def test_signal_slot():

    from functools import partial

    def return_args(returns, *args, **kwargs):
        returns.update({
            'args': args,
            'kwargs': kwargs
        })

    class View(object):
        def __init__(self):
            self.clear()

        def clear(self):
            self.args = None
            self.kwargs = None

        def set(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    # Basic structures
    signal_to_func = Signal()
    assert len(signal_to_func._slots) == 0
    assert len(signal_to_func._methods) == 0

    # Assign a slot
    returns = {}
    slot = partial(return_args, returns)
    signal_to_func.connect(slot)
    signal_to_func()
    assert len(returns) > 0
    assert len(returns['args']) == 0
    assert len(returns['kwargs']) == 0

    # Signal with arguments
    returns.clear()
    an_arg = 'an arg'
    signal_to_func(an_arg)
    assert len(returns['args']) > 0
    assert returns['args'][0] == an_arg
    signal_to_func(a_kwarg=an_arg)
    assert len(returns['kwargs']) > 0
    assert returns['kwargs']['a_kwarg'] == an_arg

    # Signal with methods
    signal_to_method = Signal()
    view = View()
    signal_to_method.connect(view.set)
    signal_to_method()
    assert len(view.args) == 0
    assert len(view.kwargs) == 0

    # Signal with methods and arguments
    view.clear()
    signal_to_method(an_arg)
    assert view.args[0] == an_arg
    view.clear()
    signal_to_method(a_kwarg=an_arg)
    assert view.kwargs['a_kwarg'] == an_arg

    # Delete some slots
    returns.clear()
    signal_to_func.disconnect(slot)
    signal_to_func(an_arg, a_kwarg=an_arg)
    assert len(returns) == 0
    view.clear()
    signal_to_method.disconnect(view.set)
    signal_to_method(an_arg, a_kwarg=an_arg)
    assert view.args is None
    assert view.kwargs is None

    # Test initialization
    a_signal = Signal(slot, view.set)
    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert returns['args'][0] == an_arg
    assert returns['kwargs']['a_kwarg'] == an_arg
    assert view.args[0] == an_arg
    assert view.kwargs['a_kwarg'] == an_arg

    # Clear a signal
    a_signal.clear()
    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert len(returns) == 0
    assert view.args is None
    assert view.kwargs is None

    # Enable/disable
    a_signal = Signal(slot, view.set)
    assert a_signal.enabled

    a_signal.enabled = False
    assert not a_signal.enabled
    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert len(returns) == 0
    assert view.args is None
    assert view.kwargs is None

    a_signal.enabled = True
    assert a_signal.enabled
    a_signal(an_arg, a_kwarg=an_arg)
    assert returns['args'][0] == an_arg
    assert returns['kwargs']['a_kwarg'] == an_arg
    assert view.args[0] == an_arg
    assert view.kwargs['a_kwarg'] == an_arg

    a_signal.set_enabled(False, push=True)
    assert not a_signal.enabled
    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert len(returns) == 0
    assert view.args is None
    assert view.kwargs is None

    a_signal.reset_enabled()
    assert a_signal.enabled
    a_signal(an_arg, a_kwarg=an_arg)
    assert returns['args'][0] == an_arg
    assert returns['kwargs']['a_kwarg'] == an_arg
    assert view.args[0] == an_arg
    assert view.kwargs['a_kwarg'] == an_arg

    # Single shots
    a_signal = Signal()
    a_signal.connect(slot, single_shot=True)
    a_signal.connect(view.set, single_shot=True)
    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert returns['args'][0] == an_arg
    assert returns['kwargs']['a_kwarg'] == an_arg
    assert view.args[0] == an_arg
    assert view.kwargs['a_kwarg'] == an_arg

    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert len(returns) == 0
    assert view.args is None
    assert view.kwargs is None

    # Clearing single shots
    a_signal = Signal()
    a_signal.connect(slot, single_shot=True)
    a_signal.connect(view.set, single_shot=True)
    a_signal.clear(single_shot=True)
    returns.clear()
    view.clear()
    a_signal(an_arg, a_kwarg=an_arg)
    assert len(returns) == 0
    assert view.args is None
    assert view.kwargs is None


if __name__ == '__main__':
    test_signal_slot()

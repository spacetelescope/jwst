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
from collections import namedtuple
from functools import partial
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

    def emit(self, *args, **kwargs):
        """Invoke slots attached to the signal

        No return of results is expected.

        Parameters
        ----------
        args: (arg[, ...])
            Positional arguments to pass to the slots.

        kwargs: {key: value[, ...]}
            Keyword arguments to pass to the slots.
        """
        for _ in self.call(*args, **kwargs):
            pass

    __call__ = emit

    def call(self, *args, **kwargs):
        """Generator returning result of each slot connected.

        Parameters
        ----------
        args: (arg[, ...])
            Positional arguments to pass to the slots.

        kwargs: {key: value[, ...]}
            Keyword arguments to pass to the slots.

        Returns
        -------
        generator
            A generator returning the result from each slot.
        """
        for slot in self.slots:
            try:
                yield slot(*args, **kwargs)
            except Exception as exception:
                logger.debug(
                    'Signal {}: Slot {} raised {}'.format(
                        self.__class__.__name_,
                        slot,
                        exception
                    )
                )

    def reduce(self, *args, **kwargs):
        """Return a reduction of all the slots

        Parameters
        ----------
        args: (arg[, ...])
            Positional arguments to pass to the slots.

        kwargs: {key: value[, ...]}
            Keyword arguments to pass to the slots.

        Returns
        -------
        result: object or (object [,...])
            The result or tuple of results. See `Notes`.


        Notes
        -----

        Each slot is given the results of the previous
        slot as a new positional argument list. As such, if multiple
        arguments are required, each slot should return a tuple that
        can then be passed as arguments to the next function.

        The keywoard arguments are simply passed to each slot unchanged.

        There is no guarantee on order which the slots are invoked.

        """
        result = None
        for slot in self.slots:
            result = slot(*args, **kwargs)
            args = result
            if not isinstance(args, tuple):
                args = (args, )
        return result

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

    @property
    def slots(self):
        """Generator returning slots"""
        if not self.enabled:
            return

        # No recursive signalling
        self.set_enabled(False, push=True)

        try:
            for slot in self._slots:
                yield slot.func
            for obj, slots in self._methods.items():
                for slot in slots:
                    yield partial(slot.func, obj)
        finally:
            # Clean out single shots
            self._slots = [
                slot
                for slot in self._slots
                if not slot.single_shot
            ]
            for obj in self._methods:
                self._methods[obj] = [
                    slot
                    for slot in self._methods[obj]
                    if not slot.single_shot
                ]
            self.reset_enabled()


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

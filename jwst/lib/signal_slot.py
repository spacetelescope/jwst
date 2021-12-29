""" A signal/slot implementation
"""
from collections import namedtuple
import inspect
import logging


__all__ = ['Signal',
           'Signals',
           'SignalsNotAClass']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

"""Slot data structure
"""
Slot = namedtuple('Slot', ['func', 'single_shot'])


class Signal():
    """Signal

    A Signal, when triggered, call the connected slots.

    Parameters
    ----------
    funcs: func[, ...]
        Remaining arguments will be functions to connect
        to this signal.

    Attributes
    ----------
    enabled: bool
        If True, the slots are called. Otherwise, nothing
        happens when triggered.

    """

    def __init__(self, *funcs):
        self._slots = list()
        self._enabled = True
        self._states = list()

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

        The keyword arguments are simply passed to each slot unchanged.

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
        slot = Slot(
            func=func,
            single_shot=single_shot
        )
        self._slots.append(slot)

    def disconnect(self, func):
        self._slots = [
            slot
            for slot in self._slots
            if slot.func != func
        ]

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
        else:
            self._slots = [
                slot
                for slot in self._slots
                if not slot.single_shot
            ]

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
        finally:
            # Clean out single shots
            self._slots = [
                slot
                for slot in self._slots
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

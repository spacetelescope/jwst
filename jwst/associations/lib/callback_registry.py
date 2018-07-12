"""Callback registry"""

from ...lib.signal_slot import Signal

__all__ = ['CallbackRegistry']


class CallbackRegistry():
    """Callback registry"""

    def __init__(self):
        self.registry = dict()

    def add(self, event, callback):
        """Add a callback to an event"""
        try:
            signal = self.registry[event]
        except KeyError:
            signal = Signal()
        signal.connect(callback)
        self.registry[event] = signal

    def reduce(self, event, *args):
        """Peform a reduction on the event args

        Parameters
        ----------
        args: [arg[,...]]
            The args to filter

        Returns
        -------
        The reduced results.
        If no results can be determined,
        such as if no callbacks were registered,
        `None` is returned.

        Notes
        -----
        Each function is given the results of the previous
        function. As such, if the data has more than one
        object, the return of each function should be a tuple that can
        then be passed as arguments to the next function.

        There is no guarantee on order which the registered
        callbacks are made. Currently, the callbacks are in a list.
        Hence, the callbacks will be called in the order registered.

        """
        result = self.registry[event].reduce(*args)
        return result

    def add_decorator(self, event):
        """Add callbacks by decoration

        Parameters
        ----------
        event: str
            The name of event to attach the object to.
        """
        def decorator(func):
            self.add(event, func)
            return func
        return decorator

    __call__ = add_decorator


# ##########
# Test suite
# ##########

class TestSuite():

    @staticmethod
    def fn_single_1(arg):
        return arg + 1

    @staticmethod
    def fn_single_2(arg):
        return arg + arg

    @staticmethod
    def fn_double_1(arg1, arg2):
        return (arg1 + 1, arg2 - 1)

    @staticmethod
    def fn_double_2(*args):
        return (args[0] + args[0], args[1] + args[1])

    def test_basics(self):
        cbr = CallbackRegistry()
        cbr.add('event1', TestSuite.fn_single_1)
        assert cbr.reduce('event1', 1) == 2

    def test_reduce(self):
        cbr = CallbackRegistry()
        cbr.add('event1', TestSuite.fn_single_1)
        cbr.add('event1', TestSuite.fn_single_2)
        assert cbr.reduce('event1', 1) == 4

    def test_double(self):
        cbr = CallbackRegistry()
        cbr.add('event1', TestSuite.fn_double_1)
        assert cbr.reduce('event1', 1, 2) == (2, 1)

        cbr.add('event2', TestSuite.fn_double_2)
        assert cbr.reduce('event2', 2, 1) == (4, 2)

    def test_double_reduce(self):
        cbr = CallbackRegistry()
        cbr.add('event1', TestSuite.fn_double_1)
        cbr.add('event1', TestSuite.fn_double_2)
        result = cbr.reduce('event1', 1, 2)
        assert result in [(4, 2), (3, 3)]

    def test_decorator(self):
        cbr = CallbackRegistry()

        @cbr('event1')
        def fn_single_1(arg):
            return arg + 1

        @cbr('event1')
        def fn_single_2(arg):
            return arg + arg

        assert cbr.reduce('event1', 1) in (3, 4)

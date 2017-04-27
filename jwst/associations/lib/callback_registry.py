"""Callback registry"""

from collections import defaultdict

__all__ = ['CallbackRegistry']


class CallbackRegistry(object):
    """Callback registry"""

    def __init__(self):
        self.registry = defaultdict(list)

    def add(self, event, callback):
        """Add a callback to an event"""
        self.registry[event].append(callback)

    def reduce(self, event, *args):
        """Peform a reduction on the event args

        Parameters
        ----------
        args: [arg[,...]]
            The args to reduce

        Returns
        -------
        The reduction
        """
        for callback in self.registry[event]:
            result = callback(*args)
            args = result
            if not isinstance(args, tuple):
                args = (args, )
        return result

    def __call__(self, event):
        """Add callback by calling instance

        Allows instances to be used as a decorator
        """
        def wrapper(callback):
            self.add(event, callback)
        return wrapper


# ##########
# Test suite
# ##########

class TestSuite(object):

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

    def test_double_reduce(self):
        cbr = CallbackRegistry()
        cbr.add('event1', TestSuite.fn_double_1)
        cbr.add('event1', TestSuite.fn_double_2)
        assert cbr.reduce('event1', 1, 2) == (4, 2)

    def test_decorator(self):
        cbr = CallbackRegistry()

        @cbr('event1')
        def fn_single_1(arg):
            return arg + 1

        @cbr('event1')
        def fn_single_2(arg):
            return arg + arg

        assert cbr.reduce('event1', 1) == 4

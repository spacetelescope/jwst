"""
Callback registry.

This registry stores events which may be triggered on
an AssociationRegistry given some signal; for instance,
finalization of associations.
"""

from jwst.lib.signal_slot import Signal

__all__ = ["CallbackRegistry"]


class CallbackRegistry:
    """Callback registry."""

    def __init__(self):
        self.registry = {}

    def add(self, event, callback):
        """Add a callback to an event."""
        try:
            signal = self.registry[event]
        except KeyError:
            signal = Signal()
        signal.connect(callback)
        self.registry[event] = signal

    def reduce(self, event, *args):
        """
        Perform a reduction on the event args.

        Parameters
        ----------
        *args : [arg[,...]]
            The args to filter

        Returns
        -------
        list, dict or None
            The reduced results; if no results can be determined,
            such as if no callbacks were registered, `None` is returned.

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
        """
        Add callbacks by decoration.

        Parameters
        ----------
        event : str
            The name of event to attach the object to.

        Returns
        -------
        decorator
            The decorator attached to the provided event name.
        """

        def decorator(func):
            self.add(event, func)
            return func

        return decorator

    __call__ = add_decorator

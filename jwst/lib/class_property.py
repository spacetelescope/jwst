"""Allow property access to both instance and class attributes"""

class ClassProperty:
    """Allow class-level property descriptors

    This property allows definition of a non-data descriptor to
    be defined that is accessible from both the class and instances
    of the class.

    Parameters
    ----------
    fget: function
        The function to use as the `__get__` method.

    Notes
    -----
    On class access, only the `__get__` method of a descriptor is used.
    All other descriptor methods are ignored. As such, allowing the definitions
    of `__set__` and `__delete__` are not allowed.

    The code is based on the `Property` example descriptor in the python
    documentation.

    Examples
    --------
    Usage is the same as the Python `property` descriptor. However, access
    is available from both the class and instances.

    >>> class MyClass:
    ...     @ClassProperty
    ...     def baz(obj):
    ...         return obj._baz
    ...
    ...     _baz = 'class'

    >>> MyClass.baz
    'class'

    >>> mc = MyClass()
    >>> mc.baz
    'class'

    >>> mc._baz = 'instance'
    >>> mc.baz
    'instance'

    >>> MyClass.baz
    'class'
    """
    def __init__(self, fget=None, doc=None):
        self.fget = fget
        if doc is None and fget is not None:
            doc = fget.__doc__
        self.__doc__ = doc

    def __get__(self, obj, objtype=None):
        if obj is None:
            obj = objtype
        if self.fget is None:
            raise AttributeError("unreadable attribute")
        return self.fget(obj)

    def getter(self, fget):
        return type(self)(fget, self.__doc__)

"""Constraints
"""
from copy import deepcopy


class SimpleConstraint:
    """A basic constraint

    Parameters
    ----------
    init: dict
        dict where the key:value pairs define
        the following parameters

    value: object or None
        Value that must be matched.

    force_unique: bool
        If the constraint is satisfied, reset `value`
        to the value of the source.

    test: function
        The test function for the constraint.
        Takes two arguments:
            - constraint
            - object to compare against.
        Returns a boolean.
        Default is `SimpleConstraint.eq`

    Attributes
    ----------
    All `Parameters` are also `Attributes`

    Examples
    --------

    Create a constraint where the attribute `attr` of an object
    matches the value `my_value`:

    >>> from jwst.associations.lib.constraint import SimpleConstraint
    >>> c = SimpleConstraint(value='my_value')
    >>> print(c)
    SimpleConstraint({'value': 'my_value' })

    To check a constraint, call `check_and_set`. A successful match
    will return a `SimpleConstraint` and a reprocess list.
    >>> item = 'my_value'
    >>> new_c, reprocess = c.check_and_set(item)
    SimpleConstraint, []

    If it doesn't match, `False` will be returned.
    >>> bad_item = 'not_my_value'
    >>> c.check_and_set(bad_item)
    False, []

    A `SimpleConstraint` can also be initialized by a `dict`
    of the relevant parameters:
    >>> init = {'value': 'my_value'}
    >>> c = SimpleConstraint(init)
    >>> print(c)
    SimpleConstraint({'value': 'my_value'})

    If the value to check is `None`, the `SimpleConstraint` will
    succesfully match whatever object given. However, a new `SimpleConstraint`
    will be returned where the `value` is now set to whatever the attribute
    was of the object.
    >>> c = SimpleConstraint(value=None, sources=['attr'])
    >>> new_c, reprocess = c.check_and_set(item)
    >>> print(result)
    SimpleConstraint({'value': 'my_value'})

    This behavior can be overriden by the `force_unique` paramter:
    >>> c = SimpleConstraint(value=None, sources=['attr'], force_unique=False)
    >>> result, reprocess = c.check_and_set(item)
    >>> print(result)
    SimpleConstraint({'value': None})
    """

    def __init__(self, init=None, **kwargs):

        # Defined attributes
        self.value = None
        self.force_unique = True
        self.test = self.eq

        # The reprocess list. For this class, an empty list
        # is always returned.
        self._reprocess = []

        if init is not None:
            self.__dict__.update(init)
        else:
            self.__dict__.update(kwargs)

    def check_and_set(self, item):
        """Check and set the constraint

        Returns
        -------
        success: SimpleConstraint or False
            If successful, a copy of the constraint
            is returned with modified value.
        """
        if self.value is not None:
            satisfied = self.test(item)
            if satisfied:
                return self, self._reprocess
            else:
                return False, self._reprocess

        updated = self
        if self.force_unique:
            updated = deepcopy(self)
            updated.value = item
        return updated, self._reprocess

    def eq(self, item):
        """True if constraint.value and item are equal."""
        return self.value == item

    def __str__(self):
        result = '{}({})'.format(
            self.__class__.__name__,
            str(self.__dict__)
        )
        return result


class Constraint:
    """Constraint that is made up of SimpleConstraint

    Parameters
    ----------
    init: object or [object[,...]]
        A single object or list of objects where the
        objects are as follows.
        - SimpleConstraint or subclass
        - Constraint

    reduce: function
        A reduction function with signature `x(iterable)`
        where `iterable` is the `components` list. Returns
        boolean indicating state of the components.
        Default value is `Constraint.all`

    Attributes
    ----------
    constraints: [Constraint[,...]]
        `Constraint`s or `SimpleConstaint`s that
        make this constraint.

    reduce: function
        A reduction function with signature `x(iterable)`
        where `iterable` is the `components` list. Returns
        boolean indicating state of the components.
        Predefined functions are:
        - `all`: True if all components return True
        - `any`: True if any component returns True
    """
    def __init__(self, init=None, reduce=None):
        self.reduce = reduce
        if self.reduce is None:
            self.reduce = self.all
        self.constraints = []
        if init is None:
            pass
        elif isinstance(init, list):
            self.constraints = init
        elif isinstance(init, Constraint):
            self.reduce = init.reduce
            self.constraints = deepcopy(init.constraints)
        elif isinstance(init, SimpleConstraint):
            self.constraints = [init]
        else:
            raise TypeError(
                'Invalid initialization value type {}.'
                '\nValid types are `SimpleConstaint`, `Constraint`,'
                '\nor subclass.'.format(type(init))
            )

    def check_and_set(self, item):
        """Check and set the constraint

        Returns
        -------
        2-tuple of (`Constraint`, reprocess)
        """
        # Do we have positive?
        results = [
            constraint.check_and_set(item)
            for constraint in self.constraints
        ]
        result = self.reduce(results)

        # If a positive, replace positive returning
        # constraints in the list.
        all_reprocess = []
        new_constraint = False
        if result:
            new_constraint = Constraint(reduce=self.reduce)
            for idx, (constraint, reprocess) in enumerate(results):
                all_reprocess.extend(reprocess)
                if constraint:
                    new_constraint.constraints.append(constraint)
                else:
                    new_constraint.constraints.append(self.constraints[idx])

        return new_constraint, all_reprocess

    @staticmethod
    def all(results):
        constraints, reprocess = zip(*results)
        return all(constraints)

    @staticmethod
    def any(results):
        constraints, reprocess = zip(*results)
        return any(constraints)

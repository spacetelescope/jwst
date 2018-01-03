"""Constraints
"""
import abc
from copy import deepcopy
from itertools import chain
import re

from .utilities import (
    evaluate,
    getattr_from_list,
    is_iterable
)

__all__ = [
    'AttrConstraint',
    'Constraint',
    'ConstraintTrue',
    'SimpleConstraint',
]


class SimpleConstraintABC(abc.ABC):
    """Simple Constraint ABC

    Parameters
    ----------
    init: dict
        dict where the key:value pairs define
        the following parameters

    value: object or None
        Value that must be matched.

    name: str or None
        Option name for constraint

    **kwargs: key:value pairs
        Other initialization parameters
    """

    def __init__(self, init=None, value=None, name=None, **kwargs):

        # Defined attributes
        self.value = value
        self.name = name

        if init is not None:
            self.__dict__.update(init)
        else:
            self.__dict__.update(kwargs)

    @abc.abstractmethod
    def check_and_set(self, item):
        """Check and set the constraint

        Returns
        -------
        success: obj or False
            If constraint is satisfied, some type of object, may be `True`,
            is returned.
        """
        return self, []

    # Make iterable to work with `Constraint`.
    # Since this is a leaf, simple return ourselves.
    def __iter__(self):
        yield self

    def __str__(self):
        result = '{}({})'.format(
            self.__class__.__name__,
            str(self.__dict__)
        )
        return result


class ConstraintTrue(SimpleConstraintABC):
    """Always return True"""

    def check_and_set(self, item):
        return super(ConstraintTrue, self).check_and_set(item)


class SimpleConstraint(SimpleConstraintABC):
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

    name: str or None
        Option name for constraint

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

    def __init__(self, init=None, force_unique=True, test=None, **kwargs):

        # Defined attributes
        super(SimpleConstraint, self).__init__(init=init, **kwargs)
        self.force_unique = force_unique
        if test is None:
            self.test = self.eq
        else:
            self.test = test

        # The reprocess list. For this class, an empty list
        # is always returned.
        self._reprocess = []

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


class AttrConstraint(SimpleConstraintABC):
    """Test attribute of an item

    Parameters
    ----------
    sources: [str[,...]]
        List of attributes to query

    value: str or None
        The value to check for. If None and
        `force_unique`, any value in the first
        available source will become the value.

    force_unique: bool
        If the initial value of `value` is None,
        `value` will be set to the first source.
        Otherwise, this will be left as None.

    required: bool
        One of the sources must exist. Otherwise,
        return as a matched constraint.

    onlyif: function
        Boolean function that takes `item` as argument.
        If True, the rest of the condition is checked. Otherwise
        return as a matched condition

    force_reprocess: bool
        Add item back onto the reprocess list.

    evaluate: bool
        Evaluate the item's value before checking condition.

    found_values: set(str[,...])
        Set of actual found values for this condition.

    name: str or None
        Name of the constraint.

    invalid_values: [str[,...]]
        List of values that are invalid in an item.
        Will cause a non-match.
    """

    def __init__(self,
                 init=None,
                 evaluate=False,
                 force_undefined=False,
                 force_unique=True,
                 invalid_values=None,
                 onlyif=None,
                 required=True,
                 **kwargs):

        # Attributes
        super(AttrConstraint, self).__init__(init=init, **kwargs)
        self.evaluate = evaluate
        self.force_unique = force_unique
        self.force_undefined = force_undefined
        if invalid_values is None:
            self.invalid_values = []
        else:
            self.invalid_values = invalid_values
        if onlyif is None:
            self.onlyif = lambda item: True
        else:
            self.onlyif = onlyif
        self.required = True

        self.found_values = set()

    def check_and_set(self, item):
        """Check and set constraints based on item

        Parameters
        ----------
        item: dict
            The item to check on.

        Returns
        -------
        matching_constraint, reprocess: AttrConstraint or False
            A 2-tuple consisting of:
            - matching_constraint if a successful match
            - List of items that need to be checked again.
        """
        reprocess = []
        match = deepcopy(self)

        # Only perform check on specified `onlyif` condition
        if not match.onlyif(item):
            if match.force_reprocess:
                reprocess.append(item)
            return (match, reprocess)

        # Get the condition information.
        try:
            source, value = getattr_from_list(
                item,
                match.sources,
                invalid_values=match.invalid_values
            )
        except KeyError:
            if match.required and not match.force_undefined:
                return False, reprocess
            else:
                return match, reprocess
        else:
            if match.force_undefined:
                return False, reprocess

        # If the value is a list, build the reprocess list
        if match.evaluate:
            evaled = evaluate(value)
            if is_iterable(evaled):
                for avalue in evaled:
                    new_item = deepcopy(item)
                    new_item[source] = str(avalue)
                    reprocess.append(new_item)
                return False, reprocess
            value = str(evaled)

        # Check condition
        if match.value is not None:
            if not meets_conditions(
                    value, match.value
            ):
                return False, reprocess

        # At this point, the constraint has passed.
        # Fix the conditions.
        escaped_value = re.escape(value)
        match.found_values.add(escaped_value)
        if match.value is None or match.force_unique:
            match.value = escaped_value
            match.sources = [source]
            match.force_unique = False

        # That's all folks
        return match, reprocess


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

    name: str or None
        Optional name for constraint.

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

    Notes
    -----
    Named constraints can be accessed directly through indexing:

    >>> c = Constraint(SimpleConstaint(name='simple', value='a_value'))
    >>> c['simple']
    SimpleConstraint('value': 'a_value')
    """
    def __init__(self, init=None, reduce=None, name=None):
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
        matches, reprocess = zip(*results)
        match = self.reduce(matches)

        all_reprocess = chain(*reprocess)

        # If a positive, replace positive returning
        # constraints in the list.
        new_constraint = False
        if match:
            new_constraint = Constraint(self)
            for idx, constraint in enumerate(matches):
                if constraint:
                    new_constraint.constraints[idx] = constraint

        return new_constraint, all_reprocess

    @staticmethod
    def all(matches):
        return all(matches)

    @staticmethod
    def any(matches):
        return any(matches)

    # Make iterable
    def __iter__(self):
        for constraint in chain(*map(iter, self.constraints)):
            yield constraint

    # Index implementaion
    def __getitem__(self, key):
        """Retrieve a named constraint"""
        for constraint in self.constraints:
            name = getattr(constraint, 'name', None)
            if name is not None and name == key:
                return constraint
            try:
                found = constraint[key]
            except (KeyError, TypeError):
                pass
            else:
                return found
        raise KeyError('Constraint {} not found'.format(key))

    def __setitem__(self, key, value):
        """Not implemented"""
        raise NotImplemented('Cannot set constraints by index.')

    def __delitem__(self, key):
        """Not implemented"""
        raise NotImplementedError('Cannot delete a constraint by index.')

    def __str__(self):
        """String representation"""
        result = 'Constraint {} with {} on:\n{}'.format(
            getattr(self, 'name', 'Unnamed'),
            str(self.reduce),
            ''.join([
                str(constraint)
                for constraint in self.constraints
            ])
        )

        return result + '\n'


# ---------
# Utilities
# ---------
def meets_conditions(value, conditions):
    """Check whether value meets any of the provided conditions

    Parameters
    ----------
    values: str
        The value to be check with.

    condition: regex,
        Regular expressions to match against.

    Returns
    -------
    True if any condition is meant.
    """

    if not is_iterable(conditions):
        conditions = [conditions]
    for condition in conditions:
        condition = ''.join([
            '^',
            condition,
            '$'
        ])
        match = re.match(condition, value, flags=re.IGNORECASE)
        if match:
            return True
    return False

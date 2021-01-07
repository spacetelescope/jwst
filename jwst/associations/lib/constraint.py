"""Constraints
"""
import abc
import collections
from copy import deepcopy
from itertools import chain
import logging
import re
import typing

from .process_list import ProcessList
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

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class SimpleConstraintABC(abc.ABC):
    """Simple Constraint ABC

    Parameters
    ----------
    init : dict
        dict where the key:value pairs define
        the following parameters

    value : object or None
        Value that must be matched.

    name : str or None
        Option name for constraint

    **kwargs : key:value pairs
        Other initialization parameters

    Attributes
    ----------
    matched : bool
        Last call to `check_and_set`
    """

    # Attributes to show in the string representation.
    _str_attrs = ('name', 'value')

    def __init__(self, init=None, value=None, name=None, **kwargs):

        # Defined attributes
        self.value = value
        self.name = name
        self.matched = False

        if init is not None:
            self.__dict__.update(init)
        else:
            self.__dict__.update(kwargs)

    @abc.abstractmethod
    def check_and_set(self, item):
        """Check and set the constraint

        Returns
        -------
        success, reprocess : bool, [ProcessList[,...]]
            Returns 2-tuple of

                - True if check is successful.
                - List of `ProcessList`.
        """
        self.matched = True
        return self.matched, []

    def copy(self):
        """Copy ourselves"""
        return deepcopy(self)

    @property
    def dup_names(self): #  -> dict[str, list[typing.Union[SimpleConstraint, Constraint]]]
        """Return dictionary of constraints with duplicate names

        This method is meant to be overridden by classes
        that need to traverse a list of constraints.

        Returns
        -------
        dups : {str: [constraint[,...]][,...]}
            Returns a mapping between the duplicated name
            and all the constraints that define that name.
        """
        return {}

    def get_all_attr(self, attribute: str): # -> list[tuple[SimpleConstraint, typing.Any]]:
        """Return the specified attribute

        This method is meant to be overridden by classes
        that need to traverse a list of constraints.

        Parameters
        ----------
        attribute : str
            The attribute to retrieve

        Returns
        -------
        [(self, value)] : [(SimpleConstraint, object)]
            The value of the attribute in a tuple. If there is no attribute,
            an empty tuple is returned.
        """
        value = getattr(self, attribute)
        if value is not None:
            return [(self, value)]
        return []

    # Make iterable to work with `Constraint`.
    # Since this is a leaf, simple return ourselves.
    def __iter__(self):
        yield self

    def __repr__(self):
        result = '{}({})'.format(
            self.__class__.__name__,
            str(self.__dict__)
        )
        return result

    def __str__(self):
        result = '{}({})'.format(
            self.__class__.__name__,
            {
                str_attr: getattr(self, str_attr)
                for str_attr in self._str_attrs
            }
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
    init : dict
        dict where the key:value pairs define
        the following parameters

    value : object or None
        Value that must be matched.
        If None, any retrieved value will match.

    sources : func(item) or None
        Function taking `item` as argument used to
        retrieve a value to check against.
        If None, the item itself is used as the value.

    force_unique : bool
        If the constraint is satisfied, reset `value`
        to the value of the source.

    test : function
        The test function for the constraint.
        Takes two arguments:

            - constraint
            - object to compare against.

        Returns a boolean.
        Default is `SimpleConstraint.eq`

    name : str or None
        Option name for constraint

    reprocess_on_match : bool
        Reprocess the item if the constraint is satisfied.

    reprocess_on_fail : bool
        Reprocess the item if the constraint is not satisfied.

    work_over : ProcessList.[BOTH, EXISTING, RULES]
        The condition on which this constraint should operate.

    reprocess_rules : [rule[,..]] or None
        List of rules to be applied to.
        If None, calling function will determine the ruleset.
        If empty, [], all rules will be used.

    Attributes
    ----------
    All `Parameters` are also `Attributes`

    Examples
    --------

    Create a constraint where the attribute `attr` of an object
    matches the value `my_value`:

    >>> c = SimpleConstraint(value='my_value')
    >>> print(c)
    SimpleConstraint({'name': None, 'value': 'my_value'})

    To check a constraint, call `check_and_set`. A successful match
    will return a tuple of `True` and a reprocess list.
    >>> item = 'my_value'
    >>> c.check_and_set(item)
    (True, [])

    If it doesn't match, `False` will be returned.
    >>> bad_item = 'not_my_value'
    >>> c.check_and_set(bad_item)
    (False, [])

    A `SimpleConstraint` can also be initialized by a `dict`
    of the relevant parameters:
    >>> init = {'value': 'my_value'}
    >>> c = SimpleConstraint(init)
    >>> print(c)
    SimpleConstraint({'name': None, 'value': 'my_value'})

    If the value to check is `None`, the `SimpleConstraint` will
    succesfully match whatever object given. However, a new `SimpleConstraint`
    will be returned where the `value` is now set to whatever the attribute
    was of the object.
    >>> c = SimpleConstraint(value=None)
    >>> matched, reprocess = c.check_and_set(item)
    >>> print(c)
    SimpleConstraint({'name': None, 'value': 'my_value'})

    This behavior can be overriden by the `force_unique` paramter:
    >>> c = SimpleConstraint(value=None, force_unique=False)
    >>> matched, reprocess = c.check_and_set(item)
    >>> print(c)
    SimpleConstraint({'name': None, 'value': None})

    """

    def __init__(
            self,
            init=None,
            sources=None,
            force_unique=True,
            test=None,
            reprocess_on_match=False,
            reprocess_on_fail=False,
            work_over=ProcessList.BOTH,
            reprocess_rules=None,
            **kwargs
    ):

        # Defined attributes
        self.sources = sources
        self.force_unique = force_unique
        self.test = test
        self.reprocess_on_match = reprocess_on_match
        self.reprocess_on_fail = reprocess_on_fail
        self.work_over = work_over
        self.reprocess_rules = reprocess_rules
        super(SimpleConstraint, self).__init__(init=init, **kwargs)

        # Give defaults some real meaning.
        if self.sources is None:
            self.sources = lambda item: item
        if test is None:
            self.test = self.eq

    def check_and_set(self, item):
        """Check and set the constraint

        Returns
        -------
        success, reprocess : bool, [ProcessList[,...]]
            Returns 2-tuple of

                - True if check is successful.
                - List of `ProcessList`.
        """
        source_value = self.sources(item)

        satisfied = True
        if self.value is not None:
            satisfied = self.test(self.value, source_value)
        self.matched = satisfied

        if self.matched:
            if self.force_unique:
                self.value = source_value

        # Determine reprocessing
        reprocess = []
        if (self.matched and self.reprocess_on_match) or \
           (not self.matched and self.reprocess_on_fail):
            reprocess.append(ProcessList(
                items=[item],
                work_over=self.work_over,
                rules=self.reprocess_rules
            ))

        return self.matched, reprocess

    def eq(self, value1, value2):
        """True if constraint.value and item are equal."""
        return value1 == value2


class AttrConstraint(SimpleConstraintABC):
    """Test attribute of an item

    Parameters
    ----------
    sources : [str[,...]]
        List of attributes to query

    value : str, function or None
        The value to check for. If None and
        `force_unique`, any value in the first
        available source will become the value.
        If function, the function takes no arguments
        and returns a string.

    evaluate : bool
        Evaluate the item's value before checking condition.

    force_reprocess : ProcessList.state or False
        Add item back onto the reprocess list using
        the specified `ProcessList` work over state.

    force_unique : bool
        If the initial value is `None` or a list of possible values,
        the constraint will be modified to be the value first matched.

    invalid_values : [str[,...]]
        List of values that are invalid in an item.
        Will cause a non-match.

    name : str or None
        Name of the constraint.

    only_on_match : bool
        If `force_reprocess`, only do the reprocess
        if the entire constraint is satisfied.

    onlyif : function
        Boolean function that takes `item` as argument.
        If True, the rest of the condition is checked. Otherwise
        return as a matched condition

    required : bool
        One of the sources must exist. Otherwise,
        return as a matched constraint.

    Attributes
    ----------
    found_values : set(str[,...])
        Set of actual found values for this condition.

    matched : bool
        Last result of `check_and_set`

    """

    # Attributes to show in the string representation.
    _str_attrs = ('name', 'sources', 'value')

    def __init__(self,
                 init=None,
                 sources=None,
                 evaluate=False,
                 force_reprocess=False,
                 force_undefined=False,
                 force_unique=True,
                 invalid_values=None,
                 only_on_match=False,
                 onlyif=None,
                 required=True,
                 **kwargs):

        # Attributes
        self.sources = sources
        self.evaluate = evaluate
        self.force_reprocess = force_reprocess
        self.force_undefined = force_undefined
        self.force_unique = force_unique
        self.invalid_values = invalid_values
        self.only_on_match = only_on_match
        self.onlyif = onlyif
        self.required = required
        super(AttrConstraint, self).__init__(init=init, **kwargs)

        # Give some defaults real meaning.
        if invalid_values is None:
            self.invalid_values = []
        if onlyif is None:
            self.onlyif = lambda item: True

        # Haven't actually matched anything yet.
        self.found_values = set()
        self.matched = False

    def check_and_set(self, item):
        """Check and set constraints based on item

        Parameters
        ----------
        item : dict
            The item to check on.

        Returns
        -------
        success, reprocess : bool, [ProcessList[,...]]
            Returns 2-tuple of

                - True if check is successful.
                - List of `ProcessList`.
        """
        reprocess = []

        # Only perform check on specified `onlyif` condition
        if not self.onlyif(item):
            if self.force_reprocess:
                reprocess.append(
                    ProcessList(
                        items=[item],
                        work_over=self.force_reprocess,
                        only_on_match=self.only_on_match,
                    )
                )
            self.matched = True
            return (self.matched, reprocess)

        # Get the condition information.
        try:
            source, value = getattr_from_list(
                item,
                self.sources,
                invalid_values=self.invalid_values
            )
        except KeyError:
            if self.required and not self.force_undefined:
                self.matched = False
                return self.matched, reprocess
            else:
                self.matched = True
                return self.matched, reprocess
        else:
            if self.force_undefined:
                self.matched = False
                return self.matched, reprocess

        # If the value is a list, build the reprocess list
        if self.evaluate:
            evaled = evaluate(value)
            if is_iterable(evaled):
                reprocess_items = []
                for avalue in evaled:
                    new_item = deepcopy(item)
                    new_item[source] = str(avalue)
                    reprocess_items.append(new_item)
                reprocess.append(ProcessList(
                    items=reprocess_items,
                ))
                self.matched = False
                return self.matched, reprocess
            value = str(evaled)

        # Check condition
        if self.value is not None:
            if callable(self.value):
                match_value = self.value()
            else:
                match_value = self.value
            if not meets_conditions(
                    value, match_value
            ):
                self.matched = False
                return self.matched, reprocess

        # At this point, the constraint has passed.
        # Fix the conditions.
        escaped_value = re.escape(value)
        self.found_values.add(escaped_value)
        if self.force_unique:
            self.value = escaped_value
            self.sources = [source]
            self.force_unique = False

        # If required to reprocess, add to the reprocess list.
        if self.force_reprocess:
            reprocess.append(
                ProcessList(
                    items=[item],
                    work_over=self.force_reprocess,
                    only_on_match=self.only_on_match
                )
            )

        # That's all folks
        self.matched = True
        return self.matched, reprocess


class Constraint:
    """Constraint that is made up of SimpleConstraints

    Parameters
    ----------
    init : object or [object[,...]]
        A single object or list of objects where the
        objects are as follows.
        - SimpleConstraint or subclass
        - Constraint

    reduce : function
        A reduction function with signature `x(iterable)`
        where `iterable` is the `components` list. Returns
        boolean indicating state of the components.
        Default value is `Constraint.all`

    name : str or None
        Optional name for constraint.

    reprocess_on_match : bool
        Reprocess the item if the constraint is satisfied.

    reprocess_on_fail : bool
        Reprocess the item if the constraint is not satisfied.

    work_over : ProcessList.[BOTH, EXISTING, RULES]
        The condition on which this constraint should operate.

    reprocess_rules : [rule[,..]] or None
        List of rules to be applied to.
        If None, calling function will determine the ruleset.
        If empty, [], all rules will be used.

    Attributes
    ----------
    constraints : [Constraint[,...]]
        List of `Constraint` or `SimpleConstaint` that
        make this constraint.

    matched : bool
        Result of the last `check_and_set`

    reduce : function
        A reduction function with signature `x(iterable)`
        where `iterable` is the `components` list. Returns
        boolean indicating state of the components.
        Predefined functions are:
        - `all`: True if all components return True
        - `any`: True if any component returns True

    Notes
    -----
    Named constraints can be accessed directly through indexing:

    >>> c = Constraint(SimpleConstraint(name='simple', value='a_value'))
    >>> c['simple']  # doctest: +SKIP
    SimpleConstraint({'sources': <function SimpleConstraint.__init__.<locals>.<lambda> at 0x7f8be05f5730>,
                      'force_unique': True,
                      'test': <bound method SimpleConstraint.eq of SimpleConstraint({...})>,
                      'reprocess_on_match': False,
                      'reprocess_on_fail': False,
                      'work_over': 1,
                      'reprocess_rules': None,
                      'value': 'a_value',
                      'name': 'simple',
                      'matched': False})
    """
    def __init__(
            self,
            init=None,
            reduce=None,
            name=None,
            reprocess_on_match=False,
            reprocess_on_fail=False,
            work_over=ProcessList.BOTH,
            reprocess_rules=None
    ):
        self.constraints = []

        # Initialize from named parameters
        self.reduce = reduce
        self.name = name
        self.reprocess_on_match = reprocess_on_match
        self.reprocess_on_fail = reprocess_on_fail
        self.work_over = work_over
        self.reprocess_rules = reprocess_rules

        # Initialize from a structure.
        if init is None:
            pass
        elif isinstance(init, list):
            self.constraints = init
        elif isinstance(init, Constraint):
            self.reduce = init.reduce
            self.name = init.name
            self.reprocess_on_match = init.reprocess_on_match
            self.reprocess_on_fail = init.reprocess_on_fail
            self.work_over = init.work_over
            self.reprocess_rules = init.reprocess_rules
            self.constraints = deepcopy(init.constraints)
        elif isinstance(init, SimpleConstraintABC):
            self.constraints = [init]
        else:
            raise TypeError(
                'Invalid initialization value type {}.'
                '\nValid types are `SimpleConstaint`, `Constraint`,'
                '\nor subclass.'.format(type(init))
            )

        # Give some defaults real meaning.
        self.matched = False
        if self.reduce is None:
            self.reduce = self.all

    def check_and_set(self, item, work_over=ProcessList.BOTH):
        """Check and set the constraint

        Returns
        -------
        success, reprocess : bool, [ProcessList[,...]]
            Returns 2-tuple of

                - success : True if check is successful.
                - List of `ProcessList`.
        """
        if work_over not in (self.work_over, ProcessList.BOTH):
            return False, []

        # Do we have positive?
        self.matched, reprocess = self.reduce(item, self.constraints)

        # Determine reprocessing
        if (self.matched and self.reprocess_on_match) or \
           (not self.matched and self.reprocess_on_fail):
            reprocess.append([ProcessList(
                items=[item],
                work_over=self.work_over,
                rules=self.reprocess_rules
            )])

        return self.matched, list(chain(*reprocess))

    def append(self, constraint):
        """Append a new constraint"""
        self.constraints.append(constraint)

    def copy(self):
        """Copy ourselves"""
        return deepcopy(self)

    @staticmethod
    def all(item, constraints):
        """Return positive only if all results are positive."""

        # If there are no constraints, there is nothing to match.
        # Result is false.
        if len(constraints) == 0:
            return False, []

        # Find all negatives. Note first negative
        # that requires reprocessing and how many
        # negatives do not.
        all_match = True
        negative_reprocess = None
        to_reprocess = []
        for constraint in constraints:
            match, reprocess = constraint.check_and_set(item)

            if match:
                if all_match:
                    to_reprocess.append(reprocess)
            else:
                all_match = False

                # If not match and no reprocessing, then fail
                # completely. However, if there is reprocessing, take
                # the first one. Continue to check to ensure
                # there is no further complete fail.
                if len(reprocess) == 0:
                    negative_reprocess = None
                    break
                elif negative_reprocess is None:
                    negative_reprocess = [reprocess]

        if not all_match:
            if negative_reprocess is not None:
                to_reprocess = negative_reprocess
            else:
                to_reprocess = []

        return all_match, to_reprocess

    @staticmethod
    def any(item, constraints):
        """Return the first successful constraint."""
        # If there are no constraints, there is nothing to match.
        # Result is false.
        if len(constraints) == 0:
            return False, []

        to_reprocess = []
        for constraint in constraints:
            match, reprocess = constraint.check_and_set(item)
            if match:
                to_reprocess = [reprocess]
                break
            to_reprocess.append(reprocess)
        return match, to_reprocess

    @staticmethod
    def notany(item, constraints):
        """True if none of the constraints match"""
        match, to_reprocess = Constraint.any(item, constraints)
        return not match, to_reprocess

    @property
    def dup_names(self): # -> dict[str, list[typing.Union[SimpleConstraint, Constraint]]]:
        """Return dictionary of constraints with duplicate names

        This method is meant to be overridden by classes
        that need to traverse a list of constraints.

        Returns
        -------
        dups : {str: [constraint[,...]][,...]}
            Returns a mapping between the duplicated name
            and all the constraints that define that name.
        """
        attrs = self.get_all_attr('name')
        constraints, names = zip(*attrs)
        dups = [name for name, count in collections.Counter(names).items() if count > 1]
        result = collections.defaultdict(list)
        for name, constraint in zip(names, constraints):
            if name in dups:
                result[name].append(constraint)

        # Turn off the defaultdict factory.
        result.default_factory = None
        return result

    def get_all_attr(self, attribute: str): # -> list[tuple[typing.Union[SimpleConstraint, Constraint], typing.Any]]:
        """Return the specified attribute

        This method is meant to be overridden by classes
        that need to traverse a list of constraints.

        Parameters
        ----------
        attribute : str
            The attribute to retrieve

        Returns
        -------
        result : [(SimpleConstraint or Constraint, object)[,...]]
            The list of values of the attribute in a tuple. If there is no attribute,
            an empty tuple is returned.

        Raises
        ------
        AttributeError
            If the attribute is not found.
        """
        result = []
        value = getattr(self, attribute)
        if value is not None:
            result = [(self, value)]
        for constraint in self.constraints:
            result.extend(constraint.get_all_attr(attribute))

        return result

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
        raise NotImplementedError('Cannot set constraints by index.')

    def __delitem__(self, key):
        """Not implemented"""
        raise NotImplementedError('Cannot delete a constraint by index.')

    def __repr__(self):
        result = '{}(name={}).{}([{}])'.format(
            self.__class__.__name__,
            str(getattr(self, 'name', None)),
            str(self.reduce.__name__),
            ''.join([
                repr(constraint)
                for constraint in self.constraints
            ])
        )
        return result

    def __str__(self):
        result = '\n'.join([
            str(constraint)
            for constraint in self
            if constraint.name is not None
        ])
        return result


# ---------
# Utilities
# ---------
def meets_conditions(value, conditions):
    """Check whether value meets any of the provided conditions

    Parameters
    ----------
    values : str
        The value to be check with.

    condition : regex,
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

"""Constraints - use these to define the rules governing association candidate types."""

import abc
import collections
from copy import deepcopy
from itertools import chain
import logging
import re

from .process_list import ListCategory, ProcessList
from .utilities import evaluate, getattr_from_list, is_iterable
from jwst.associations.pool import PoolRow

__all__ = [
    "AttrConstraint",
    "Constraint",
    "ConstraintTrue",
    "SimpleConstraint",
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class SimpleConstraintABC(abc.ABC):
    """
    Simple Constraint ABC.

    Parameters
    ----------
    init : dict
        Dictionary where the key:value pairs define
        the following parameters

    value : object or None
        Value that must be matched.

    name : str or None
        Option name for constraint

    **kwargs : key:value pairs
        Other initialization parameters

    Attributes
    ----------
    found_values : set(str[,...])
        Set of actual found values for this condition. True SimpleConstraints
        do not normally set this; the value is not different than `value`.

    matched : bool
        Last call to `check_and_set`
    """

    # Attributes to show in the string representation.
    _str_attrs: tuple = ("name", "value")

    def __new__(cls, *args, **kwargs):  # noqa: ARG003
        """
        Force creation of the constraint attribute dict before anything else.

        Returns
        -------
        ~jwst.associations.lib.constraint.SimpleConstraintABC
            New instance of class.
        """
        obj = super().__new__(cls)
        obj._ca_history = collections.deque()  # noqa: SLF001
        obj._constraint_attributes = {}  # noqa: SLF001
        return obj

    def __init__(self, init=None, value=None, name=None, **kwargs):
        # Defined attributes
        self.value = value
        self.name = name
        self.matched = False
        self.found_values = set()

        if init is not None:
            self._constraint_attributes.update(init)
        else:
            self._constraint_attributes.update(kwargs)

    def __getattr__(self, name):
        """
        Retrieve user defined attribute.

        Returns
        -------
        any
            Attribute corresponding to provided name.
        """
        if name.startswith("_"):
            return super().__getattribute__(name)
        if name in self._constraint_attributes:
            return self._constraint_attributes[name]
        raise AttributeError(f"No such attribute {name}")

    def __setattr__(self, name, value):
        """Store all attributes in the user dictionary."""
        if not name.startswith("_"):
            self._constraint_attributes[name] = value
        else:
            object.__setattr__(self, name, value)

    @abc.abstractmethod
    def check_and_set(self, item):
        """
        Check and set the constraint.

        Returns
        -------
        success, reprocess : bool, [~jwst.associations.ProcessList[,...]]
            Returns 2-tuple of

                - True if check is successful.
                - List of `~jwst.associations.ProcessList`.
        """
        self.matched = True
        self.found_values.add(self.value)
        return self.matched, []

    @property
    def dup_names(self):
        """
        Return dictionary of constraints with duplicate names.

        This method is meant to be overridden by classes
        that need to traverse a list of constraints.

        Returns
        -------
        dups : {str: [constraint[,...]][,...]}
            Returns a mapping between the duplicated name
            and all the constraints that define that name.
        """
        return {}

    @property
    def id(self):
        """
        Return identifier for the constraint.

        Returns
        -------
        id : str
            The identifier
        """
        return f"{self.__class__.__name__}:{self.name}"

    def copy(self):
        """
        Copy self.

        Returns
        -------
        object
            Deepcopy of self.
        """
        return deepcopy(self)

    def get_all_attr(self, attribute, name=None):
        """
        Return the specified attribute.

        This method exists solely to support `Constraint.get_all_attr`.
        This obviates the need for class/method checking.

        Parameters
        ----------
        attribute : str
            The attribute to retrieve

        name : str or None
            Only return attribute if the name of the current constraint
            matches the requested named constraints. If None, always
            return value.

        Returns
        -------
        [(self, value)] : [(SimpleConstraint, object)]
            The value of the attribute in a tuple. If there is no attribute,
            an empty tuple is returned.
        """
        if name is None or name == self.name:
            value = getattr(self, attribute, None)
            if value is not None:
                if not isinstance(value, (list, set)) or len(value):
                    return [(self, value)]
        return []

    def restore(self):
        """Restore constraint state."""
        try:
            self._constraint_attributes = self._ca_history.pop()
        except IndexError:
            logger.debug("No more attribute history to restore from. restore is a NOOP")

    def preserve(self):
        """Save the current state of the constraints."""
        ca_copy = self._constraint_attributes.copy()
        ca_copy["found_values"] = self._constraint_attributes["found_values"].copy()
        self._ca_history.append(ca_copy)

    # Make iterable to work with `Constraint`.
    # Since this is a leaf, simple return ourselves.
    def __iter__(self):
        yield self

    def __repr__(self):
        result = f"{self.__class__.__name__}({str(self._constraint_attributes)})"
        return result

    def __str__(self):
        result = (
            f"{self.__class__.__name__}("
            f"{ ({str_attr: getattr(self, str_attr) for str_attr in self._str_attrs}) })"
        )
        return result


class ConstraintTrue(SimpleConstraintABC):
    """Always return True."""

    def check_and_set(self, item):
        """
        Empty method to return True from parent abstract method.

        Returns
        -------
        bool
            True from parent abstract method.
        """
        return super(ConstraintTrue, self).check_and_set(item)


class SimpleConstraint(SimpleConstraintABC):
    """
    A basic constraint.

    Examples
    --------
    Create a constraint where the attribute `attr` of an object
    matches the value `my_value`:

    >>> c = SimpleConstraint(value="my_value")
    >>> print(c)
    SimpleConstraint({'name': None, 'value': 'my_value'})

    To check a constraint, call `check_and_set`. A successful match
    will return a tuple of `True` and a reprocess list.
    >>> item = "my_value"
    >>> c.check_and_set(item)
    (True, [])

    If it doesn't match, `False` will be returned.
    >>> bad_item = "not_my_value"
    >>> c.check_and_set(bad_item)
    (False, [])

    A `SimpleConstraint` can also be initialized by a `dict`
    of the relevant parameters:
    >>> init = {"value": "my_value"}
    >>> c = SimpleConstraint(init)
    >>> print(c)
    SimpleConstraint({'name': None, 'value': 'my_value'})

    If the value to check is `None`, the `SimpleConstraint` will
    successfully match whatever object given. However, a new `SimpleConstraint`
    will be returned where the `value` is now set to whatever the attribute
    was of the object.
    >>> c = SimpleConstraint(value=None)
    >>> matched, reprocess = c.check_and_set(item)
    >>> print(c)
    SimpleConstraint({'name': None, 'value': 'my_value'})

    This behavior can be overridden by the `force_unique` parameter:
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
        work_over=ListCategory.BOTH,
        reprocess_rules=None,
        **kwargs,
    ):
        """
        Initialize a new SimpleConstraint.

        Parameters
        ----------
        init : dict
            Dictionary where the key:value pairs define
            the following parameters.

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

        reprocess_on_match : bool
            Reprocess the item if the constraint is satisfied.

        reprocess_on_fail : bool
            Reprocess the item if the constraint is not satisfied.

        work_over : ListCategory.[BOTH, EXISTING, RULES]
            The condition on which this constraint should operate.

        reprocess_rules : [rule[,..]] or None
            List of rules to be applied to.
            If None, calling function will determine the ruleset.
            If empty, [], all rules will be used.
        """
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
        """
        Check and set the constraint.

        Returns
        -------
        success, reprocess : bool, [~jwst.associations.ProcessList[,...]]
            Returns 2-tuple of

                - True if check is successful.
                - List of `~jwst.associations.ProcessList`.
        """
        source_value = self.sources(item)

        satisfied = True
        if self.value is not None:
            satisfied = self.test(self.value, source_value)
        self.matched = satisfied

        if self.matched:
            if self.force_unique:
                self.value = source_value
            self.found_values.add(self.value)

        # Determine reprocessing
        reprocess = []
        if (self.matched and self.reprocess_on_match) or (
            not self.matched and self.reprocess_on_fail
        ):
            reprocess.append(
                ProcessList(
                    items=[item],
                    work_over=self.work_over,
                    rules=self.reprocess_rules,
                    trigger_constraints=[self.id],
                )
            )

        return self.matched, reprocess

    def eq(self, value1, value2):
        """
        Check if constraint.value and item are equal.

        Parameters
        ----------
        value1 : any
            The first value to compare.
        value2 : any
            The second value to compare.

        Returns
        -------
        bool
            True if the two values are deemed equal.
        """
        return value1 == value2


class AttrConstraint(SimpleConstraintABC):
    """
    Test attribute of an item.

    Attributes
    ----------
    found_values : set(str[,...])
        Set of actual found values for this condition.
    matched : bool
        Last result of `check_and_set`
    """

    # Attributes to show in the string representation.
    _str_attrs = ("name", "sources", "value")

    def __init__(
        self,
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
        **kwargs,
    ):
        """
        Initialize a new AttrConstraint.

        Parameters
        ----------
        sources : [str[,...]]
            List of attributes to query
        evaluate : bool
            Evaluate the item's value before checking condition.
        force_reprocess : ListCategory.state or False
            Add item back onto the reprocess list using
            the specified `~jwst.associations.ProcessList` work over state.
        force_unique : bool
            If the initial value is `None` or a list of possible values,
            the constraint will be modified to be the value first matched.
        invalid_values : [str[,...]]
            List of values that are invalid in an item.
            Will cause a non-match.
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
        """
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
        super().__init__(init=init, **kwargs)

        # Give some defaults real meaning.
        if invalid_values is None:
            self.invalid_values = []
        if onlyif is None:
            self.onlyif = lambda _item: True

        # Haven't actually matched anything yet.
        self.found_values = set()
        self.matched = False

    def check_and_set(self, item):
        """
        Check and set constraints based on item.

        Parameters
        ----------
        item : dict
            The item to check on.

        Returns
        -------
        success, reprocess : bool, [~jwst.associations.ProcessList[,...]]
            Returns 2-tuple of

                - True if check is successful.
                - List of `~jwst.associations.ProcessList`.
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
                        trigger_constraints=[self.id],
                    )
                )
            self.matched = True
            return self.matched, reprocess

        # Get the condition information.
        try:
            source, value = getattr_from_list(
                item, self.sources, invalid_values=self.invalid_values
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

        evaled = value
        if self.evaluate:
            evaled = evaluate(value)

        # If the constraint has no value to check against, and given
        # value evaluates to a list, the item must be duplicated,
        # with each value from its list, and all the new items reprocessed.
        # Otherwise, the value is the value to set the constraint by.
        if self.value is None:
            if is_iterable(evaled):
                reprocess.append(reprocess_multivalue(item, source, evaled, self))
                self.matched = False
                return self.matched, reprocess
            value = str(evaled)

        # Else, the constraint does have a value. Check against it.
        else:
            if callable(self.value):
                match_value = self.value()
            else:
                match_value = self.value
            if not is_iterable(evaled):
                evaled = [evaled]
            for evaled_item in evaled:
                value = str(evaled_item)
                if meets_conditions(value, match_value):
                    break
            else:
                # The condition is not matched, leave now.
                self.matched = False
                return self.matched, reprocess

            # A match was found. If there is a list of potential values,
            # set them up for reprocessing.
            next_evaleds = [next_evaled for next_evaled in evaled if next_evaled != evaled_item]
            if next_evaleds:
                reprocess.append(reprocess_multivalue(item, source, next_evaleds, self))

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
                    only_on_match=self.only_on_match,
                    trigger_constraints=[self.id],
                )
            )

        # That's all folks
        self.matched = True
        return self.matched, reprocess


class Constraint:
    """
    Constraint that is made up of SimpleConstraints.

    Attributes
    ----------
    constraints : [Constraint[,...]]
        List of `Constraint` or `SimpleConstraint` that
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

    >>> c = Constraint(SimpleConstraint(name="simple", value="a_value"))
    >>> c["simple"]  # doctest: +SKIP
    SimpleConstraint({'sources': <function SimpleConstraint.__init__.<locals>.<lambda>,
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
        work_over=ListCategory.BOTH,
        reprocess_rules=None,
    ):
        """
        Initialize a new Constraint.

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

        work_over : ListCategory.[BOTH, EXISTING, RULES]
            The condition on which this constraint should operate.

        reprocess_rules : [rule[,..]] or None
            List of rules to be applied to.
            If None, calling function will determine the ruleset.
            If empty, [], all rules will be used.
        """
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
                f"Invalid initialization value type {type(init)}."
                "\nValid types are `SimpleConstraint`, `Constraint`,"
                "\nor subclass."
            )

        # Give some defaults real meaning.
        self.matched = False
        if self.reduce is None:
            self.reduce = self.all

    @property
    def dup_names(self):
        """
        Return dictionary of constraints with duplicate names.

        This method is meant to be overridden by classes
        that need to traverse a list of constraints.

        Returns
        -------
        dups : {str: [constraint[,...]][,...]}
            Returns a mapping between the duplicated name
            and all the constraints that define that name.
        """
        attrs = self.get_all_attr("name")
        constraints, names = zip(*attrs, strict=True)
        dups = [name for name, count in collections.Counter(names).items() if count > 1]
        result = collections.defaultdict(list)
        for name, constraint in zip(names, constraints, strict=True):
            if name in dups:
                result[name].append(constraint)

        # Turn off the defaultdict factory.
        result.default_factory = None
        return result

    @property
    def id(self):
        """
        Return identifier for the constraint.

        Returns
        -------
        id : str
            The identifier
        """
        return f"{self.__class__.__name__}:{self.name}"

    def append(self, constraint):
        """Append a new constraint."""
        self.constraints.append(constraint)

    def check_and_set(self, item, work_over=ListCategory.BOTH):
        """
        Check and set the constraint.

        Returns
        -------
        success, reprocess : bool, [~jwst.associations.ProcessList[,...]]
            Returns 2-tuple of

                - success : True if check is successful.
                - List of `~jwst.associations.ProcessList`.
        """
        if work_over not in (self.work_over, ListCategory.BOTH):
            return False, []

        # Do we have positive?
        self.matched, reprocess = self.reduce(item, self.constraints)

        # Determine reprocessing
        if (self.matched and self.reprocess_on_match) or (
            not self.matched and self.reprocess_on_fail
        ):
            reprocess.append(
                [
                    ProcessList(
                        items=[item],
                        work_over=self.work_over,
                        rules=self.reprocess_rules,
                        trigger_constraints=[self.id],
                    )
                ]
            )

        return self.matched, list(chain(*reprocess))

    def copy(self):
        """
        Copy ourselves.

        Returns
        -------
        object
            Deepcopy of self.
        """
        return deepcopy(self)

    def get_all_attr(self, attribute, name=None):
        """
        Return the specified attribute for specified constraints.

        Parameters
        ----------
        attribute : str
            The attribute to retrieve

        name : str or None
            Only return attribute if the name of the current constraint
            matches the requested named constraints. If None, always
            return value.

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
        if name is None or name == self.name:
            value = getattr(self, attribute, None)
            if value is not None:
                result = [(self, value)]
        for constraint in self.constraints:
            result.extend(constraint.get_all_attr(attribute, name=name))

        return result

    def preserve(self):
        """Preserve all constraint states."""
        for constraint in self.constraints:
            constraint.preserve()

    def restore(self):
        """Restore all constraint states."""
        for constraint in self.constraints:
            constraint.restore()

    @staticmethod
    def all(item, constraints):
        """
        Return positive only if all results are positive.

        Parameters
        ----------
        item : ACID
            The candidate.
        constraints : list[Constraint, ...]
            The list of constraints to check.

        Returns
        -------
        bool, list(Constraint, ...) or None
            True if all constraints positive, with empty list.
            If no constraints, False and empty list. Otherwise
            False with list of constraints to reprocess.
        """
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
        """
        Return the first successful constraint.

        Parameters
        ----------
        item : ACID
            The candidate.
        constraints : list[Constraint, ...]
            The list of constraints to check.

        Returns
        -------
        bool, list(Constraint, ...) or None
            False, [] if no match or constraints to reprocess.
            True, list(Constraints) if match found, and any constraints
            to reprocess listed.
        """
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
        """
        Check if none of the constraints match; true if none do.

        Parameters
        ----------
        item : ACID
            The candidate.
        constraints : list[Constraint, ...]
            The list of constraints to check.

        Returns
        -------
        bool
            True if none of the constraints match.
        """
        match, to_reprocess = Constraint.any(item, constraints)
        return not match, to_reprocess

    @staticmethod
    def notall(item, constraints):
        """
        Check if not all of the constraints match; true if not all do.

        Parameters
        ----------
        item : ACID
            The candidate.
        constraints : list[Constraint, ...]
            The list of constraints to check.

        Returns
        -------
        bool
            True if not all constraints match.
        """
        match, to_reprocess = Constraint.all(item, constraints)
        return not match, to_reprocess

    def __delitem__(self, key):
        """Not implemented."""
        raise NotImplementedError("Cannot delete a constraint by index.")

    # Make iterable
    def __iter__(self):
        yield from chain(*map(iter, self.constraints))

    # Index implementation
    def __getitem__(self, key):
        """
        Retrieve a named constraint.

        Parameters
        ----------
        key : str
            The key to retrieve a value with.

        Returns
        -------
        jwst.associations.lib.constraint.Constraint
            The constraint to be retrieved.
        """
        for constraint in self.constraints:
            name = getattr(constraint, "name", None)
            if name is not None and name == key:
                return constraint
            try:
                found = constraint[key]
            except (KeyError, TypeError):
                pass
            else:
                return found
        raise KeyError(f"Constraint {key} not found")

    def __repr__(self):
        result = "{}(name={}).{}([{}])".format(
            self.__class__.__name__,
            str(getattr(self, "name", None)),
            str(self.reduce.__name__),
            "".join([repr(constraint) for constraint in self.constraints]),
        )
        return result

    def __setitem__(self, key, value):
        """Not implemented."""
        raise NotImplementedError("Cannot set constraints by index.")

    def __str__(self):
        result = "\n".join([str(constraint) for constraint in self if constraint.name is not None])
        return result


# Utilities
def meets_conditions(value, conditions):
    """
    Check whether value meets any of the provided conditions.

    Parameters
    ----------
    value : str
        The value to be check with.
    conditions : regex,
        Regular expressions to match against.

    Returns
    -------
    bool
        True if any condition is meant.
    """
    if not is_iterable(conditions):
        conditions = [conditions]
    for condition in conditions:
        condition = "".join(["^", condition, "$"])
        match = re.match(condition, value, flags=re.IGNORECASE)
        if match:
            return True
    return False


def reprocess_multivalue(item, source, values, constraint):
    """
    Complete reprocessing of items that have a list of values.

    Parameters
    ----------
    item : dict
        The item.
    source : str
        The attribute which has the multi-values.
    values : list
        The list of values
    constraint : Constraint
        The constraint which is triggering the reprocessing.

    Returns
    -------
    process_list : ProcessList
        The process list to put on the reprocess queue
    """
    reprocess_items = []
    for value in values:
        new_item = PoolRow(item)
        new_item[source] = str(value)
        reprocess_items.append(new_item)
    process_list = ProcessList(items=reprocess_items, trigger_constraints=[constraint.id])
    return process_list

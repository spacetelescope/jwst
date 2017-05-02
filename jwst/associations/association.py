from ast import literal_eval
from collections import MutableMapping
from copy import deepcopy
from datetime import datetime
import json
import jsonschema
import logging
from nose.tools import nottest
import re

from astropy.extern import six
from numpy.ma import masked

from . import __version__
from .exceptions import (
    AssociationNotValidError
)
from .lib.ioregistry import IORegistry

__all__ = ['Association']


# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Timestamp template
_TIMESTAMP_TEMPLATE = '%Y%m%dt%H%M%S'


class Association(MutableMapping):
    """Association Base Class

    Parameters
    ----------
    version_id: str or None
        Version_Id to use in the name of this association.
        If None, nothing is added.

    Raises
    ------
    AssociationError
        If a member doesn't match any of the registered associations.

    Attributes
    ----------
    instance: dict-like
        The instance is the association data structure.
        See `data` below

    meta: dict
        Information about the association.

    data: dict
        The association. The format of this data structure
        is determined by the individual assocations and, if
        defined, valided against their specified schema.

    schema_file: str
        The name of the output schema that an association
        must adhere to.

    registry: AssociationRegistry
        The registry this association came from.

    asn_name: str
        The suggested file name of association

    asn_rule: str
        The name of the rule
    """

    # Assume no registry
    registry = None

    # Default force a constraint to use first value.
    DEFAULT_FORCE_UNIQUE = False

    # Default require that the constraint exists or otherwise
    # can be explicitly checked.
    DEFAULT_REQUIRE_CONSTRAINT = True

    # Global constraints
    GLOBAL_CONSTRAINTS = {}

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = None

    # Initialize a global IO registry
    ioregistry = IORegistry()

    def __init__(
            self,
            version_id=None,
    ):

        self.data = dict()
        self.add_constraints(deepcopy(self.GLOBAL_CONSTRAINTS))
        self.run_init_hook = True
        self.meta = {}

        self.version_id = version_id

        self.data.update({
            'asn_type': 'None',
            'asn_rule': self.asn_rule,
            'version_id': self.version_id,
            'code_version': __version__,
        })

    @classmethod
    def create(cls, member, version_id=None):
        """Create association if member belongs

        Parameters
        ----------
        member: dict
            The member to initialize the association with.

        version_id: str or None
            Version_Id to use in the name of this association.
            If None, nothing is added.

        Returns
        -------
        (association, reprocess_list)
            2-tuple consisting of:
            - association: The association or, if the member does not
                this rule, None
            - [ProcessList[, ...]]: List of members to process again.
        """
        asn = cls(version_id=version_id)
        matches, reprocess = asn.add(member)
        if not matches:
            return None, reprocess
        return asn, reprocess

    @property
    def asn_name(self):
        return 'unamed_association'

    @classmethod
    def _asn_rule(cls):
        return cls.__name__

    @property
    def asn_rule(self):
        return self._asn_rule()

    @classmethod
    def validate(cls, asn):
        """Validate an association against this rule

        Parameters
        ----------
        asn: Association or association-like
            The association structure to examine

        Returns
        -------
        valid: bool
            True if valid. Otherwise the `AssociationNotValidError` is raised

        Raises
        ------
        AssociationNotValidError
            If there is some reason validation failed.

        Notes
        -----
        The base method checks against the rule class' schema
        If the rule class does not define a schema, a warning is issued
        but the routine will return True.
        """
        if not hasattr(cls, 'schema_file'):
            logger.warning(
                'Cannot validate: {} has no schema. Presuming OK.'.format(cls)
            )
            return True

        if isinstance(asn, cls):
            asn_data = asn.data
        else:
            asn_data = asn

        with open(cls.schema_file, 'r') as schema_file:
            asn_schema = json.load(schema_file)

        try:
            jsonschema.validate(asn_data, asn_schema)
        except (AttributeError, jsonschema.ValidationError) as err:
            raise AssociationNotValidError('Validation failed')
        return True

    def dump(self, format='json', **kwargs):
        """Serialize the association

        Parameters
        ----------
        format: str
            The format to use to dump the association into.

        kwargs: dict
            List of arguments to pass to the registered
            routines for the current association type.

        Returns
        -------
        (name, serialized):
            Tuple where the first item is the suggested
            base name for the file.
            Second item is the serialization.

        Raises
        ------
        AssociationError
            If the operation cannot be done

        AssociationNotValidError
            If the given association does not validate.
        """
        if self.is_valid:
            return self.ioregistry[format].dump(self, **kwargs)
        else:
            raise AssociationNotValidError(
                'Association {} is not valid'.format(self)
            )

    @classmethod
    def load(
            cls,
            serialized,
            format=None,
            validate=True,
            **kwargs
    ):
        """Marshall a previously serialized association

        Parameters
        ----------
        serialized: object
            The serialized form of the association.

        format: str or None
            The format to force. If None, try all available.

        validate: bool
            Validate against the class' defined schema, if any.

        kwargs: dict
            Other arguments to pass to the `load` method

        Returns
        -------
        The Association object

        Raises
        ------
        AssociationNotValidError
            Cannot create or validate the association.

        Notes
        -----
        The `serialized` object can be in any format
        supported by the registered I/O routines. For example, for
        `json` and `yaml` formats, the input can be either a string or
        a file object containing the string.
        """
        if format is None:
            formats = [
                format
                for format_name, format in cls.ioregistry.items()
            ]
        else:
            formats = [cls.ioregistry[format]]

        for format in formats:
            try:
                asn = format.load(
                    cls, serialized, **kwargs
                )
            except AssociationNotValidError:
                continue
            else:
                break
        else:
            raise AssociationNotValidError(
                'Cannot translate "{}" to an association'.format(serialized)
            )

        # Validate
        if validate:
            cls.validate(asn)

        return asn

    @property
    def is_valid(self):
        """Check if association is valid"""
        try:
            self.__class__.validate(self)
        except AssociationNotValidError:
            return False
        return True

    def add(self, member, check_constraints=True):
        """Add the member to the association

        Parameters
        ----------
        member: dict
            The member to add.

        check_constraints: bool
            If True, see if the member should belong to this association.
            If False, just add it.

        Returns
        -------
        (matches, reprocess_list)
            2-tuple consisting of:
            - bool: True if the all constraints are satisfied
            - [ProcessList[, ...]]: List of members to process again.
        """
        if check_constraints:
            matches, reprocess = self.test_and_set_constraints(member)

        if matches:
            if self.run_init_hook:
                self._init_hook(member)
            self._add(member)
            self.run_init_hook = False

        return matches, reprocess

    @nottest
    def test_and_set_constraints(self, member):
        """Test whether the given dictionaries match parameters for
        for this association

        Parameters
        ----------
        member: dict
            The parameters to check/set for this association.
            This can be a list of dictionaries.

        Returns
        -------
        (matches, reprocess_list)
            2-tuple consisting of:
            - bool: True if the all constraints are satisfied
            - [ProcessMember[, ...]]: List of members to process again.

        Notes
        -----
        If a constraint is present, but does not have a value,
        that constraint is set, and, by definition, matches.
        """

        # If there are no conditions, the default is a
        # successful test
        matches = True

        reprocess = []
        constraints = deepcopy(self.constraints)
        for constraint, conditions in constraints.items():
            test = conditions.get('test', self.match_member)
            matches, new_reprocess = test(member, constraint, conditions)
            reprocess.extend(new_reprocess)
            if not matches:
                break
        else:
            self.constraints = constraints
        return matches, reprocess

    def match_member(self, member, constraint, conditions):
        """Use member info to match to the conditions

        Parameters
        ----------
        member: dict
            The member to retrieve the values from

        constraint: str
            The name of the constraint

        conditions: dict
            The conditions structure

        Returns
        -------
        (matches, reprocess_list)
            2-tuple consisting of:
            - bool: True if the all constraints are satisfied
            - [ProcessList[, ...]]: List of members to process again.
        """
        reprocess = []

        # Get the condition information.
        try:
            input, value = getattr_from_list(
                member,
                conditions['inputs'],
                invalid_values=self.INVALID_VALUES
            )
        except KeyError:
            if not conditions.get('force_undefined', False) and \
               conditions.get(
                   'required',
                   self.DEFAULT_REQUIRE_CONSTRAINT
               ):
                return False, reprocess
            else:
                return True, reprocess
        else:
            if conditions.get('force_undefined', False):
                return False, reprocess

        # If the value is a list, build the reprocess list
        evaled = evaluate(value)
        current_value = None
        if is_iterable(evaled):
            process_members = []
            for avalue in evaled:
                if not current_value:
                    current_value = avalue
                    continue
                new_member = deepcopy(member)
                new_member[input] = str(avalue)
                process_members.append(new_member)
            if len(process_members) > 0:
                reprocess.append(ProcessList(
                    process_members,
                    [type(self)]
                ))
            evaled = current_value

        evaled_str = str(evaled)
        if conditions['value'] is not None:
            if not meets_conditions(
                    evaled_str, conditions['value']
            ):
                logging.debug(
                    'Constraint {}'
                    ' does not match association.'.format(constraint)
                )
                return False, reprocess

        # At this point, the constraint has passed.
        # Fix the conditions.
        if conditions['value'] is None or \
           conditions.get('force_unique', self.DEFAULT_FORCE_UNIQUE):
            conditions['value'] = re.escape(evaled_str)
            conditions['inputs'] = [input]
            conditions['force_unique'] = False

        # That's all folks
        return True, reprocess

    def match_constraint(self, member, constraint, conditions):
        """Generic constraint checking

        Parameters
        ----------
        member: dict
            The member to retrieve the values from

        constraint: str
            The name of the constraint

        conditions: dict
            The conditions structure

        Returns
        -------
        (matches, reprocess_list)
            2-tuple consisting of:
            - bool: True if the all constraints are satisfied
            - [ProcessList[, ...]]: List of members to process again.
        """
        reprocess = []
        evaled_str = conditions['inputs'](member)
        if conditions['value'] is not None:
            if not meets_conditions(
                    evaled_str, conditions['value']
            ):
                logging.debug(
                    'Constraint {}'
                    ' does not match association.'.format(constraint)
                )
                return False, reprocess

        # At this point, the constraint has passed.
        # Fix the conditions.
        if conditions['value'] is None or \
           conditions.get('force_unique', self.DEFAULT_FORCE_UNIQUE):
            conditions['value'] = re.escape(evaled_str)
            conditions['force_unique'] = False

        # That's all folks
        return True, reprocess

    def add_constraints(self, new_constraints):
        """Add a set of constraints to the current constraints."""

        try:
            constraints = self.constraints
        except AttributeError:
            constraints = {}
            self.constraints = constraints
        for constraint, value in six.iteritems(new_constraints):
            constraints[constraint] = constraints.get(constraint, value)

    def constraints_to_text(self):
        yield 'Constraints:'
        for name, conditions in self.constraints.items():
            if conditions.get('force_undefined', False):
                yield '    {:s}: Is Invalid'.format(name)
            else:
                yield '    {:s}: {}'.format(name, conditions['value'])

    def _init_hook(self, member):
        """Post-check and pre-member-adding initialization."""
        pass

    def _add(self, member):
        """Add a member, association-specific"""
        raise NotImplementedError(
            'Association._add must be implemented by a specific assocation rule.'
        )

    def _add_items(self, items, **kwargs):
        """ Force adding items to the association

        Parameters
        ----------
        items: [object[, ...]]
            A list of items to make members of the association.

        Notes
        -----
        This is a low-level shortcut into adding members, such as file names,
        to an association. All defined shortcuts and other initializations are
        by-passed, resulting in a potentially unusable association.
        """
        try:
            self['members'].update(items)
        except KeyError:
            self['members'] = items

    # #################################################
    # Methods required for implementing MutableMapping
    # #################################################
    def __getitem__(self, key):
        return self.data[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.data[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.data[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __keytransform__(self, key):
        return key

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()

    def values(self):
        return self.data.values()


class ProcessList(object):
    """A Process list

    Parameters
    ----------
    members: [member[, ...]]
        The list of members to process

    rules: [Association[, ...]]
        List of rules to process the members against.
    """

    def __init__(self, members, rules):
        self.members = members
        self.rules = rules


# #########
# Utilities
# #########
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


def make_timestamp():
    timestamp = datetime.utcnow().strftime(
        _TIMESTAMP_TEMPLATE
    )
    return timestamp


def getattr_from_list(adict, attributes, invalid_values=None):
    """Retrieve value from dict using a list of attributes

    Parameters
    ----------
    adict: dict
        dict to retrieve from

    attributes: list
        List of attributes

    invalid_values: set
        A set of values that essentially mean the
        attribute does not exist.

    Returns
    -------
    (attribute, value)
        Returns the value and the attribute from
        which the value was taken.

    Raises
    ------
    KeyError
        None of the attributes are found in the dict.
    """
    if invalid_values is None:
        invalid_values = set()

    for attribute in attributes:
        try:
            result = adict[attribute]
        except KeyError:
            continue
        else:
            if result is masked:
                continue
            if result not in invalid_values:
                return attribute, result
            else:
                continue
    else:
        raise KeyError


def evaluate(value):
    """Evaluate a value

    Parameters
    ----------
    value: str
        The string to evaluate.

    Returns
    -------
    type or str
        The evaluation. If the value cannot be
        evaluated, the value is simply returned
    """
    try:
        evaled = literal_eval(value)
    except (ValueError, SyntaxError):
        evaled = value
    return evaled


def is_iterable(obj):
    return not isinstance(obj, six.string_types) and \
        not isinstance(obj, tuple) and \
        hasattr(obj, '__iter__')

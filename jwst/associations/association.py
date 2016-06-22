import logging
from collections import namedtuple
from copy import deepcopy
from datetime import datetime
from inspect import getfile, getmembers, isclass, ismodule
from itertools import count
import json
import jsonschema
from nose.tools import nottest
from os.path import abspath, basename, dirname, expanduser, expandvars, join
import re
import sys

from astropy.extern import six
from numpy.ma import masked


__all__ = [
    'Association',
    'AssociationRegistry',
    'AssociationError',
    'SERIALIZATION_PROTOCOLS',
]


# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Library files
_ASN_RULE = 'association_rules.py'

# User-level association definitions start with...
USER_ASN = 'Asn_'

# Timestamp template
_TIMESTAMP_TEMPLATE = '%Y%m%dT%H%M%S'


class AssociationError(Exception):
    """Basic failure of an association"""


class AssociationRegistry(dict):
    """The available assocations

    Parameters
    ----------
    definition_files: [str,]
        The files to find the association definitions in.

    include_default: bool
        True to include the default definitions.

    global_constraints: dict
        Constraints to be added to each rule.
    """

    def __init__(self,
                 definition_files=None,
                 include_default=True,
                 global_constraints=None):
        super(AssociationRegistry, self).__init__()

        # Setup constraints that are to be applied
        # to every rule.
        if global_constraints is None:
            global_constraints = {}

        if definition_files is None:
            definition_files = []
        if include_default:
            definition_files.insert(0, libpath(_ASN_RULE))
        if len(definition_files) <= 0:
            raise AssociationError('No rule definition files specified.')
        Utility = type('Utility', (object,), {})
        for fname in definition_files:
            logger.debug('Import rules files "{}"'.format(fname))
            module = import_from_file(fname)
            logger.debug('Module="{}"'.format(module))
            for class_name, class_object in get_classes(module):
                logger.debug('class_name="{}"'.format(class_name))
                if class_name.startswith(USER_ASN):
                    class_object.GLOBAL_CONSTRAINTS = global_constraints
                    self.__setitem__(class_name, class_object)
                if class_name == 'Utility':
                    Utility = type('Utility', (class_object, Utility), {})
        self.Utility = Utility

    def match(self, member, timestamp=None, ignore=None):
        """See if member belongs to any of the association defined.

        Parameters
        ----------
        member: dict
            A member, like from a Pool, to find assocations for.

        timestamp: str
            If specified, a string appened to association names.
            Generated if not specified.

        ignore: list
            A list of associations to ignore when looking for a match.
            Intended to ensure that already created associations
            are not re-created.

        Returns
        -------
        [association,]
            A list of associations this member belongs to.
        """
        logger.debug('Starting...')
        associations = []
        for name, rule in self.items():
            if rule not in ignore:
                try:
                    associations.append(rule(member, timestamp))
                except AssociationError as error:
                    logger.debug('Rule "{}" not matched'.format(name))
                    logger.debug('Error="{}"'.format(error))
                    continue
        if len(associations) == 0:
            raise AssociationError('Member does not match any rules.')
        return associations


class Association(object):
    """An Association

    Parameters
    ----------
    member: dict
        The member to initialize the association with.

    timestamp: str
        Timestamp to use in the name of this association. Should conform
        to the datetime.strftime format '%Y%m%dT%M%H%S'. If None, class
        instantiation will create this string using current time.

    Raises
    ------
    AssociationError
        If a member doesn't match any of the registered associations.

    Attributes
    ----------
    meta: dict
        Information about the association.

    data: dict
        The association. The format of this data structure
        is determined by the individual assocations and, if
        defined, valided against their specified schema.

    schema_file: str
        The name of the output schema that an association
        must adhere to.
    """

    # Default force a constraint to use first value.
    DEFAULT_FORCE_UNIQUE = False

    # Default require that the constraint exists or otherwise
    # can be explicitly checked.
    DEFAULT_REQUIRE_CONSTRAINT = True

    # Global constraints
    GLOBAL_CONSTRAINTS = {}

    # Associations of the same type are sequenced.
    _sequence = count(1)

    def __init__(self, member, timestamp=None):

        self.add_constraints(deepcopy(self.GLOBAL_CONSTRAINTS))
        self.test_and_set_constraints(member)

        # Member belongs to us!
        # Continue initialization.
        self.sequence = six.advance_iterator(self._sequence)
        if timestamp is not None:
            self.timestamp = timestamp
        else:
            self.timestamp = make_timestamp()
        self.meta = {}
        self.data = {
            'asn_type': 'None',
            'asn_rule': self.__class__.__name__,
            'creation_time': self.timestamp
        }

        # Peform further initializations before actually
        # adding the member to this association.
        self._init_hook(member)

        self.add(member, check_constraints=False)

    @property
    def asn_name(self):
        return 'unamed_association'

    def dump(self, protocol='json'):
        """Serialize the association

        Parameters
        ----------
        protocol: ('json',)
            The format to use for serialization.

        Returns
        -------
        (name, serialized):
            Tuple where the first item is the suggested
            base name for the file.
            Second item is the serialization.
        """
        return SERIALIZATION_PROTOCOLS[protocol].serialize(self)

    @classmethod
    def load(cls, serialized):
        """Marshall a previously serialized association

        Parameters
        ----------
        serialized: object
            The serialized form of the association.

        Returns
        -------
        The Association object
        """
        asn = None
        for protocol in SERIALIZATION_PROTOCOLS:
            try:
                asn = SERIALIZATION_PROTOCOLS[protocol].unserialize(serialized)
                break
            except AssociationError:
                continue
        else:
            raise AssociationError('Cannot translate "{}" to an association'.format(serialized))

        return asn

    def to_json(self):
        """Create JSON representation.

        Returns
        -------
        (name, str):
            Tuple where the first item is the suggested
            base name for the JSON file.
            Second item is the string containing the JSON serialization.
        """

        # Validate
        schema_path = libpath(self.schema_file)
        with open(schema_path, 'r') as schema_file:
            adb_schema = json.load(schema_file)
        jsonschema.validate(self.data, adb_schema)

        return (
            self.asn_name,
            json.dumps(self.data, indent=4, separators=(',', ': '))
        )

    @classmethod
    def from_json(cls, serialized):
        """Unserialize an assocation from JSON

        Parameters
        ----------
        serialized: str or file object
            The JSON to read

        Returns
        -------
        association: dict
            The association
        """
        try:
            asn = json.loads(serialized)
        except TypeError:
            try:
                asn = json.load(serialized)
            except IOError:
                raise AssociationError('Containter is not JSON: "{}"'.format(serialized))

        return asn

    def add(self, member, check_constraints=True):
        """Add the member to the association

        Parameters
        ----------
        member: dict
            The member to add.

        check_constraints: bool
            If True, see if the member should belong to this association.
            If False, just add it.
        """
        if check_constraints:
            self.test_and_set_constraints(member)

        self._add(member)

    @nottest
    def test_and_set_constraints(self, member):
        """Test whether the given dictionaries match parameters for
        for this association

        Parameters
        ----------
        member: dict
            The parameters to check/set for this association.
            This can be a list of dictionaries.

        Raises
        ------
        AssociationError
            If a match fails.

        Notes
        -----
        If a constraint is present, but does not have a value,
        that constraint is set, and, by definition, matches.
        """
        for constraint, conditions in self.constraints.items():
            try:
                input, value = getattr_from_list(member, conditions['inputs'])
            except KeyError:
                if conditions.get('required', self.DEFAULT_REQUIRE_CONSTRAINT):
                    raise AssociationError(
                        'Constraint {} not present in member.'.format(constraint)
                    )
                else:
                    conditions['value'] = 'Constraint not present and ignored'
                    continue
            if conditions['value'] is not None:
                if not meets_conditions(
                        value, conditions['value']
                ):
                    raise AssociationError(
                        'Constraint {} does not match association.'.format(constraint)
                    )
            if conditions['value'] is None or \
               conditions.get('force_unique', self.DEFAULT_FORCE_UNIQUE):
                logger.debug('Input="{}" Value="{}"'.format(input, value))
                conditions['value'] = re.escape(value)
                conditions['input'] = [input]
                conditions['force_unique'] = False

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
        for c, p in self.constraints.items():
            yield '    {:s}: {}'.format(c, p['value'])

    @classmethod
    def reset_sequence(cls):
        cls._sequence = count(1)

    def _init_hook(self, member):
        """Post-check and pre-member-adding initialization."""
        pass

    def _add(self, member):
        """Add a member, association-specific"""
        raise NotImplementedError('Association._add must be implemented by a specific assocation rule.')


# Utilities
def import_from_file(filename):
    path = expandvars(expanduser(filename))
    module_name = basename(path).split('.')[0]
    folder = dirname(path)
    sys.path.insert(0, folder)
    module = __import__(module_name)
    sys.path.pop(0)
    return module


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

    if isinstance(conditions, six.string_types):
        conditions = [conditions]
    for condition in conditions:
        match = re.match(condition, value, flags=re.IGNORECASE)
        if match:
            return True
    return False


def libpath(filepath):
    '''Return the full path to the module library.'''

    return join(dirname(abspath(getfile(Association))),
                'lib',
                filepath)


def make_timestamp():
    timestamp = datetime.utcnow().strftime(
        _TIMESTAMP_TEMPLATE
    )
    return timestamp


def getattr_from_list(adict, attributes):
    """Retrieve value from dict using a list of attributes

    Parameters
    ----------
    adict: dict
        dict to retrieve from

    attributes: list
        List of attributes

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
    for attribute in attributes:
        try:
            result = adict[attribute]
            if result is masked:
                raise KeyError
            return attribute, result
        except KeyError:
            continue
    else:
        raise KeyError


def get_classes(module):
    """Recursively get all classes in the module"""
    logger.debug('Called.')
    for class_name, class_object in getmembers(
            module,
            lambda o: isclass(o) or ismodule(o)
    ):
        logger.debug('name="{}" object="{}"'.format(class_name, class_object))
        if ismodule(class_object) and class_name.startswith('asn_'):
            for sub_name, sub_class in get_classes(class_object):
                yield sub_name, sub_class
        elif isclass(class_object):
            yield class_name, class_object


# Available serialization protocols
ProtocolFuncs = namedtuple('ProtocolFuncs', ['serialize', 'unserialize'])
SERIALIZATION_PROTOCOLS = {
    'json': ProtocolFuncs(
        serialize=Association.to_json,
        unserialize=Association.from_json)
}

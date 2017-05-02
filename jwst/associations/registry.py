"""Association Registry"""
from importlib import import_module
from inspect import (
    getmembers,
    isclass,
    ismodule
)
import logging
from os.path import (
    basename,
    dirname,
    expanduser,
    expandvars,
)
import sys

from . import libpath
from .exceptions import (
    AssociationError,
    AssociationNotValidError,
    AssociationProcessMembers,
)
from .lib.callback_registry import CallbackRegistry

__all__ = ['AssociationRegistry']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Library files
_ASN_RULE = 'association_rules.py'

# User-level association definitions start with...
USER_ASN = 'Asn_'


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

    name: str
        An identifying string, used to prefix rule names.

    include_bases: bool
        If True, include base classes not considered
        rules.
    """

    # Callback registry
    callback = CallbackRegistry()

    def __init__(self,
                 definition_files=None,
                 include_default=True,
                 global_constraints=None,
                 name=None,
                 include_bases=False):
        super(AssociationRegistry, self).__init__()

        # Generate a UUID for this instance. Used to modify rule
        # names.
        self.name = name

        # Precache the set of rules
        self._rule_set = set()

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

        self.schemas = []
        Utility = type('Utility', (object,), {})
        for fname in definition_files:
            logger.debug('Import rules files "{}"'.format(fname))
            module = import_from_file(fname)
            logger.debug('Module="{}"'.format(module))
            self.schemas += [
                schema
                for schema in find_member(module, 'ASN_SCHEMA')
            ]
            for class_name, class_object in get_classes(module):
                logger.debug('class_name="{}"'.format(class_name))
                if include_bases or class_name.startswith(USER_ASN):
                    try:
                        rule_name = '_'.join([self.name, class_name])
                    except TypeError:
                        rule_name = class_name
                    rule = type(rule_name, (class_object,), {})
                    rule.GLOBAL_CONSTRAINTS = global_constraints
                    rule.registry = self
                    self.__setitem__(rule_name, rule)
                    self._rule_set.add(rule)
                if class_name == 'Utility':
                    Utility = type('Utility', (class_object, Utility), {})
        self.Utility = Utility

    @property
    def rule_set(self):
        return self._rule_set

    def match(self, member, version_id=None, allow=None, ignore=None):
        """See if member belongs to any of the associations defined.

        Parameters
        ----------
        member: dict
            A member, like from a Pool, to find assocations for.

        version_id: str
            If specified, a string appened to association names.
            If None, nothing is used.

        allow: [type(Association), ...]
            List of rules to allow to be matched. If None, all
            available rules will be used.

        ignore: list
            A list of associations to ignore when looking for a match.
            Intended to ensure that already created associations
            are not re-created.

        Returns
        -------
        (associations, reprocess_list): 2-tuple
            associations: [association,...]
                List of associations member belongs to. Empty if none match
            reprocess_list: [AssociationReprocess, ...]
                List of reprocess events.
        """
        logger.debug('Starting...')
        if allow is None:
            allow = self.rule_set
        if ignore is None:
            ignore = []
        associations = []
        process_list = []
        for name, rule in self.items():
            if rule not in ignore and rule in allow:
                logger.debug('Checking membership for rule "{}"'.format(rule))
                try:
                    associations.append(rule(member, version_id))
                except AssociationError as error:
                    logger.debug('Rule "{}" not matched'.format(name))
                    logger.debug('Reason="{}"'.format(error))
                except AssociationProcessMembers as process_event:
                    logger.debug('Process event "{}"'.format(process_event))
                    process_list.append(process_event)
                else:
                    logger.debug('Member belongs to rule "{}"'.format(rule))
        return associations, process_list

    def validate(self, association):
        """Validate a given association against schema

        Parameters
        ----------
        association: association-like
            The data to validate

        Returns
        -------
        rules: list
            List of rules that validated

        Raises
        ------
        AssociationNotValidError
            Association did not validate
        """

        # Change rule validation from an exception
        # to a boolean
        def is_valid(rule, association):
            try:
                rule.validate(association)
            except AssociationNotValidError:
                return False
            else:
                return True

        results = [
            rule
            for rule_name, rule in self.items()
            if is_valid(rule, association)
        ]

        if len(results) == 0:
            raise AssociationNotValidError(
                'Structure did not validate: "{}"'.format(association)
            )
        return results

    def load(
            self,
            serialized,
            format=None,
            validate=True,
            first=True,
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

        first: bool
            A serialization potentially matches many rules.
            Only return the first succesful load.

        kwargs: dict
            Other arguments to pass to the `load` method

        Returns
        -------
        The Association object, or the list of association objects.

        Raises
        ------
        AssociationError
            Cannot create or validate the association.
        """
        results = []
        for rule_name, rule in self.items():
            try:
                results.append(
                    rule.load(
                        serialized,
                        format=format,
                        validate=validate,
                        **kwargs
                    )
                )
            except (AssociationError, AttributeError) as err:
                lasterr = err
                continue
            if first:
                break
        if len(results) == 0:
            logger.debug('Data did not validate against any rule.')
            raise lasterr
        if first:
            return results[0]
        else:
            return results

    def finalize(self, associations):
        """Finalize newly generated associations

        Parameters
        ----------
        assocations: [association[, ...]]
            The list of associations
        """
        finalized = self.callback.filter('finalize', associations)
        return finalized

# Utilities
def import_from_file(filename):
    """Import a file as a module

    Parameters
    ---------
    filename: str
        The file to import

    Returns
    -------
    module: python module
        The imported module
    """
    path = expandvars(expanduser(filename))
    module_name = basename(path).split('.')[0]
    folder = dirname(path)
    sys.path.insert(0, folder)
    try:
        module = import_module(module_name)
    finally:
        sys.path.pop(0)
    return module


def find_member(module, member):
    """Find all instances of member in module or sub-modules

    Parameters
    ----------
    module: module
        The module to recursively search through.

    member: str
        The member to find.

    Returns
    -------
    values: iterator
        Iterator that returns all values of the member
    """
    for name, value in getmembers(module):
        if ismodule(value) and name.startswith('asn_'):
            for inner_value in find_member(value, member):
                yield inner_value
        elif name == member:
            yield value


def get_classes(module):
    """Recursively get all classes in the module

    Parameters
    ----------
    module: python module
        The module to examine

    Returns
    -------
    class members: generator
        A generator that will yield all class members in the module.
    """
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


# ##########
# Unit Tests
# ##########
def test_import_from_file():
    from copy import deepcopy
    from pytest import raises as pytest_raises
    from tempfile import NamedTemporaryFile

    current_path = deepcopy(sys.path)
    with NamedTemporaryFile() as junk_fh:
        junk_path = junk_fh.name
        with pytest_raises(ImportError):
            module = import_from_file(junk_path)
        assert current_path == sys.path

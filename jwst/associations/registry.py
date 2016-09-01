"""Association Registry"""
from importlib import import_module
from inspect import (
    getmembers,
    isclass,
    ismodule
)
import json
import jsonschema
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
    AssociationNotValidError
)

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
                if class_name.startswith(USER_ASN):
                    class_object.GLOBAL_CONSTRAINTS = global_constraints
                    self.__setitem__(class_name, class_object)
                if class_name == 'Utility':
                    Utility = type('Utility', (class_object, Utility), {})
        self.Utility = Utility

    def match(self, member, timestamp=None, ignore=None, reprocess_cb=None):
        """See if member belongs to any of the associations defined.

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

        reprocess_cb: func(member)
            If a member should be reprocessed,
            this function is called.

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
                    associations.append(rule(member, timestamp, reprocess_cb))
                except AssociationError as error:
                    logger.debug('Rule "{}" not matched'.format(name))
                    logger.debug('Error="{}"'.format(error))
                    continue
        if len(associations) == 0:
            raise AssociationError('Member does not match any rules.')
        return associations

    def validate(self, association):
        """Validate a given association against schema

        Parameters
        ----------
        association: dict
            The data to validate

        Returns
        -------
        schemas: list
            List of schemas which validated

        Raises
        ------
        AssociationNotValidError
            Association did not validate
        """
        results = []
        for schema_file in self.schemas:
            with open(schema_file, 'r') as handle:
                schema = json.load(handle)
            try:
                jsonschema.validate(association, schema)
            except jsonschema.ValidationError:
                continue
            else:
                results.append(schema)

        if len(results) == 0:
            raise AssociationNotValidError(
                'Structure did not valid: "{}"'.format(association)
            )
        return results


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

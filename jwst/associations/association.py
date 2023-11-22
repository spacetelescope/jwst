from collections.abc import MutableMapping
from copy import deepcopy
from datetime import datetime, timezone
import json
import jsonschema
import logging
import re
import os
import warnings

from . import __version__
from .exceptions import (
    AssociationNotValidError
)
from .lib.constraint import (
    Constraint,
    meets_conditions
)
from stpipe.format_template import FormatTemplate
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
    version_id : str or None
        Version ID to use in the name of this association.
        If None, nothing is added.

    Raises
    ------
    AssociationError
        If an item doesn't match.

    Attributes
    ----------
    instance : dict-like
        The instance is the association data structure.
        See `data` below

    meta : dict
        Information about the association.

    data : dict
        The association. The format of this data structure
        is determined by the individual associations and, if
        defined, validated against their specified schema.

    schema_file : str
        The name of the output schema that an association
        must adhere to.
    """

    registry = None
    """Registry this rule has been placed in."""

    DEFAULT_FORCE_UNIQUE = False
    """Default whether to force constraints to use unique values."""

    DEFAULT_REQUIRE_CONSTRAINT = True
    """Default require that the constraint exists or otherwise
    can be explicitly checked.
    """

    DEFAULT_EVALUATE = False
    """Default do not evaluate input values"""

    GLOBAL_CONSTRAINT = None
    """Global constraints"""

    INVALID_VALUES = None
    """Attribute values that indicate the
    attribute is not specified.
    """

    ioregistry = IORegistry()
    """The association IO registry"""

    def __init__(
            self,
            version_id=None,
    ):

        self.data = dict()
        self.run_init_hook = True
        self.meta = {}

        self.version_id = version_id

        self.data.update({
            'asn_type': 'None',
            'asn_rule': self.asn_rule,
            'version_id': self.version_id,
            'code_version': __version__,
        })

        # Setup constraints
        # These may be predefined by a rule.
        try:
            constraints = self.constraints
        except AttributeError:
            constraints = Constraint()
        if self.GLOBAL_CONSTRAINT is not None:
            constraints.append(self.GLOBAL_CONSTRAINT.copy())
        self.constraints = constraints

    @classmethod
    def create(cls, item, version_id=None):
        """Create association if item belongs

        Parameters
        ----------
        item : dict
            The item to initialize the association with.

        version_id : str or None
            Version ID to use in the name of this association.
            If None, nothing is added.

        Returns
        -------
        (association, reprocess_list)
            2-tuple consisting of:
                - association or None: The association or, if the item does not
                  match this rule, None
                - [ProcessList[, ...]]: List of items to process again.
        """
        asn = cls(version_id=version_id)
        matches, reprocess = asn.add(item)
        if not matches:
            return None, reprocess
        return asn, reprocess

    @property
    def asn_name(self):
        """Suggest filename for the association"""
        return 'unnamed_association'

    @classmethod
    def _asn_rule(cls):
        return cls.__name__

    @property
    def asn_rule(self):
        """Name of the rule"""
        return self._asn_rule()

    @classmethod
    def validate(cls, asn):
        """Validate an association against this rule

        Parameters
        ----------
        asn : Association or association-like
            The association structure to examine

        Returns
        -------
        valid : bool
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
            logger.debug('Validation failed:\n%s', err)
            raise AssociationNotValidError('Validation failed') from err

        # Validate no path data for expnames
        for product in asn_data["products"]:
            members = product['members']
            for member in members:
                fpath, fname = os.path.split(member["expname"])
                if len(fpath) > 0:
                    err_str = "'expname' contains path, but should only be a filename."
                    err_str += "  All input files should be in a single directory"
                    err_str += ", so no path is needed."
                    logger.debug(err_str)
                    warnings.warn(err_str, UserWarning)
        return True

    def dump(self, format='json', **kwargs):
        """Serialize the association

        Parameters
        ----------
        format : str
            The format to use to dump the association into.

        kwargs : dict
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
        serialized : object
            The serialized form of the association.

        format : str or None
            The format to force. If None, try all available.

        validate : bool
            Validate against the class' defined schema, if any.

        kwargs : dict
            Other arguments to pass to the `load` method

        Returns
        -------
        association : Association
            The association.

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
                format_func
                for format_name, format_func in cls.ioregistry.items()
            ]
        else:
            formats = [cls.ioregistry[format]]

        for format_func in formats:
            try:
                asn = format_func.load(
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

    def add(self, item, check_constraints=True):
        """Add the item to the association

        Parameters
        ----------
        item : dict
            The item to add.

        check_constraints : bool
            If True, see if the item should belong to this association.
            If False, just add it.

        Returns
        -------
        (match, reprocess_list)
            2-tuple consisting of:
                - bool : True if match
                - [ProcessList[, ...]]: List of items to process again.
        """
        if self.is_item_member(item):
            return True, []

        match = not check_constraints
        if check_constraints:
            match, reprocess = self.check_and_set_constraints(item)

        if match:
            if self.run_init_hook:
                self._init_hook(item)
                self.run_init_hook = False
            self._add(item)

        # If a constraint `force_match` exists, set the `match`
        # result to the value of the constraint.
        try:
            force_match = self.constraints['force_match'].value
        except (KeyError, TypeError):
            pass
        else:
            if force_match is not None:
                match = force_match

        return match, reprocess

    def check_and_set_constraints(self, item):
        """Check whether the given dictionaries match parameters for
        for this association

        Parameters
        ----------
        item : dict
            The parameters to check/set for this association.
            This can be a list of dictionaries.

        Returns
        -------
        (match, reprocess)
            2-tuple consisting of:
                - bool : Did constraint match?
                - [ProcessItem[, ...]]: List of items to process again.

        """
        cached_constraints = deepcopy(self.constraints)
        match, reprocess = cached_constraints.check_and_set(item)
        if match:
            self.constraints = cached_constraints

        # Set the association type for all reprocessed items.
        for process_list in reprocess:
            process_list.trigger_rules.update([type(self)])
            if process_list.rules is None:
                process_list.rules = [type(self)]

        return match, reprocess

    def match_constraint(self, item, constraint, conditions):
        """Generic constraint checking

        Parameters
        ----------
        item : dict
            The item to retrieve the values from

        constraint : str
            The name of the constraint

        conditions : dict
            The conditions structure

        Returns
        -------
        (matches, reprocess_list)
            2-tuple consisting of:
                - bool : True if the all constraints are satisfied
                - [ProcessList[, ...]]: List of items to process again.
        """
        reprocess = []
        evaled_str = conditions['inputs'](item)
        if conditions['value'] is not None:
            if not meets_conditions(
                    evaled_str, conditions['value']
            ):
                return False, reprocess

        # At this point, the constraint has passed.
        # Fix the conditions.
        escaped_value = re.escape(evaled_str)
        conditions['found_values'].add(escaped_value)
        if conditions['value'] is None or \
           conditions.get('force_unique', self.DEFAULT_FORCE_UNIQUE):
            conditions['value'] = escaped_value
            conditions['force_unique'] = False

        # That's all folks
        return True, reprocess

    def finalize(self):
        """Finalize association

        Finalize or close-off this association. Perform validations,
        modifications, etc. to ensure that the association is
        complete.

        Returns
        -------
        associations : [association[, ...]] or None
            List of fully-qualified associations that this association
            represents.
            `None` if a complete association cannot be produced.

        """
        if self.is_valid:
            return [self]
        else:
            return None

    def is_item_member(self, item):
        """Check if item is already a member of this association

        Parameters
        ----------
        item : dict
            The item to add.

        Returns
        -------
        is_item_member : bool
            True if item is a member.
        """
        raise NotImplementedError(
            'Association.is_item_member must be implemented by a specific association rule.'
        )

    def _init_hook(self, item):
        """Post-check and pre-item-adding initialization."""
        pass

    def _add(self, item):
        """Add a item, association-specific"""
        raise NotImplementedError(
            'Association._add must be implemented by a specific association rule.'
        )

    def _add_items(self, items, **kwargs):
        """ Force adding items to the association

        Parameters
        ----------
        items : [object[, ...]]
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


# #########
# Utilities
# #########
def finalize(asns):
    """Finalize associations by calling their `finalize_hook` method

    Notes
    -----
    This is a functioning example of a finalize callback, and can be used
    as the generic callback. Suggested usage is as follows:

    .. code-block:: python

       from jwst.associations.association import finalize as generic_finalize
       RegistryMarker.callback('finalize')(generic_finalize)
    """
    finalized_asns = list(filter(
        lambda asn: asn is not None,
        map(lambda asn: asn.finalize(), asns)
    ))
    return finalized_asns


def make_timestamp():
    timestamp = datetime.now(timezone.utc).strftime(
        _TIMESTAMP_TEMPLATE
    )
    return timestamp


# Define default product name filling
format_product = FormatTemplate()

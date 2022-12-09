"""Association Registry"""
import importlib.util
from inspect import (
    getmembers,
    isclass,
    isfunction,
    ismethod,
    ismodule
)
import logging
from os.path import (
    basename,
    expanduser,
    expandvars,
)
from enum import EnumMeta

from . import libpath
from .exceptions import (
    AssociationError,
    AssociationNotValidError
)
from .lib.callback_registry import CallbackRegistry

__all__ = [
    'AssociationRegistry',
    'RegistryMarker'
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Library files
_ASN_RULE = 'association_rules.py'


class AssociationRegistry(dict):
    """The available associations

    Parameters
    ----------
    definition_files : [str,]
        The files to find the association definitions in.

    include_default : bool
        True to include the default definitions.

    global_constraints : Constraint
        Constraints to be added to each rule.

    name : str
        An identifying string, used to prefix rule names.

    include_bases : bool
        If True, include base classes not considered
        rules.

    Notes
    -----
    The general workflow is as follows:

        * Create the registry
            >>> from jwst.associations.registry import AssociationRegistry
            >>> registry = AssociationRegistry()

        * Create associations from an item
            >>> associations, reprocess = registry.match(item) # doctest: +SKIP

        * Finalize the associations
            >>> final_asns = registry.callback.reduce('finalize', associations) # doctest: +SKIP

    In practice, this is one step in a larger loop over all items to
    be associated. This does not account for adding items to already
    existing associations. See :py:func:`~jwst.associations.generate` for more information.
    """

    def __init__(self,
                 definition_files=None,
                 include_default=True,
                 global_constraints=None,
                 name=None,
                 include_bases=False):
        super().__init__()

        # Generate a UUID for this instance. Used to modify rule
        # names.
        self.name = name

        # Callback registry
        self.callback = CallbackRegistry()

        # Precache the set of rules
        self._rule_set = set()

        if definition_files is None:
            definition_files = []
        if include_default:
            definition_files.insert(0, libpath(_ASN_RULE))
        if len(definition_files) <= 0:
            raise AssociationError('No rule definition files specified.')

        self.schemas = []
        self.Utility = type('Utility', (object,), {})
        for fname in definition_files:
            module = import_from_file(fname)
            self.populate(
                module,
                global_constraints=global_constraints,
                include_bases=include_bases
            )

    @property
    def rule_set(self):
        """Rules within the Registry"""
        return self._rule_set

    def match(self, item, version_id=None, allow=None, ignore=None):
        """See if item belongs to any of the associations defined.

        Parameters
        ----------
        item : dict
            An item, like from a Pool, to find associations for.

        version_id : str
            If specified, a string appended to association names.
            If None, nothing is used.

        allow : [type(Association), ...]
            List of rules to allow to be matched. If None, all
            available rules will be used.

        ignore : list
            A list of associations to ignore when looking for a match.
            Intended to ensure that already created associations
            are not re-created.

        Returns
        -------
        (associations, reprocess_list): 2-tuple
            associations : [association,...]
                List of associations item belongs to. Empty if none match
            reprocess_list : [AssociationReprocess, ...]
                List of reprocess events.
        """
        if allow is None:
            allow = self.rule_set
        if ignore is None:
            ignore = []
        associations = []
        process_list = []
        for name, rule in self.items():
            if rule not in ignore and rule in allow:
                asn, reprocess = rule.create(item, version_id)
                process_list.extend(reprocess)
                if asn is not None:
                    associations.append(asn)
        return associations, process_list

    def validate(self, association):
        """Validate a given association

        Parameters
        ----------
        association : association-like
            The data to validate

        Returns
        -------
        rules : list
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
        """Load a previously serialized association

        Parameters
        ----------
        serialized : object
            The serialized form of the association.

        format : str or None
            The format to force. If None, try all available.

        validate : bool
            Validate against the class' defined schema, if any.

        first : bool
            A serialization potentially matches many rules.
            Only return the first succesful load.

        kwargs : dict
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
            raise lasterr
        if first:
            return results[0]
        else:
            return results

    def populate(self,
                 module,
                 global_constraints=None,
                 include_bases=None
                 ):
        """Parse out all rules and callbacks in a module and add them to the registry

        Parameters
        ----------
        module : module
            The module, and all submodules, to be parsed.
        """
        for name, obj in get_marked(module, include_bases=include_bases):

            # Add rules.
            if include_bases or obj._asnreg_role == 'rule':
                try:
                    self.add_rule(name, obj, global_constraints=global_constraints)
                except TypeError:
                    logger.debug(f'Could not add object {obj} as a rule due to TypeError')
                continue

            # Add callbacks
            if obj._asnreg_role == 'callback':
                for event in obj._asnreg_events:
                    self.callback.add(event, obj)
                    continue

            # Add schema
            if obj._asnreg_role == 'schema':
                self.schemas.append(obj._asnreg_schema)
                continue

            # Add utility classes
            if obj._asnreg_role == 'utility':
                self.Utility = type(
                    'Utility',
                    (obj, self.Utility),
                    {}
                )

    def add_rule(self, name, obj, global_constraints=None):
        """Add object as rule to registry

        Parameters
        ----------
        name : str
            Name of the object

        obj : object
            The object to be considered a rule

        global_constraints : dict
            The global constraints to attach to the rule.
        """
        try:
            rule_name = '_'.join([self.name, name])
        except TypeError:
            rule_name = name
        if not valid_class(obj):
            raise TypeError(f'Object cannot be used as rule: {obj}')
        rule = type(rule_name, (obj,), {})
        rule.GLOBAL_CONSTRAINT = global_constraints
        rule.registry = self
        self.__setitem__(rule_name, rule)
        self._rule_set.add(rule)


class RegistryMarker:
    """Mark rules, callbacks, and modules for inclusion into a registry"""

    class Schema:
        def __init__(self, obj):
            self._asnreg_role = 'schema'
            self._asnreg_schema = obj
            RegistryMarker.mark(self)

        @property
        def schema(self):
            return self._asnreg_schema

    @staticmethod
    def mark(obj):
        """Mark that an object should be part of the registry

        Parameters
        ----------
        obj : object
            The object to mark

        Returns
        -------
        obj
            Object that has been marked. Returned to enable
            use as a decorator.

        Notes
        -----
        The following attributes are added to the object:

        - _asnreg_mark : True
              Attribute added to object and is set to True

        - _asnreg_role : str or None
              If not already assigned, the role is left
              unspecified using None.

        """
        obj._asnreg_marked = True
        obj._asnreg_role = getattr(obj, '_asnreg_role', None)
        return obj

    @staticmethod
    def rule(obj):
        """Mark object as rule

        Parameters
        ----------
        obj : object
            The object that should be treated as a rule

        Returns
        -------
        obj : object
            Return object to enable use as a decorator.

        Notes
        -----
        The following attributes are added to the object:

        - _asnreg_role : 'rule'
              Attributed added to object and set to `rule`

        - _asnreg_mark : True
              Attributed added to object and set to True
        """
        obj._asnreg_role = 'rule'
        RegistryMarker.mark(obj)
        return obj

    @staticmethod
    def callback(event):
        """Mark object as a callback for an event

        Parameters
        ----------
        event : str
            Event this is a callback for.

        obj : func
            Function, or any callable, to be called
            when the corresponding event is triggered.

        Returns
        -------
        func
            Function to use as a decorator for the object to be marked.

        Notes
        -----
        The following attributes are added to the object:

        - _asnreg_role : 'callback'
              The role the object as been assigned.
        - _asnreg_events : [event[, ...]]
              The events this callable object is a callback for.
        - _asnreg_mark : True
              Indicated that the object has been marked.
        """
        def decorator(func):
            try:
                events = func._asnreg_events
            except AttributeError:
                events = list()
            events.append(event)
            RegistryMarker.mark(func)
            func._asnreg_role = 'callback'
            func._asnreg_events = events
            return func
        return decorator

    @staticmethod
    def schema(filename):
        """Mark a file as a schema source"""
        schema = RegistryMarker.Schema(filename)
        return schema

    @staticmethod
    def utility(class_obj):
        """Mark the class as a Utility class"""
        class_obj._asnreg_role = 'utility'
        RegistryMarker.mark(class_obj)
        return class_obj

    @staticmethod
    def is_marked(obj):
        """Has an objected been marked?"""
        return hasattr(obj, '_asnreg_marked')


# Utilities
def import_from_file(filename):
    """Import a file as a module

    Parameters
    ---------
    filename : str
        The file to import

    Returns
    -------
    module : python module
        The imported module
    """
    path = expandvars(expanduser(filename))
    module_name = basename(path).split('.')[0]
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def get_marked(module, predicate=None, include_bases=False):
    """Recursively get all executable objects

    Parameters
    ----------
    module : python module
        The module to examine

    predicate : bool func(object)
        Determinant of what gets returned.
        If None, all object types are examined

    include_bases : bool
        If True, include base classes not considered
        rules.

    Returns
    -------
    class object : generator
        A generator that will yield all class members in the module.
    """
    def is_method(obj):
        return isfunction(obj) or ismethod(obj)

    for name, obj in getmembers(module, predicate):
        if isclass(obj):
            for sub_name, sub_obj in get_marked(obj, predicate=is_method):
                yield sub_name, sub_obj
            if RegistryMarker.is_marked(obj) or include_bases:
                yield name, obj
        elif RegistryMarker.is_marked(obj):
            if ismodule(obj):
                for sub_name, sub_obj in get_marked(
                        obj, predicate=predicate, include_bases=include_bases
                ):
                    yield sub_name, sub_obj
            else:
                yield name, obj


def valid_class(obj):
    """Verify if a given object could be used as a rule class

    Parameters
    ------------
    obj : obj
        Object to check

    Returns
    --------
    is_valid : bool
        True if the object could be considered a rule class
    """
    is_valid = type(obj) is not EnumMeta and isclass(obj)
    return is_valid

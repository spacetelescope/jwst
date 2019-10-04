.. currentmodule:: jwst.associations.association
                   
.. _reference-association:

Association Rules
=================

Association definitions, or ``rules``, are Python classes, all based on
:py:class:`~jwst.associations.Association`. The base class provides only a
framework, much like an abstract base class; all functionality must be
implemented in sub-classes.

Any subclass that is intended to produce an association is referred to
as a ``rule``. Any rule subclass must have a name that begins with the
string ``Asn_``. This is to ensure that any other classes involved in
defining the definition of the rule classes do not get used as rules
themselves, such as the :py:class:`~jwst.associations.Association` itself.

Association Dynamic Definition
------------------------------

Associations are created by matching members to rules. However, an
important concept to remember is that an association is defined by
both the rule matched, and by the initial member that matched it.
The following example will illustrate this concept.

For JWST :ref:`Level 3<asn-level3-techspecs>`, many associations created must
have members that all share the same filter. To avoid writing rules
for each filter, the rules have a condition that states that it
doesn't matter what filter is specified, as long as the association
contains all the same filter.

To accomplish this, the association defines a constraint where filter
must have a valid value, but can be any valid value. When the
association is first attempted to be instantiated with a member, and
that member has a valid filter, the association is created. However,
the constraint on filter value in the newly created association is
modified to match exactly the filter value that the first member had.
Now, when other members are attempted to be added to the association,
the filter of the new members must match exactly with what the
association is expecting.

This dynamic definition allows rules to be written where each value of
a specific attribute of a member does not have to be explicitly
stated. This provides for very robust, yet concise, set of rule definitions.

User-level API
--------------

Core Keys
^^^^^^^^^

To be repetitive, the basic association is simply a dict (default) or
list. The structure of the dict is completely determined by the rules.
However, the base class defines the following keys:

    ``asn_type``
        The type of the association.

    ``asn_rule``
        The name of the rule.

    ``version_id``
        A version number for any associations created by this rule.

    ``code_version``
        The version of the generator library in use.

These keys are accessed in the same way any dict key is accessed::

  asn = Asn_MyAssociation()
  print(asn['asn_rule'])
  
  #--> MyAssociation

.. _ref-asn-core-methods:

Core Methods
^^^^^^^^^^^^

These are the methods of an association rule deal with creation or returning the
created association. A rule may define other methods, but the
following are required to be implemented.

    :meth:`create() <Association.create>`
          Create an association.

    :meth:`add() <Association.add>`
          Add a member to the current association.

    :meth:`dump() <Association.dump>`
          Return the string serialization of the association.

    :meth:`load() <Association.load>`
          Return the association from its serialization.


Creation
^^^^^^^^

To create an association based on a member, the ``create`` method of the
rule is called::

  (association, reprocess_list) = Asn_SomeRule.create(member)

``create`` returns a 2-tuple: The first element is the association and the
second element is a list of ``reprocess`` instances.

If the member matches the conditions for the rule, an association is
returned. If the member does not belong, ``None`` is returned for the
association.

Whether an association is created or not, it is possible a list of
``reprocess`` instances may be returned. This list represents the
expansion of the pool in :ref:`member-with-lists`

Addition
^^^^^^^^

To add members to an existing association, one uses the :py:meth:`Association.add
<jwst.associations.Association>` method::

  (matches, reprocess_list) = association.add(new_member)

If the association accepts the member, the ``matches`` element of the
2-tuple will be ``True``.

Typically, one does not deal with a single rule, but a collection of
rules. For association creation, one typically uses an
:py:class:`~jwst.associations.AssociationRegistry` to collect all the rules a pool will be
compared against. Association registries provide extra functionality to
deal with a large and varied set of association rules.

Saving and Loading
^^^^^^^^^^^^^^^^^^

Once created, an association can be serialized using its
:meth:`Association.dump <jwst.associations.association.Association.dump>` method.
Serialization creates a string representation of the association which
can then be saved as one wishes. Some code that does a basic save
looks like::

  file_name, serialized = association.dump()
  with open(file_name, 'w') as file_handle:
      file_handle.write(serialized)

Note that ``dump`` returns a 2-tuple. The first element is the suggested
file name to use to save the association. The second element is the
serialization.

To retrieve an association, one uses the :meth:`Association.load
<jwst.associations.association.Association.load>` method::

  with open(file_name, 'r') as file_handle:
      association = Association.load(file_handle)

:meth:`Association.load
<jwst.associations.association.Association.load>` will only validate
the incoming data against whatever schema or other validation checks
the particular subclass calls for. The generally preferred method for
loading an association is through the
:func:`jwst.associations.load_asn` function.

Defining New Associations
-------------------------

All association rules are based on the
:class:`~jwst.associations.association.Association` base class. This
class will not create associations on its own; subclasses must be
defined. What an association is and how it is later used is completely
left to the subclasses. The base class itself only defines the
framework required to create associations. The rest of this section
will discuss the minimum functionality that a subclass needs to
implement in order to create an association.

.. _class-naming:

Class Naming
^^^^^^^^^^^^

The :py:class:`~jwst.associations.AssociationRegistry` is used to store
the association rules. Since rules are defined by Python classes, a
way of indicating what the final rule classes are is needed. By
definition, rule classes are classes that begin with the string ``Asn_``.
Only these classes are used to produce associations.

Core Attributes
^^^^^^^^^^^^^^^

Since rule classes will potentially have a large number of attributes
and methods, the base :class:`Association
<jwst.associations.association.Association>` class defines two
attributes: ``data``, which contains the actual association, and ``meta``,
the structure that holds auxiliary information needed for association
creation. Subclasses may redefine these attributes as they see fit.
However, it is suggested that they be used as conceptually defined here.

``data`` Attribute
""""""""""""""""""

``data`` contains the association itself. Currently, the base class
predefines ``data`` as a dict. The base class itself is a subclass of
`~py3:collections.abc.MutableMapping`. Any instance behaves as a dict. The contents of that
dict is the contents of the ``data`` attribute. For example::

  asn = Asn_MyAssociation()
  asn.data['value'] = 'a value'
  
  assert asn['value'] == 'a value'
  # True

  asn['value'] = 'another value'
  assert asn.data['value'] == 'another value'
  # True

Instantiation
^^^^^^^^^^^^^

Instantiating a rule, in and of itself, does nothing more than setup
the constraints that define the rule, and basic structure
initialization.

Implementing ``create()``
^^^^^^^^^^^^^^^^^^^^^^^^^

The base class function performs the following steps:

    - Instantiates an instance of the rule
    - Calls ``add()`` to attempt to add the member to the instance

If ``add()`` returns ``matches==False``, then ``create`` returns ``None`` as the
new association.

Any override of this method is expected to first call ``super``. On
success, any further initialization may be performed.

Implementing ``add()``
^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~jwst.associations.association.Association.add` method adds
members to an association.

If a member does belong to the association, the following events
occur:

Constraint Modification
    Any wildcard constraints are modified so that any further matching
    must match exactly the value provided by the current member.

``self._init_hook()`` is executed
    If a new association is being created, the
    rule's ``_init_hook`` method is executed, if defined. This allows a
    rule to do further initialization before the member is
    officially added to the association.

``self._add()`` is executed
    The rule class must define ``_add()``. This method officially adds
    the member to the association.

Implementing ``dump()`` and ``load()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The base ``Association`` class defines the
:meth:`~jwst.associations.association.Association.dump` and
:meth:`~jwst.associations.association.Association.load` methods to
serialize the data structure pointing to by the ``data`` attribute. If
the new rule uses the ``data`` attribute for storing the association
information, no further overriding of these methods is necessary.

However, if the new rule does not define ``data``, then these methods
will need be overridden.

Rule Registration
=================

In order for a rule to be used by ``generate``, the rule must be loaded
into an ``AssociationRegistry``.  Since a rule is just a class that is
defined as part of a, most likely, larger module, the registry needs
to know what classes are rules. Classes to be used as rules are marked
with the ``RegistryMarker.rule`` decorator as follows::

  # myrules.py
  from jwst.associations import (Association, RegistryMarker)

  @RegistryMarker.rule
  class MyRule(Association):
      ...

Then, when the rule file is used to create an ``AssociationRegistry``,
the class ``MyRule`` will be included as one of the available rules::

.. doctest-skip::
   
  >>> from jwst.associations import AssociationRegistry
  >>> registry = AssociationRegistry('myrules.py', include_default=False)
  >>> print(registry)
      {'MyRule': <class 'abc.MyRule'>}


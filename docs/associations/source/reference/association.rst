.. currentmodule:: jwst.associations.association
                   
.. _reference-association:

*****************
Association Rules
*****************

Association definitions, or `rules`, are Python classes, all based on
the :ref:`association`. The base class provides only a framework, much
like an abstract base class; all functionality must be implemented in
subclasses.

Any subclass that is intended to produce an association is referred to
as a `rule`. Any rule subclass must have a name that begins with the
string `Asn_`. This is to ensure that any other classes involved in
defining the definition of the rule classes do not get used as rules
themselves, such as the :ref:`association` itself.

Association Dynamic Definition
==============================

Associations are created by matching members to rules. However, an
important concept to remember is that an association is defined by
both the rule matched, and by the initial member that matched it.
Below is an example of what this means.

For JWST :ref:`level3-associations`, many associations created must
have members that all share the same filter. To avoid writing rules
for each filter, the rules have a condition that states that it
doesn't matter what filter is specified, as long as the association
contains all the same filter.

To accomplish this, the association defines a constraint where filter
must have a valid value, but can be any valid value. When the
association is first attempted to be instantiated with a member, and
that member has a valid filter, the association is created. However,
the constraint on filter value in the newly created association is
modified to match exacly the filter value that the first member had.
Now, when other members are attempted to be added to the association,
the filter of the new members must match exactly with what the
association is expecting.

This dynamic definition allows rules to be written where each value of
a specific attribute of a member does not have to be explicitly
stated. This provides for very robust, yet concise, rules.

User-level API
--------------

Core Keys
=========

To be repetitive, the basic association is simply a dict (default) or
list. The structure of the dict is compeletly determined by the rules.
However, the base class defines the following keys:

    `asn_type`
        The type of the association.

    `asn_rule`
        The name of the rule.

    `version_id`
        A version number for any associations created by this rule.

    `code_version`
        The version of the generator library in use.

These keys are accessed in the same way any dict key is accessed::

  asn = Asn_MyAssociation()
  print(asn['asn_rule'])
  
  #--> MyAssociation

.. _ref-asn-core-methods:

Core Methods
============

All methods of an association rule deal with creation or returning the
created association.

    :meth:`Association`
          Create an association.

    :meth:`add() <Association.add>`
          Add a member to the current association.

    :meth:`dump() <Association.dump>`
          Return the string serialization of the association.

    :meth:`load() <Association.load>`
          Return the association from its serialization.

There are other methods that deal with the details of association
creation and rule definition.

Creation
========

To create an association based on a member, one simply attempts to
instantiate the rule using the member as an argument::

  association = Asn_SomeRule(member)

If, for reasons determined within the class instantiation code itself,
the member should not belong the association, an `AssociationError` is
raised. Otherwise, the association is successfully created, with the
member as its first member.

To add members to an existing association, one uses the :meth:`Association.add
<jwst.associations.association.Association.add>` method::

  association.add(new_member)

If the association accepts the member, no further action occurs. If
the member does not belong in the association, again, an
`AssociationError` is raised.

Typically, one does not deal with a single rule, but a collection of
rules. For association creation, one typically uses an
:ref:`AssociationRegistry` to collect all the rules a pool will be
compared against. Associaton registries provide extra functionality to
deal with a large and varied set of association rules.

Saving and Loading
==================
Once created, an association can be serialized using its
:meth:`Association.dump <jwst.associations.association.Association.dump>` method.
Serialization creates a string representation of the association which
can then be saved as one wishes. Some code that does a basic save
looks like::

  file_name, serialized = association.dump()
  with open(file_name, 'w') as file_handle:
      file_handle.write(serialized)

Note that `dump` returns a 2-tuple. The first element is the suggested
file name to use to save the association. The second element is the
serialization.

To retrieve an association, one uses the :meth:`Association.load
<jwst.associations.association.Association.load>` method::

  with open(file_name, 'r') as file_handle:
      association = Association.load(file_handle)

Creating versus Using Associations
==================================


Defining New Associations
-------------------------

All association rules are based on the
:class:`Association <jwst.associations.association.Association>` base class. This
class will not create associations on its own; subclasses must be
defined. What an association is and how it is later used is completely
left to the subclasses. The base class itself only defines the
framework required to create associations. The rest of this section
will discuss the minimum functionality that a subclass needs to
implement in order to create an association.

Class Naming
============

The :ref:`generate` machinery uses :ref:`AssociationRegistry` to store
the association rules. Since rules are defined by Python classes, a
way of indicating what the final rule classes are is needed. By
definition, rule classes are classes that begin with the string `Asn_`.
Only these classes are used to produce associations.

Core Attributes
===============

Since rule classes will potentially have a large number of attributes
and methods, the base :class:`Association
<jwst.associations.association.Association>` class defines two
attributes: `data`, which contains the actual association, and `meta`,
the structure that holds auxiliary information needed for association
creation. Subclasses may redfine these attributes as they see fit.
However, it is suggested that they be used as conceptually defined here.

`data` Attribute
~~~~~~~~~~~~~~~~

`data` contains the association itself. Currently, the base class
predefines `data` as a dict. The base class itself is a subclass of
`MutableMapping`. Any instance behaves as a dict. The contents of that
dict is the contents of the `data` attribute. For example::

  asn = Asn_MyAssociation()
  asn.dict['value'] = 'a value'
  
  assert asn['value'] == 'a value'
  # True

  asn['value'] = 'another value'
  assert asn.data['value'] == 'another value'
  # True

Initialization
==============

When overriding the `__init__()` function, the base classes ultimately
call `add()` to confirm that the given member actually belongs to the
current rule. Therefore, any initialization required to ensure that
the member can actually be added should be done before calling
`super().__init__()`. Any post-member-add operations can be done
afterwards, noting that if the `add()` fails, an `AssociationError` is raised.

Implementing `add()`
====================

*TBD*

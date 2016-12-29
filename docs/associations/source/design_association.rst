.. currentmodule:: jwst.associations.association

.. Links
.. _schema: http://json-schema.org/

.. _design-association:

*******************
Design: Association
*******************

At the end of the day, an association is simply a list or dict
containing the members that belong to that association.

Not to be confused with associations, the association definitions, or
`rules`, are Python classes, all based on the :ref:`association`.
Association rules are discussed in :ref:`design-rules`.

Usage
=====

Associations are read from an association file, sometimes referred to
as an associaion table, using the :meth:`Association.load
<jwst.associations.association.Association.load>` method::

  with open('a_level3_association.json', 'r') as file_handle:
      a_level3_association = Association.load(file_handle)

As stated above, an association is either a simple list of members, or
a dict where one of the attributes is the member list. The contents of
an association is completely determined by the association rules that
created the association. However, the default association rules create
a dict with the following keys:

    `asn_type`
        The type of the association.

    `asn_rule`
        The name of the rule.

    `version_id`
        A version number for any associations created by this rule.

    `code_version`
        The version of the generator library in use.

These keys are accessed in the same way any dict key is accessed::

  with open('a_level3_image_association.json', 'r') as file_handle:
      a_level3_image_association = Association.load(file_handle)

  print(a_level3_image_association['asn_rule'])
  
  #--> Asn_Image
  
For JWST's :ref:`level3-associations`, associations are dicts where
the primary attribute is `products`, the list of Level3 products the
association is meant to produce. Each product has a `members`
attribute which is the list of exposures that go into that product. An
example of accessing the member list would be::

  members = a_level3_association['products][0]['members']

Creation and Modification
=========================

Associations are initially created given a :ref:`pool <design-pool>`
and :ref:`rule set <design-rules>` and run through the :ref:`generator
<design-generator>`. The associations are saved as files, often as
`JSON` or `YAML` files. Since these files are simply text files, to
add or remove members from an association, they that can be edited
with any text editor. Similarly, to create new associations, a copy
can be made of an existing association and then edited.

Often, associations have a schema_, or definition, associated with
them. If this is the case, any editing or create of association files
must match the association releated to said schema. For an example,
see the discussion of JWST's :ref:`level3-associations`.

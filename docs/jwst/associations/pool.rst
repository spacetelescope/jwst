.. _design-pool:

================
Association Pool
================
  
Association pools are tables. Pools are instantiated using the
:class:`~jwst.associations.AssociationPool` class. This class is a subclass of
:ref:`astropy Table <astropy:astropy-table>`.
As such, any file that is supported by :ref:`astropy:table_io` can be
used as an association pool.

Each row of a pool defines a ``member``, and the columns define the
attributes of that member. It is these attributes that the generator
uses to determine which members go into which associations.

Regardless of any implied or explicit typing of data by a table file,
internally all data are converted to lowercase strings. It is left up to the
individual association definitions on how they will use these
attributes.

For JWST Stage 2/Stage 3 associations, there is a special case: If an
attribute has a value that is equivalent to a Python list::

    [element, ...]

the list will be expanded by the Stage 2/Stage 3 associations. This
expansion is explained in :ref:`member-with-lists`.

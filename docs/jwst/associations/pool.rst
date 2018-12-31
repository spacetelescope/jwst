.. _design-pool:

==================
 Association Pool
==================

.. _astropy Table:
   http://docs.astropy.org/en/stable/table/index.html
   
Association pools are simply tables. Pools are instantiated using the
:py:class:`~jwst.associations.AssociationPool`. This class is simply a subclass of `astropy
Table`_. As such, any file that is supported by  astropy I/O can be
used as an association pool.

Each row of a pool defines a ``member``, and the columns define the
attributes of that member. It is these attributes that the generator
uses to determine which members go into which associations.

Regardless of any implied or explicit typing of data by a table file,
internally all data are converted to lowercase strings. It is left up to the
individual association definitions on how they will use these
attributes.

For JWST Level2/Level3 associations, there is a special case. If an
attribute has a value that is equivalent to a Python list::

  [element, ...]

the list will be expanded by the Level2/Level3 associations. This
expansion is explained in :ref:`member-with-lists`

.. _metadata:


Metadata
========

Metadata information associated with a data model is accessed through
its `meta` member.  For example, to access the date that an
observation was made::

    print(model.meta.observation.date)

Metadata values are automatically type-checked against the schema when
they are set. Therefore, setting a keyword which expects a number to a
string will raise an exception::

    >>> from jwst.datamodels import ImageModel
    >>> model = ImageModel()
    >>> model.meta.target.ra = "foo"    # doctest: +SKIP
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "site-packages/jwst.datamodels/schema.py", line 672, in __setattr__
        object.__setattr__(self, attr, val)
      File "site-packages/jwst.datamodels/schema.py", line 490, in __set__
        val = self.to_basic_type(val)
      File "site-packages/jwst.datamodels/schema.py", line 422, in to_basic_type
        raise ValueError(e.message)
    ValueError: 'foo' is not of type u'number'

The set of available metadata elements is defined in a YAML Schema
that ships with `jwst.datamodels`.

There is also a utility method for finding elements in the metadata
schema.  `search_schema` will search the schema for the given
substring in metadata names as well as their documentation.  The
search is case-insensitive::

    >>> from jwst.datamodels import ImageModel
    >>> # Create a model of the desired type
    >>> model = ImageModel()
    >>> # Call `search_schema` on it to find possibly related elements.
    >>> model.search_schema('target')
    meta.target
    <BLANKLINE>
    meta.target.catalog_name
    <BLANKLINE>
    meta.target.dec
    <BLANKLINE>
    meta.target.dec_uncertainty
    <BLANKLINE>
    meta.target.proper_motion_dec
    <BLANKLINE>
    meta.target.proper_motion_epoch
    <BLANKLINE>
    meta.target.proper_motion_ra
    <BLANKLINE>
    meta.target.proposer_dec
    <BLANKLINE>
    meta.target.proposer_name
    <BLANKLINE>
    meta.target.proposer_ra
    <BLANKLINE>
    meta.target.ra
    <BLANKLINE>
    meta.target.ra_uncertainty
    <BLANKLINE>
    meta.target.source_type
    <BLANKLINE>
    meta.target.source_type_apt
    <BLANKLINE>
    meta.target.type
    <BLANKLINE>
    meta.visit.internal_target
    <BLANKLINE>


An alternative method to get and set metadata values is to use a
dot-separated name as a dictionary lookup.  This is useful for
databases, such as CRDS, where the path to the metadata element is
most conveniently stored as a string.  The following two lines are
equivalent::

    print(model['meta.observation.date'])
    print(model.meta.observation.date)

Working with lists
==================

Unlike ordinary Python lists, lists in the schema may be restricted to
only accept a certain set of values.  Items may be added to lists in
two ways: by passing a dictionary containing the desired key/value
pairs for the object, or using the lists special method `item` to
create a metadata object and then assigning that to the list.

For example, suppose the metadata element `meta.transformations` is a
list of transformation objects, each of which has a `type` (string)
and a `coeff` (number) member.  We can assign elements to the list in
the following equivalent ways::

.. doctest-skip::

    >>> trans = model.meta.transformations.item()
    >>> trans.type = 'SIN'
    >>> trans.coeff = 42.0
    >>> model.meta.transformations.append(trans)
    >>> model.meta.transformations.append({'type': 'SIN', 'coeff': 42.0})

When accessing the items of the list, the result is a normal metadata
object where the attributes are type-checked::

.. doctest-skip::
  
    >>> trans = model.meta.transformations[0]
    >>> print(trans)
    <jwst.datamodels.schema.Transformations object at 0x123a810>
    >>> print(trans.type)
    SIN
    >>> trans.type = 42.0
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "site-packages/jwst.datamodels/schema.py", line 672, in __setattr__
         object.__setattr__(self, attr, val)
      File "site-packages/jwst.datamodels/schema.py", line 490, in __set__
         val = self.to_basic_type(val)
      File "site-packages/jwst.datamodels/schema.py", line 422, in to_basic_type
         raise ValueError(e.message)
    ValueError: 42.0 is not of type u'string'

JSON Schema
===========

The `jwst.datamodels` library defines its metadata using `Draft 4 of
the JSON Schema specification
<http://tools.ietf.org/html/draft-zyp-json-schema-04>`_, but
jwst.datamodels uses YAML for the syntax.  A good resource for
learning about JSON schema is the book `Understanding JSON Schema
<http://spacetelescope.github.com/understanding-json-schema>`_.  The
mapping from Javascript to Python concepts (such as Javascript “array”
== Python “list”) is added where applicable.

In addition to the standard JSON Schema keywords, ``jwst.datamodels``
also supports the following additional keywords.

Arrays
''''''

The following keywords have to do with validating n-dimensional arrays:

- ``ndim``: The number of dimensions of the array.

- ``max_ndim``: The maximum number of dimensions of the array.

- ``datatype``: For defining an array, ``datatype`` should be a string.
  For defining a table, it should be a list.

- **array**: ``datatype`` should be one of the following strings,
  representing fixed-length datatypes:

  bool8, int8, int16, int32, int64, uint8, uint16, uint32,
  uint64, float16, float32, float64, float128, complex64,
  complex128, complex256

Or, for fixed-length strings, an array ``[ascii, XX]`` where
``XX`` is the maximum length of the string.

(Datatypes whose size depend on the platform are not supported
since this would make files less portable).

- **table**: ``datatype`` should be a list of dictionaries.  Each
  element in the list defines a column and has the following keys:

  - ``datatype``: A string to select the type of the column.
    This is the same as the ``datatype`` for an array (as
    described above).

  - ``name`` (optional): An optional name for the column.

  - ``shape`` (optional): The shape of the data in the column.
    May be either an integer (for a single-dimensional shape),
    or a list of integers.

FITS-specific Schema Attributes
'''''''''''''''''''''''''''''''

`jwst.datamodels` also adds some new keys to the schema language in
order to handle reading and writing FITS files.  These attributes all
have the prefix ``fits_``.

- ``fits_keyword``: Specifies the FITS keyword to store the value in.
  Must be a string with a maximum length of 8 characters.

- ``fits_hdu``: Specifies the FITS HDU to store the value in.  May be
  a number (to specify the nth HDU) or a name (to specify the
  extension with the given ``EXTNAME``).  By default this is set to 0,
  and therefore refers to the primary HDU.

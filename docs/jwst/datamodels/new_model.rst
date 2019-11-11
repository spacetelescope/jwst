.. -*- coding: utf-8 -*-


Creating a new model
====================

This tutorial describes the steps necessary to define a new model type
using `jwst.datamodels`.

For further reading and details, see the reference materials in
:ref:`metadata`.

In this tutorial, we'll go through the process of creating a new type
of model for a file format used for storing the bad pixel mask for
JWST's MIRI instrument.  This file format has a 2D array containing a
bit field for each of the pixels, and a table describing what each of
the bits in the array means.

.. note::

  While an attempt is made to present a real-world example here, it
  may not reflect the actual final format of this file type, which is
  still subject to change at the time of this writing.

This example will be built as a third-party Python package, i.e. not
part of `jwst.datamodels` itself.  Doing so adds a few extra wrinkles
to the process, and it's most helpful to show what those wrinkles are.
To skip ahead and just see the example in its entirety, see the
``examples/custom_model`` directory within the `jwst.datamodels` source
tree.

Directory layout
----------------

The bare minimum directory layout for a Python package that creates a
custom model is as below::

  .
  |-- lib
  |   |--- __init__.py
  |   |--- bad_pixel_mask.py
  |   |--- schemas
  |   |--- bad_pixel_mask.schema.yaml
  |   |--- tests
  |       |--- __init__.py
  |       |--- test_bad_pixel_mask.py
  |       |--- data
  |       |--- bad_pixel_mask.fits
  |--- setup.py

The main pieces are the new schema in ``bad_pixel_mask.schema.yaml``,
the custom model class in ``bad_pixel_mask.py``, a distutils-based
`setup.py` file to install the package, and some unit tests and
associated data.  Normally, you would also have some code that *uses*
the custom model included in the package, but that isn't included in
this minimal example.

The schema file
----------------

Let's start with the schema file, ``bad_pixel_mask.schema.yaml``.
There are a few things it needs to do:

   1) It should contain all of the core metadata from the core schema
      that ships with `jwst.datamodels`.  In JSON Schema parlance, this
      schema "extends" the core schema.  In object-oriented
      programming terminology, this could be said that our schema
      "inherits from" the core schema.  It's all the same thing.

   2) Define the pixel array containing the information about each of
      the bad pixels.  This will be an integer for each pixel where
      each bit is ascribed a particular meaning.

   3) Define a table describing what each of the bit fields in the
      pixel array means.  This will have three columns: one for the
      bit field's number (a power of 2), one for a name token to
      identify it, and one with a human-readable description.

At the top level, every JSON schema must be a mapping (dictionary) of
type "object", and should include the core schema:

.. code-block:: yaml

  allOf:
     - $ref: "http://jwst.stsci.edu/schemas/core.schema.yaml"
     - type: object
       properties:
          ...

There's a lot going on in this one item.  ``$ref`` declares the schema
fragment that we want to include (the "base class" schema).  Here, the
``$ref`` mapping causes the system to go out and fetch the content at
the given URL, and then replace the mapping with that content.

The ``$ref`` URL can be a relative URL, in which case it is relative
to the schema file where ``$ref`` is used.  In our case, however, it's
an absolute URL.  Before you visit that URL to see what's there, I'll
save you the trouble: there is nothing at that HTTP address.  The host
``jwst.stsci.edu`` is recognized as a "special" address by the
system that causes the schema to be looked up alongside installed
Python code.  For example, to refer to a (hypothetical)
``my_instrument`` schema that ships with a Python package called
``astroboy``, use the following URL::

  http://jwst.stsci.edu/schemas/astroboy/my_instrument.schema.yaml

The "package" portion may be omitted to refer to schemas in the
`jwst.datamodels` core, which is how we arrive at the URL we're using
here::

  http://jwst.stsci.edu/schemas/core.schema.yaml

.. note::

   At some time in the future, we will actually be hosting schemas at
   a URL similar to the one above.  This will allow schemas to be
   shared with tools built in languages other than Python.  Until we
   have that hosting established, this works quite well and does not
   require any coordination among Python packages that define new
   models.  Keep an eye out if you use this feature, though -- the
   precise URL used may change.

The next part of the file describes the array data, that is, things
that are Numpy arrays on the Python side and images or tables on the
FITS side.

First, we describe the main ``"dq"`` array.  It's declared to be
2-dimensional, and each element is an unsigned 32-bit integer:

.. code-block:: yaml

    properties:
      dq:
        title: Bad pixel mask
        fits_hdu: DQ
        default: 0
        ndim: 2
        datatype: uint16

The next entry describes a table that will store the mapping between
bit fields and their meanings.  This table has four columns:

   - ``BIT``: The value of the bit field (a power of 2)

   - ``VALUE``: The value resulting when raising 2 to the BIT power

   - ``NAME``: The name used to refer to the bit field

   - ``DESCRIPTION``: A longer, human-readable description of the bit field

.. code-block:: yaml

        dq_def:
          title: DQ flag definitions
          fits_hdu: DQ_DEF
          dtype:
            - name: BIT
              datatype: uint32
            - name: VALUE
              datatype: uint32
            - name: NAME
              datatype: [ascii, 40]
            - name: DESCRIPTION
              datatype: [ascii, 80]


And finally, we add a metadata element that is specific to this
format.  To avoid recomputing it repeatedly, we'd like to store a sum
of all of the "bad" (i.e. non-zero) pixels stored in the bad pixel
mask array.  In the model, we want to refer to this value as
``model.meta.bad_pixel_count``.  In the FITS file, lets store this in
the primary header in a keyword named ``BPCOUNT``:

.. code-block:: yaml

        meta:
          properties:
            bad_pixel_count:
              type: integer
              title: Total count of all bad pixels
              fits_keyword: BPCOUNT


That's all there is to the schema file, and that's the hardest part.

The model class
----------------

Now, let's see how this schema is tied in with a new Python class for
the model.

First, we need to import the `DataModel` class, which is the base
class for all models::

  from jwst.datamodels import DataModel

Then we create a new Python class that inherits from `DataModel`, and
set its `schema_url` class member to point to the schema that we just
defined above::

  class MiriBadPixelMaskModel(DataModel):
      schema_url = "bad_pixel_mask.schema.yaml"

Here, the `schema_url` has all of the "magical" URL abilities
described above when we used the ``$ref`` feature.  However, here we
are using a relative URL.  In this case, it is relative to the file in
which this class is defined, with a small twist to avoid intermingling
Python code and schema files: It looks for the given file in a
directory called ``schemas`` inside the directory containing the
Python module in which the class is defined.

As an alternative, we could just as easily have said that we want to
use the ``image`` schema from the core without defining any extra
elements, by setting `schema_url` to::

  schema_url = "http://jwst.stsci.edu/schemas/image.schema.yaml"

.. note::

  At this point you may be wondering why both the schema and the class
  have to inherit from base classes.  Certainly, it would have been
  more convenient to have the inheritance on the Python side
  automatically create the inheritance on the schema side (or vice
  versa).  The reason we can't is that the schema files are designed
  to be language-agnostic: it is possible to use them from an entirely
  different implementation of the `jwst.datamodels` framework possibly
  even written in a language other than Python.  So the schemas need
  to "stand alone" from the Python classes.  It's certainly possible
  to have the schema inherit from one thing and the Python class
  inherit from another, and the `jwst.datamodels` framework won't and
  can't really complain, but doing that is only going to lead to
  confusion, so just don't do it.

Within this class, we'll define a constructor.  All model constructors
must take the highly polymorphic ``init`` value as the first argument.
This can be a file, another model, or all kinds of other things.  See
the docstring of `jwst.datamodels.DataModel.__init__` for more
information.  But we're going to let the base class handle that
anyway.

The rest of the arguments are up to you, but generally it's handy to
add a couple of keyword arguments so the user can data arrays when
creating a model from scratch.  If you don't need to do that, then
technically writing a new constructor for the model is optional:

.. code-block:: python

    def __init__(self, init=None, dq=None, dq_def=None, **kwargs):
        """
        A data model to represent MIRI bad pixel masks.

        Parameters
        ----------
        init : any
            Any of the initializers supported by `~jwst.datamodels.DataModel`.

        dq : numpy array
            The data quality array.

        dq_def : numpy array
            The data quality definitions table.
        """
        super(MiriBadPixelMaskModel, self).__init__(init=init, **kwargs)

        if dq is not None:
            self.dq = dq

        if dq_def is not None:
            self.dq_def = dq_def


The ``super..`` line is just the standard Python way of calling the
constructor of the base class.  The rest of the constructor sets the
arrays on the object if any were provided.

The other methods of your class may provide additional conveniences on
top of the underlying file format.  This is completely optional and if
your file format is supported well enough by the underlying schema
alone, it may not be necessary to define any extra methods.

In the case of our example, it would be nice to have a function that,
given the name of a bit field, would return a new array that is `True`
wherever that bit field is true in the main mask array.  Since the
order and content of the bit fields are defined in the `dq_def`
table, the function should use it in order to do this work:

.. code-block:: python

    def get_mask_for_field(self, name):
        """
        Returns an array that is `True` everywhere a given bitfield is
        True in the mask.

        Parameters
        ----------
        name : str
            The name of the bit field to retrieve

        Returns
        -------
        array : boolean numpy array
            `True` everywhere the requested bitfield is `True`.  This
            is the same shape as the mask array.  This array is a copy
            and changes to it will not affect the underlying model.
        """
        # Find the field value that corresponds to the given name
        field_value = None
        for value, field_name, title in self.dq_def:
            if field_name == name:
                field_value = value
                break
        if field_value is None:
            raise ValueError("Field name {0} not found".format(name))

        # Create an array that is `True` only for the requested
        # bit field
        return self.dq & field_value

One thing to note here: this array is semantically a "copy" of the
underlying data.  Most Numpy arrays in the model framework are
mutable, and we expect that changing their values will update the
model itself, and be saved out by subsequent saves to disk.  Since the
array we are returning here has no connection back to the model's main
data array (``mask``), it's helpful to remind the user of that in the
docstring, and not present it as a member or property, but as a getter
function.

.. note::

   Since handling bit fields like this is such a commonly useful
   thing, it's possible that this functionality will become a part of
   `jwst.datamodels` itself in the future.  However, this still stands
   as a good example of something someone may want to do in a custom
   model class.

Lastly, remember the ``meta.bad_pixel_count`` element we defined
above?  We need some way to make sure that whenever the file is
written out that it has the correct value.  The model may have been
loaded and modified.  For this, `DataModel` has the `on_save` method
hook, which may be overridden by the subclass to add anything that
should happen just before saving:

.. code-block:: python

    def on_save(self, path):
        super(MiriBadPixelMaskModel, self).on_save(path)

        self.meta.bad_pixel_count = np.sum(self.mask != 0)

Note that here, like in the constructor, it is important to "chain up"
to the base class so that any things that the base class wants to do
right before saving also happen.

The `setup.py` script
---------------------

Writing a distutils `setup.py` script is beyond the scope of this
tutorial but it's worth noting one thing.  Since the schema files are
not Python files, they are not automatically picked up by distutils,
and must be included in the ``package_data`` option.  A complete, yet
minimal, ``setup.py`` is presented below:

.. code-block:: python

  #!/usr/bin/env python

  from distutils.core import setup

  setup(
      name='custom_model',
      description='Custom model example for jwst.datamodels',
      packages=['custom_model', 'custom_model.tests'],
      package_dir={'custom_model': 'lib'},
      package_data={'custom_model': ['schemas/*.schema.yaml'],
                    'custom_model.tests' : ['data/*.fits']}
      )

Using the new model
-------------------

The new model can now be used.  For example, to get the locations of
all of the "hot" pixels::

   from custom_model.bad_pixel_mask import MiriBadPixelMaskModel

   with MiriBadPixelMaskModel("bad_pixel_mask.fits") as dm:
       hot_pixels = dm.get_mask_for_field('HOT')

A table-based model
-------------------

In addition to n-dimensional data arrays, models can also contain tabular
data. For example, the photometric correction reference file used in the
JWST calibration pipeline consists of a table with several columns. The schema
file for one of these models looks like this:

.. code-block:: yaml

    title: NIRISS SOSS photometric flux conversion data model
    allOf:
      - $ref: "referencefile.schema.yaml"
      - $ref: "keyword_exptype.schema.yaml"
      - $ref: "keyword_pexptype.schema.yaml"
      - $ref: "keyword_pixelarea.schema.yaml"
      - type: object
        properties:
          phot_table:
            title: Photometric flux conversion factors table
            fits_hdu: PHOTOM
            datatype:
              - name: filter
                datatype: [ascii, 12]
              - name: pupil
                datatype: [ascii, 15]
              - name: order
                datatype: int16
              - name: photmj
                datatype: float32
              - name: uncertainty
                datatype: float32
              - name: nelem
                datatype: int16
              - name: wavelength
                datatype: float32
                ndim: 1
              - name: relresponse
                datatype: float32
                ndim: 1
              - name: reluncertainty
                datatype: float32
                ndim: 1

In this particular table the first 6 columns contain scalar entries of types
string, float, and integer. The entries in the final 3 columns, on the other
hand, contain 1-D float arrays (vectors). The "ndim" attribute is used to
specify the number of dimensions the arrays are allowed to have.

The corresponding python module containing the data model class is quite
simple:

.. code-block:: python

    class NisSossPhotomModel(ReferenceFileModel):
        """
        A data model for NIRISS SOSS photom reference files.
        """
        schema_url = "nissoss_photom.schema"


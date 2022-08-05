.. _datamodels:

About models
============

The purpose of the data model is to abstract away the peculiarities of
the underlying file format.  The same data model may be used for data
created from scratch in memory, or loaded from FITS or ASDF files or
some future file format.


Hierarchy of models
-------------------

There are different data model classes for different kinds of data.

One model instance, many arrays
-------------------------------

Each model instance generally has many arrays that are associated with
it.  For example, the `ImageModel` class has the following arrays
associated with it:

    - `data`: The science data
    - `dq`: The data quality array
    - `err`: The error array

The shape of these arrays must be broadcast-compatible.  If you try to
assign an array to one of these members that is not
broadcast-compatible with the data array, an exception is raised.

Working with models
===================

Creating a data model from scratch
----------------------------------

To create a new `ImageModel`, just call its constructor.  To create a
new model where all of the arrays will have default values, simply
provide a shape as the first argument::

    from jwst.datamodels import ImageModel
    with ImageModel((1024, 1024)) as im:
        ...

In this form, the memory for the arrays will not be allocated until
the arrays are accessed.  This is useful if, for example, you don’t
need a data quality array -- the memory for such an array will not be
consumed::

  # Print out the data array.  It is allocated here on first access
  # and defaults to being filled with zeros.
  print(im.data)

If you already have data in a numpy array, you can also create a model
using that array by passing it in as a data keyword argument::

    data = np.empty((50, 50))
    dq = np.empty((50, 50))
    with ImageModel(data=data, dq=dq) as im:
        ...

Creating a data model from a file
---------------------------------

The `jwst.datamodels.open` function is a convenient way to create a
model from a file on disk.  It may be passed any of the following:

    - a path to a FITS file
    - a path to an ASDF file
    - a `astropy.io.fits.HDUList` object
    - a readable file-like object

The file will be opened, and based on the nature of the data in the
file, the correct data model class will be returned.  For example, if
the file contains 2-dimensional data, an `ImageModel` instance will be
returned.  You will generally want to instantiate a model using a
`with` statement so that the file will be closed automatically when
exiting the `with` block.

::

    from jwst import datamodels
    with datamodels.open("myimage.fits") as im:
        assert isinstance(im, datamodels.ImageModel)

If you know the type of data stored in the file, or you want to ensure
that what is being loaded is of a particular type, use the constructor
of the desired concrete class.  For example, if you want to ensure
that the file being opened contains 2-dimensional image data::

    from jwst.datamodels import ImageModel
    with ImageModel("myimage.fits") as im:
        # raises exception if myimage.fits is not an image file
        pass

This will raise an exception if the file contains data of the wrong
shape.

Saving a data model to a file
-----------------------------

Simply call the `save` method on the model instance.  The format to
save into will either be deduced from the filename (if provided) or
the `format` keyword argument::

    im.save("myimage.fits")

.. note::

   Unlike ``astropy.io.fits``, `save` always clobbers the output file.


Copying a model
---------------

To create a new model based on another model, simply use its `copy`
method.  This will perform a deep-copy: that is, no changes to the
original model will propagate to the new model::

    new_model = old_model.copy()

It is also possible to copy all of the known metadata from one
model into a new one using the update method::

    new_model.update(old_model)

History information
-------------------

Models contain a list of history records, accessed through the
`history` attribute.  This is just an ordered list of strings --
nothing more sophisticated.

To get to the history::

    entries = model.history
    for entry in entries:
      pass

To add an entry to the history, first create the entry by calling
`stdatamodels.util.create_history_entry` and appending the entry to the model
history::

    import stdatamodels
    entry = stdatamodels.util.create_history_entry("Processed through the frobulator step")
    model.history.append(entry)

These history entries are stored in ``HISTORY`` keywords when saving
to FITS format. As an option, history entries can contain a dictionary
with a description of the software used. The dictionary must have the
following keys:

  ``name``: The name of the software
  ``author``: The author or institution that produced the software
  ``homepage``: A URI to the homepage of the software
  ``version``: The version of the software

The calling sequence to create  a history entry with the software
description is::

  entry =  stdatamodels.util.create_history_entry(description, software=software_dict)

where the second argument is the dictionary with the keywords
mentioned.

Looking at the contents of a model
----------------------------------

Use ``model.info()`` to look at the contents of a data model. It renders
the underlying ASDF tree starting at the root or a specified ``node``.
The number of displayed rows is controlled by the ``max_row`` argument::

  im.info()
  root.tree (AsdfObject)
  ├─asdf_library (Software)
  │ ├─author (str): Space Telescope Science Institute
  │ ├─homepage (str): http://github.com/spacetelescope/asdf
  │ ├─name (str): asdf
  │ └─version (str): 2.5.2a1.dev12+g12aa460
  ├─history (dict)
  │ └─extensions (list) ...
  ├─data (ndarray): shape=(2048, 2048), dtype=float32
  ├─dq (ndarray): shape=(2048, 2048), dtype=uint32
  ├─err (ndarray): shape=(2048, 2048), dtype=float32
  ├─meta (dict)
  │ ├─aperture (dict) ...
  │ ├─bunit_data (str): DN/s
  │ ├─bunit_err (str): DN/s
  │ ├─cal_step (dict) ...
  │ ├─calibration_software_revision (str): 3bfd782b
  │ ├─calibration_software_version (str): 0.14.3a1.dev133+g3bfd782b.d20200216
  │ ├─coordinates (dict) ...
  │ └─28 not shown
  ├─var_poisson (ndarray): shape=(2048, 2048), dtype=float32
  ├─var_rnoise (ndarray): shape=(2048, 2048), dtype=float32
  └─extra_fits (dict) ...
  Some nodes not shown.


Searching a model
-----------------

``model.search()`` can be used to search the ASDF tree by ``key`` or
``value``::

  im.search(key='filter')

  root.tree (AsdfObject)
  └─meta (dict)
  ├─instrument (dict)
  │ └─filter (str): F170LP
  └─ref_file (dict)
    └─filteroffset (dict)



Converting from ``astropy.io.fits``
===================================

This section describes how to port code that uses ``astropy.io.fits``
to use `jwst.datamodels`.

.. _datamodels-open:

Opening a file
--------------

Instead of::

    astropy.io.fits.open("myfile.fits")

use::

    from jwst.datamodels import ImageModel
    with ImageModel("myfile.fits") as model:
        ...

In place of `ImageModel`, use the type of data one expects to find in
the file.  For example, if spectrographic data is expected, use
`SpecModel`.  If it doesn't matter (perhaps the application is only
sorting FITS files into categories) use the base class `DataModel`.

An alternative is to use::

    from jwst import datamodels
    with datamodels.open("myfile.fits") as model:
        ...

The `datamodels.open()` method checks if the `DATAMODL` FITS keyword has
been set, which records the DataModel that was used to create the file.
If the keyword is not set, then `datamodels.open()` does its best to
guess the best DataModel to use.

Accessing data
--------------

Data should be accessed through one of the pre-defined data members on
the model (`data`, `dq`, `err`).  There is no longer a need to hunt
through the HDU list to find the data.

Instead of::

    hdulist['SCI'].data

use::

    model.data

Accessing keywords
------------------

The data model hides direct access to FITS header keywords.  Instead,
use the :ref:`metadata` tree.

There is a convenience method, `find_fits_keyword` to find where a
FITS keyword is used in the metadata tree::

    >>> from jwst.datamodels import DataModel
    >>> # First, create a model of the desired type
    >>> model = DataModel()
    >>> model.find_fits_keyword('DATE-OBS')
    [u'meta.observation.date']

This information shows that instead of::

    print(hdulist[0].header['DATE-OBS'])

use::

    print(model.meta.observation.date)

Extra FITS keywords
-------------------

When loading arbitrary FITS files, there may be keywords that are not
listed in the schema for that data model.  These "extra" FITS keywords
are put under the model in the `_extra_fits` namespace.

Under the `_extra_fits` namespace is a section for each header data
unit, and under those are the extra FITS keywords.  For example, if
the FITS file contains a keyword `FOO` in the primary header, its
value can be obtained using::

    model._extra_fits.PRIMARY.FOO

This feature is useful to retain any extra keywords from input files
to output products.

To get a list of everything in `_extra_fits`::

    model._extra_fits._instance

returns a dictionary of of the instance at the model._extra_fits node.

`_instance` can be used at any node in the tree to return a dictionary
of rest of the tree structure at that node.

Environment Variables
---------------------

There are a number of environment variables that affect how models are read.

PASS_INVALID_VALUES
  Used by `~jwst.datamodels.DataModel` when instantiating
  a model from a file. If ``True``, values that do not validate the schema will
  still be added to the metadata. If ``False``, they will be set to ``None``.
  Default is ``False``.

STRICT_VALIDATION
  Used by `~jwst.datamodels.DataModel` when instantiating a model from a file.
  If ``True``, schema validation errors will generate an exception.
  If ``False``, they will generate a warning.
  Default is ``False``.

SKIP_FITS_UPDATE
  Used by `~jwst.datamodels.DataModel` when instantiating a
  model from a FITS file. When ``False``, models opened from FITS files will
  proceed and load the FITS header values into the model. When ``True`` and the
  FITS file has an ASDF extension, the loading/validation of the FITS header
  will be skipped, loading the model only from the ASDF extension. If not
  defined, the instantiation routines will determine whether the loading/validation
  of the FITS header can be skipped or not.

DMODEL_ALLOWED_MEMORY
  Implemented by the utility function
  `jwst.datamodels.util.check_memory_allocation` and used by
  `~jwst.outlier_detection.OutlierDetectionStep` and
  `~jwst.resample.ResampleStep`. When defined, determines how much of currently
  available memory should be used to instantiated an output resampled image. If
  not defined, no check is made.

  Examples would be: ``1.0`` would allow all available memory to be used. ``0.5``
  would allow only half the available memory to be used.

For flag or boolean variables, any value in ``('true', 't', 'yes', 'y')`` or a
non-zero number, will evaluate as ``True``. Any value in ``('false', 'f', 'no',
'n', '0')`` will evaluate as ``False``. The values are case-insensitive.

All of the environment variables have equivalent function arguments in the API
for the relevant code. The environment variables are used only if explicit
values had not been used in a script. In other words, values in code override
environment variables.

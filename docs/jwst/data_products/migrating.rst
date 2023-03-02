Migrating deprecated products
-----------------------------

On rare occasion, the model schemas are changed in such a way as to
break compatibility with data products produced by earlier versions
of this package.  When these older files are opened the software
will report validation errors:

.. doctest-skip::

  >>> from stdatamodels.jwst import datamodels
  >>> datamodels.open("jw95115001001_02102_00001_nrs1_x1d.fits")
  ...
  ValueError: Column names don't match schema...

In some cases it will be possible to update the file to the
new format using the ``migrate_data`` tool included with this package:
::

    $ migrate_data jw95115001001_02102_00001_nrs1_x1d.fits --in-place

It can also be run on multiple files:
::

    $ migrate_data *_x1d.fits --in-place

Or configured to write updated files to a separate output directory:
::

    $ migrate_data *_x1d.fits --output-dir some/other/directory

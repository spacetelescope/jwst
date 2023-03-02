.. _blender_handbook:

Role of Model Blender
======================

One problem with combining data from multiple exposures stems from not being able
to keep track of what kind of data was used to create the final product.  The
final product only reports one value for each of the metadata attributes from the
model schema used to describe the science data, where each of multiple inputs may
have had different values for each attribute.  The model_blender package solves
this problem by allowing the user to define rules that can be used to determine a
final single value from all the input values for each attribute, using separate
rules for each attribute as appropriate.  This package also creates a FITS binary
table that records the input attribute values for all the input models used to
create the final product, allowing the user to select what attributes to keep in
this table.

This code works by

  - reading in all the input datamodels (either already in-memory or from FITS files)
  - evaluating the rules for each attribute as defined in the model's schema
  - determining from definitions in the input model's schema what attributes to keep in the table
  - applying each attributes rule to the set of input values to determine the final output value
  - updating the output model's metadata with the new values
  - generating the output table with one row for each input model's values


Using model_blender
===================
The model blender package requires

  - all the input models be available
  - the output product has already been generated

Both the input models and output product could be provided as either a datamodel
instance in memory or as the name of a FITS file on disk.  The primary advantage
to working strictly in-memory with datamodel instances comes from minimizing the
amount of disk I/O needed for this operation which can result in significantly
more efficient (read that: faster) processing.

.. note::

  The generated output model will be considered to contain a default
  (or perhaps even empty) set of :ref:`metadata` based on some
  model defined in 
  .. comment out until stdatamodels is released
  .. ref DataModels<datamodels>.
  DataModels.
  This metadata will be replaced
  **in-place** when running :ref:`blender_api`.

The simplest way to run model blender only requires calling a single interface:

.. code-block:: python

  from jwst.model_blender import blendmeta
  blendmeta.blendmodels(product, inputs=input_list)

where

  - `product`: the datamodel (or FITS filename) for the already combined product
  - `input_list`: list of input datamodels or FITS filenames for all inputs used
    to create the `product`


The output product will end up with new metadata attribute values and a new HDRTAB
FITS binary table extension in the FITS file when the product gets saved to disk.


Customizing the behavior
========================
By default, `blendmodels` will not write out the updated product model to disk.
This allows the user or calling program to revise or apply data-specific logic
to redefine the output value for any of the output product's metadata attributes.
For example, when combining multiple images, the WCS information does not represent
any combination of the input WCS attributes.  Instead, the user can have
their own processing code replace the *blended* WCS attributes with one that was
computed separately using a complex, accurate algorithm.  This is, in fact, what
the resample step does to create the final resampled output product whenever it is
called by steps in the JWST pipeline.

Additional control over the behavior of model_blender comes from editing the
schema for the input datamodels where the rules for each attribute are defined.
A sample definition from the core schema demonstrates the basic syntax used for
any model blending definitions::

          time_end:
            title: UTC time at end of exposure
            type: string
            fits_keyword: TIME-END
            blend_rule: last
            blend_table: True

Any attribute without an entry for `blend_rule` will use the default rule of
`first` which selects the first value from all inputs in the order provided as the
final output value.  Any attribute with a `blend_table` rule will insure that
the specific attribute will be included in the output HDRTAB binary table appended
to the product model when it gets written out to disk as a FITS file.

The full set of rules included in the package are described in
:ref:`blender_rules` and include common list/array operations such as
(but not limited to):

  - minimum
  - maximum
  - first
  - last
  - mean
  - zero

These can then be used to customize the output value for any given attribute
should the rule provided by default with the schema installed with the
JWST environment not be correct for the user's input data.  The user can simply
edit the schema definition installed in their JWST environment to apply custom
rules for blending the data being processed.

.. _blender_handbook:

Role of Model Blender
======================

Model Blender performs 2 operations:

  - Accumulate, blend and set metadata when combining several models (
    as is done during resampling) using combination rules defined in the datamodel
    schemas.
  - Record the metadata of input models in a table added to the
    combined model.

One problem with combining data from multiple exposures stems from not being able
to keep track of what kind of data was used to create the final product.  The
final product only reports one value for each of the metadata attributes from the
model schema used to describe the science data, where each of multiple inputs may
have had different values for each attribute.  The model_blender package solves
this problem by allowing the user to define rules that can be used to determine a
final single value from all the input values for each attribute, using separate
rules for each attribute as appropriate.  This package also creates a
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


Using blendmodels
-----------------

The simplest way to run model blender only requires calling a single interface:

.. code-block:: python

  from jwst.model_blender import blendmeta
  blendmeta.blendmodels(product, input_list)

where

  - `product`: the datamodel for the already combined product
  - `input_list`: list of input datamodels used to create the `product`


The output product will end up with new metadata attribute values and a new ``hdrtab``
attribute that will produce a ``HDRTAB`` table extension when the product is saved
to disk.

Using ModelBlender
------------------

Alternatively, metadata blending can be done incrementally by using `ModelBlender`
to accumulate metadata from several models and then blending the result into an
output model.

.. code-block:: python

   from jwst.model_blender import ModelBlender
   blender = ModelBlender()
   for model in input_models:
       blender.accumulate(model)
   blender.finalize_model(product)

This produces a ``product`` identical to a call to `blendmodels` described above.


Customizing the behavior
========================


Combined model metadata
-----------------------

The metadata assigned to the combined (output) model will be determined from:

  - the first combined model metadata
  - any metadata with a ``blend_rule``

Looking at this sample from the core schema::

          start_time:
            title: "[d] exposure start time in MJD"
            type: number
            fits_keyword: EXPSTART
            blend_table: True
            blend_rule: min

Since ``start_time`` uses ``blend_rule: min`` the ``start_time`` for the combined
model will be the minimum of all input models. If ``blend_rule`` was not defined, the
``start_time`` for the combined model would be the ``start_time`` of the first model.

.. warning::

   All blended models must be of the same type.

.. warning::

   Since model blender copies the metadata of the first model, all metadata
   from that model will be added to combined model. This is true even for
   metadata that is not defined in the schema.

.. warning::

   Metadata nested in arrays will not be blended. For example,
   ``meta.background.polynomial_info`` contains a list/array of
   dictionaries. Model blender will not blend this attribute.

Input model metadata Table
--------------------------

As seen in the above example, ``start_time`` defines ``blend_table: True``.
This tells model blender to include this attribute in the table added to the
combined model that lists attributes from input models. If the attribute
has a ``fits_keyword``, the resulting column will use the keyword for the
column name (if not defined, the attribute name will be used). So for
the above example, the resulting table will contain an ``EXPSTART`` column
(and not a ``meta.exposure.start_time`` column) with a row listing
the ``start_time`` for each input model.

If an input model is missing an attribute a ``nan`` will be stored in the
corresponding table cell.

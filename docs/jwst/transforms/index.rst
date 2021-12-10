==========
Transforms
==========

.. toctree::
   :maxdepth: 2


ASDF Schema Definitions
-----------------------

This package defines `ASDF <http://asdf-standard.readthedocs.io>`__ schemas
that are used for validating the representation of transforms in the ASDF
format. These schemas contain useful documentation about the associated types,
and can also be used by other implementations that wish to interoperate with
these transform definitions.

.. asdf-autoschemas::
   :schema_root: ../jwst/transforms/resources/schemas
   :standard_prefix: stsci.edu/jwst_pipeline

   grating_equation-1.0.0
   gwa_to_slit-1.0.0
   logical-1.0.0
   miri_ab2slice-1.0.0
   nircam_grism_dispersion-1.0.0
   niriss_grism_dispersion-1.0.0
   niriss_soss-1.0.0
   refraction_index_from_prism-1.0.0
   rotation_sequence-1.0.0
   slit_to_msa-1.0.0
   snell-1.0.0


.. automodapi:: jwst.transforms
.. automodapi:: jwst.transforms.models

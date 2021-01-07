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
   :schema_root: ../jwst/transforms/schemas
   :standard_prefix: stsci.edu/jwst_pipeline

   coords-0.7.0
   grating_equation-0.7.0
   gwa_to_slit-0.7.0
   logical-0.7.0
   miri_ab2slice-0.7.0
   nircam_grism_dispersion-0.7.0
   niriss_grism_dispersion-0.7.0
   niriss_soss-0.7.0
   refraction_index_from_prism-0.7.0
   rotation_sequence-0.7.0
   slit_to_msa-0.7.0
   snell-0.7.0
   v23tosky-0.7.0


.. automodapi:: jwst.transforms
.. automodapi:: jwst.transforms.models

Description
===========

:Class: `jwst.nsclean.NSCleanStep`
:Alias: nsclean

.. note::

   This step, which is intended to be called in the
   :ref:`calwebb_spec2 <calwebb_spec2>` pipeline on NIRSpec rate data,
   is now deprecated.  In future builds, it will be replaced by
   the :ref:`clean_flicker_noise <clean_flicker_noise_step>`
   step, called in the :ref:`calwebb_detector1 <calwebb_detector1>`
   pipeline on ramp data.

Overview
========
This step currently runs a version of the
:ref:`clean_flicker_noise <clean_flicker_noise_step>` algorithm,
with slightly different parameters and default values, intended
to be backwards-compatible with the previous implementation of
the step developed specifically for NIRSpec data. The only
significant difference is that in this step, the cleaning is
performed directly on the rate data, rather than group images.

See the :ref:`clean_flicker_noise <clean_flicker_noise_step>`
documentation for more information, or the
:ref:`nsclean step arguments <nsclean_arguments>` for the default
values used for this step.

If this step is run as part of the 
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline when processing an
association that includes background or imprint images, these
images will be processed using the same ``nsclean`` parameters 
as the science image.  If a user supplied mask is provided, it
will be used for all images in an association. Otherwise, separate
masks will be calculated for each image using the same input
parameters.

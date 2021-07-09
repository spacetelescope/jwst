.. _wfss_contam_step_args:

Step Arguments
==============
The ``wfss_contam`` step uses the following optional arguments.

``--save_simulated_image``
  A boolean indicating whether the full-frame simulated grism image containing all
  simulated spectra within the field-of-view should be saved to a file. The file
  name uses a product type suffix of "simul".
  Defaults to ``False``.

``--save_contam_images``
  A boolean indicating whether the estimated contamination images for each source
  cutout should be saved to a file. The file name uses a product type suffix of "contam".
  The resulting file has one SCI extension for each source contained in the input
  grism image.
  Defaults to ``False``.

``--maximum_cores``
  The fraction of available cores that will be
  used for multi-processing in this step. The default value is 'none' which does not use
  multi-processing. The other options are 'quarter', 'half', and 'all'. Note that these
  fractions refer to the total available cores and on most CPUs these include physical
  and virtual cores.

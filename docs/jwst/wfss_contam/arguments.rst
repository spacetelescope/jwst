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

``--orders``
  A list indicating which grism orders to simulate. The default value is None, which
  means all orders that are defined in the wavelength and specwcs reference files
  for that instrument will be simulated.
  To specify a single order from the command line, use e.g. "0," or "1,"
  (the comma allows the code to identify this as a list).

``--magnitude_limit``
  A float indicating the magnitude limit for sources to be included in the simulation.
  The magnitude limit is taken from the input source catalog's isophotal AB magnitude column.
  The limit is scaled according to the relative sensitivity of each spectral order, such that
  fewer sources are included in orders with lower sensitivity.
  The default value is None, which means no magnitude limit is applied to any of the orders
  and all sources are included.

``--wl_oversample``
  Indicates the oversampling factor for the wavelength grid used in the
  simulation of the dispersed spectra. Defaults to 2.

``--max_pixels_per_chunk``
  Sets the maximum number of direct image pixels to run through the grism transforms at once.
  Decreasing this value will reduce the memory usage of the step, but will
  increase the runtime. Note that if a single source has more pixels than the value of 
  ``max_pixels_per_chunk``, that source will be skipped. Defaults to 50,000.

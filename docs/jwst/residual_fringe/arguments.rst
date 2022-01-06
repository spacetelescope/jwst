.. _resid_fringe_arguments:

Step Arguments
==============

The residual fringe step has two step arguments that can be used to specify wavelength regions in which no correction will be determined.
The two arguments give lists of minimum and maximum wavelength values, respectively, for the regions to be ignored.
The two lists must contain an equal number of elements.

``ignore_region_min`` [float, default = None]
  The minimum wavelengths for the region(s) to be ignored, given as a comma-separated list.

``ignore_region_max`` [float, default = None]
  The maximum wavelengths for the region(s) to be ignored, given as a comma-separated list.


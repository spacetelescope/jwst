.. _resid_fringe_arguments:

Step Arguments
==============

The residual fringe step has two step arguments that can be used to ignore determining a residual fringe
correction in a wavelength range. These two parameters must have the same number of elements. 

``ignore_region_min [list]``
  This is the minumum wavelength of the region to ignore.

  ``ignore_region_max [list]``
  This is the maximum  wavelength of the region to ignore.

Regions are defined by corresponding minimum and maximum values. For example ``ignore_region_min``="4.79," and
``ignore_region_max``="4.82", would ignore the wavelengths between 4.79 to 4.82 and not determine a residual fringe correction for data in this wavelength range. 


.. _white_light_arguments:

Step Arguments
==============
The ``white_light`` step has two step-specific arguments to allow
wavelength limits during the flux summation. One or both may be specified.

``--min_wavelength`` (float, default=None)
  If specified, the ``white_light`` step will sum
  from the specified wavelength to either a specified ``max_wavelength``
  or the end of the flux array.

``--max_wavelength`` (float, default=None)
  If specified, the ``white_light`` step will sum
  from either a specified ``min_wavelength`` or the beginning of the
  flux array to the specified wavelength.

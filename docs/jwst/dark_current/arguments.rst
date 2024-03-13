Step Arguments
==============

The dark current step has these step-specific arguments:

*  ``--dark_output`` (string, default=None)

If the ``dark_output`` argument is given with a filename for its value,
the frame-averaged dark data that are created within the step will be
saved to that file.

*  ``--average_dark_current`` (float, default=None)

The average dark current for this detector in units of electrons per second.
This will be used to calculate the Poisson noise contribution due to the dark
current. This parameter is be a scalar quantity; if a 2D array is desired to
describe the dark current pixel-to-pixel, this must be specified in the dark
reference file.

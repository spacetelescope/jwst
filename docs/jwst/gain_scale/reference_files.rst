Reference File
==============

The GAIN reference file is used by both the gain_scale step and ramp fitting.
For gain_scale, the only purpose of the reference file is to retrieve the
`GAINFACT` keyword value from its header.  If the `ramp_fit` step, which also
uses the gain reference file, succeeded in finding the `GAINFACT` keyword in
this reference file, it will store the value in the `GAINFACT` keyword in the
science data, in which case the `gain_scale` step will not reload the gain
reference file.

.. include:: ../includes/gain_format.rst


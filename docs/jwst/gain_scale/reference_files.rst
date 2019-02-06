Reference File
==============

The ``gain_scale`` step uses the GAIN reference file.
It requires this reference file only to
get the value of the "GAINFACT" keyword in the header of the file.
This is the value used to rescale the science data. The ``ramp_fit``
step also uses the GAIN reference file and if it succeeded in finding
the "GAINFACT" keyword when it was executed, it will have already
stored the keyword value in the science data, for later use by the
``gain_scale`` step. In this case the ``gain_scale`` step will not read
the GAIN reference file again when it runs.

.. include:: ../references_general/gain_reffile.inc

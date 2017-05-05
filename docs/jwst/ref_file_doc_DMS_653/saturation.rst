==========
Saturation
==========

The saturation level of each pixel (in units of DN) is stored as a 2-D image
in the reference file. For each group in the science data file, the pipeline
compares each pixel's DN value with its saturation level.  If the pixel exceeds
the saturation level, then the SATURATED flag is set for that pixel in the
corresponding plane of the GROUPDQ array – and in all subsequent planes.
No saturation check is performed on pixels for which the flag NO_SAT_CHECK is set.

.. toctree::
   :maxdepth: 4

   saturation_reference_files

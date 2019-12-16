Reference Files
===============
The ``straylight`` step uses the :ref:`REGIONS <regions_reffile>` reference
file, which stores locations of the slice regions on the detectors. This
reference file provides 2-D detector images in which each pixel is set to
the number of the corresponding slice (or 0 for inter-slice pixels) at each of 9
different throughput levels ranging from 10% - 90%.  While ``assign_wcs``
uses a fairly exclusive slice mask (selecting only pixels with high
throughput for science analysis), the ``straylight`` step uses a very
inclusive 20% threshhold to define slice pixels and obtain a cleaner
sample of inter-slice pixels.

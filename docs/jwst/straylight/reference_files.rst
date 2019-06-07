Reference Files
===============

:ref:`REGIONS: <regions_reffile>` Stores location of the regions on the detector

The ``straylight`` step uses the :ref:`REGIONS <regions_reffile>` reference
file to determine the location of the MRS slices on the detector. This
reference file provides 2-D detector images in which each pixel is set to
the number of the corresponding slice (or 0 for inter-slice pixels) at each of 9
different throughput levels ranging from 10% - 90%.  While ``assign_wcs``
uses a fairly exclusive slice mask (selecting only pixels with high
throughput for science analysis), the ``straylight`` step uses a very
inclusive 20% threshhold to define slice pixels and obtain a cleaner
sample of inter-slice pixels.

By default the ``straylight`` step does not use a reference file.
It instead uses meta information added by the ``assign_wcs`` step that
maps each pixel to an IFU slice number or the regions between slices
within the 2-D slice image.

If the step option ``method`` is set to "Nearest", a more simplistic
algorithm is used to compute the correction, which requires the use of
the ``STRAYMASK`` reference file. It contains a simple mask indicating
which pixels lie within the inter-slice regions of the image.

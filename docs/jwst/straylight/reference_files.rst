Reference Files
===============

By default the ``straylight`` step does not use a reference file.
It instead uses meta information added by the ``assign_wcs`` step that
maps each pixel to an IFU slice number or the regions between slices
within the 2-D slice image.

If the step option ``method`` is set to "Nearest", a more simplistic
algorithm is used to compute the correction, which requires the use of
the ``STRAYMASK`` reference file. It contains a simple mask indicating
which pixels lie within the inter-slice regions of the image.

.. include:: ../references_general/straymask_reffile.rst

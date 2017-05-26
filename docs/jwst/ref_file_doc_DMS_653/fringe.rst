======
Fringe
======

This step applies a fringe correction to the SCI data of an input data set by dividing the SCI and ERR arrays by a fringe
reference image.  In particular, the SCI array from the fringe reference file is divided into the SCI and ERR arrays of
the science data set. Only pixels that have valid values in the SCI array of the reference file will be corrected.
This correction is applied only to MIRI MRS (IFU) mode exposures, which are always single full-frame 2-D images.

.. toctree::
   :maxdepth: 4

   fringe_reference_files

Description
============

Overview
--------
The ``extract_2d`` step extracts 2D arrays from spectral images. The extractions
are performed within all of the SCI, ERR, and DQ arrays of the input image
model. The extracted SCI, ERR, and DQ arrays are stored as one or more ``slit``
objects in an output MultiSlitModel.

Assumptions
-----------
This step uses region information stored in the WCS object of the data model,
which is populated by the ``assign_wcs`` step. Hence the ``assign_wcs`` step
must be applied to the science exposure before running this step.

Algorithm
---------
The algorithm currently works only for NIRSpec Fixed Slit and NIRISS Single
Object Slitless Spectroscopy (SOSS) mode exposures.

* For a NIRSPEC Fixed Slit exposure using either full frame or ALLSLITS
  subarray readout, one or more of the regions containing the 5 fixed slits
  will be extracted. If the step parameter ``which_subarray`` is left
  unspecified, the default behavior is to extract regions for all 5 of the
  available slits on detector NRS2 and the 4 slits on detector NRS1.
  The region for just one of the slits may be extracted by specifying the slit
  name with the ``which_subarray`` argument, using one of the following
  accepted names: ``S1600A1``, ``S200A1``, ``S200A2``, ``S200B1``, ``S400A1``.
  Note that the ``S200B1`` slit can only be extracted from exposures on the
  NRS2 detector.

* For a NIRISS full-frame SOSS observation, a single slit will be extracted,
  corresponding to the same region in the image as the NIRISS 2048x256
  subarray, which is the default for SOSS exposures.

The corner locations of the regions to be extracted are determined from the
region definitions contained in the exposure's WCS, which has definitions for
each of the available slits for the instrument mode. The corners of each
extraction region are computed from the lower and upper x/y limits of the
region's vertices in the WCS.

The output MultiSlit data model will have the meta data associated with each
slit region populated with the name and region characteristic for the slits,
corresponding to the FITS keywords ``SLTNAME``, ``SLTSTRT1``, ``SLTSIZE1``,
``SLTSTRT2``, and ``SLTSIZE2``.

Step Arguments
==============
The extract_2d step has one optional argument:

* ``--which_subarray``: name (string value) of a specific slit region to
  extract. The default value of None will cause all known slits for the
  instrument mode to be extracted. Currently only used for NIRspec fixed slit
  exposures.

Reference Files
===============
This step does not require any reference files.


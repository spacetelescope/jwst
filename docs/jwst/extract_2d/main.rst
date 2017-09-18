
Description
============

Overview
--------
The ``extract_2d`` step extracts 2D arrays from spectral images. The extractions
are performed within all of the SCI, ERR, and DQ arrays of the input image
model. It also computes an array of wavelengths. The SCI, ERR, DQ and WAVELENGTH
arrays are stored as one or more ``slit`` objects in an output MultiSlitModel
and saved as separate extensions in the output FITS file.

Assumptions
-----------
This step uses the ``bounding_box`` attribute of the WCS stored in the data model,
which is populated by the ``assign_wcs`` step. Hence the ``assign_wcs`` step
must be applied to the science exposure before running this step.

Algorithm
---------
The step is currently applied only to NIRSpec Fixed Slit and NIRSPEC MSA observations.

If the step parameter ``which_subarray`` is left unspecified, the default behavior is
to extract all slits which fall within a detector. Only one slit may be extracted by
specifying the slit name with the ``which_subarray`` argument, using one of the following
accepted names: ``S1600A1``, ``S200A1``, ``S200A2``, ``S200B1``, ``S400A1`` or ``S200B1``
in the case of NIRSPEC FS exposure or any of the slitlet names in the case of the MSA.

To find out what slits are available for extraction:

  >>> from jwst.assign_wcs import nirspec
  >>> nirspec.get_open_slits(input_model)


The corner locations of the regions to be extracted are determined from the
``bounding_box`` contained in the exposure's WCS, which defines the range of valid inputs
along each axis. The input coordinates are in the image frame, i.e. subarray shifts
are accounted for.

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

* ``--apply_wavecorr``: bool (default is True). Flag indicating whether to apply the
   Nirspec wavelength zero-point correction.

Reference Files
===============
To apply the Nirspec wavelength zero-point correction, this step uses the
WAVECORR reference file. The zero-point correction is applied to observations
with EXP_TYPE of "NRS_FIXEDSLT", "NRS_BRIGHTOBJ" or "NRS_MSASPEC". This is an optional
correction (on by default). It can be turned off by specifying ``apply_wavecorr=False``
when running the step.

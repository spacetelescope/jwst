
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

For WFSS modes in NIRCAM and NIRSS, no ``bounding_box`` has been attached
to the datamodel. This is to keep the WCS flexible enough to be used with any
source catalog that may be associated with the dispersed image. Instead, there
is a helper method that is used to calculate the bounding boxes that contain
the dispersed spectra for each object. One box is made for each order. ``extract2d``
uses the source catalog referenced in the input models meta information to create
the list of objects and their corresponding bounding box, this list is used to make
the 2D cutouts from the dispersed image.

Algorithm
---------
The step is currently applied only to NIRSpec Fixed Slit, NIRSPEC MSA, NIRSPEC TSO,
NIRCAM WFSS and NIRISS WFSS observations.

For NIRSPEC:

If the step parameter ``slit_name`` is left unspecified, the default behavior is
to extract all slits which project on the detector. Only one slit may be extracted by
specifying the slit name with the ``slit_name`` argument, using one of the following
accepted names: ``S1600A1``, ``S200A1``, ``S200A2``, ``S200B1`` or ``S400A1``
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


For NIRCAM WFSS and NIRISS WFSS :

If the step parameter ``grism_objects`` is left unspecified, the default behavior
is to use the source catalog that is specified in the input model's meta information,
``input_model.meta.source_catalog.filename``. Otherwise, a user can submit of list of
``GrismObjects`` that contains information about the objects that should be extracted.
The ``GrismObject`` list can be created automatically by using the method in
``jwst.assign_wcs.utils.create_grism_bbox``. This method also uses the name of the source
catalog saved in the input model's meta information. If it's better to construct a list
of ``GrismObjects`` outside of these, the ``GrismObject`` itself can be imported from
``jwst.transforms.models``.


Step Arguments
==============
The extract_2d step has two optional arguments for NIRSPEC observations:

* ``--slit_name``: name (string value) of a specific slit region to
  extract. The default value of None will cause all known slits for the
  instrument mode to be extracted. Currently only used for NIRspec fixed slit
  exposures.

* ``--apply_wavecorr``: bool (default is True). Flag indicating whether to apply the
   Nirspec wavelength zero-point correction.


For NIRCAM and NIRISS, the extract_2d step has only one optional argument:

* ``--grism_objects``: list (default is empty). A list of ``jwst.transforms.models.GrismObject``.


Reference Files
===============
To apply the Nirspec wavelength zero-point correction, this step uses the
WAVECORR reference file. The zero-point correction is applied to observations
with EXP_TYPE of "NRS_FIXEDSLT", "NRS_BRIGHTOBJ" or "NRS_MSASPEC". This is an optional
correction (on by default). It can be turned off by specifying ``apply_wavecorr=False``
when running the step.

NIRCAM WFSS and NIRISS WFSS observations use the wavelengthrange reference file in order
to construct the bounding boxes around each objects orders. If a list of ``GrismObject``
is supplied, then no reference file is neccessary.

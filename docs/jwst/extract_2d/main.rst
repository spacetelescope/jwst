Description
===========

Overview
--------
The ``extract_2d`` step extracts 2D arrays from spectral images. The extractions
are performed within all of the SCI, ERR, and DQ arrays of the input image
model, as well as any variance arrays that may be present. It also computes an
array of wavelengths to attach to the extracted data. The extracted arrays
are stored as one or more ``slit`` objects in an output MultiSlitModel
and saved as separate tuples of extensions in the output FITS file.

Assumptions
-----------
This step uses the ``bounding_box`` attribute of the WCS stored in the data model,
which is populated by the ``assign_wcs`` step. Hence the ``assign_wcs`` step
must be applied to the science exposure before running this step.

For NIRCam and NIRISS WFSS modes, no ``bounding_box`` has been attached
to the data model. This is to keep the WCS flexible enough to be used with any
source catalog that may be associated with the dispersed image. Instead, there
is a helper method that is used to calculate the bounding boxes that contain
the dispersed spectra for each object. One box is made for each spectral order of
each object. The ``extract2d`` step uses the source catalog referenced in the input
model's meta information to create the list of objects and their corresponding
bounding box. This list is used to make the 2D cutouts from the dispersed image.

NIRCam TSGRISM exposures do not use a source catalog, so instead it relies on the
assumption that the source of interest is located at the aperture reference point.
More details are given below.

Algorithm
---------
This step is currently applied only to NIRSpec Fixed-Slit, NIRSpec MOS, NIRSpec TSO
(BrightObj), NIRCam and NIRISS WFSS, and NIRCam TSGRISM observations.

NIRSpec
+++++++

If the step parameter ``slit_name`` is left unspecified, the default behavior is
to extract all slits that project onto the detector. A single slit may be extracted by
specifying the slit name with the ``slit_name`` argument. In the case of a NIRSpec
fixed-slit exposure the allowed slit names are: "S1600A1", "S200A1", "S200A2", "S200B1"
and "S400A1". For NIRSpec MOS exposures, the slit name is the slitlet number from the
MSA metadata file, corresponding to the value of the "SLTNAME" keyword in FITS products.

To find out what slits are available for extraction:

  >>> from jwst.assign_wcs import nirspec
  >>> nirspec.get_open_slits(input_model) # doctest: +SKIP


The corner locations of the regions to be extracted are determined from the
``bounding_box`` contained in the exposure's WCS, which defines the range of valid inputs
along each axis. The input coordinates are in the image frame, i.e. subarray shifts
are accounted for.

The output MultiSlit data model will have the meta data associated with each
slit region populated with the name and region characteristic for the slits,
corresponding to the FITS keywords "SLTNAME", "SLTSTRT1", "SLTSIZE1",
"SLTSTRT2", and "SLTSIZE2."  Keyword "DISPAXIS" (dispersion direction)
will be copied from the input file to each of the output cutout images.


NIRCam WFSS and NIRISS WFSS
+++++++++++++++++++++++++++

If the step parameter ``grism_objects`` is left unspecified, the default behavior
is to use the source catalog that is specified in the input model's meta information,
``input_model.meta.source_catalog.filename``. Otherwise, a user can submit of list of
``GrismObjects`` that contains information about the objects that should be extracted.
The ``GrismObject`` list can be created automatically by using the method in
``jwst.assign_wcs.utils.create_grism_bbox``. This method also uses the name of the source
catalog saved in the input model's meta information. If it's better to construct a list
of ``GrismObjects`` outside of these, the ``GrismObject`` itself can be imported from
``jwst.transforms.models``.

The dispersion direction will be documented by copying keyword "DISPAXIS"
(1 = horizontal, 2 = vertical) from the input file to the output cutout.


NIRCam TSGRISM
++++++++++++++

There is no source catalog created for TSO grism observations, because no associated
direct images are obtained from which to derive such a catalog. So the ``extract_2d``
step relies on the fact that the source of interest is placed at the aperture reference
point to determine the source location. The aperture reference location, in units of
image x and y pixels, is read from the keywords "XREF_SCI" and "YREF_SCI" in the SCI
extension header of the input image. These values are used to set the source location
for all computations involving the extent of the spectral trace and pixel wavelength
assignments.

NIRCam subarrays used for TSGRISM observations always have their "bottom" edge located
at the physical bottom edge of the detector and vary in size vertically.
The source spectrum trace will always be centered somewhere near row 34 in the vertical
direction (dispersion running parallel to rows) of the dispersed image.
So the larger subarrays just result in a larger region of sky above the spectrum.

For TSGRISM, ``extract_2d`` always produces a cutout that is 64 pixels in height
(cross-dispersion direction), regardless of whether the original image is full-frame
or subarray.
This cutout height is equal to the height of the smallest available subarray (2048 x 64).
This is to allow area within the cutout for sampling the background in later steps,
such as ``extract_1d``. The slit height is a parameter that a user can set
(during reprocessing) to tailor their results, but the entire extent of the image in
the dispersion direction (along the image x-axis) is always included in the cutout.

The dispersion direction is horizontal for this mode, and it will be
documented by copying the keyword "DISPAXIS" (with value 1) from the input file
to the output cutout.


Step Arguments
==============
The ``extract_2d`` step has various optional arguments that apply to certain observation
modes. For NIRSpec observations there is one applicable argument:

``--slit_name``
  name [string value] of a specific slit region to extract. The default value of None
  will cause all known slits for the instrument mode to be extracted.

For NIRCam and NIRISS WFSS, the ``extract_2d`` step has three optional arguments:

``--grism_objects``
  list (default is empty). A list of ``jwst.transforms.models.GrismObject``.

``--mmag_extract``
  float (default is 99.) The minimum magnitude object to extract.

``--extract_orders``
  list. The list of orders to extract. The default is taken from the
  ``wavelengthrange`` reference file.


For NIRCam TSGRISM, the ``extract_2d`` step has two optional arguments:

``--extract_orders``
  list. The list of orders to extract. The default is taken from the ``wavelengthrange``
  reference file.

``--extract_height``
  int. The cross-dispersion size (in units of pixels) to extract.

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
``input_model.meta.source_catalog.filename``. Otherwise, a user can submit a list of
``GrismObjects`` that contains information about the objects that should be extracted.
The ``GrismObject`` list can be created automatically by using the method in
``jwst.assign_wcs.utils.create_grism_bbox``. This method also uses the name of the source
catalog saved in the input model's meta information. If it's better to construct a list
of ``GrismObjects`` outside of these, the ``GrismObject`` itself can be imported from
``jwst.transforms.models``.

The dispersion direction will be documented by copying keyword "DISPAXIS"
(1 = horizontal, 2 = vertical) from the input file to the output cutout.

Following are examples of how to customize the list of grism objects used in extract_2d.
The input file must have a WCS assigned to it by running ``assign_wcs``. The default values
for  wavelength range and the spectral orders are stored in the ``wavelengthrange``
reference file which can be retrieved from CRDS. A user may supply a different
wavelength range by passing `None` to ``reference_files``. In this case the spectral
orders to be extracted and their corresponding wavelength range will be taken
from the ``wavelength_range`` parameter which is a dictionary ``{spectral_order: (lam_min, lam_max)}``.

.. doctest-skip::

  >>> from jwst.datamodels import ImageModel
  >>> input_model = ImageModel("nircam_wfss_assign_wcs.fits")

Retrieve the wavelengthrange file specific for this mode:

.. doctest-skip::

  >>> from jwst.extract_2d import Extract2dStep
  >>> step = Extract2dStep()
  >>> refs = {}
  >>> for ref_type in step.reference_file_types:
  ...     refs[ref_type] = step.get_reference_file(input_model, ref_type)
  >>> print(refs)
  {'wavelengthrange': '/crds/references/jwst/niriss/jwst_niriss_wavelengthrange_0002.asdf'}

Create a list of grism objects for a specified spectral order with a limited
minimum magnitude, and a specified half height of the extraction box in
cross-dispersion direction. The ``wfss_extract_half_height`` parameter applies only to
point sources.

.. doctest-skip::

  >>> from jwst.assign_wcs.util import create_grism_bbox
  >>> grism_objects = create_grism_bbox(im, refs, mmag_extract=17,
  ... extract_orders=[1], wfss_extract_half_height=10)
  >>> print(len(grism_objects))
  6
  >>> print(grism_objects[0])
  GrismObject(sid=12, order_bounding={1: ((246, 266), (1367, 1581))}, sky_centroid=<SkyCoord (ICRS): (ra, dec) in deg
  (85.19582803, -69.53656873)>, partial_order={1: False}, waverange={1: (1.29, 1.71)}, sky_bbox_ll=<SkyCoord (ICRS): (ra, dec) in deg
  (85.19917182, -69.53721616)>, sky_bbox_lr=<SkyCoord (ICRS): (ra, dec) in deg
  (85.19270524, -69.53718398)>, sky_bbox_ur=<SkyCoord (ICRS): (ra, dec) in deg
  (85.19276186, -69.53579839)>, sky_bbox_ul=<SkyCoord (ICRS): (ra, dec) in deg
  (85.19922801, -69.53583056)>, xcentroid=1574.0825945473498, ycentroid=254.2556654610221)

Create a list of grism objects for a specified spectral order and wavelength range.
Use the source ID, ``sid`` to modify the extraction limits for specific objects.
The computed extraction limits are in the ``order_bounding`` attribute ordered ``(y, x)``.

.. doctest-skip::

  >>> from jwst.assign_wcs.util import create_grism_bbox
  >>> grism_objects = create_grism_bbox(im, mmag_extract=17, wavelength_range={1: (3.01, 4.26)})
  >>> print([obj.sid for obj in grism_objects])
  [12, 26, 31, 37, 57]
  >>> print(grism_objects[-1])
  id: 57
  order_bounding {1: ((995, 1114), (-18, 407))}
  sky_centroid: <SkyCoord (ICRS): (ra, dec) in deg
      (85.23831544, -69.52207261)>
  sky_bbox_ll: <SkyCoord (ICRS): (ra, dec) in deg
      (85.24337262, -69.5231152)>
  sky_bbox_lr: <SkyCoord (ICRS): (ra, dec) in deg
      (85.2351383, -69.52307624)>
  sky_bbox_ur: <SkyCoord (ICRS): (ra, dec) in deg
      (85.23522188, -69.5209249)>
  sky_bbox_ul:<SkyCoord (ICRS): (ra, dec) in deg
      (85.24345537, -69.52096386)>
  xcentroid: 767.278551509201
  ycentroid: 1053.7806251513593
  partial_order: {1: True}
  waverange: {1: (3.01, 4.26)}
  >>> grism_object[-1].order_bounding[1] = ((1000, 1110), (0, 450))
  >>> print(grism_object[-1].order_bounding
  {1: ((1000, 1110), (0, 450))})

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


For NIRCam TSGRISM, the ``extract_2d`` step has one optional argument:

``--tsgrism_extract_height``
  int The cross-dispersion size (in units of pixels) to extract.


For NIRCam and NIRISS WFSS mode, the ``extract_2d`` step has one optional argument:

  ``--wfss_extract_half_height``
    int. The cross-dispersion half size, in pixels.

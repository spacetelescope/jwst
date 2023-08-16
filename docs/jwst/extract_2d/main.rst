Description
===========

:Class: `jwst.extract_2d.Extract2dStep`
:Alias: extract_2d

Overview
--------
The ``extract_2d`` step extracts 2D arrays from spectral images. The extractions
are performed within all of the SCI, ERR, and DQ arrays of the input image
model, as well as any variance arrays that may be present. It also computes an
array of wavelengths to attach to the extracted data. The extracted arrays
are stored as one or more ``slit`` objects in an output ``MultiSlitModel``
and saved as separate tuples of extensions in the output FITS file.

Assumptions
-----------
This step uses the ``bounding_box`` attribute of the WCS stored in the data model,
which is populated by the ``assign_wcs`` step. Hence the ``assign_wcs`` step
must be applied to the science exposure before running this step.

For NIRCam and NIRISS WFSS modes, no ``bounding_box`` is attached to the data
model by the ``assign_wcs`` step.
This is to keep the WCS flexible enough to be used with any
source catalog or similar list of objects that may be associated with the dispersed image.
Instead, there
is a helper method that is used to calculate the bounding boxes that contain
the spectral trace for each object. One box is made for each spectral order of
each object. In regular, automated processing, the ``extract_2d`` step uses the
source catalog specified in the input
model's meta information to create the list of objects and their corresponding
bounding box. This list is used to make the 2D cutouts from the dispersed image.

NIRCam TSGRISM exposures do not use a source catalog, so the step instead relies on the
assumption that the source of interest is located at the aperture reference point
and centers the extraction there.
More details are given below.

Algorithm
---------
This step is currently applied only to NIRSpec Fixed-Slit, NIRSpec MOS, NIRSpec TSO
(BrightObj), NIRCam and NIRISS WFSS, and NIRCam TSGRISM observations.

NIRSpec Fixed Slit and MOS
++++++++++++++++++++++++++

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

The output ``MultiSlitModel`` data model will have the meta data associated with each
slit region populated with the name and region characteristic for the slits,
corresponding to the FITS keywords "SLTNAME", "SLTSTRT1", "SLTSIZE1",
"SLTSTRT2", and "SLTSIZE2."  Keyword "DISPAXIS" (dispersion direction)
will be copied from the input file to each of the output cutout images.


NIRCam and NIRISS WFSS
++++++++++++++++++++++

During normal, automated processing of WFSS grism images, the 
step parameter ``grism_objects`` is left unspecified, in which case the ``extract_2d``
step uses the source catalog that is specified in the input model's meta information,
``input_model.meta.source_catalog.filename`` ("SCATFILE" keyword) to define the
list of objects to be extracted.
Otherwise, a user can submit a list of ``GrismObjects`` that contains information
about the objects that they want extracted.
The ``GrismObject`` list can be created automatically by using the method in
``jwst.assign_wcs.utils.create_grism_bbox``. This method also uses the name of the source
catalog saved in the input model's meta information. If it's better to construct a list
of ``GrismObjects`` outside of these, the ``GrismObject`` itself can be imported from
``jwst.transforms.models``.

The dispersion direction will be documented by copying keyword "DISPAXIS"
(1 = horizontal, 2 = vertical) from the input file to the output cutout.

The ``wfss_mmag_extract`` and ``wfss_nbright`` parameters both affect which objects
from a source catalog will be retained for extraction. The rejection or retention of
objects proceeds as follows:

1. As each object is read from the source catalog, they are immediately rejected if 
   their isophotal_abmag > ``wfss_mmag_extract``, meaning that only objects brighter than
   ``wfss_mmag_extract`` will be retained. The default ``wfss_mmag_extract`` value of
   ``None`` retains all objects.

2. If the computed footprint (bounding box) of the spectral trace of an object lies
   completely outside the field of view of the grism image, it is rejected.

3. The list of objects retained after the above two filtering steps have been applied is
   sorted based on ``isophotal_abmag`` (listed for each source in the source catalog) and
   only the brightest ``wfss_nbright`` objects are retained. The default value of
   ``wfss_nbright`` is currently 1000.

All remaining objects are then extracted from the grism image.

WFSS Examples
^^^^^^^^^^^^^
The extraction of sources from WFSS grism images is a multi-step process, as outlined above.
Here we show detailed examples of how to customize the list of WFSS grism objects to be
extracted, in order to better explain the various steps.
First, the input file (or data model) must aleady have a WCS object assigned to it by running
the :ref:`assign_wcs <assign_wcs_step>` step. The default values
for the wavelength range of each spectral order to be extracted are also required;
they are stored in the ``wavelengthrange`` reference file, which can be retrieved from CRDS.

Load the grism image, which is assumed to have already been processed through ``assign_wcs``,
into an `ImageModel` data model (used for all 2-D "images", regardless of whether
they actually contain imaging data or dispersed spectra):

.. doctest-skip::

  >>> from stdatamodels.jwst.datamodels import ImageModel
  >>> input_model = ImageModel("jw12345001001_03101_00001_nis_assign_wcs.fits")

Load the ``extract_2d`` step and retrieve the ``wavelengthrange`` reference file
specific for this mode:

.. doctest-skip::

  >>> from jwst.extract_2d import Extract2dStep
  >>> step = Extract2dStep()
  >>> refs = {}
  >>> reftype = 'wavelengthrange'
  >>> refs[reftype] = step.get_reference_file(input_model, reftype)
  >>> print(refs)
  {'wavelengthrange': '/crds/jwst/references/jwst_niriss_wavelengthrange_0002.asdf'}

Create a list of grism objects for a specified spectral order with a limited
minimum magnitude and a specified half-height of the extraction box in the
cross-dispersion direction via the ``wfss_extract_half_height`` parameter.
Note that the half-height parameter only applies to point sources.

.. doctest-skip::

  >>> from jwst.assign_wcs.util import create_grism_bbox
  >>> grism_objects = create_grism_bbox(input_model, refs, mmag_extract=17,
  ... extract_orders=[1], wfss_extract_half_height=10)
  >>> print(len(grism_objects))
  6
  >>> print(grism_objects[0])
  id: 432
  order_bounding {1: ((array(1113), array(1471)), (array(1389), array(1609)))}
  sky_centroid: <SkyCoord (ICRS): (ra, dec) in deg
      (3.59204081, -30.40553435)>
  sky_bbox_ll: <SkyCoord (ICRS): (ra, dec) in deg
      (3.59375611, -30.40286617)>
  sky_bbox_lr: <SkyCoord (ICRS): (ra, dec) in deg
      (3.59520565, -30.40665425)>
  sky_bbox_ur: <SkyCoord (ICRS): (ra, dec) in deg
      (3.58950974, -30.4082754)>
  sky_bbox_ul:<SkyCoord (ICRS): (ra, dec) in deg
      (3.5880604, -30.40448726)>
  xcentroid: 1503.6541213285695
  ycentroid: 1298.2882813663837
  partial_order: {1: False}
  waverange: {1: (0.97, 1.32)}
  is_extended: True
  isophotal_abmag: 16.185488680084294

Create a list of grism objects for a specified spectral order and wavelength range.
In this case we don't use the default wavelength range limits from the ``wavelengthrange``
reference file, but instead designate custom limits via the ``wavelength_range`` parameter
passed to the ``create_grism_bbox`` function, which is a dictionary of the form
``{spectral_order: (wave_min, wave_max)}``. 
Use the source ID, ``sid``, to identify the object(s) to be modified.
The computed extraction limits are stored in the ``order_bounding`` attribute,
which is ordered ``(y, x)``.

.. doctest-skip::

  >>> from jwst.assign_wcs.util import create_grism_bbox
  >>> grism_objects = create_grism_bbox(input_model, mmag_extract=18,
  ... wavelength_range={1: (3.01, 4.26)})
  >>> print([obj.sid for obj in grism_objects])
  [12, 26, 31, 37, 104]
  >>> print(grism_objects[-1])
  id: 104
  order_bounding {1: ((array(1165), array(1566)), (array(1458), array(1577)))}
  sky_centroid: <SkyCoord (ICRS): (ra, dec) in deg
      (3.57958792, -30.40926139)>
  sky_bbox_ll: <SkyCoord (ICRS): (ra, dec) in deg
      (3.58060118, -30.40800999)>
  sky_bbox_lr: <SkyCoord (ICRS): (ra, dec) in deg
      (3.58136873, -30.41001654)>
  sky_bbox_ur: <SkyCoord (ICRS): (ra, dec) in deg
      (3.57866098, -30.4107869)>
  sky_bbox_ul:<SkyCoord (ICRS): (ra, dec) in deg
      (3.57789348, -30.40878033)>
  xcentroid: 1513.4964315117466
  ycentroid: 1920.6251490007467
  partial_order: {1: False}
  waverange: {1: (3.01, 4.26)}
  is_extended: True
  isophotal_abmag: 17.88278103874113
  >>> grism_object[-1].order_bounding[1] = ((1200, 1500), (1480, 1520))
  >>> print(grism_object[-1].order_bounding
  {1: ((1200, 1500), (1480,1520))}

The ``grism_objects`` list created in the above examples can now be used
as input to the ``extract_2d`` step in order to extract the particular objects
defined in that list:

.. doctest-skip::

  >>> result = step.call(input_model, grism_objects=grism_objects)

``result`` is a ``MultiSlitModel`` data model, containing one ``SlitModel``
instance for each extracted object, which includes meta data that identify
each object and the actual extracted data arrays, e.g.:

.. doctest-skip::

  >>> print(len(result.slits))
  8
  >>> result.slits[0].source_id
  104
  >>> result.slits[0].data
  array([[..., ..., ...]])


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

There are several arguments available for Wide-Field Slitless Spectroscopy (WFSS) and
Time-Series (TSO) grism spectroscopy:

``--tsgrism_extract_height``
  int. The cross-dispersion extraction size, in units of pixels. Only applies to TSO
  mode.

``--wfss_extract_half_height``
  int. The cross-dispersion half size of the extraction region, in pixels, applied to
  point sources. Only applies to WFSS mode.

``--wfss_mmag_extract``
  float (default is ``None``). The minimum (faintest) magnitude object to extract, based on
  the value of `isophotal_abmag` in the source catalog. Only applies to WFSS mode.

``--wfss_nbright``
  int (default is 1000). The number of brightest source catalog objects to extract.
  Can be used in conjunction with ``wfss_mmag_extract``. Only applies to WFSS mode.

``--extract_orders``
  list. The list of spectral orders to extract. The default is taken from the
  ``wavelengthrange`` reference file. Applies to both WFSS and TSO modes.

``--grism_objects``
  list (default is empty). A list of ``jwst.transforms.models.GrismObject``.

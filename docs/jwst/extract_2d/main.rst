
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
the list of objects and their corresponding bounding box. This list is used to make
the 2D cutouts from the dispersed image.

Algorithm
---------
The step is currently applied only to NIRSpec Fixed Slit, NIRSPEC MSA, NIRSPEC TSO,
NIRCAM and NIRISS WFSS, and NIRCAM TSGRISM observations.

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


For NIRCAM TSGRISM:

There is no source catalog created for TSO observations because the source is always
placed on the same pixel, the user can only vary the size of the subarray. All of the
subarrays have their "bottom" edge located at the physical bottom edge of the detector
and grow in size vertically. The source spectrum trace will always be centered
somewhere near row 34 in the vertical direction (dispersion running parallel to rows).
So the larger subarrays will just result in larger amount of sky above the spectrum.

`extract_2d` will always produce a cutout that is 64 pixels in height
(cross-dispersion direction) for all subarrays and full frame exposures,
which is equal to the height of the smallest available subarray (2048 x 64).
This is to allow area within the cutout for sampling the background in later steps,
such as `extract_1d`. The slit height is a parameter that a user can set
(during reprocessing) to tailor their results. 


Step Arguments
==============
The `extract_2d` step has two optional arguments for NIRSPEC observations:

* ``--slit_name``: name (string value) of a specific slit region to
  extract. The default value of None will cause all known slits for the
  instrument mode to be extracted. Currently only used for NIRspec fixed slit
  exposures.

* ``--apply_wavecorr``: bool (default is True). Flag indicating whether to apply the Nirspec wavelength zero-point correction.


For NIRCAM and NIRISS WFSS, the `extract_2d` step has three optional arguments:

* ``--grism_objects``: list (default is empty). A list of ``jwst.transforms.models.GrismObject``.

* ``--mmag_extract``: float. (default is 99.) the minimum magnitude object to extract

* ``--extract_orders``: list. The list of orders to extract. The default is taken from the `wavelengthrange` reference file.


For NIRCAM TSGRISM, the `extract_2d` step has two optional arguments:

* ``--extract_orders``: list. The list of orders to extract. The default is taken from the `wavelengthrange` reference file.

* ``--extract_height``: int. The cross-dispersion size to extract


Reference Files
===============
To apply the Nirspec wavelength zero-point correction, this step uses the
WAVECORR reference file. The zero-point correction is applied to observations
with EXP_TYPE of "NRS_FIXEDSLT", "NRS_BRIGHTOBJ" or "NRS_MSASPEC". This is an optional
correction (on by default). It can be turned off by specifying ``apply_wavecorr=False``
when running the step.

NIRCAM WFSS and NIRISS WFSS observations use the `wavelengthrange` reference file in order
to construct the bounding boxes around each objects orders. If a list of ``GrismObject``s
is supplied, then no reference file is neccessary.

For NIRCAM WFSS and TSGRIM modes and NIRISS WFSS mode the `wavelengthrange` file contains the wavelength limits to use when caluclating the minimum and maximum dispersion extents on the detector. It also contains the default list of orders that should be extracted for each filter.
To be consistent with other modes, and for conveniance,it also lists the orders and filters that are valid  with the file.

:order: A list of orders this file covers
:wavelengthrange: A list containing the list of [order, filter, wavelength min, wavelength max] 
:waverange_selector: The list of FILTER names available
:extract_orders: A list containing the list of orders to extract for each filter


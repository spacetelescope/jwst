.. _extract1d_reffile:

EXTRACT1D Reference File
------------------------
The EXTRACT1D reference file contains information needed to guide
the 1D extraction process. It is a text file, with the information
in JSON format.

.. include:: ../references_general/extract1d_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for EXTRACT1D
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in EXTRACT1D reference files,
because they are used as CRDS selectors
(see :ref:`extract1d_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
EXTRACT1D reference files are text files, with the information stored in
JSON format.
All the information is specified in a list with key ``apertures``.  Each
element of this list is a dictionary, one for each aperture (e.g. a slit)
that is supported by the given reference file.  The particular dictionary
to use is found by matching the slit name in the science data with the
value of key ``id``.  Key ``spectral_order`` is optional, but if it is
present, it must match the expected spectral order number.

The following keys are supported for non-IFU data (see below for IFU keys).
Key ``id`` is the primary criterion for selecting which element of
the ``apertures`` list to use.  The slit name (except for a full-frame
input image) is compared with the values of ``id`` in the ``apertures``
list to select the appropriate aperture.
In order to allow the possibility of multiple
spectral orders for the same slit name, there may be more than one element
of ``apertures`` with the same value for key ``id``.  These should then be
distinguished by using the secondary selection criterion ``spectral_order``.
In this case, the various spectral orders would likely have different
extraction locations within the image, so different elements of ``apertures``
are needed in order to specify those locations.
If key ``dispaxis`` is specified, its value will be used to set the
dispersion direction within the image.  If ``dispaxis`` is
not specified, the dispersion direction will be taken to be the axis
along which the wavelengths change more rapidly.
Key ``region_type`` can be omitted, but if it is specified, its value must
be "target".  The remaining keys specify the characteristics of the source
and background extraction regions.

* id: the slit name, e.g. "S200A1" (string)
* spectral_order: the spectral order number (optional); this can be either
  positive or negative, but it should not be zero (int)
* dispaxis: dispersion direction, 1 for X, 2 for Y (int)
* xstart: first pixel in the horizontal direction, X (int)
* xstop: last pixel in the horizontal direction, X (int)
* ystart: first pixel in the vertical direction, Y (int)
* ystop: last pixel in the vertical direction, Y (int)
* src_coeff: this takes priority for specifying the source extraction region
  (list of lists of float)
* bkg_coeff: for specifying background subraction regions
  (list of lists of float)
* independent_var: "wavelength" or "pixel" (string)
* smoothing_length: width of boxcar for smoothing background regions along
  the dispersion direction (odd int)
* bkg_fit: the type of background fit or computation (string)
* bkg_order: order of polynomial fit to background regions (int)
* extract_width: number of pixels in cross-dispersion direction (int)

See :ref:`extract-1d-for-slits` for more details on how these parameters
are used in the extraction process.

For IFU cube data, the following keys are used instead of those above.

* id: the slit name, but this can be "ANY" (string)
* x_center: X pixel coordinate of the target (pixels, float, the default
  is the center of the image along the X axis)
* y_center: Y pixel coordinate of the target (pixels, float, the default
  is the center of the image along the Y axis)
* radius: the radius of the circular extraction aperture (pixels, float,
  default is one quarter of the smaller of the image axis lengths)
* subtract_background: if true, subtract a background determined from an
  annulus with inner and outer radii given by ``inner_bkg`` and
  ``outer_bkg`` (boolean)
* inner_bkg: of the inner edge of the background annulus (pixels, float,
  default = ``radius``)
* outer_bkg: of the outer edge of the background annulus (pixels, float,
  default = ``inner_bkg * sqrt(2)``)
* method: one of "exact", "subpixel", or "center", the method
  used by photutils for computing the overlap between apertures and pixels
  (string, default is "exact")
* subpixels: if ``method`` is "subpixel", pixels will be resampled by this
  factor in each dimension (int, the default is 5)

See :ref:`extract-1d-for-ifu` for more details on how these parameters
are used in the extraction process.

Example EXTRACT1D Reference File
--------------------------------
The following JSON was taken as an example from reference file
jwst_niriss_extract1d_0003.json::

  {
      "REFTYPE": "EXTRACT1D",
      "INSTRUME": "NIRISS",
      "TELESCOP": "JWST",
      "DETECTOR": "NIS",
      "EXP_TYPE": "NIS_SOSS",
      "PEDIGREE": "GROUND",
      "DESCRIP": "NIRISS SOSS extraction params for ground testing",
      "AUTHOR": "M.Wolfe, H.Bushouse",
      "HISTORY": "Build 7.1 of the JWST Calibraton pipeline. The regions are rectagular and do not follow the trace.",
      "USEAFTER": "2015-11-01T00:00:00",
      "apertures": [
        {
        "id": "FULL",
        "region_type": "target",
        "bkg_coeff": [[2014.5],[2043.5]],
        "xstart":    4,
        "xstop":  2044,
        "ystart": 1792,
        "ystop":  1972,
        "dispaxis": 1,
        "extract_width": 181
        },

        {
        "id": "SUBSTRIP256",
        "region_type": "target",
        "bkg_coeff": [[221.5],[251.5]],
        "xstart":   4,
        "xstop": 2044,
        "ystart": 20,
        "ystop":  220,
        "dispaxis": 1,
        "extract_width": 201
        },

        {
        "id": "SUBSTRIP96",
        "region_type": "target",
        "bkg_coeff": [[1.5],[8.5],[92.5],[94.5]],
        "xstart":   4,
        "xstop": 2044,
        "ystart":  10,
        "ystop":   92,
        "dispaxis": 1,
        "extract_width": 83
        }]
  }


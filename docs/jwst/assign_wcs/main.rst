
Description
===========

:Class: `jwst.assign_wcs.AssignWcsStep`
:Alias: assign_wcs


The ``assign_wcs`` step is run at the beginning of the Stage 2 pipelines - :ref:`calwebb_image2 <calwebb_image2>`
and :ref:`calwebb_spec2 <calwebb_spec2>` - for both imaging and spectroscopic exposures.
It associates a WCS object with each science exposure. The WCS object transforms
positions in the detector frame to positions in a world coordinate frame - ICRS and wavelength.
In general there may be intermediate coordinate frames depending on the instrument.
The WCS is saved in the ASDF extension of the FITS file. It can be accessed as an attribute of
the meta object when the FITS file is opened as a data model (``model.meta.wcs``).

The forward direction of the transforms is from detector to world coordinates
and the input pixel coordinates are 0-indexed.

The ``assign_wcs`` step expects to find the basic WCS keywords in the
"SCI" extension header of the input FITS file. Distortion and spectral models are stored in reference files in the
`ASDF <http://asdf-standard.readthedocs.org/en/latest/>`__  format.

For each observing mode, determined by the value of ``EXP_TYPE`` in the science header,
``assign_wcs`` retrieves reference files from CRDS and creates a pipeline of transforms from
input frame ``detector`` to a frame ``v2v3``. This part of the WCS pipeline may include
intermediate coordinate frames. The basic WCS keywords are used to create
the transform from frame ``v2v3`` to frame ``world``. All of this information is used to
create and populate the WCS object for the exposure.

For display of imaging data with software like DS9 that relies on specific WCS information,
a SIP-based approximation to the WCS is fit. The results are FITS keywords stored in
``model.meta.wcsinfo``. This is not an exact fit, but is accurate to ~0.01 pixel by default,
and is sufficient for display purposes. This step is performed by default for imaging modes,
but can be switched off, and parameters controlling the SIP fit can
also be adjusted.  Note that if these parameters are changed, the equivalent parameters
for the ``tweakreg`` step should be adjusted to match.

The ``assign_wcs`` step can accept as input either a ``rate`` product, which is the result of
averaging over all integrations in an exposure, or a ``rateints`` product, which is a 3D cube of
per-integration images.

The ``assign_wcs`` WCS implementation is based on `gwcs <https://gwcs.readthedocs.io/en/latest/>`__ and
uses `asdf <http://asdf.readthedocs.io/en/latest/>`__ to define and store reference files and transforms.

.. Note:: In addition to CRDS reference files, applying ``assign_wcs`` to NIRSpec MOS
   exposures depends critically on an MSA metadata file to provide information
   for MOS slitlets in use and their constituent shutters. See :ref:`msa_metadata <msa_metadata>`
   for detailed information about the MSA metadata files and their contents.

Basic WCS keywords and the transform from ``v2v3`` to ``world``
---------------------------------------------------------------

All JWST instruments use the following FITS header keywords to
define the transform from ``v2v3`` to ``world``:

``RA_REF``, ``DEC_REF`` - a fiducial point on the sky, ICRS [deg]

``V2_REF``, ``V3_REF`` - a point in the V2V3 system that maps to ``RA_REF``, ``DEC_REF`` [arcsec]

``ROLL_REF`` - local roll angle associated with each aperture [deg]

``RADESYS`` - standard coordinate system [ICRS]

These quantities are used to create a 3D Euler angle rotation between the V2V3 spherical system,
associated with the telescope, and a standard celestial system.

For spectroscopic data, ``assign_wcs`` populates the keyword ``DISPAXIS``
with an integer value that indicates whether the dispersion direction is
oriented more nearly along the horizontal (DISPAXIS = 1) or vertical
(DISPAXIS = 2) direction in the image frame.


Using the WCS interactively
---------------------------

Once a FITS file is opened as a `DataModel` the WCS can be accessed as an attribute
of the meta object. Calling it as a function with detector positions as inputs returns the
corresponding world coordinates. Using MIRI LRS fixed slit as an example:

.. doctest-skip::

  >>> from stdatamodels.jwst.datamodels import ImageModel
  >>> exp = ImageModel('miri_fixedslit_assign_wcs.fits')
  >>> ra, dec, lam = exp.meta.wcs(x, y)
  >>> print(ra, dec, lam)
  (329.97260532549336, 372.0242999250267, 5.4176100046836675)

The WFSS modes for NIRCam and NIRISS have a slightly different calling structure.
In addition to the (x, y) coordinates, they need to know other information about the
spectrum or source object. In the JWST backward direction (going from the sky to
the detector) the WCS model also looks for the wavelength and order and returns
the (x,y) location of that wavelength+order on the dispersed image and the original
source pixel location, as entered, along with the order that was specified:

.. doctest-skip::

  >>> from stdatamodels.jwst.datamodels import ImageModel
  >>> exp = ImageModel('nircam_wfss_assign_wcs.fits')
  >>> x, y, x0, y0, order = exp.meta.wcs(x0, y0, wavelength, order)
  >>> print(x0, y0, wavelength, order)
  (365.523884327, 11.6539963919, 2.557881113, 2)
  >>> print(x, y, x0, y0, order)
  (1539.5898464615102, 11.6539963919, 365.523884327, 11.6539963919, 2)

Similarly, for all slit-like NIRSpec spectroscopic modes (MOS, FS, BOTS),
the assigned WCS needs to know the slit ID in order to return valid coordinates.
For example, to retrieve world coordinates for a pixel in slit 12 of a MOS observation
at pixel (x, y) = (804, 522):

.. doctest-skip::

  >>> from stdatamodels.jwst.datamodels import ImageModel
  >>> exp = ImageModel('nirspec_mos_assign_wcs.fits')
  >>> ra, dec, lam, slit_id = exp.meta.wcs(804, 522, 12)
  >>> print(ra, dec, lam, slit_id)
  (321.15970971929175, -16.549348214686127, 3.235814824179365, 12.0)

For MOS observations, the slit ID is the name of the slit, as specified by the
"slitlet_id" field in the :ref:`MSA metadata file<msa_metadata>`.
For fixed slit observations, the slit ID is a fixed integer for each slit as shown in the table below.
These values can also be retrieved by slit name using the :func:`jwst.assign_wcs.nrs_fs_slit_id` function.

.. list-table:: NIRSpec Fixed Slit IDs
   :header-rows: 1

   * - Slit Name
     - Slit ID
   * - S200A1
     - -101
   * - S200A2
     - -102
   * - S400A1
     - -103
   * - S1600A1
     - -104
   * - S200B1
     - -105

For these modes, when they are processed through the :ref:`extract_2d <extract_2d_step>`
step, a new WCS is assigned to each extracted slit that fixes the slit
ID to a specific value, so it is no longer required on input and not reported on output.
For example, for a NIRSpec fixed slits exposure, which has been processed through the
extract_2d step, coordinates for the first slit image at (x, y) = (56, 15) can be calculated
with:

.. doctest-skip::

  >>> exp = datamodels.MultiSlitModel('nrs1_fixed_assign_wcs_extract_2d.fits')
  >>> ra, dec, lam = exp.slits[0].meta.wcs(56, 15)
  >>> print(ra, dec, lam)
  (46.25382856669748 46.279084130418504 0.9024513743123106)

The WCS also provides access to intermediate coordinate frames
and transforms between any two frames in the WCS pipeline in the forward or
backward directions. For this same fixed slit exposure:

.. doctest-skip::

  >>> exp.slits[0].meta.wcs.available_frames
  ['detector', 'sca', 'gwa', 'slit_frame', 'msa_frame', 'oteip', 'v2v3', 'v2v3vacorr', 'world']
  >>> detector2msa = exp.slits[0].meta.wcs.get_transform('detector', 'msa_frame')
  >>> detector2msa(56, 15)
  (0.02697267383337021, -0.0025054994862709653, 9.024513743123106e-07)
  >>> msa2detector = exp.slits[0].meta.wcs.get_transform('msa_frame', 'detector')
  >>> msa2detector(0.027, -0.0025, 9.02e-07)
  (55.28220392304502, 14.868098132739078)

WCS of NIRSpec IFU exposures
----------------------------

For NIRSpec IFU data, the assigned WCS may be either a slice-based WCS which uses slice
numbers (zero-indexed) to identify the "slits" containing valid coordinates, or it may be a
coordinate-based WCS that does not require the slice ID to perform coordinate transformations
for any pixel in the IFU image. The coordinate-based WCS is the default choice; if the slice-based
WCS is needed for diagnostic purposes, it can be produced by specifying ``nrs_ifu_slice_wcs = True``
in the parameters for the ``assign_wcs`` step.

With the coordinate-based WCS, it is possible to
retrieve coordinates for the entire image at once, for example:

.. doctest-skip::

  >>> import numpy as np
  >>> from stdatamodels.jwst.datamodels import IFUImageModel
  >>> exp = datamodels.IFUImageModel('nirspec_ifu_assign_wcs.fits')
  >>> y, x = np.mgrid[:exp.data.shape[0], :exp.data.shape[1]]
  >>> ra, dec, lam = exp.meta.wcs(x, y)

Coordinates for pixels outside of slice regions will be NaN.

For reference, the slice ID assigned to each pixel is available in an image
attached to the model at the end of the ``assign_wcs`` step, in the ``regions``
attribute (saved in the FITS extension REGIONS).  This image can be used to calculate
coordinates for a particular IFU slice, as needed.  For example, using the same model
and (x, y) grid as above, to get the coordinates for slice 1:

.. doctest-skip::

  >>> slice_idx = 1
  >>> in_slice = (exp.regions == slice_idx)
  >>> slice_ra, slice_dec, slice_lam = exp.meta.wcs(x[in_slice], y[in_slice])

Slice IDs in the regions image are one-indexed, with values between 1 and 30.
Zero indicates a pixel outside any valid slice region.

WCS of slitless grism exposures
-------------------------------

The WCS forward transforms for slitless grism exposures (``NIS_WFSS``, ``NRC_WFSS``, ``NRC_TSGRISM``)
take as input the ``x, y`` coordinates on the dispersed image, the ``x0, y0`` coordinate of
the center of the object in the direct image and ``spectral order``. They return the ``x0, y0`` coordinate of the center
of the object in the direct image, ``wavelength`` and ``spectral order``.

For NIRISS WFSS data the reference files contain a reference value for the filter wheel
position angle. The trace is rotated about an angle which is the difference between
the reference and actual angles.

For WFSS modes (``NIS_WFSS``, ``NRC_WFSS``), an approximation of the GWCS object
associated with a direct image with the same instrument configuration as the grism image
is saved as FITS WCS in the headers of grism images.

Corrections Due to Spacecraft Motion
------------------------------------

The WCS transforms contain two corrections due to motion of the observatory.

Absolute velocity aberration is calculated onboard when acquiring the guide star, but
differential velocity aberration effects are calculated during the ``assign_wcs`` step.
This introduces corrections in the conversion from sky coordinates to observatory
V2/V3 coordinates, and is stored in the WCS under the ``v2v3vacorr`` frame.

For spectroscopic data, a relativistic Doppler correction is applied to all wavelengths to place
observations into the barycentric reference frame. This correction factor is applied to the WCS
wavelength solution created during the ``assign_wcs`` step, such that extracted spectral products
will have wavelength arrays in the barycentric frame.

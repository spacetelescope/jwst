Description
===========

:Classes: `jwst.resample.resample_spec_step.ResampleSpecStep`
:Alias: resample_spec

This routine will resample each input 2D spectral image based on the WCS and
distortion information, and will combine multiple resampled images
into a single rectified product.  The distortion information should have
been incorporated into the image using the
:ref:`assign_wcs <assign_wcs_step>` step.

The ``resample`` step can take as input either:

#. a single 2D input image, or
#. an association table (in JSON format)

The underlying resampling algorithm is the same as is used in the
:ref:`resample step <resample_step>`: 2D spectral images are resampled into
a rectified spatial/spectral output WCS via the drizzle algorithm.  See
the documentation for that step for more information on error propagation
and output extensions, as well as further references for the drizzle
algorithm.

The parameters for the drizzle operation are provided via
:ref:`resample_spec step parameters <resample_spec_step_args>`.

Output WCS
----------

The output WCS for the resampled spectral image is created from a combination of
the WCS information of all input images.  The first input image is taken as the
reference for the output WCS and expanded to include the spatial and spectral
range of all inputs.

The full spatial range is linearly sampled along the slit at the specified output
pixel scale.  The wavelength range is not rebinned to a linear sampling: the wavelengths
for each output dispersion element are assigned from the median value across spatial
elements in the input images.

GWCS
^^^^

The full celestial and spectral WCS information is stored in the GWCS object
in the metadata for each output slit image (``meta.wcs``).  The WCS is derived from
that assigned in the :ref:`assign_wcs <assign_wcs_step>` step, but drops several
intermediate frames that are no longer used after resampling.  For NIRSpec fixed slit and MOS,
the output frames available are:

* ``detector``: image x, y
* ``slit_frame``: slit x, slit y, wavelength
* ``world``: RA, Dec, wavelength

For MIRI LRS, the output frames available are:

* ``detector``: image x, y
* ``world``: RA, Dec, wavelength

FITS WCS
^^^^^^^^

For user convenience and FITS display purposes, a simplified spatial-spectral WCS is
also stored in the FITS metadata for the output products.

The spatial coordinates are set to a simple linear scaling by the output pixel scale
in arcsec. The scaling value is stored in the standard FITS WCS keyword (``CDELTi``). The
origin for the spatial scale (``CRPIXi``) is set to the planned source location if available;
the center of the cross-dispersion axis if not.  This encodes spatial scaling information
but not celestial position: for full spatial information, the GWCS object should be used.

The wavelengths for the output image are nonlinear so they cannot be represented by simple
FITS header keywords, but they do have the same value for all spatial elements.
The values are extracted from the ``world`` frame in the GWCS at the central spatial index and
stored in a new ``WCS-TABLE`` extension. Wavelengths for each slit in the output product are
stored in separate columns in the table, named for the slit (e.g., ``wave_slit_S200A1`` for
NIRSpec fixed slit S200A1), or just ``wavelength`` if no slit name is available.
The ``WCS-TABLE`` extension follows the FITS standard described in
*Representations of spectral coordinates in FITS*, Greisen, et al., 2006, **A & A**, 446, 747-771.

Reference Files
---------------
The ``resample_spec`` step does not use any reference files.

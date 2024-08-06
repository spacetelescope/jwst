Description
===========

:Classes: `jwst.resample.ResampleSpecStep`
:Alias: resample_spec

This routine will resample each input 2D spectral image based on the WCS and
distortion information, and will combine multiple resampled images
into a single rectified product.  The distortion information should have
been incorporated into the image using the
:ref:`assign_wcs <assign_wcs_step>` step.

The ``resample`` step can take as input either:

#. a single 2D input image
#. an association table (in json format)

The underlying resampling algorithm is the same as is used in the
:ref:`resample step<resample_step>`: 2D spectral images are resampled into
a rectified spatial/spectral output WCS via the drizzle algorithm.  See
the documentation for that step for more information on error propagation
and output extensions, as well as further references for the drizzle
algorithm.

The output WCS is created from a combination of the WCS information of
all input images.  The first input image is taken as the reference
for the output WCS and expanded to include the spatial and spectral
range of all inputs.

The parameters for the drizzle operation are provided via ``resample_spec``
step parameters, which may be overridden by a step parameter reference
file from CRDS.

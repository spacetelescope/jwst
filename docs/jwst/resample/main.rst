Description
===========

This routine will resample each input 2D image based on the WCS and
distortion information, and will combine multiple resampled images
into a single undistorted product.  The distortion information should have
been incorporated into the image using the
:ref:`assign_wcs <assign_wcs_step>` step.

The ``resample`` step can take as input either:

  * a single 2D input image
  * an association table (in json format)

The defined parameters for the drizzle operation itself get
provided by the DRIZPARS reference file (from CRDS).  The exact values
used depends on the number of input images being combined and the filter
being used. Other information may be added as selection criteria later,
but for now, only basic information is used.

The output product gets defined using the WCS information of all inputs,
even if it is just a single input image. The output WCS defines a
field-of-view that encompasses the undistorted footprints on the sky
of all the input images with the same orientation and plate scale
as the first listed input image.

This step uses the interface to the C-based cdriz routine to do the
resampling via the drizzle method.  The input-to-output pixel
mapping is determined via a mapping function derived from the
WCS of each input image and the WCS of the defined output product.
This mapping function gets passed to cdriz to drive the actual
drizzling to create the output product.

Use the ``resample_spec`` step for spectroscopic data.  The dispersion
direction is needed for this case, and this is obtained from the
DISPAXIS keyword.

A full description of the drizzling algorithm, and parameters for
drizzling, can be found in the
`DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_.

Description
===========

:Class: `jwst.targ_centroid.TargCentroidStep`
:Alias: targ_centroid

The ``targ_centroid`` step determines the location of a point source in the detector
frame using the target acquisition (TA) verification image. The measured source
position is used by downstream calibration steps, such as
:ref:`pathloss <pathloss_step>`, to apply corrections that account
for source positioning within the aperture.

This step is applicable to MIRI Low-Resolution Spectrometer (LRS) observations
where a TA verification image is available. Pre-launch assumptions that telescope
pointing information would be sufficient for accurate pathloss corrections have
been superseded by in-flight experience, which demonstrates the need to use
actual TA verification images when available.

Upon successful completion of this step, the status keyword S_TACNTR is
set to "COMPLETE", the position of the source in the TA verification is stored in
the model's ``ta_xpos`` and ``ta_ypos`` attributes (FITS keywords TA_XPOS and TA_YPOS),
and the expected source position in the science exposure, applying any
necessary filter offset, is stored in the model's
``source_xpos`` and ``source_ypos`` attributes (FITS keywords SRCXPOS and SRCYPOS).

When run as part of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, the
step is executed after the :ref:`srctype <srctype_step>` step.

Step Inputs
-----------

The ``targ_centroid`` step can accept either:

* A single science exposure (ImageModel) along with a separate TA verification
  image file specified via the ``ta_file`` parameter, or
* An association file that includes both the science exposure and the TA
  verification image as separate members.

If no TA verification image is found in either location, the step will be skipped.
If the ``ta_file`` parameter is provided and an association with a TA verification
image is also provided, the file specified by ``ta_file`` will take precedence.

Algorithm
---------

The ``targ_centroid`` step performs the following operations:

#. **Load reference files**: The step retrieves the following reference files from CRDS:

   * FILTEROFFSET: Provides column and row offsets specific to the observation filter

#. **Assign WCS to TA verification image**: If the TA verification image does not already have
   a WCS assigned, the step invokes the :ref:`assign_wcs <assign_wcs_step>` step
   to compute and attach the WCS to the TA verification image data model.

#. **Determine reference position**: The step computes the expected source location
   based on the TA verification image WCS and metadata.

#. **Extract cutout**: A cutout of the TA verification image is extracted,
   centered on the expected source position.

#. **Find the centroid**: The centroid of the source within the cutout is determined
   using the :func:`~photutils.centroids.centroid_2dg` method,
   which fits the source with a 2D Gaussian.
   This yields the fitted x and y coordinates of the source in
   cutout pixel coordinates.

#. **Store TA position**: The measured source position in the TA verification
   image is stored in the input data model's ``ta_xpos`` and ``ta_ypos``
   attributes for reference.

#. **Transform coordinates**: The fitted position is transformed from TA verification image cutout
   coordinates to science image coordinates, accounting for any subarray
   offsets and/or dither offsets.

#. **Apply filter offsets**: The FILTEROFFSET reference file corrections are
   applied to the final x and y position.

#. **Store source position**: The corrected source position is stored in the output
   data model's ``source_xpos`` and ``source_ypos`` attributes for use by
   subsequent calibration steps.

Step Outputs
------------

The input science exposure is returned unmodified, except with new metadata attributes:

* ``ta_xpos``: The measured x-coordinate of the source
  in the TA verification image (0-indexed pixels)
* ``ta_ypos``: The measured y-coordinate of the source
  in the TA verification image (0-indexed pixels)
* ``source_xpos``: The x-coordinate of the source in the science detector
  coordinate system (0-indexed pixels).
* ``source_ypos``: The y-coordinate of the source in the science detector
  coordinate system (0-indexed pixels).
* ``meta.cal_step.targ_centroid`` keyword set to "COMPLETE"

Description
===========

:Class: `jwst.targ_centroid.TargCentroidStep`
:Alias: targ_centroid

The ``targ_centroid`` step determines the location of a point source in the detector
frame using the target acquisition (TA) verification image. The measured source
position is used by downstream calibration steps, such as
:ref:`pathloss <pathloss_step>`, to apply corrections that account
for source positioning within the aperture.

This step is applicable to MIRI LRS (Low Resolution Spectrometer) observations
where a TA verification image is available. Pre-launch assumptions that telescope
pointing information would be sufficient for accurate pathloss corrections have
been superseded by in-flight experience, which demonstrates the need to use
actual TA verification images when available.

Upon successful completion of this step, the status keyword S_TACNTR will be
set to "COMPLETE" and the source position will be stored in the output data
model's ``source_xpos`` and ``source_ypos`` attributes.

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

#. **Load reference files**: The step retrieves three reference files from CRDS:

   * SPECWCS: Contains the reference position (x_ref, y_ref) for the slit and
     slitless configurations
   * FILTEROFFSET: Provides column and row offsets specific to the observation filter
   * PATHLOSS: Used for slit observations to weight the model point source by the pathloss
     response during fitting

#. **Determine reference position**: Based on the SPECWCS reference file metadata,
   the step identifies the expected source location. If the TA verification image
   uses a subarray, appropriate coordinate transformations are applied.

#. **Extract cutout**: A cutout of the TA verification image is extracted,
   centered on the reference position.

#. **Build source model**: An initial source model is constructed:

   * The model is an Airy disk at approximately diffraction-limited resolution
   * The initial position is set to the flux-weighted centroid of the cutout region
   * Free parameters x position, y position, flux amplitude, and Airy
     disk radius are all allowed to vary during fitting

#. **Apply pathloss weighting** (slit mode only): For slit observations, the
   source model is multiplied by the pathloss correction at the central wavelength
   of the observation filter. The pathloss correction ensures the fit accounts for
   the aperture response, which is relevant at longer wavelengths where the PSF size
   approaches the slit width.

#. **Fit the model**: The model is fit to the data in the cutout using a
   least-squares optimizer to determine the best-fit x and y position.

#. **Transform coordinates**: The fitted position is transformed from cutout
   coordinates to full-frame detector coordinates, accounting for any subarray
   offsets.

#. **Apply filter offsets**: The FILTEROFFSET reference file corrections are
   applied to the final x and y position.

#. **Store results**: The corrected source position is stored in the output
   data model's ``source_xpos`` and ``source_ypos`` attributes for use by
   subsequent calibration steps.

Step Outputs
------------

The input science exposure is returned unmodified, except with three new attributes:

* ``source_xpos``: The x-coordinate of the source in the full-frame detector
  coordinate system (0-indexed pixels)
* ``source_ypos``: The y-coordinate of the source in the full-frame detector
  coordinate system (0-indexed pixels)
* ``meta.cal_step.targ_centroid`` keyword set to "COMPLETE"

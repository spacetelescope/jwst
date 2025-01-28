Description
===========

:Class: `jwst.clean_flicker_noise.CleanFlickerNoiseStep`
:Alias: clean_flicker_noise

Overview
--------
The ``clean_flicker_noise`` step removes flicker noise from calibrated
ramp images, after the :ref:`jump <jump_step>` step and prior to
performing the :ref:`ramp_fitting <ramp_fitting_step>` step.

For NIR detectors, the noise addressed by this step is 1/f noise, which
appears as faint banding along the detector slow axis.  For NIRSpec and
NIRISS, 1/f noise looks like vertical stripes; for NIRCam, it appears
as horizontal stripes.

MIRI images also often show a faint vertical banding, similar to 1/f noise
but from a different physical source.  This type of flicker noise can be
corrected by similar methods, so it is also addressed by this step.

To correct for flicker noise, the algorithm requires that the noise
generated between one group readout and the next be isolated as much
as possible from the astronomical signal. The residual noise may then
be fit in frequency space, via an FFT, or may be characterized by a
median along the rows or columns, as appropriate for the detector.
The modeled noise is then directly subtracted from each group readout.

The correction is available for any type of NIRSpec, NIRISS, or NIRCam
exposure. For MIRI, it is available only for imaging exposures.

Creation of a scene mask
------------------------
One of the key components of the correction is knowing which pixels can
be used to fit the background noise.  The step builds a scene mask
on the fly from a draft rate image, generated from the input ramp data,
which is used to mark usable and unusable pixels. The mask is a 2D
Boolean array, having the same size as the image, with
pixels set to True interpreted as being OK to use.

The process of building the mask varies somewhat depending on the
observing mode of the image being processed. Some features are common
to all modes, while others are mode-specific. The following sections
describe each type of masking that can be applied. At the end, there
is a summary of the types applied to each instrument mode.

The user-settable step parameter `save_mask` can be used to save the
scene mask to a file, if desired (see the
:ref:`step arguments <clean_flicker_noise_arguments>`).

Note that a user may also supply their own mask image as input to the step,
in which case the process of creating a mask is skipped. The step parameter
`user_mask` is used to specify an input mask.  If specified, the input
mask must contain a datamodel matching the shape of a 'rate' or 'rateints'
file generated for the input science data (either `~jwst.datamodels.ImageModel`
or `~jwst.datamodels.CubeModel`).  To generate a custom mask, it may be
easiest to save the mask output from a default run of the ``clean_flicker_noise``
step on the input data, then modify it as needed.

NIRSpec IFU Slices
^^^^^^^^^^^^^^^^^^
For IFU images, the majority of the mask is based on knowing which
pixels are contained within the footprints of the IFU slices. To do
this, the image's World Coordinate System (WCS) object is queried in
order to determine which pixels are contained within each of the 30
slices. Pixels within each slice are set to False (do not use) in the
mask.

NIRSpec MOS/FS Slits
^^^^^^^^^^^^^^^^^^^^
The footprints of each open MOS slitlet or fixed slit are flagged in
a similar way as IFU slices. For MOS and FS images, the WCS object is
queried to determine which pixels are contained within each open
slit/slitlet and they are set to False in the mask.

NIRSpec MSA Failed Open Shutters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pixels affected by stuck open MSA shutters are masked, because they
may contain signal. This is accomplished by setting all pixels flagged by the
:ref:`msaflagopen <msaflagopen_step>` step with DQ value "MSA_FAILED_OPEN"
to False in the mask.

NIRSpec Fixed-Slit Region Pixels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Full-frame MOS and IFU images may contain signal from the always open
fixed slits, which appear in a fixed region in the middle of each image.
The entire region containing the fixed slits is masked out when
processing MOS and IFU images. The masked region is currently hardwired
in the step to image indexes [1:2048, 923:1116], where the indexes are
in x, y order and in 1-indexed values.

Note, however, that it is possible to plan one or more fixed slit targets
alongside MSA slitlets in MOS observations. In this situation, the fixed
slit region is not automatically masked.

MIRI Imaging Metering Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The MIRI imager has a metering structure covering a large part of the
detector area. These regions must be masked in order to fit and
remove the background data in the science areas of the detector.
Regions marked with DQ value "DO_NOT_USE" by the
:ref:`flat_field <flatfield_step>` step are set to False in the
scene mask.

Missing Data
^^^^^^^^^^^^
Any pixel in the draft rate image that has a value of NaN or exactly zero
is flagged as False in the mask. This typically includes any reference
pixels that are present in the exposure.

Outliers
^^^^^^^^
For imaging modes, bright, compact sources must be distinguished
from the background and masked in order to fit and remove the
smooth background level.

For spectral modes that already have significant masking applied,
pixels in the unilluminated areas of the region can still contain anomalous
signal, due to uncaught cosmic rays, hot pixels, etc.

For both modes, a sigma-clipping routine is employed to find such outliers
within the draft rate image and set them to False in the mask. All pixels with
values greater than :math:`median+n\_sigma*sigma` are assumed to contain
signal and are set to False in the scene mask. In addition, all pixels
with values less than :math:`median-3.0*sigma` are assumed to be bad pixels,
and are also set to False in the scene mask.

Mode-Specific Masking Steps
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following table indicates which flavors of masking are applied to
images from each instrument and observing mode.

.. |c| unicode:: U+2713 .. checkmark

+--------------------------+-----+-----+-----+-------+--------+--------+
|                          |    NIRSpec      | MIRI  | NIRCam | NIRISS |
+--------------------------+-----+-----+-----+-------+--------+--------+
|                          | IFU | MOS |  FS | Image | All    | All    |
+==========================+=====+=====+=====+=======+========+========+
| IFU Slices\ :sup:`1`     | |c| |     |     |       |        |        |
+--------------------------+-----+-----+-----+-------+--------+--------+
| Slits/Slitlets\ :sup:`1` |     | |c| | |c| |       |        |        |
+--------------------------+-----+-----+-----+-------+--------+--------+
| MSA_FAILED_OPEN\ :sup:`1`| |c| | |c| |     |       |        |        |
+--------------------------+-----+-----+-----+-------+--------+--------+
| Non-science\ :sup:`1`    |     |     |     | |c|   |        |        |
+--------------------------+-----+-----+-----+-------+--------+--------+
| FS Region\ :sup:`1`      | |c| | |c| |     |       |        |        |
+--------------------------+-----+-----+-----+-------+--------+--------+
| Missing Data             | |c| | |c| | |c| | |c|   | |c|    | |c|    |
+--------------------------+-----+-----+-----+-------+--------+--------+
| Outliers                 | |c| | |c| | |c| | |c|   | |c|    | |c|    |
+--------------------------+-----+-----+-----+-------+--------+--------+

:sup:`1`\ These steps are only applied if the
:ref:`step parameter <clean_flicker_noise_arguments>`
`mask_science_regions` is set to True.

Correction Algorithm
--------------------

The detailed process for fitting and removing flicker noise is as follows.
See the :ref:`step arguments <clean_flicker_noise_arguments>` for more
information on all referenced parameters.

#. From the calibrated ramp input, make a draft rate (`single_mask` = True)
   or rateints (`single_mask` = False) file.

#. Create a scene mask from the rate data.

   #. If `mask_science_regions` is set and the input is NIRSpec data,
      run :ref:`assign_wcs <assign_wcs_step>` and
      :ref:`msaflagopen <msaflagopen_step>` on the draft rate data,
      then mask any known science areas or failed-open MSA shutters.

      This will mask out regions that are likely to contain significant
      astronomical signal.

   #. If `mask_science_regions` is set and the input is MIRI imaging data,
      run :ref:`flat_field <flatfield_step>` on the draft rate data,
      and extract just the DQ plane from the output. Pixels flagged
      as 'DO_NOT_USE' by the flat fielding process are masked.

      This will mask out regions of the detector under the metering
      structure.

   #. If `apply_flat_field` is set and a flat file is available, divide the
      draft rate data by the flat image.

   #. Iteratively sigma clip the data to get a center value (mean or median)
      and sigma value (standard deviation).

   #. If `fit_histogram` is set, compute a histogram from 4-sigma clipped
      values and fit a Gaussian to it to refine the center and sigma values.

   #. Mask data more than 3 * sigma below the center as bad values.

   #. Mask data more than `n_sigma` * sigma above the center as signal
      (not background).

#. Iterate over each integration and group in the data, to fit and correct
   for noise.

   #. Make a diff image (current group – previous group) to correct.

   #. If `apply_flat_field` is set and a flat file is available, divide the
      diff image by the flat image.

   #. Fit and remove a background level, using the scene mask to identify
      background pixels.

      #. Clip the background data in the diff image to remove more outliers.

      #. If `background_method` = 'median', the background value is a simple
         median of the remaining values.

      #. If `background_method` = 'model', the background data is fit with
         a low-resolution model via the photutils
         `Background2D <https://photutils.readthedocs.io/en/latest/api/photutils.background.Background2D.html>`_
         utility. The resolution box size is set by `background_box_size`.

      #. Subtract the background level from the diff image and clip again
         to `n_sigma` * sigma, with sigma recomputed from the
         background-subtracted data in the remaining background pixels.

   #. Fit and remove the residual noise in the background-subtracted image.

      #. If `fit_method` = 'fft', the ``nsclean`` library is called to fit
         and remove the noise in frequency space.

      #. If `fit_method` = 'median', the noise is fit with a simple median
         along the appropriate detector axis and subtracted from the
         background-subtracted image.

         If `fit_by_channel` = True, and the data is a NIR full-frame exposure,
         the median value is computed and subtracted independently for each
         detector channel.

   #. Restore the background level to the cleaned, background-subtracted
      diff image.  Also restore the flat structure if needed by multiplying the
      cleaned diff by the flat image.

   #. Add the cleaned diff back to a cleaned version of the previous
      group image.

References
==========

The FFT cleaning algorithm implementation is based on NSClean,
developed by Bernard Rauscher. Details on the source of the correlated
noise and the algorithm used by the ``nsclean`` library to fit and
remove it can be found in
`Rauscher 2024 <https://ui.adsabs.harvard.edu/abs/2023arXiv230603250R/abstract>`__.

The background fitting and median cleaning algorithm are based on
the image1overf algorithm, developed by Chris Willott, and available
on GitHub at `chriswillott/jwst <https://github.com/chriswillott/jwst>`__.
The algorithm was adapted to the `clean_flicker_noise` step and is released
under the BSD license for the JWST calibration pipeline by permission
of the author.

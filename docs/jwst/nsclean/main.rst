Description
===========

:Class: `jwst.nsclean.NSCleanStep`
:Alias: nsclean

Overview
========
The ``nsclean`` step applies an algorithm for removing correlated read
noise from NIRSpec images. The noise often appears as faint vertical
banding and so-called "picture frame noise." The algorithm uses dark
(unilluminated) areas of an image to fit a background model in Fourier
space. When the fit is subtracted, it removes nearly all correlated noise.
Compared to simpler strategies, like subtracting a rolling median, this
algorithm is more thorough and uniform. It is also computationally
undemanding, typically requiring only a few seconds to clean a full-frame
image.

The correction can be applied to any type of NIRSpec exposure, including
IFU, MOS, fixed slit, and Bright Object Time Series (BOTS), in both full-frame
and subarray readouts. Time series (3D) data are corrected one integration
at a time.

.. note::

   The step is currently not capable of processing images taken using the
   "ALLSLITS" subarray. Other subarray types are allowed.

Details on the source of the correlated noise and the algorithm used
in the ``nsclean`` step to fit and remove it can be found in
`Rauscher 2023 <https://ui.adsabs.harvard.edu/abs/2023arXiv230603250R/abstract>`_.

Upon completion of the step, the step status keyword "S_NSCLEN" gets set
to "COMPLETE" in the output science data.

Assumptions
===========
As described below, the creation of a pixel mask depends on the presence
of a World Coordinate System (WCS) object for the image, which is
constructed by the :ref:`assign_wcs <assign_wcs_step>` step.
In addition, creating a mask for IFU and MOS images depends on
the presence of DQ flags assigned by the
:ref:`msaflagopen <msaflagopen_step>` step.
It is therefore required that those steps be run before attempting to
apply ``nsclean``.

Creation of an image mask
=========================
One of the key components of the correction is knowing which pixels are
unilluminated and hence can be used in fitting the background noise.
The step builds a mask on the fly for each image, which is used to mark
useable and unuseable pixels. The mask is a 2D boolean array, having the same
size as the image, with pixels set to True interpreted as being OK to use.
The process of building the mask varies somewhat depending on the
observing mode of the image being processed. Some features are common
to all modes, while others are mode-specific. The following sections
describe each type of masking that can be applied and at the end there
is a summary of the types applied to each image mode.

The user-settable step parameter `save_mask` can be used to save the
mask to a file, if desired (see :ref:`nsclean step arguments <nsclean_arguments>`).

Note that a user may supply their own mask image as input to the step,
in which case the process of creating a mask is skipped. The step parameter
`user_mask` is used to specify an input mask.

IFU Slices
----------
For IFU images the majority of the mask is based on knowing which
pixels are contained within the footprints of the IFU slices. To do
this, the image's World Coordinate System (WCS) object is queried in
order to determine which pixels are contained within each of the 30
slices. Pixels within each slice are set to False (do not use) in the
mask.

MOS/FS Slits
------------
The footprints of each open MOS slitlet or fixed slit are flagged in
a similar way as IFU slices. For MOS and FS images, the WCS object is
queried to determine which pixels are contained within each open
slit/slitlet and they are set to False in the mask.

MSA Failed Open Shutters
------------------------
Pixels affected by stuck open MSA shutters are masked, because they
may contain signal. This is accomplished by setting all pixels flagged by the
:ref:`msaflagopen <msaflagopen_step>` step with DQ value "MSA_FAILED_OPEN"
to False in the mask.

NaN Pixels
----------
Any pixel in the input image that has a value of NaN is temporarily reset
to zero for input to the fitting routine and flagged as False in the mask.
Upon completion of the noise subtraction, this population of pixels is
set back to NaN again in the output (corrected) image.

Fixed-Slit Region Pixels
------------------------
Full-frame MOS and IFU images may contain signal from the always open
fixed slits, which appear in fixed region in the middle of each image.
The entire region containing the fixed slits is masked out when
processing MOS and IFU images. The masked region is currently hardwired
in the step to image indexes [1:2048, 923:1116], where the indexes are
in x, y order and in 1-indexed values.

Left/Right Reference Pixel Columns
----------------------------------
Full-frame images contain 4 columns of reference pixels on the left and
right edges of the image. These are not to be used in the fitting
algorithm and hence are set to False in the mask.

Outliers
--------
Pixels in the unilluminated regions of the region can contain anomalous
signal, due to uncaught Cosmic Rays, hot pixels, etc. A sigma-clipping
routine is employed to find such outliers within the input image and set
them to False in the mask. All pixels with values greater than
:math:`median+n_sigma*sigma` are set to False in the mask.
Here `median` and `sigma` are computed
from the image using the astropy.stats `sigma_clipped_stats` routine,
using the image mask to exclude pixels that have already been flagged
and a clipping level of 5 sigma. `n_sigma` is a user-settable step
parameter, with a default value of 5.0
(see :ref:`nsclean step arguments <nsclean_arguments>`).

Mode-Specific Masking Steps
---------------------------
The following table indicates which flavors of masking are applied to
images from each type of observing mode.

.. |c| unicode:: U+2713 .. checkmark

+--------------------------+-----+-----+-----+
|                          |     | Mode|     |
+--------------------------+-----+-----+-----+
| Masking                  | IFU | MOS |  FS |
+==========================+=====+=====+=====+
| IFU Slices\ :sup:`1`     | |c| |     |     |
+--------------------------+-----+-----+-----+
| Slits/Slitlets\ :sup:`1` |     | |c| | |c| |
+--------------------------+-----+-----+-----+
| MSA_FAILED_OPEN          | |c| | |c| | |c| |
+--------------------------+-----+-----+-----+
| NaN Pixels               | |c| | |c| | |c| |
+--------------------------+-----+-----+-----+
| FS Region                | |c| | |c| |     |
+--------------------------+-----+-----+-----+
| Reference Pix            | |c| | |c| | |c| |
+--------------------------+-----+-----+-----+
| Outliers                 | |c| | |c| | |c| |
+--------------------------+-----+-----+-----+

:sup:`1`\ The application of these steps can be turned on and off via
the step parameter `mask_spectral_regions`. This parameter controls
whether the "IFU Slices" and "Slits/Slitlets" portions of the masking
are applied.

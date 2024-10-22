.. _outlier-detection-ifu:

Integral Field Unit (IFU) Data
==============================

This module serves as the interface for applying ``outlier_detection`` to IFU
observations, like those taken with NIRSpec and MIRI. A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a :py:class:`~jwst.datamodels.ModelContainer` object,
serves as the basic format for all processing performed by this step.

After launch it was discovered that the bad pixels on the MIRI detectors vary with time.
The pixels varied from usable to unusable, and at times, back to usable  on a time frame that was too short
(sometimes as short as 2 days)  to fold into the bad pixel mask applied in the 
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline. At this time it is believed that NIRSpec IFU data
also have bad pixels that vary with time, though the time variation is still under study.
The outlier detection step is designed to flag these pixels as outliers, in addition
to cosmic ray hits that were not flagged by the jump step.

The basis of the outlier detection flagging for IFU data  is to look for pixels on the detector
that are regularly discrepant from their neighbors, with a sharper division than could be explained
by the detector PSF.
This routine performs the following operations:

#. Extract parameter settings for the input ModelContainer and merge them with any user-provided values.
   See :ref:`outlier detection arguments <outlier_detection_step_args>` for the full list of parameters.

#. Loop over cal files, computing nearest-neighbor differences for each pixel
   in the along-dispersion direction.
   For MIRI, with the dispersion axis along the y axis, the neighbors that are used to
   to find the differences are to the left and right of each pixel being examined.
   For NIRSpec, with the dispersion along the x axis, the neighbors that are used to
   find the differences are above and below the pixel being examined.
   The smaller of the two (left, right or up, down) differences is stored as the difference value for each
   pixel. This avoids artifacts from bright edges.

#. Compare the nearest-neighbor differences across science exposures to find the minimum
   neighbor difference at each detector pixel.

#. Determine a local spatial median of the minimum difference array using a median filter with a kernel size
   set by the user according to the ``kernel_size`` parameter.

#. Normalize the minimum difference array by the local median.

#. Select outliers by flagging those normalized minimum values larger than the ``threshold_percent``
   parameter.

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

.. _outlier-detection-ifu:

Outlier Detection for IFU Data
==============================

This module serves as the interface for applying ``outlier_detection`` to IFU
observations, like those taken with NIRSpec and MIRI.  The code implements the
basic outlier detection algorithm searching for pixels that are consistent outliers
in the calibrated images created by the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.
After launch it was discovered the
bad pixels on the MIRI detectors vary with time. The pixels varied from usable to unusable, and at times,
back to usable  on a time frame that was too short (sometimes as short as 2 days)  to fold into the bad pixel mask applied
in the :ref:`calwebb_detector1 <calwebb_detector1>` pipeline. 
At this time it is believed that NIRSpec IFU data also have bad pixels that vary with time, though the time variation is
still under study.

The basis of the outlier detection flagging for IFU data  is to look for pixels on the detector
that are regularly discrepant from their neighbors, with a sharper division than could be explained
by the detector PSF. The algorithm flags pixels that are  outliers when compared to their neighbors for a set of
input files contained in an association. The neighbor pixel differences are the neighbors in the spatial direction.
For MIRI data neighbor differences are found to the left and right of every
science pixel, while for NIRSpec data the neighbor differences are
determined from pixels above and below every science pixel. The difference between the MIRI MRS and NIRPSpec algorithm for finding the
spatial pixel differences is due to the opposite dispersion directions between the two instruments. 
The  pixel differences for each input model 
in the association are determined and stored in a stack of pixel differences. 
For each pixel the minimum difference
through this stack is determined and normalized. The normalization uses a local median of the difference array
(set by the kernel size). A pixel is flagged as an outlier if  this normalized minimum difference
is greater than the input threshold percentage.  Pixels that are found to be outliers are flaged in in the DQ array.



.. automodapi:: jwst.outlier_detection.outlier_detection_ifu

.. _outlier-detection-ifu:

Outlier Detection for IFU Data
==============================

This module serves as the interface for applying ``outlier_detection`` to IFU
observations, like those taken with NIRSpec and MIRI.  The code implements the
basic outlier detection algorithm searching for pixels that are consistent outliers
on the data after the :ref:`calwebb_spec2 <calwebb_spec2>`. After launch it was discovered the
bad pixels on the MIRI detectors vary with time. The pixels varied from usable to unsable, and at times,
back to usable  on a time frame that was too short (sometimes as short as 2 days)  to fold into the bad pixel mask applied
in the :ref:`calwebb_detector1 <calwebb_detector1>`. 
At this time it is believed that NIRSpec IFU data also have bad pixels that vary with time, though the time variation is
still under study.

An algorithm was developed to flag pixels that are  outliers when compared to their neighbors for a set of
input files contained in association. The neighbor pixel differences are defined by the dispersion axis. For MIRI data,
the dispersion axis is along the y-axis and differences are found to the left and right of every
science pixel. For NIRSpec data, the dispersion axis is along the x-axis and differences are
found to between the pixels above and below every science pixel. The pixel differences for each input model 
in the assications is determined and this stored in a stack of pixel differences. 
For each pixel the minimum difference
through this stack is determined and normalized. The normalization uses a local median of the difference array
(set by the kernel size). The input model DQ array is flagged as an outlier if this normalized minimum difference
is greater than the input threshold percentage. 


.. automodapi:: jwst.outlier_detection.outlier_detection_ifu

.. _outlier-detection-tso:

Time-Series Observations (TSO) Data
===================================

This module serves as the interface for applying ``outlier_detection`` to time
series observations.

Normal imaging data benefit from combining all integrations into a
single image. TSO data's value, however, comes from looking for variations from one
integration to the next.  The outlier detection algorithm, therefore, gets run with 
a few variations to accommodate the nature of these 3D data.
A :py:class:`~jwst.datamodels.CubeModel` object serves as the basic format for all
processing performed by this step. This routine performs the following operations:

#. Convert input data into a CubeModel (3D data array) if a ModelContainer
   of 2D data arrays is provided.

#. Do not attempt resampling; data are assumed to be aligned and have an identical WCS.
   This is true automatically for a CubeModel.

#. Apply a bad pixel mask to the input data based on the input DQ arrays and the ``good_bits``
   parameter.

#. Compute a median cube by combining all planes in the CubeModel pixel-by-pixel using a
   rolling-median algorithm, in order to flag outliers integration-by-integration but
   preserve real time variability. The ``rolling_window_width`` parameter specifies the
   number of integrations over which to compute the median.

#. If the ``save_intermediate_results`` parameter is set to True, write the rolling-median
   CubeModel to disk with the suffix ``_median.fits``.

#. Perform a statistical comparison frame-by-frame between the rolling-median cube and 
   the input data. The formula used is the same as for imaging data without resampling:
   
   .. math:: | image\_input - image\_median | > SNR * input\_err

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

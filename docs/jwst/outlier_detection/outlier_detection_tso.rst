.. _outlier-detection-tso:

Outlier Detection for TSO Data
==============================

Time-series observations (TSO) result in input data stored as a 3D CubeModel
where each plane in the cube represents a separate integration without changing the
pointing.  Normal imaging data benefit from combining all integrations into a
single image. TSO data's value, however, comes from looking for variations from one
integration to the next.  The outlier detection algorithm, therefore, gets run with 
a few variations to accomodate the nature of these 3D data.

#. Input data is converted into a CubeModel (3D data array) if a ModelContainer
   of 2D data arrays is provided.

#. The median image is created without resampling the input data

   * The median image is created by combining all planes in the 
     CubeModel pixel-by-pixel using a rolling-median algorithm, in order
     to flag outliers integration-by-integration but preserve real time variability.
   * The ``rolling_window_width`` parameter specifies the number of integrations over
     which to compute the median. The default is 25. If the number of integrations
     is less than or equal to ``rolling_window_width``, a simple median is used instead.
   * The ``maskpt`` parameter sets the percentage of the weight image values to
     use, and any pixel with a weight below this value gets flagged as "bad" and
     ignored when the median is taken.
   * The rolling-median CubeModel (3D data array) is written out to disk as `_<asn_id>_median.fits`
     if the ``save_intermediate_results`` parameter is set to True.
   * All integrations are aligned already, so no resampling or shifting needs to be performed

#. A matched median gets created by combining the single median frame with the 
   noise model for each input integration.

#. A statistical comparison is performed between the matched median with each input integration.  

#. The input data model DQ arrays are updated with the mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection_tso

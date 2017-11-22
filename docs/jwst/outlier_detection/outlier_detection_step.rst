.. _outlier_B7_design_:

Python Step Design: OutlierDetectionStep
========================================

This module provides the sole interface to all methods of performing outlier
detection on JWST observations.  The ``OutlierDetectionStep`` supports multiple
algorithms and determines the appropriate algorithm for the type of observation
being processed.  This step supports:


This step uses the following logic to apply the appropriate algorithm to the
input data:

* Interpret inputs (ASN table, ModelContainer or CubeModel)
  to identify all input observations to be processed

* Read in type of exposures in input by interpreting ``meta.exposure.type`` from inputs

* Read in parameters set by user.

* Select outlier detection algorithm based on exposure type

  - **Images**: like those taken with NIRCam, will use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection`
  - **coronagraphic observations**: use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` with resampling turned off
  - **time-series observations(TSO)**: both imaging and spectroscopic modes, will use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` with resampling turned off
  - **NIRISS and MIRI IFU observations**: use :py:class:`~jwst.outlier_detection.outlier_detection_ifu.OutlierDetectionIFU`
  - **long-slit spectroscopic observations**: use :py:class:`~jwst.outlier_detection.outlier_detection_spec.OutlierDetectionSpec`

* Read in reference files used by outlier detection, currently:

  - **readnoise**: based on :py:class:`~jwst.datamodels.ReadnoiseModel`
  - **gain**: based on :py:class:`~jwst.datamodels.GainModel`

* Instantiate and run outlier detection class determined for the exposure type
  using reference files and parameter values interpreted from inputs.

* Return input_models with DQ arrays updated with flags for identified outliers


.. automodapi:: jwst.outlier_detection.outlier_detection_step

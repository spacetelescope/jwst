.. _outlier_B7_design_:

Python Step Design: OutlierDetectionStep
========================================

This module provides the sole interface to all methods of performing outlier
detection on JWST observations.  The ``OutlierDetectionStep`` supports multiple
algorithms and determines the appropriate algorithm for the type of observation
being processed.  This step supports:

* **Image modes**: 'NRC_IMAGE', 'MIR_IMAGE', 'NIS_IMAGE', 'FGS_IMAGE'
* **Spectroscopic modes**: 'NRC_GRISM', 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NIS_WFSS'
* **Time-Series-Observation(TSO) Spectroscopic modes**: 'NIS_SOSS', 'MIR_LRS-SLITLESS', 'NRC_TSGRISM', 'NRS_BRIGHTOBJ'
* **IFU Spectroscopic modes**: 'NRS_IFU', 'MIR_MRS'
* **TSO Image modes**:'NRC_TSIMAGE'
* **Coronagraphic Image modes**: 'NRC_CORON', 'MIR_LYOT', 'MIR_4QPM'


This step uses the following logic to apply the appropriate algorithm to the
input data:

* Interpret inputs (ASN table, ModelContainer or CubeModel)
  to identify all input observations to be processed

* Read in type of exposures in input by interpreting ``meta.exposure.type`` from inputs

* Read in parameters set by user.

* Select outlier detection algorithm based on exposure type

  - **Images**: like those taken with NIRCam, will use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection`
  - **coronagraphic observations**: use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` with resampling turned off
  - **Time-Series Observations(TSO)**: both imaging and spectroscopic modes, will use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` with resampling turned off
  - **NIRSpec and MIRI IFU observations**: use :py:class:`~jwst.outlier_detection.outlier_detection_ifu.OutlierDetectionIFU`
  - **Long-slit spectroscopic observations**: use :py:class:`~jwst.outlier_detection.outlier_detection_spec.OutlierDetectionSpec`

* Read in reference files used by outlier detection, currently:

  - **readnoise**: based on :py:class:`~jwst.datamodels.ReadnoiseModel`
  - **gain**: based on :py:class:`~jwst.datamodels.GainModel`

* Instantiate and run outlier detection class determined for the exposure type
  using reference files and parameter values interpreted from inputs.

* Return input_models with DQ arrays updated with flags for identified outliers


.. automodapi:: jwst.outlier_detection.outlier_detection_step

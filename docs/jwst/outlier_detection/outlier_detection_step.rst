.. _outlier_design:

Python Step Design: OutlierDetectionStep
-----------------------------------------

This module provides the sole interface to all methods of performing outlier
detection on JWST observations.  The ``outlier_detection`` step supports multiple
algorithms and determines the appropriate algorithm for the type of observation
being processed.  This step supports:

* **Image modes**: 'FGS_IMAGE', 'MIR_IMAGE', 'NRC_IMAGE', 'NIS_IMAGE'
* **Spectroscopic modes**: 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC'
* **Time-Series-Observation(TSO) Spectroscopic modes**: 'MIR_LRS-SLITLESS', 'NRC_TSGRISM', 'NIS_SOSS', 'NRS_BRIGHTOBJ'
* **IFU Spectroscopic modes**: 'MIR_MRS', 'NRS_IFU'
* **TSO Image modes**: 'NRC_TSIMAGE'
* **Coronagraphic Image modes**: 'MIR_LYOT', 'MIR_4QPM', 'NRC_CORON'


This step uses the following logic to apply the appropriate algorithm to the
input data:

* Interpret inputs (ASN table, ModelContainer or CubeModel)
  to identify all input observations to be processed

* Read in type of exposures in input by interpreting ``meta.exposure.type`` from inputs

* Read in parameters set by user

* Select outlier detection algorithm based on exposure type

  - **Images**: like those taken with NIRCam, will use
    :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` as described
    in :ref:`outlier-detection-imaging`
  - **Coronagraphic observations**:
    use :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` with
    resampling turned off as described in :ref:`outlier-detection-imaging`
  - **Time-Series Observations(TSO)**: both imaging and spectroscopic modes, use
    :py:class:`~jwst.outlier_detection.outlier_detection.OutlierDetection` with resampling
    turned off as described in :ref:`outlier-detection-imaging`
  - **IFU observations**: use
    :py:class:`~jwst.outlier_detection.outlier_detection_ifu.OutlierDetectionIFU` as
    described in :ref:`outlier-detection-ifu`
  - **Long-slit spectroscopic observations**: use
    :py:class:`~jwst.outlier_detection.outlier_detection_spec.OutlierDetectionSpec` as
    described in :ref:`outlier-detection-spec`

* Instantiate and run outlier detection class determined for the exposure type
  using parameter values interpreted from inputs.

* Return input models with DQ arrays updated with flags for identified outliers


.. automodapi:: jwst.outlier_detection.outlier_detection_step

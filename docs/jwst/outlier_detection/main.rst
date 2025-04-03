.. _outlier_design:

Overview
========

This module provides the sole interface to all methods of performing outlier
detection on JWST observations.

Processing multiple datasets together allows for the identification of bad pixels
or cosmic rays that remain in each of the input images, often at levels which
were not detectable by the :ref:`jump <jump_step>` step.
The ``outlier_detection`` step supports multiple
algorithms and determines the appropriate algorithm for the type of observation
being processed.  This step supports:

* **Image modes**: 'FGS_IMAGE', 'MIR_IMAGE', 'NRC_IMAGE', 'NIS_IMAGE'
   - See :ref:`outlier-detection-imaging` for algorithm details
* **Slit-like Spectroscopic modes**: 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC'
   - See :ref:`outlier-detection-spec` for algorithm details
* **Time-Series-Observation (TSO) modes**: 'MIR_LRS-SLITLESS', 'NRC_TSGRISM', 'NIS_SOSS', 'NRS_BRIGHTOBJ', 'NRC_TSIMAGE', as well as TSOs obtained with the MIRI imager ('MIR_IMAGE' with TSOVISIT=True).
   - See :ref:`outlier-detection-tso` for algorithm details
* **IFU Spectroscopic modes**: 'MIR_MRS', 'NRS_IFU'
   - See :ref:`outlier-detection-ifu` for algorithm details
* **Coronagraphic Image modes**: 'MIR_LYOT', 'MIR_4QPM', 'NRC_CORON'
   - See :ref:`outlier-detection-coron` for algorithm details

This step uses the following logic to apply the appropriate algorithm to the
input data:

#. Interpret inputs (Association, ModelContainer, ModelLibrary, or CubeModel)
   to identify all input observations to be processed

#. Read in parameters set by user. See :ref:`outlier_detection_step_args` for the full list
   of parameters.

#. Select outlier detection algorithm based on exposure type in input model ``meta.exposure.type``.

#. Instantiate and run outlier detection class determined for the exposure type
   using parameter values interpreted from inputs.

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

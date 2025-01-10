Outlier Detection Transition Guide
==================================

As of jwst release 1.18, the outlier detection step has been replaced by
separate steps for the various observing modes. This page serves as a transition
guide for migrating to use the new steps.

The new steps and the exposure types processed through those steps are:

.. list-table:: Outlier Detection Steps
   :header-rows: 1
   
   * - Step
     - Exposure Types
   * - `OutlierDetectionCoronStep`
     - 'MIR_LYOT', 'MIR_4QPM', 'NRC_CORON'
   * - `OutlierDetectionIFUStep`
     - 'MIR_MRS', 'NRS_IFU'
   * - `OutlierDetectionImagingStep`
     - 'FGS_IMAGE', 'MIR_IMAGE', 'NRC_IMAGE', 'NIS_IMAGE'
   * - `OutlierDetectionSpecStep`
     - 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC'
   * - `OutlierDetectionTSOStep`
     - 'MIR_LRS-SLITLESS', 'NRC_TSGRISM', 'NIS_SOSS', 'NRS_BRIGHTOBJ', 'NRC_TSIMAGE'

Command-line syntax
-------------------

To run the individual step on its own, use the following command line syntax:

::

   strun outlier_detection_coron input.fits

instead of the previous syntax:

::
    
   strun outlier_detection input.fits # NO LONGER WORKS


The new steps have the same optional arguments as the original step, but only
the ones relevant to that observing mode.
For example, the `OutlierDetectionCoronStep` has the
`--save_intermediate_results`, `--good_bits`, `--snr`, and `--maskpt` arguments,
but no longer has the `--kernel` argument, which didn't have any effect
for coronagraphic observations in the first place.

Running the new steps within a pipeline works the same as before. To set options
for the new steps within a pipeline, use the `--steps.<step_name>` syntax, e.g.

::

   strun calwebb_image3 --steps.outlier_detection_coron.save_intermediate_results=True input.fits

Note that for `calwebb_spec3`, the outlier detection step is different for IFU data than
for slit-like data. For IFU data, use `--steps.outlier_detection_ifu`, and for fixed slit
and MSA data, use `--steps.outlier_detection_spec`.

Python syntax
-------------
The new steps can be run in Python using the following syntax:

::
  
   from jwst.outlier_detection_coron import OutlierDetectionCoronStep
   OutlierDetectionCoronStep.call(input.fits, save_intermediate_results=True)

The new pipelines can be run in Python using, e.g.,:

::

   from jwst.pipeline import Image3Pipeline
   Image3Pipeline.call(input.fits, steps={'outlier_detection_coron': {'save_intermediate_results': True}})

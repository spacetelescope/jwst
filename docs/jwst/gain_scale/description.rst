Description
============

The ``gain_scale`` step rescales pixel values in JWST countrate
science data products in order to correct for the effect of using
a non-standard detector gain setting. The countrate data are
rescaled to make them appear as if they had been obtained using
the standard gain setting.

This currently only applies to NIRSpec exposures that are read out
using a subarray pattern, in which case a gain setting of 2 is used
instead of the standard setting of 1. Note that this only applies
to NIRSpec subarray data obtained after April 2017, which is when
the change was made in the instrument flight software to use gain=2.
NIRSpec subarray data obtained previous to that time used the
standard gain=1 setting.

The `gain_scale` step is applied at the end of the
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline, after the
:ref:`ramp_fit <ramp_fitting_step>` step has been applied. It is applied
to both the "rate" and "rateints" products from
:ref:`ramp_fit <ramp_fitting_step>`, if both
types of products were created. The science (SCI) and error (ERR)
arrays are multiplied by the gain factor, and the Poisson
variance (VAR_POISSON) and read noise variance (VAR_RNOISE) arrays
are multiplied by the square of the gain factor.

The scaling factor is obtained from the "GAINFACT" keyword in the
header of the gain reference file. Normally the
:ref:`ramp_fit <ramp_fitting_step>` step
reads that keyword value during its execution and stores the value in
the science data "GAINFACT" keyword, so that the gain reference file
does not have to be loaded again by the ``gain_scale`` step. If, however,
the step does not find that keyword populated in the science data, it
loads the gain reference file to retrieve it. If all attempts to
find the scaling factor fail, the step is skipped.

Gain reference files for instruments or modes that use the standard
gain setting will typically not have the "GAINFACT" keyword in their
header, which causes the ``gain_scale`` step to be skipped. Alternatively,
gain reference files for modes that use the standard gain can have
GAINFACT=1.0, in which case the correction is benign.

Upon successful completion of the step, the "S_GANSCL" keyword in the
science data is set to "COMPLETE".

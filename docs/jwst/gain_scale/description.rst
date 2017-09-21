Description
============

The `gain_scale` step rescales pixel values in JWST countrate
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

The `gain_scale` step is applied at the end of the `calwebb_detector1`
pipeline, after the `ramp_fit` step has been applied. It is applied
to both the `rate` and `rateints` products from `ramp_fit`, if both
types of products were created. The science (`SCI`) and error (`ERR`)
arrays are both rescaled.

The scaling factor is obtained from the `GAINFACT` keyword in the
header of the gain reference file. Normally the `ramp_fit` step will
read that keyword value during its execution and store the value in
the science data keyword `GAINFACT`, so that the gain reference file
does not have to be loaded again by the `gain_scale` step. If, however,
the step does not find that keyword populated in the science data, it
will load the gain reference file to retreive it. If all attempts to
find the scaling factor fail, the step will be skipped.

Gain reference files for instruments or modes that use the standard
gain setting will typically not have the `GAINFACT` keyword in their
header, which will cause the `gain_scale` step to be skipped. Alternatively,
gain reference files for modes that use the standard gain can have
`GAINFACT=1.0`, in which case the correction will be benign.

Upon successful completion of the step, the `S_GANSCL` keyword in the
science data will be set to "COMPLETE."

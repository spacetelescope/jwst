Description
===========

:Class: `jwst.rscd.RscdStep`
:Alias: rscd

Assumptions
-----------
This correction is currently only implemented for MIRI data.
It is assumed this step occurs before the dark subtraction, but after linearity
correction.

Background
__________

The MIRI Focal Plane System (FPS) consists of the detectors and the electronics to control them.
There are a number of non-ideal detector and readout effects that produce reset offsets,
nonlinearities at the start of an integration, non-linear ramps with increasing signal,
latent images, and drifts in the slopes.

The MIRI readout electronics use field effect transistors (FETs) in their operation,
which have been shown to be the source of the ramp offsets, the nonlinearities at the start
of an integration, and the overall changes in the slopes. The FET acts as a switch to allow
charge to build up and to also initialize (clear) the charge.
However, the reset FETs do not instantaneously reset the level. Instead, the exponential
adjustment of the FET after a reset causes the initial frames in an integration to be offset
from their expected values. Between exposures, the MIRI detectors are continually reset;
however, for a multiple integration exposure there is a single reset between integrations.
The effects of this decay are reduced in the first integration because a number
of resets have occurred from the last exposure and the effect has partially decayed  by the time
it takes to read out the last exposure, set up the next exposure, and begin exposing.
Because of these physical transients, the JWST pipeline includes a dedicated RSCD step
that automatically applies these group_skip flags to ensure that the subsequent Jump Detection
and Ramp Fitting steps only use the linear portion of the integration

The Reset Switch Charge Decay (RSCD) step corrects for these effects by simply
flagging the first N groups as DO_NOT_USE. 


Algorithm
_________

The RSCD (Reset Switch Charge Decay) step identifies and flags groups at the beginning of
MIRI integrations that are affected by non-linear transients. These transients are caused
by the exponential settling of the detector’s Field Effect Transistor (FET) switches
immediately following a reset.

This step flags the N groups at the beginning of all integrations
as bad (the "DO_NOT_USE" bit is set in the
GROUPDQ flag array). The number of groups to skip is depends on the readout pattern,
subarray size and integration number. To maintain the statistical viability of the ramp, the step
only applies flags if the integration contains at least three more groups than the required
skip number (Groups > `n_skip` + 3). If this condition is not met, the step is bypassed to allow
later pipeline stages enough data points to perform a linear fit.

Standard RSCD correction flags the first N groups as DO_NOT_USE. However, for very bright sources,
the pixel might saturate immediately after those skipped groups. If the algorithm blindly skips
the RSCD groups, it might leave the pixel with zero or one valid group, making it impossible to
calculate a flux (slope). Instead the algorithm "backs off" the number of skipped groups for
specific pixels that are at risk of losing all their unsaturated data.
Because reducing the RSCD skip introduces some non-linear FET transient data back into the fit,
these pixels are flagged The algorithm flags them in the PIXELDQ array as FLUX_ESTIMATED to warn
the users that the flux value may be slightly biased by the RSCD effect.
If only one group is left valid, the algorithm records information header (more information given
the table below). This allows the :ref:`ramp_fitting <ramp_fitting_step>` to still derive a flux value (provided the user has enabled suppress_one_group = False).


This step results in the data contained in the the first `n_skip` groups
being excluded from subsequent steps, such as :ref:`jump detection <jump_step>`
and :ref:`ramp_fitting <ramp_fitting_step>`.

Only the GROUPDQ array is modified. The SCI and ERR arrays remain unchanged. The PIXELDQ arrays are
only updated in the case of bright saturating data when the RSCD skip count is lowered
to preserve valid groups. In this case,  the FLUX_ESTIMATED flag isadded  to indicate a potential
bias from the FET transient.

RSCD Keywords   Meaning

=============   ======================================================================================================
INT1SKIP        Number of groups skipped in integration 1
INT2SKIP        Number of groups skipped in integration 2 and higher
INT1UGP1        Number of pixels where 1st the group is kept for integration 1
INT2UGP2        Number of pixels where 1st the group is kept for integration 2 and higher
INT1BORS        Number of pixels where RSCD reduced the groups skipped due to saturation for integration 1
INT2BORS        Number of pixels where RSCD reduced the groups skipped due to saturation for integration 2 and higher

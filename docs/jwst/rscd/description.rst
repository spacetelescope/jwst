Description
===========

:Class: `jwst.rscd.RscdStep`
:Alias: rscd

Assumptions
-----------
This correction is currently only implemented for MIRI data and is only applied
to integrations after the first integration (i.e. this step does not correct the
first integration).
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
The effects of this decay are not measurable in the first integration because a number
of resets have occurred from the last exposure and the effect has decayed away by the time
it takes to read out the last exposure, set up the next exposure, and begin exposing.
There are low level reset effects in the first integration that are related to the strength of the dark
current and can be removed with an integration-dependent dark.

The Reset Switch Charge Decay (RSCD) step corrects for these effects by simply
flagging the first N groups as DO_NOT_USE.  An actual correction algorithm allowing for the first N groups to be
used is under development.

Algorithm
_________

This correction is only applied to integrations > 1.
This step flags the N groups at the beginning of all 2nd and higher integrations
as bad (the "DO_NOT_USE" bit is set in the
GROUPDQ flag array), but only if the total number of groups in each
integration is greater than N+3.
This results in the data contained in the the first N groups
being excluded from subsequent steps, such as :ref:`jump detection <jump_step>`
and :ref:`ramp_fitting <ramp_fitting_step>`.
No flags are added if NGROUPS <= N+3, because doing so would leave too few good
groups to work with in later steps.

Only the GROUPDQ array is modified. The SCI, ERR, and PIXELDQ arrays are unchanged.


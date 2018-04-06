
Description
===========

Assumptions
-----------
This correction is currently only implemented for MIRI data and is only for integrations
after the first integration (i.e. this step does not correct the first integration).
It is assumed this step occurs before the dark subtraction, but after linearity.

Background
__________

The MIRI Focal Plane System (FPS) consists of the detectors and the electronics to control them.
There are a number of non-ideal detector and readout effects which produce reset offsets,
nonlinearities at the start of an integration, non-linear ramps with increasing signal,
latent images and drifts in the slopes. 

The manner in which the MIRI readout electronics operate have been
shown to be the source of the reset offsets and nonlinearities at the start of the integration.
Basically the MIRI reset electronis use field effect transisitors (FETs) in their operation.  The FET acts as switch
to allow charge to build up and to also initialize (clear) the charge. However, the reset FETS do not instanteously
reset the level, instead the expontenial adjustment of the  FET after a reset causes the initial frames in an integration
to be offset from their expected values.  The Reset Switch Charge Decay (RSCD) step corrects for the slow adjustment of the
FET output to its asymptotic level after a reset. This correction is made for integrations > 1 and is based on the signal
level in the last frame of the previous integration in the exposure. Between exposures the MIRI detectors
are conintually reset; however for a multiple integration exposure there is a single reset between integrations.
The reset switch charge decay has an e-folding time scale ~ 1.3 * frame time so the affects of this decay are
not measureable in the first integration  because a number of resets have occurred from the last exposure and
the effect has decayed away by the time it takes to  readout out the last exposure, set up the next exposure and begin
exposuring. There are low level reset effects in the first integration that are related to the strength of the dark
current and can be removed with an integration dependent dark.


For MIRI multiple integration data, the reset switch decay causes the
the initial groups  in  integrations after the first one  to be offset from
their expected  linear accumulation of signal.
The most significant deviations ocurr in groups 1 and 2. The amplitude of the different between the expected value
and the measured value varies for even and odd rows and is related to the signal in the last frame of the last integration.

The MIRI reset electronics also cause an offset, zero-point offset in multiple integration data. Subsequent integrations after
the first integration start a lower DN level. The ampitude of this offset is proporational
to the singnal level in the previous integration. Fortunately this offset is constant for all the groups in the integration,
thus has no impact on the slopes determined for each integration.




Algorithm
_________
The rscd correction step applies the reset switch charge decay reference file. Based on READOUT pattern
(FAST or SLOW) and  Subarray type (FULL or one of MIRI defined subarray types) the reference file contains
the scale factor and decay time (tau)  for even and odd rows to corrected for the reset effects. The
accumulated DN of the pixel  from the previous integration is estimated by extrapolating the ramp using the second to last 
and third to last groups. For each pixel the group values are corrected according the formula:

    corrected value = input vaule + dn_accumulated * scale * exp(-T / tau),

where T is the time since the last frame in the last integration.

The correction aglorithm is slightly modified if the previous integration saturated. In this case the scale factor 
in the above equation is calculated using an estimation of the what the last frame in the previous integration
would of been if saturation did not exist. This estimate last frame is from a linear fit of the non-saturating
groups in the ramp. 
 

Subarrays
----------

Currently the rscd correction for subarray data is the same as it is for full array data. However,
we anticipate a seperate set of correction coefficients in the future.

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

The manner in which the MIRI readout electronics operate have been
shown to be the source of the ramp offsets, nonlinearities at the start of the integration, and overall changes in slopes.
Basically the MIRI reset electronics use field effect transistors (FETs) in their operation.  The FET acts as a switch
to allow charge to build up and to also initialize (clear) the charge. However, the reset FETS do not instantaneously
reset the level, instead the exponential adjustment of the  FET after a reset causes the initial frames in an integration
to be offset from their expected values. Between exposures the MIRI detectors
are continually reset; however for a multiple integration exposure there is a single reset between integrations.
The effects of this decay are
not measurable in the first integration  because a number of resets have occurred from the last exposure and
the effect has decayed away by the time it takes to read out the last exposure, set up the next exposure and begin
exposing. There are low level reset effects in the first integration that are related to the strength of the dark
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

..
    This text refers to an earlier version of the enhanced RSCD correction.
    It needs updating to the latest version of this correction once that has been
    decided and the code updated.

    The step applies an exponential decay correction based on coefficients in the "RSCD"
    reference file. The reference files are selected based on readout pattern
    (READPATT=FAST or SLOW) and subarray type (FULL or one of the MIRI defined subarray types).
    The reference file contains the information necessary to derive the scale factor and decay time
    to correct for the reset effects. The correction differs for even and odd row numbers.

    The correction to be added to the input data has the form::

        corrected data = input data data + dn_accumulated * scale * exp(-T / tau)  (Equation 1)

    where T is the time since the last group in the previous integration, tau is the exponential time constant and
    dn_accumulated is the DN level that was accumulated for the pixel from the previous integration.
    Because of the last frame effect the value of the last group in an integration is not measured accurately. Therefore,
    the accumulated DN of the pixel from the previous integration (last group value)  is estimated by extrapolating
    the ramp using the second to last  and third to last groups.

    In the case where the previous integration does not saturate the :math:`scale` term in Equation 1  is determined as follows:

     :math:`scale = b{1}* [Counts{2}^{b{2}} * [1/exp(Counts{2}/b{3}) -1] \; \; Equation \;  2`

    The terms :math:`b{2}` and :math:`b{3}` are read in from the RSCD reference file.
    The following two additional equations are needed to calculate the :math:`b{1}` and :math:`Counts{2}` terms:

    	  :math:`b{1} = ascale * (illum_{zpt} + illum_{slope}*N + illum2* N^2) \; \; (Equation \; 2.1)`
    	  :math:`Counts{2} = Final \, DN \, in \, the \,  last \, group \, in \; the \, last \, integration
    	  \, - Crossover \, Point \; \; (Equation \; 2.2)`


    In equation 2.1, N is the number of groups per integration and :math:`ascale`, :math:`illum_{zpt}`,
    :math:`illum_{slope}`, and :math:`illum2` are read in from the RSCD reference file. The :math:`Crossover \, Point`
    in equation 2.2 is also read in from the RSCD reference file.

    If the previous integration saturates, the  :math:`scale` term in Equation 1 is found in the  following manner:

       :math:`scale_\text{sat} = slope * Counts{3} + sat_\text{mzp} \; \; (Equation \; 3)`

    where :math:`Counts{3}` is an  estimate of what the last group in the previous integration would have been if
    saturation did not exist. The :math:`slope` in equation 3  is calculated according to the formula:

       :math:`slope = sat_{zp} + sat_{slope} * N + sat_2*N^2 + evenrow_{corrections} \; \; (Equation 3.1)`.

    The terms :math:`sat_\text{mzp}`, :math:`sat_{zp}`, :math:`sat_2`, :math:`evenrow_{corrections}`
    are read in from the RSCD reference file.

    All fourteen  parameters :math:`tau`, :math:`b{1}`, :math:`b{2}`, :math:`b{3}`, :math:`illum_{zpt}`,
    :math:`illum_{slope}`, :math:`illum2`, :math:`Crossover Point`, :math:`sat_{zp}`, :math:`sat_{slope}`, :math:`sat_2`,
    :math:`sat_{scale}`, :math:`sat_\text{mzp}`, and :math:`evenrow_{corrections}` are found in the RSCD reference files.
    There is a separate set for even and odd rows for each readout (READPATT) mode and subarray type.

    Subarrays
    ----------

    Currently the RSCD correction for subarray data is the same as it is for full array data. However,
    we anticipate a separate set of correction coefficients in the future.

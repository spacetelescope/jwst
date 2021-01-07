Description
===========
Based on a model, this step computes the number of traps that are
expected to have captured or released a charge during an exposure.
The released charge is proportional to the persistence signal, and
this will be subtracted (group by group) from the science data.  An
image of the number of filled traps at the end of the exposure will
be written as an output file, in order to be used as input for
correcting the persistence of a subsequent exposure.

There may be an input traps-filled file (defaults to 0), giving the number
of traps that are filled in each pixel.  There is one plane of this 3-D image
for each "trap family," sets of traps having similar capture and decay
parameters.  The traps-filled file is therefore coupled with the trappars
reference table, which gives parameters family-by-family.  There are currently
three trap families.

If an input traps-filled file was specified, the contents of that file will
be updated (decreased) to account for trap decays from the EXPEND of the
traps-filled file to the EXPSTART of the current science file before starting
the processing of the science data.

When processing a science image, the traps-filled file is the basis for
computing the number of trap decays, which are computed group-by-group.  On
the other hand, the trap-density file is the basis for predicting trap
captures, which are computed at the end of each integration.  The
traps-filled file will be updated (decreased by the number of traps that
released a charge) after processing each group of the science image.  The
traps-filled file will then be increased by the number of traps that were
predicted to have captured a charge by the end of each integration.

There is often a reset at the beginning of each integration, and if so,
that time (a frame time) will be included in the trap capture for each
integration, and it will be included for the tray decay for the first
group of each integration.

The number of trap decays in a given time interval is computed as follows:

.. math::
    n\_decays = trapsfilled \cdot (1 - exp(-\Delta t / \tau))

where trapsfilled is the number of filled traps, i.e. the value of the
traps-filled image at the
beginning of the time interval, for the current trap family and at the
current pixel; :math:`\Delta t` is the time interval (seconds) over which
the decay is computed; and :math:`\tau` is the reciprocal of the absolute
value of the decay parameter (column name "decay_param") for the current
trap family.  Since this is called for each group, the value of the
traps-filled image must be updated at the end of each group.

For each pixel, the persistence in a group is the sum of the trap decays
over all trap families.  This persistence is subtracted from the science
data for the current group. Pixels that have large persistence values
subtracted from them are flagged in the DQ array, as information to the
user (see the Step Arguments section).

Trap capture is more involved than is trap decay.  The computation of trap
capture is different for an impulse (e.g. a cosmic-ray event) than for a
ramp, and saturation also affects capture.  Computing trap capture needs
an estimate of the ramp slope, and it needs the locations (pixel number and
group number) of cosmic-ray jumps.  At the time of writing, the ``persistence``
step is run before the ``jump`` step, so the GROUPDQ array in the input to
``persistence`` does not contain the information that is required to account
for cosmic-ray events.

Because the ``persistence`` step is run before ``ramp_fit``, the persistence step does
not have the value of the slope, so the step must compute its own estimate
of the slope.  The algorithm is as follows.  First of all, the slope must be
computed before the loop over groups in which trap decay is computed and
persistence is corrected, since that correction will in general change the
slope.  Within an integration, the difference is taken between groups of the
ramp.  The difference is set to a very large value if a group is saturated.
(The "very large value" is the larger of :math:`10^5` and twice the maximum
difference between groups.)  The difference array is then sorted.  All the
differences affected by saturation will be at the high end.  Cosmic-ray
affected differences should be just below, except for jumps that are smaller
than some of the noise.  We can then ignore saturated values and jumps by
knowing how many of them there are (which we know from the GROUPDQ array).
The average of the remaining differences is the slope.  The slope is needed
with two different units.  The `grp_slope` is the slope in units of DN
(data numbers) per group.  The `slope` is in units of
(DN / persistence saturation limit) / second, where "persistence saturation
limit" is the (pixel-dependent) value (in DN) from the PERSAT reference file.

The number of traps that capture charge is computed at the end of each
integration.  The number of captures is computed in three phases:  the
portion of the ramp that is increasing smoothly from group to group;
the saturated portion (if any) of the ramp; the contribution from
cosmic-ray events.

For the smoothly increasing portion of the ramp, the time interval over
which traps capture charge is
nominally :math:`nresets \cdot tframe + ngroups \cdot tgroup`
where nresets is the number of resets at the beginning of the integration,
tframe is the frame time, and tgroup is the group time.
However, this time must be reduced by the group time multiplied by the
number of groups for which the data value exceeds the persistence saturation
limit.  This reduced value is :math:`Delta t` in the expression below.

The number of captures in each pixel during the integration is:

.. math::
    trapsfilled = 2 \cdot &(trapdensity \cdot slope^2 \\
                      &\cdot (\Delta t^2 \cdot (par0 + par2) / 2
                       + par0 \cdot (\Delta t \cdot \tau + \tau^2) \\
                       &\cdot exp(-\Delta t / \tau) - par0 \cdot \tau^2))

where par0 and par2 are the values from columns "capture0" and "capture2"
respectively, from the trappars reference table, and :math:`\tau` is the
reciprocal of the absolute value from column "capture1", for the row
corresponding to the current trap family.  trapdensity is the
relative density of traps, normalized to a median of 1.  :math:`\Delta t`
is the time interval in seconds over which
the charge capture is to be computed, as described above.  slope is the
ramp slope (computed before the loop over groups), in units of fraction
of the persistence saturation limit per second.  This returns the number
of traps that were predicted to be filled during the integration, due to
the smoothly increasing portion of the ramp.  This is passed as input to
the function that computes the additional traps that were filled due to
the saturated portion of the ramp.

"Saturation" in this context means that the data value in a group exceeds
the persistence saturation limit, i.e. the value in the PERSAT reference
file.  filled_during_integration is (initially) the array of the number of
pixels that were filled, as returned by the function for the smoothly
increasing portion of the ramp.  In the function for computing decays
for the saturated part of the ramp, for pixels that are saturated in the
first group, filled_during_integration
is set to :math:`trapdensity \cdot par2` (column "capture2").  This accounts
for "instantaneous" traps, ones that fill over a negligible time scale.

The number of "exponential" traps (as opposed to instantaneous) is:

.. math::
    exp\_filled\_traps = filled\_during\_integration - trapdensity \cdot par2

and the number of traps that were empty and could be filled is:

.. math::
    empty\_traps = trapdensity \cdot par0 - exp\_filled\_traps

so the traps that are filled depending on the exponential component is:

.. math::
    new\_filled\_traps = empty\_traps \cdot (1 - exp(-sattime / \tau))

where sattime is the duration in seconds over which the pixel was saturated.

Therefore, the total number of traps filled during the current integration is:

.. math::
    filled\_traps = filled\_during\_integration + new\_filled\_traps

This value is passed to the function that computes the additional traps
that were filled due to cosmic-ray events.

The number of traps that will be filled due to a cosmic-ray event depends
on the amount of time from the CR event to the end of the integration.  Thus,
we must first find (via the flags in the GROUPDQ extension) which groups and
which pixels were affected by CR hits.  This is handled by looping over
group number, starting with the second group (since we currently don't flag
CRs in the first group), and selecting all pixels with a jump.  For these
pixels, the amplitude of the jump is computed to be the difference between
the current and previous groups minus grp_slope (the slope in DN per group).
If a jump is negative, it will be set to zero.

If there was a cosmic-ray hit in group number k, then

.. math::
    \Delta t = (ngroups - k - 0.5) \cdot tgroup

is the time from the CR-affected group to the end of the integration, with
the approximation that the CR event was in the middle (timewise) of the group.
The number of traps filled as a result of this CR hit is:

.. math::
    crfilled = 2 \cdot trapdensity \cdot jump
                \cdot (par0 \cdot (1 - exp(-\Delta t / \tau)) + par2)

and the number of filled traps for the current pixel will be incremented
by that amount.

Input
=====
The input science file is a RampModel.

A trapsfilled file (TrapsFilledModel) may optionally be passed as input
as well.  This normally would be specified unless the previous exposure
with the current detector was taken more than several hours previously,
that is, so long ago that persistence from that exposure could be ignored.
If none is provided, an array filled with 0 will be used as the starting
point for computing new traps-filled information.

Output
======
The output science file is a RampModel, a persistence-corrected copy of
the input data.

A second output file will be written, with suffix "_trapsfilled".  This
is a TrapsFilledModel, the number of filled traps at each pixel at the end
of the exposure.  This takes into account the capture of charge by traps
due to the current science exposure, as well as the release of charge
from traps given in the input trapsfilled file, if one was specified.  Note
that this file will always be written, even if no input_trapsfilled file
was specified.  This file should be passed as input to the next run of the
persistence step for data that used the same detector as the current run.
Pass this file using the input_trapsfilled argument.

If the user specifies ``save_persistence=True``, a third output file will
be written, with suffix "_output_pers".  This is a RampModel matching the
output science file, but this gives the persistence that was subtracted
from each group in each integration.

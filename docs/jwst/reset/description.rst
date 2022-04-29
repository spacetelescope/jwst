Description
===========

:Class: `jwst.reset.ResetStep`
:Alias: reset

The reset correction is a MIRI step that attempts to correct
for the reset anomaly effect. This effect is caused by the non-ideal behavior of the FET upon resetting in the dark
causing the initial frames in an integration to be offset from their expected values. Another MIRI effect caused by
resetting the detectors is the RSCD effect (see :ref:`rscd <rscd_step>`). 


Assumptions
-----------
The reset correction is a MIRI-specific correction. It will not be applied to data from  other instruments. 



Background
__________

For MIRI exposures, the initial groups in each integration suffer from two
effects related to the resetting of the detectors. The first effect is that the
first few groups after a reset do not fall
on the expected linear accumulation of signal.
The most significant deviations occur in groups 1 and 2.
This behavior is relatively uniform detector-wide. The second effect,
on the other hand, is the appearance of
significant extra spatial structure in these initial
groups, before fading out in later groups.

The reset anomaly effect fades out by ~group 15 for full array data. It takes a few more groups
for the effect to fade away on subarray data. The time constant of the effect seems to be closely
related to the group number and not time since reset.

For multiple integration data, the reset anomaly also varies in amplitude
for the first few integrations before settling down to a relatively
constant correction for integrations greater than the second integration for full array
data. Because of the shorter readout time, the subarray data requires a few
more integrations before the effect is relatively stable from integration
to integration.

Algorithm
_________
The reset correction step applies the reset reference file.
The reset reference file contains an integration dependent
correction for the first N groups, where N is defined by the reset
correction reference file.

The format of the reset reference file is NCols X NRows X NGroups X NInts.
For full frame data, the current implementation uses a reset anomaly reference file,
which contains a correction for the first 15 groups for
all integrations.  The reference file contains two corrections: one for the first integration
and a second one for all other integrations. The correction 
was determined so that the correction is forced to be zero on group 15.  For each integration in the input science data,
the reset corrections are subtracted, group-by-group, integration-by-
integration. If the input science data contains more groups than the
reset correction, then correction for those groups is zero. If the
input science data contains more integrations than the reset correction
then the correction corresponding to the last integration in the reset file
is used.

There is a single, NCols X NRowss, DQ flag image for all the integrations.
The reset DQ flag array  are combined with the science PIXELDQ array using
numpy's bitwise_or function. The ERR arrays of the science data are
currently not modified at all.

Subarrays
----------

The reset correction is  subarray-dependent, therefore this
step makes no attempt to extract subarrays from the reset reference file to
match input subarrays. It instead relies on the presence of matching subarray
reset reference files in the CRDS. In addition, the number of NGROUPS and NINTS
for subarray data varies from the full array data as well as from each other.

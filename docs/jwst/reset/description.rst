Description
===========

.. note:: The reset step is not currently applied to MIRI exposures, but
          will likely be reinstated in a future build.

Assumptions
-----------
The reset correction is a MIRI-specific correction. It is
assumed that the input science data have *NOT* had the zero group (or bias)
subtracted. The reset correction should not remove the
bias signal from the science exposure, therefore the reset correction
for the first group is defined to be zero.

Background
__________

For MIRI exposures, the initial groups in each integration suffer from two
effects related to the resetting of the detectors. The first effect is that the
first few groups after a reset do not fall
on the expected linear accumulation of signal.
The most significant deviations ocurr in groups 1 and 2.
This behavior is relatively uniform detector-wide. The second effect,
on the other hand, is the appearance of
significant extra spatial structure in these initial
groups, before fading out in later groups.

The time constant associated with the reset anomaly is
roughly a minute so for full array data the effect has faded out
by ~group 20. On subarray data, where the read time  depends on
the size of the subarray, the reset anomaly affects more
groups in an integration.

For multiple integration data, the reset anomaly also varies in amplitude
for the first set of integrations before settling down to a relatively
constant correction for integrations greater than four for full array
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
The current implementation uses a reset anomaly reference file for
full array data  containing a correction for the first 30 groups for
integrations 1-4. The reference file
was determined so that the correction is forced to be zero on the last
group for each integration.  For each integration in the input science data,
the reset corrections are subtracted, group-by-group, integration-by-
integration. If the input science data contains more groups than the
reset correction, then correction for those groups is zero. If the
input science data contains more integrations than the reset correction
then the correction corresponding to the last intergration in the reset file
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

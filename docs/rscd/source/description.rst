Description
===========

Assumptions
-----------
This correction is currently only implemented for MIRI data and is only for integrations
after the first integration (i.e. this step does not correct the first integration).
It is assumed this step occurs before the dark subtraction.

Background
__________

The MIRI Focal Plane System (FPS) consists of the detectors and the electronics to control them. 
There are a number of non-ideal detector and readout effects which produce reset offsets, 
nonlinearities at the start of an integration, non-linear ramps with increasing signal, 
latent images and drifts in the slopes. The manner in which the MIRI readout electronics operate have been
shown to be the source of the reset offsets and nonlinearities at the start of the integration. 
The reset offset, also described as the zero-point offset, are easily seen in multiple integration data, where the
first integration starts at a lower DN level than subsequent integrations. The ampitude of this offset is proporational
to the singnal level in the previous integration. Fortunately this offset is constant for all the groups in the integration, 
thus have no impact on the slopes determined for each integration.  

The readout electronics have also been shown to be the source of the nonlinearities at the start of the integration.
Basically the MIRI reset electronis use field effect transisitors (FETs) in their operation.  The FET acts a switch to allow
charge to build up and to also initialize (clear) the charge. However, the reset FETS do not instanteously reset the level, instead
the expontenail adjustment of the  FET after a reset causes the initial frames in an integration to be offset from their expected values. 


The MIRI reset electroncis use a number of field effect transistors (FETs) in their operation. It has been shown that the 
slow adjustment of the FET output 
The Reset Switch Charge Decay (RSCD) step corrects for the slow adjustment of the FET output to its asymptotic level after a reset. This correction is made for integration > 1 and is based on the signal level in the last frame of the previous integration in the exposure.
Currently this step is only implemented for MIRI data. For MIRI data
the initial groups  in an integration suffer from effects related 
to the resetting of the detectors. The first effect is that the
first few samples starting an integration after a reset do not fall
on the expected linear accumulation of signal. 
The most significant deviations ocurr in groups 1 and 2. 
This behavior is relatively uniform detector-wide. The second effect, 
on the other hand, is the appearance of 
significant extra spatial structure that appears on in these initial
groups, before fading out by later groups.  

The time constant associated with the reset anomaly is
roughly a minute so for full array data the effect has faded out
by ~group 20. On subarray data, where the read time  depends on
the size of the subarray, the reset anomaly affects more 
groups in an integration.
 
For multiple integration data the reset anomaly also varies in amplitude
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


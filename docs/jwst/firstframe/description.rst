Description
===========

:Class: `jwst.firstframe.FirstFrameStep`
:Alias: firstframe

The MIRI first frame correction step flags the first group in every integration
as bad (the "DO_NOT_USE" data quality flag is added to the GROUPDQ array), but
only if the total number of groups per integration is greater than 3.
This results in the data contained in the first group
being excluded from subsequent steps, such as jump detection and ramp fitting.
No flags are added if NGROUPS <= 3, because doing so would leave too few good
groups to work with in later steps.

Only the GROUPDQ array is modified. The SCI, ERR, and PIXELDQ arrays are unchanged.

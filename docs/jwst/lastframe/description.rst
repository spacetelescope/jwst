Description
===========

:Class: `jwst.lastframe.lastframe_step.LastFrameStep`
:Alias: lastframe

The last frame correction step is only applied to MIRI data and flags the
final group in each integration as bad (the "DO_NOT_USE" bit is set in the
GROUPDQ flag array), but only if the total number of groups in each
integration is greater than 2.
This results in the data contained in the last group
being excluded from subsequent steps, such as jump detection and ramp fitting.
No flags are added if ``NGROUPS <= 2``, because doing so would leave too few good
groups to work with in later steps.

Only the GROUPDQ array is modified. The SCI, ERR, and PIXELDQ arrays are unchanged.

Step Arguments
--------------
The ``lastframe`` step does not have any step-specific arguments.

Reference Files
---------------
The ``lastframe`` step does not use any reference files.

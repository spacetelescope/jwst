Description
===========
Based on a model, this step computes the number of traps that are
expected to have captured or released a charge during an exposure.
The released charge is proportional to the persistence signal, and
this will be subtracted (group by group) from the science data.  An
image of the number of filled traps at the end of the exposure will
be written as an output file, in order to be used as input for
correcting the persistence of a subsequent exposure.

Input
=====
The input science file is a RampModel.

A trapsfilled file (TrapsFilledModel) may optionally be passed as input
as well.  This normally would be specified unless the previous exposure
with the current detector was taken more than several hours previously,
that is, so long ago that persistence from that exposure could be ignored.

Output
======
The output science file is a RampModel, a persistence-corrected copy of
the input data.

A second output file will be written, with suffix "_trapsfilled".  This
is a TrapsFilledModel, the number of filled traps at each pixel at the end
of the exposure.  This takes into account the capture of charge by traps
due to the current science exposure, as well as the release of charge
from traps shown in the input trapsfilled file, if one was specified.

If the user specified ``save_persistence=True``, a third output file will
be written, with suffix "_output_pers".  This is a RampModel matching the
output science file, but this gives the persistence that was subtracted
from each group in each integration.

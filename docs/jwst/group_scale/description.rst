Description
===========

The `group_scale` step rescales pixel values in raw JWST science
data products to correct for instances where on-board frame averaging 
did not result in the proper downlinked values.

When multiple frames are averaged together on-board into a single
group, the sum of the frames is computed and then the sum is
divided by the number of frames to compute the average. Division by
the number of frames is accomplished by simply bit-shifting the
sum by an appropriate number of bits, corresponding to the 
decimal value of the number of frames. For example, when 2 frames
are averaged into a group, the sum is shifted by 1 bit to achieve
the equivalent of dividing by 2, and for 8 frames, the sum is
shifted by 3 bits. The number of frames that are averaged into a
group is recorded in the ``NFRAMES`` header keyword in science
products and the divisor that was used is recorded in the
``FRMDIVSR`` keyword.

This method results in the correct average only when NFRAMES is a
power of 2. When NFRAMES is not a power of 2, the next largest
divisor is used to perform the averaging. For example, when
NFRAMES=5, a divisor of 8 (bit shift of 3) is used to compute the
average. This results in averaged values for each group that
are too low by the factor NFRAMES/FRMDIVSR. This step rescales the
pixel values by multiplying all groups in all integrations by the
factor FRMDIVSR/NFRAMES.

The step decides whether rescaling is necessary by comparing the
values of the NFRAMES and FRMDIVSR keywords. If they are equal,
then the on-board averaging was computed correctly and this step
is skipped. In this case, the calibration step status keyword
``S_GRPSCL`` is set to "SKIPPED." If the keyword values are not
equal, rescaling is applied and the ``S_GRPSCL`` keyword is set
to "COMPLETE".

It is assumed that this step is always applied to raw data
before any other processing is done to the pixel values and hence
rescaling is applied only to the SCI data array of the input
product. It assumes that the ERR array has not yet been populated
and hence there's no need for rescaling that array.
The input GROUPDQ and PIXELDQ arrays are not affected by this step.

MIRI FASTGRPAVG mode
^^^^^^^^^^^^^^^^^^^^

The MIRI detector readout pattern "FASTGRPAVG" results in individual
frames being averaged together into a group, but the on-board
averaging process is done differently than for other instruments.
This results in a situation where the FRMDIVSR keyword gets assigned
a value of 4, while NFRAMES still has a value of 1, despite the fact
that 4 frames were actually averaged together to produce each
downlinked group. This mismatch in keyword values would cause
the ``group_scale`` step to think that rescaling needs to be applied.

To work around this issue, the original values of the number of frames
per group and the number of groups per integration that are downlinked
from the instrument are stored in the special keywords "MIRNFRMS" and
"MIRNGRPS", respectively, so that their values are preserved. During
Stage 1 processing in the pipeline, the value of the NFRAMES keyword is
computed from MIRNFRMS * FRMDIVSR. The result is that when 4 frames
are averaged together on board, both NFRAMES and FRMDIVSR will have a
value of 4, which allows the ``group_scale`` step to correctly
determine that no rescaling of the data is necessary.

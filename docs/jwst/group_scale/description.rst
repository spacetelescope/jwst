Description
============

The `group_scale` step rescales pixel values in raw JWST science
data products in order to correct for the effect of using a value
of NFRAMES for on-board frame averaging that is not a power of 2.

When multiple frames are averaged together on-board into a single
group, the sum of the frames is computed and then the sum is
divided by the number of frames to compute an average. Division by
the number of frames is accomplished by simply bit-shifting the
sum by an appropriate number of bits, corresponding to the 
decimal value of the number of frames. For example, when 2 frames
are averaged into a group, the sum is shifted by 1 bit to achieve
the equivalent of dividing by 2, and for 8 frames, the sum is
shifted by 3 bits. The number of frames that are averaged into a
group is recorded in the `NFRAMES` header keyword in science
products and the divisor that was used is recorded in the
`FRMDIVSR` keyword.

This method only results in the correct average when NFRAMES is a
power of 2. When NFRAMES is not a power of 2, the next largest
divisor is used to perform the averaging. For example, when
NFRAMES=5, a divisor of 8 (bit shift of 3) is used to compute the
average. This results in averaged values for every group that
are too low by the factor NFRAMES/FRMDIVSR.

This step rescales raw pixel values to the correct level by
multiplying all groups in all integrations by the factor
FRMDIVSR/NFRAMES.

It is assumed that this step will always be applied to raw data
before any other processing is done to the pixel values and hence
rescaling is applied only to the SCI data array of the input
product. It assumes that the ERR array has not yet been populated
and hence there's no need for rescaling that array.

If the step detects that the values of NFRAMES and FRMDIVSR are
equal to one another, which means the data were scaled correctly
on-board, it skips processing and returns the input data unchanged.
In this case, the calibration step status keyword `S_GRPSCL` will
be set to `SKIPPED`. After successful correction of data that
needs to be rescaled, the `S_GRPSCL` keyword will be set to
`COMPLETE`.

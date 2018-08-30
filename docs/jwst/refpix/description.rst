Description
===========

Overview
--------

With a perfect detector and readout electronics, the signal in any given
readout would differ from that in the previous readout only as a result
of detected photons.  In reality, the readout electronics imposes its own
signal on top of this.  In its simplest form, the amplifiers add a constant
value to each pixel, and this constant value is different from amplifier to
amplifier in a given group, and varies from group to group for a given
amplifier.  The magnitude of this variation is of the order of a few counts.
In addition, superposed on this signal is a variation that is mainly with
row number that seems to apply to all amplifiers within a group.

The refpix step corrects for these drifts by using the reference
pixels. NIR detectors have their reference pixels in a 4-pixel wide strip
around the edge of the detectors that are completely insensitive to light,
while the MIR detectors have a 4 columns (1 for each amplifier) of reference
pixels at the left and right edges of the detector.  They also have data read
through a fifth amplifier, which is called the reference output, but these
data are not currently used in any refpix correction.

The effect is more pronounced for the NIR detectors than for the MIR
detectors.

Input details
-------------

The input file must be a ramp, and it should contain both a science
('SCI') extension and a data quality ('DQ') extension.  The latter
extension is normally added by the dq_init step, so running this
step is a prerequisite for the refpix step.

Algorithm
---------

The algorithm for the NIR and MIR detectors is different.

NIR Detector Data
-----------------

#. The data from most detectors will have been rotated and/or flipped from their detector frame in order to give them the same orientation and parity in the telescope focal plane.  The first step is to transform them back to the detector frame so that all NIR and MIR detectors can be treated equivalently.
#. It is assumed that a superbias correction has been performed.
#. For each integration, and for each group:
    #. Calculate the mean value in the top and bottom reference pixels.  The reference pixel means for each amplifier are calculated separately, and the top and bottom means are calculated separately.  Optionally, the user can choose to calculate the means of odd and even columns separately by using the ``--odd_even_columns`` runtime parameter, as evidence has been found that there is a significant odd-even column effect in some datasets.  Bad pixels (those whose DQ flag has the DO_NOT_USE bit set) are not included in the calculation of the mean.
    #. The mean is calculated as a clipped mean with a 3-sigma rejection threshold.
    #. Average the top and bottom reference pixel mean values
    #. Subtract each mean from all pixels that the mean is representative of, i.e. by amplifier and using the odd mean for the odd pixels and even mean for even pixels if this option is selected.
    #. If the ``--use_side_ref_pixels`` option is selected, use the reference pixels up the side of the A and D amplifiers to calculate a smoothed reference pixel signal as a function of row.  A running median of height set by the runtime parameter ``side_smoothing_length`` (default value 11) is calculated for the left and right side reference pixels, and the overall reference signal is obtained by averaging the left and right signals.  A multiple of this signal (set by the runtime parameter ``side_gain``, which defaults to 1.0) is subtracted from the full group on a row-by-row basis.
#. Transform the data back to the JWST focal plane, or DMS, frame.

MIR Detector Data
-----------------

#. MIR data is already in the detector frame, so no flipping/rotation is needed
#. Subtract the first group from each group within an integration.
#. For each integration, and for each group after the first:
    #. Calculate the mean value in the reference pixels for each amplifier. The left and right side reference signals are calculated separately. Optionally, the user can choose to calculate the means of odd and even rows separately using the ``--odd_even_rows`` runtime parameter, as it has been found that there is a significant odd-even row effect.  Bad pixels (those whose DQ flag has the DO_NOT_USE bit set) are not included in the calculation of the mean. The mean is calculated as a clipped mean with a 3-sigma rejection threshold.
    #. Average the left and right reference pixel mean values
    #. Subtract each mean from all pixels that the mean is representative of, i.e. by amplifier and using the odd mean for the odd row pixels and even mean for even row pixels if this option is selected.
    #. Add the first group of each integration back to each group.

At the end of the refpix step, the S_REFPIX keyword is set to 'COMPLETE'.

Subarrays
---------

Subarrays are treated slightly differently.  Once again, the data are flipped
and/or rotated to convert to the detector frame

NIR Data
--------

If the odd_even_columns flag is set to True, then the clipped means of all
reference pixels in odd-numbered columns and those in even numbered columns
are calculated separately, and subtracted from their respective data columns.
If the flag is False, then a single clipped mean is calculated from all of
the reference pixels in each group and subtracted from each pixel.

.. note::

  In subarray data, reference pixels are identified by the dq array having the
  value of REFERENCE_PIXEL (defined in datamodels/dqflags.py).  These values
  are populated when the dq_init step is run, so it is important to run this
  step before running the refpix step on subarray data.

If the science dataset has at least 1 group with no valid reference pixels,
the refpix step is skipped and the S_REFPIX header keyword is set to 'SKIPPED'.

MIR Data
--------

The refpix correction is skipped for MIRI subarray data.

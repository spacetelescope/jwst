Description
===========

:Class: `jwst.refpix.RefPixStep`
:Alias: refpix

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

The ``refpix`` step corrects for these drifts by using the reference
pixels. NIR detectors have their reference pixels in a 4-pixel wide strip
around the edge of the detectors that are completely insensitive to light,
while the MIR detectors have 4 columns (1 column for each amplifier) of reference
pixels at the left and right edges of the detector.  They also have data read
through a fifth amplifier, which is called the reference output, but these
data are not currently used in any refpix correction.

The effect is more pronounced for the NIR detectors than for the MIR
detectors.

Input details
-------------

The input file must be a 4-D ramp and it should contain both a science
(SCI) extension and a pixel data quality (PIXELDQ) extension. The PIXELDQ
extension is normally populated by the ``dq_init`` step, so running that
step is a prerequisite for the ``refpix`` step.

Algorithms
----------

The algorithms for the NIR and MIR detectors are somewhat different.
An entirely different algorithm for NIRSpec IRS2 readout mode is
described in IRS2_.

NIR Detector Data
+++++++++++++++++

#. The data from most detectors will have been rotated and/or flipped from
   their detector frame in order to give them the same orientation and parity
   in the telescope focal plane.  The first step is to transform them back to
   the detector frame so that all NIR and MIR detectors can be treated equivalently.
#. It is assumed that a superbias correction has been performed.
#. For each integration and for each group:

    #. Calculate the mean value in the top and bottom reference pixels.
       The reference pixel means for each amplifier are calculated separately,
       and the top and bottom means are calculated separately.
       Optionally, the user can choose to calculate the means of odd and even
       columns separately by using the ``--odd_even_columns`` step parameter,
       because evidence has been found that there is a significant odd-even
       column effect in some datasets.  Bad pixels (those whose DQ flag has the
       "DO_NOT_USE" bit set) are not included in the calculation of the mean.
    #. The mean is calculated as a clipped mean with a 3-sigma rejection threshold
       using the ``scipy.stats.sigmaclip`` method.
    #. Average the top and bottom reference pixel mean values
    #. Subtract each mean from all pixels that the mean is representative of,
       i.e. by amplifier and using the odd mean for the odd pixels and even mean
       for even pixels if this option is selected.
    #. If the ``--use_side_ref_pixels`` option is selected, use the reference pixels
       up the side of the A and D amplifiers to calculate a smoothed reference pixel
       signal as a function of row.  A running median of height set by the step
       parameter ``side_smoothing_length`` (default value 11) is calculated for the
       left and right side reference pixels, and the overall reference signal is
       obtained by averaging the left and right signals.  A multiple of this signal
       (set by the step parameter ``side_gain``, which defaults to 1.0) is
       subtracted from the full group on a row-by-row basis.  Note that the ``odd_even_rows``
       parameter is ignored for NIR data when the side reference pixels are processed.
    #. Transform the data back to the JWST focal plane, or DMS, frame.

MIR Detector Data
+++++++++++++++++

#. MIR data are always in the detector frame, so no flipping/rotation is needed.
#. Subtract the first group from each group within an integration.
#. For each integration, and for each group after the first:

    #. Calculate the mean value in the reference pixels for each amplifier.
       The left and right side reference signals are calculated separately.
       Optionally, the user can choose to calculate the means of odd and even
       rows separately using the ``--odd_even_rows`` step parameter, because
       it has been found that there is a significant odd-even row effect.
       Bad pixels (those whose DQ flag has the "DO_NOT_USE" bit set) are not
       included in the calculation of the mean. The mean is calculated as a
       clipped mean with a 3-sigma rejection threshold using the
       ``scipy.stats.sigmaclip`` method.  Note that the ``odd_even_columns``,
       ``use_side_ref_pixels``, ``side_smoothing_length`` and ``side_gain``
       parameters are ignored for MIRI data.
    #. Average the left and right reference pixel mean values.
    #. Subtract each mean from all pixels that the mean is representative of,
       i.e. by amplifier and using the odd mean for the odd row pixels and even
       mean for even row pixels if this option is selected.
    #. Add the first group of each integration back to each group.

At the end of the refpix step, the S_REFPIX keyword is set to "COMPLETE".

NIRCam Frame 0
--------------

If a frame zero data cube is present in the input data, the image corresponding
to each integration is corrected in the same way as the regular science data and
passed along to subsequent pipeline steps.

Subarrays
---------

Subarrays are treated slightly differently.  Once again, the data are flipped
and/or rotated to convert to the detector frame.

NIR Data
++++++++

For single amplifier readout (NOUTPUTS keyword = 1):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the ``odd_even_columns`` flag is set to True, then the clipped means of all
reference pixels in odd-numbered columns and those in even numbered columns
are calculated separately, and subtracted from their respective data columns.
If the flag is False, then a single clipped mean is calculated from all of
the reference pixels in each group and subtracted from each pixel.

.. note::

  In subarray data, reference pixels are identified by the PIXELDQ array having the
  value of "REFERENCE_PIXEL" (defined in datamodels/dqflags.py).  These values
  are populated when the ``dq_init`` step is run, so it is important to run that
  step before running the ``refpix`` step on subarray data.

If the science dataset has at least 1 group with no valid reference pixels,
the step is skipped and the S_REFPIX header keyword is set to 'SKIPPED'.

The ``use_side_ref_pixels``, ``side_smoothing_length``, ``side_gain`` and
``odd_even_rows`` parameters are ignored for these types of data.

For 4 amplifier readout (NOUTPUTS keyword = 4):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the NOUTPUTS keyword is 4 for a subarray exposure, then the data are calibrated
the same as for full-frame exposures.  The top/bottom reference values are obtained from available
reference pixel regions, and the side reference values are used if available.  If only 1 of the
top/bottom or side reference regions are available, they are used, whereas if both are available they
are averaged.  If there are no top/bottom or side reference pixels available, then that part of
the correction is omitted.  The routine will log which parameters are valid according to
whether valid reference pixels exist.

MIR Data
++++++++

The refpix correction is skipped for MIRI subarray data.

.. _IRS2:

NIRSpec IRS2 Readout Mode
+++++++++++++++++++++++++

This section describes -- in a nutshell -- the procedure for applying the
reference pixel correction for data read out using the IRS2 readout pattern.
See the JdoxIRS2_ page for for an overview, and see Rauscher2017_ for
details.

The raw data include both the science data and interspersed reference
pixel values.  The time to read out the entire detector includes not only
the time to read each pixel of science ("normal") data and some of the
reference pixels, but also time for the transition between reading normal
data and reference pixels, as well as additional overhead at the end of
each row and between frames.  For example, it takes the same length of time
to jump from reading normal pixels to reading reference pixels as it does
to read one pixel value, about ten microseconds.

Before subtracting the reference pixel and reference output values from
the science data, some processing is done on the reference values, and the
CRDS reference file factors are applied.  IRS2 readout is only used for
full-frame data, never for subarrays.  The full detector is read out
by four separate amplifiers simultaneously, and the reference output is
read at the same time.  Each of these five readouts is the same size,
640 by 2048 pixels (for IRS2).  If the CRDS reference file includes a
DQ (data quality) BINTABLE extension, interleaved reference pixel values
will be set to zero if they are flagged as bad in the DQ extension.
The next step in this processing is to
copy the science data and the reference pixel data separately to temporary
1-D arrays (both of length 712 * 2048); this is done separately for each
amp output.  The reference output is also copied to such an array, but
there is only one of these.  When copying a pixel of science or reference
pixel data to a temporary array, the elements are assigned so that the
array indexes increase with and correspond to the time at which the
pixel value was read.  That means that the change in readout direction
from one amplifier to the next is taken into account when the data are
copied, and that there will be gaps (array elements with zero values),
corresponding to the times when reference pixels were read (or science
data, depending on which is being copied), or corresponding to the
overheads mentioned in the previous paragraph.  The gaps will then be
assigned values by interpolation (cosine-weighted, then Fourier filtered).
Note that the above is done for every group.

The ``alpha`` and ``beta`` arrays that were read from the CRDS reference
file are next applied, and this is done in Fourier space.  These are
applied to the temporary 1-D arrays of reference pixel data and to the
reference output array.  ``alpha`` and ``beta`` have shape (4, 712 * 2048)
and data type Complex64 (stored as pairs of Float32 in the reference file).
The first index corresponds to the sector number for the different
output amplifiers.  ``alpha`` is read from columns 'ALPHA_0', 'ALPHA_1',
'ALPHA_2', and 'ALPHA_3'.  ``beta`` is read from columns 'BETA_0',
'BETA_1', 'BETA_2', and 'BETA_3'.

For each integration, the following is done in a loop over groups.

Let ``k`` be the output number, i.e. an index for sectors 0 through 3.
Let ``ft_refpix`` be an array of shape (4, 712 * 2048); for each output
number ``k``, ``ft_refpix[k]`` is the Fourier transform of the temporary
1-D array of reference pixel data.  Let ``ft_refout`` be the Fourier
transform of the temporary 1-D array of reference output data.  Then: ::

    for k in range(4):
        ft_refpix_corr[k] = ft_refpix[k] * beta[k] + ft_refout * alpha[k]

For each ``k``, the inverse Fourier transform of ``ft_refpix_corr[k]`` is
the processed array of reference pixel data, which is then subtracted from
the normal pixel data over the range of pixels for output ``k``.

.. _JdoxIRS2: https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-readout-modes-and-patterns/nirspec-irs2-detector-readout-mode
.. _Rauscher2017: http://adsabs.harvard.edu/abs/2017PASP..129j5003R

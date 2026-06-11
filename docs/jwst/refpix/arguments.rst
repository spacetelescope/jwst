Step Arguments
==============
The ``refpix`` step has the following step-specific arguments.

``--odd_even_columns`` (bool, default=True)
  The top/bottom reference
  signal is calculated and applied separately for even and odd-numbered
  columns.  This applies to NIR
  data only.

``--use_side_ref_pixels`` (bool, default=True)
  The side reference pixels
  are used to calculate a reference signal for each row, which is subtracted
  from the data.  This applies to NIR
  data only.

``--side_smoothing_length`` (int, default=11)
  Height of
  the window used in calculating the running median when calculating the side
  reference signal. This applies to NIR
  data only when the ``--use_side_ref_pixels`` option is selected.

``--side_gain`` (float, default=1.0)
  The factor that the side
  reference signal is multiplied by before subtracting from the group
  row-by-row.  This applies to NIR
  data only when the ``--use_side_ref_pixels`` option is selected.

``--odd_even_rows`` (bool, default=True)
  The reference signal is
  calculated and applied separately for even- and odd-numbered rows.
  This applies to MIR data only.

``--ovr_corr_mitigation_ftr`` (float, default=3.0)
  The factor to avoid overcorrection of intermittently bad reference
  pixels in the IRS2 algorithm. This is the number of sigmas away
  from the mean. This applies
  only to NIRSpec data taken with IRS2 mode.

``--preserve_irs2_refpix`` (bool, default=False)
  Interleaved reference pixels
  in IRS2 mode will be processed along with the normal pixels and preserved
  in the output.  This option is intended for calibration or diagnostic reductions
  only. For normal science operation, this argument should always be `False`,
  so that interleaved pixels are stripped before continuing processing.

``--irs2_mean_subtraction`` (bool, default=False)
  Apply or skip a mean offset
  subtraction before IRS2 correction.  Mean values are computed across reference pixels
  sorted by amplifier and detector column parity.  Setting this to `True` may help reduce
  alternating column noise in some exposures.

``--refpix_algorithm`` (str, default='median')
  This is only relevant for all NIR full-frame
  data, and can be set to ``'median'`` to use the running median or
  ``'sirs'`` to use the Simple Improved Reference Subtraction (SIRS).

``--sigreject`` (float, default=4.0)
  The number of sigmas to reject as outliers in the
  SIRS algorithm.

``--gaussmooth`` (float, default=1.0)
  The width of Gaussian smoothing kernel to use as
  a low-pass filter.

``--halfwidth`` (int, default=30)
  The half-width of convolution kernel to build.

``--siglimit`` (float, default=3.0)
  The number of standard deviations to use in
  the iterative sigma clipping algorithm that calculates the mean of the
  reference pixels.  The value is used as both
  the lower and upper bounds in the sigma clipping algorithm.

Step Arguments
==============

The reference pixel correction step has seven step-specific arguments:

*  ``--odd_even_columns``

If the ``odd_even_columns`` argument is given, the top/bottom reference
signal is calculated and applied separately for even- and odd-numbered
columns.  The default value is True, and this argument applies to NIR
data only.

*  ``--use_side_ref_pixels``

If the ``use_side_ref_pixels`` argument is given, the side reference pixels
are used to calculate a reference signal for each row, which is subtracted
from the data.  The default value is True, and this argument applies to NIR
data only.


*  ``--side_smoothing_length``

The ``side_smoothing_length`` argument is used to specify the height of
the window used in calculating the running median when calculating the side
reference signal. The default value is 11, and this argument applies to NIR
data only when the ``--use_side_ref_pixels`` option is selected.

*  ``--side_gain``

The ``side_gain`` argument is used to specify the factor that the side
reference signal is multiplied by before subtracting from the group
row-by-row.  The default value is 1.0, and this argument applies to NIR
data only when the ``--use_side_ref_pixels`` option is selected.

*  ``--odd_even_rows``

If the ``odd_even_rows`` argument is selected, the reference signal is
calculated and applied separately for even- and odd-numbered rows.  The
default value is True, and this argument applies to MIR data only.

*  ``--ovr_corr_mitigation_ftr``

This is a factor to avoid overcorrection of intermittently bad reference
pixels in the IRS2 algorithm. This factor is the number of sigmas away
from the mean. The default value is 3.0, and this argument applies
only to NIRSpec data taken with IRS2 mode.

*  ``--preserve_irs2_refpix``

If the ``preserve_irs2_refpix`` argument is set, interleaved reference pixels
in IRS2 mode will be processed along with the normal pixels and preserved
in the output.  This option is intended for calibration or diagnostic reductions
only. For normal science operation, this argument should always be False,
so that interleaved pixels are stripped before continuing processing.

*  ``--refpix_algorithm``

The ``refpix_algorithm`` argument is only relevant for all NIR full-frame
data, and can be set to 'median' (default) to use the running median or
'sirs' to use the Simple Improved Reference Subtraction (SIRS).

*  ``--sigreject``

The ``sigreject`` argument is the number of sigmas to reject as outliers in the
SIRS algorithm. The value is expected to be a float.

*  ``--gaussmooth``

The ``gaussmooth`` argument is the width of Gaussian smoothing kernel to use as
a low-pass filter. The numerical value is expected to be a float.

*  ``--halfwidth``

The ``halfwidth`` argument is the half-width of convolution kernel to build. The
numerical value is expected to be an integer.

*  ``--irs2_mean_subtraction``

The ``irs2_mean_subtraction`` argument is a boolean to apply or skip a mean offset
subtraction before IRS2 correction.  Mean values are computed across reference pixels
sorted by amplifier and detector column parity.  Setting this option to True may help reduce
alternating column noise in some exposures.


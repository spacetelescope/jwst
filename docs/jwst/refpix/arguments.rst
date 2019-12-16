Step Arguments
==============

The reference pixel correction step has five step-specific arguments:

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

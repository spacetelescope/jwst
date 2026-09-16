Step Arguments
==============
The ``align_refs`` step has the following optional arguments.

``--median_box_length`` (integer, default=3)
  The box size to use when replacing bad pixels with the median in a surrounding box.

``--bad_bits`` (string, default="DO_NOT_USE")
  The DQ bit values from the input image DQ arrays that should be considered bad
  and replaced with the median in a surrounding box. For example, setting to
  ``"OUTLIER, SATURATED"`` (or equivalently ``"16, 2"`` or ``"18"``) will treat
  all pixels flagged as OUTLIER or SATURATED as bad, while setting to ``""`` or
  ``None`` will treat all pixels as good and omit any bad pixel replacement.

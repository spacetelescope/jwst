Step Arguments
==============

The ``charge_migration`` step has the following optional arguments:

``--signal_threshold`` (float, default=-1)
  A value in units of ADU for the science value above which
  a group's DQ will be flagged as CHARGELOSS and DO_NOT_USE.
  If the value is less than zero, no flagging will be applied.

``--correction_threshold`` (float, default=-1)
  Threshold value in ADU for attempting to fit and correct for
  the brighter-fatter effect (BFE). If the value is less than zero,
  no correction will be applied.

``--mask_pad`` (integer, default=2)
    Number of pixels to fit on either side of a pixel above the ``correction_threshold``.

``--snr_limit`` (float, default=3.0)
    Minimum signal-to-noise ratio for a fit correction to consider the fit successful.
    Specified as a factor of the read noise value. Set to zero to skip the SNR check
    for the fits. Ignored if ``correction_threshold`` is negative.

``--rms_limit`` (float, default=2.2)
    Minimum RMS improvement from a fit correction to consider the fit successful.
    Specified as a reduction factor for the initial RMS. Set to zero to skip the RMS check
    for the fits. Ignored if ``correction_threshold`` is negative.

``--save_fit_info`` (boolean, default=False)
    If True, save a table of information about fit corrections.  The table is stored as
    an `~astropy.table.Table` in ASDF format.  Columns are:

    * ``x``: x index for the pixel, zero-indexed
    * ``y``: y index for the pixel, zero-indexed
    * ``fit_success``: True if the fit was successful and the correction was applied to the pixel, False otherwise
    * ``sigma_orig``: Standard deviation for the group differences prior to correction
    * ``sigma_new``: Standard deviation for the group differences after correction
    * ``diff``: Group difference data input to the fit
    * ``model``: Modeled difference values

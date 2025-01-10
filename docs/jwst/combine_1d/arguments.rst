Step Arguments
==============

The ``combine_1d`` step has two step-specific arguments:

``--exptime_key``
  This is a case-insensitive string that identifies the metadata element
  (or FITS keyword) for the weight to apply to the input data.  The default
  is "integration_time".  If the string is "effinttm" or starts with
  "integration", the integration time (FITS keyword EFFINTTM) is used
  as the weight.  If the string is "effexptm" or starts with "exposure",
  the exposure time (FITS keyword EFFEXPTM) is used as the weight.  If
  the string is "unit_weight" or "unit weight", the same weight (1) will
  be used for all input spectra.  If the string is anything else, a warning
  will be logged and unit weight will be used.

``--sigma_clip``
  Optional factor for sigma clipping outliers when combining spectra. If
  a floating point value is provided for ``sigma_clip``, this value will be
  used to set an outlier threshold for any pixels in the input spectra that 
  deviate from the median and median absolute deviation of the inputs.
  Defaults to None (such that no clipping is performed).

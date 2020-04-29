Step Arguments
==============

The ``flat_field`` step has one step-specific argument, and it is only
relevant for NIRSpec data.

``--save_interpolated_flat``
  is a boolean that indicates whether to save to a file the NIRSpec
  flat field that was constructed on-the-fly by the step.
  The default is False (do not save).

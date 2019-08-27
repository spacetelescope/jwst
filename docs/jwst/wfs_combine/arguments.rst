Step Arguments
==============

The ``wfs_combine`` step has one step-specific argument::

  --do_refine  boolean  default=False

If set to ``True``, the nominal image offsets computed from the WCS information are
refined using image cross-correlation. See the algorithm description section for details.

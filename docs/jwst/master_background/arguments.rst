.. _msb_step_args:

Step Arguments
==============
The master background subtraction step uses the following optional arguments.

``--user_background``
  The file name of a user-supplied 1-D master background spectrum. Must be in the form
  of a standard :ref:`x1d <x1d>` product containing a single 'EXTRACT1D' extension.
  When a user background spectrum is supplied, it is used for the subtraction instead of
  a computed master background, and the name of the user-supplied file is recorded in the
  MSTRBKGD keyword in the output product(s).
  Defaults to ``None``.

``--save_background``
  A boolean indicating whether the computed 1-D master background spectrum should be saved
  to a file. The file name uses a product type suffix of "masterbg".
  If a user-supplied background is specified, this argument is ignored.
  Defaults to ``False``.

``--force_subtract``
  A boolean indicating whether or not to override the step's built-in logic for determining
  if the step should be applied. By default, the step will be skipped if the
  :ref:`calwebb_spec2 <calwebb_spec2>` :ref:`background <background_step>` step has
  already been applied. If ``--force_subtract = True``, the master background will be
  applied.

``--output_use_model``
  A boolean indicating whether to use the "filename" meta attribute in the data model to
  determine the name of the output file created by the step. Defaults to ``True``.

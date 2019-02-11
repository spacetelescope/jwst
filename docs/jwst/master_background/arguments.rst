.. _msb_step_args:

Step Arguments
==============
The master spectroscopic background subtraction step has two optional arguments.

``--user_background``
  The file name of a user-supplied 1-D master background spectrum. Must be in the form
  of a standard :ref:`x1d <x1d>` product containing a single 'EXTRACT1D' extension.
  When a user background spectrum is supplied, it is used for the subtraction instead of
  a computed master background, and the name of the user-supplied file is recorded in the
  MSTRBKGD keyword in the output product(s).
  Defaults to ``None``.

``--save_background``
  A boolean indicating whether the computed master background spectrum should be saved
  to a file. If a user-supplied background is specified, this argument is ignored.
  Defaults to ``False``.


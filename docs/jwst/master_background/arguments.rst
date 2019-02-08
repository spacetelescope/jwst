.. _msb_step_args:

Step Arguments
==============
The master spectroscopic background subtraction step has two optional arguments.

``--user_background``
  The file name of a user-supplied 1-D master background spectrum. Must be in the form
  of a standard :ref:`x1d <x1d>` product containing a single 'EXTRACT1D' extension. When set,
  the creation of a master background spectrum by the step is skipped. The name of the
  supplied file is recorded in the MSTRBKGD keyword in the output product(s).
  Defaults to ``None``.

``--save_background``
  A boolean indicating whether the computed master background spectrum should be saved
  to a file. Cannot be used when ``user_background`` is specified.
  Defaults to ``False``.


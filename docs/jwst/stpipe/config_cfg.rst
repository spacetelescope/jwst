.. _config_cfg_files:

Configuration (CFG) Files
=========================

.. note::

   The ``cfg`` format can still be used but is deprecated in favor of
   :ref:`config_asdf_files`. Please convert any processes that use ``cfg`` files
   to the ``ASDF`` format.

The ``cfg`` format for configuration files uses the well-known ini-file format.

For a given step, the step's cfg file may specify parameters and their default
values, including parameters that are typically not changed between runs.
Parameters that are usually reset for each run are not included in the cfg file,
but instead specified on the command line. An example of a cfg file for the
jump detection step is::

    name = "jump"
    class = "jwst.jump.JumpStep"
    rejection_threshold = 4.0

You can list all of the parameters for this step using::

    strun jump.cfg -h

which gives the usage, the positional arguments, and the optional arguments.

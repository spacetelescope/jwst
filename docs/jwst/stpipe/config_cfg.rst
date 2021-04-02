.. _config_cfg_files:

Configuration (CFG) Files
=========================

.. note::

   The ``cfg`` format can still be used but is deprecated in favor of
   :ref:`config_asdf_files`. Please convert any processes that use ``cfg`` files
   to the ``ASDF`` format. Note also that all ``cfg`` files that are currently
   being delivered in the package and retrieved using ``collect_pipeline_cfgs``
   set no parameters; files are empty. All steps query CRDS parameter references
   for any data-dependent parameter settings, or use coded defaults.

The ``cfg`` format for configuration files uses the well-known ini-file format.

You can use the ``collect_pipeline_cfgs`` task to get copies of all the cfg
files currently in use by the jwst pipeline software. The task takes a single
argument, which is the name of the directory to which you want the cfg files
copied. Use '.' to specify the current working directory, e.g.
::

 $ collect_pipeline_cfgs .

Each step and pipeline has their own cfg file, which are used to specify
relevant parameter values. For each step in a pipeline, the pipeline cfg file
specifies either the step's arguments or the cfg file containing the step's
arguments.

For a given step, the step's cfg file specifies parameters and their default
values; it includes parameters that are typically not changed between runs.
Parameters that are usually reset for each run are not included in the cfg file,
but instead specified on the command line. An example of a cfg file for the
jump detection step is:
::

    name = "jump"
    class = "jwst.jump.JumpStep"
    rejection_threshold = 4.0

You can list all of the parameters for this step using:
::

 $ strun jump.cfg -h

which gives the usage, the positional arguments, and the optional arguments.

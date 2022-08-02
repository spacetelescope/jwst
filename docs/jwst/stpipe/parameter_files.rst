.. _parameter_files:

Parameter Files
===============

Parameter files can be used to specify parameter values when running a
pipeline or individual steps. For JWST, parameter files are retrieved from
CRDS, just as with other reference files. If there is no match between a step,
the input data, and CRDS, the coded defaults are used. These values can be
overridden either by the command line options and/or a
local parameter file. See :ref:`Parameter Precedence` for a full description of
how a parameter gets its final value.

.. note::

   Retrieval of ``Step`` parameters from CRDS can be completely disabled by
   using the ``--disable-crds-steppars`` command-line switch, or setting the
   environment variable ``STPIPE_DISABLE_CRDS_STEPPARS`` to ``true``.

A parameter file should be used when there are parameters a user wishes to
change from the default/CRDS version for a custom run of the step. To create a
parameter file add ``--save-parameters <filename.asdf>`` to the command:
::

$ strun <step.class> <required-input-files> --save-parameters <filename.asdf>

For example, to save the parameters used for a run of the ``calwebb_image2`` pipeline, use:
::

$ strun calwebb_image2 jw82500001003_02101_00001_NRCALONG_rate.fits --save-parameters my_image2.asdf

Once saved, the file can be edited, removing parameters that should be left
at their default/CRDS values, and setting the remaining parameters to the
desired values. Once modified, the new parameter file can be used:
::

$ strun my_image2.asdf jw82500001003_02101_00001_NRCALONG_rate.fits

Note that the parameter values will reflect whatever was set on the
command-line, through a specified local parameter file, and what was
retrieved from CRDS. In short, the values will be those actually used in the
running of the step.

For more information about and editing of parameter files, see
:ref:`config_asdf_files`. Note that the older :ref:`config_cfg_files` format is
still an option, understanding that this format will be deprecated.


More information on parameter files can be found in the ``stpipe`` User's
Guide at :ref:`stpipe-user-steps`.

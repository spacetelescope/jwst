.. _parameters:

==========
Parameters
==========

Parameters, which exist at both the step level and the global pipeline level,
can be set to change various aspects of processing. Parameters can be set in
a parameter file, on the command line, or passed in as an argument when running in Python.
Note that because there are multiple ways to set parameters, there is
a hierarchy involved - overrides set on a pipeline or step object will take precedence
over values in a parameter file. See :ref:`Parameter Precedence` for a full description of
how a parameter gets its final value.

If there is need to re-use a set of parameters often, parameters can be stored
in **parameter files**. See :ref:`parameter_files` for more information.

To see what parameters are available for any given
pipeline or step, use the ``-h`` option on ``strun``. Some examples are:
::

   $ strun calwebb_detector1 -h
   $ strun jwst.dq_init.DQInitStep -h



Universal Parameters
====================

The set of parameters that are common to all pipelines and steps are referred to
as **universal parameters** and are described below. When these parameters are
set at the pipeline level, they will apply to all steps within that pipeline, unless
explicitly overridden for a specific step.

.. _intro_output_directory:

Output Directory
----------------

By default, all pipeline and step outputs will drop into the current
working directory, i.e., the directory in which the process is
running. To change this, use the ``output_dir`` parameter. See .. _python_output_directory:
for instructions when running in Python, and .. _cli_output_directory: for instructions
using the command line interface.

.. _intro_output_file:

Output File
-----------

When running a pipeline, the ``stpipe`` infrastructure automatically passes the
output data model from one step to the input of the next step, without
saving any intermediate results to disk. If you want to save the results from
individual steps, you have two options:

  - Specify ``save_results``.
    This option will save the results of the step, using a filename
    created by the step.

  - Specify a file name using ``output_file <basename>``.
    This option will save the step results using the name specified.

To do this using the Python pipeline interface, see .. _python_output_file:. To do
this when using the command line interface, see .. _cli_output_file:.


Override Reference File
-----------------------

For any step that uses a calibration reference file you always have the
option to override the automatic selection of a reference file from CRDS and
specify your own file to use. Parameters for this are of the form
``--override_<ref_type>``, where ``ref_type`` is the name of the reference file
type, such as ``mask``, ``dark``, ``gain``, or ``linearity``. When in doubt as to
the correct name, just use the ``-h`` argument to ``strun`` to show you the list
of available override parameters.

To override the use of the default linearity file selection, for example,
you would use:
::

  $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.linearity.override_linearity='my_lin.fits'

Skip
----

Another parameter available to all steps in a pipeline is ``skip``. If
``skip=True`` is set for any step, that step will be skipped, with the output of
the previous step being automatically passed directly to the input of the step
following the one that was skipped. For example, if you want to skip the
linearity correction step, one can specify the ``skip`` parameter for the
``strun`` command:
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.linearity.skip=True

Alternatively, if using a :ref:`parameter file<parameter_files>`, edit the
file to add the following snippet:
::

  steps:
  - class: jwst.linearity.linearity_step.LinearityStep
    parameters:
      skip: true

Pipeline/Step Parameters
========================

All pipelines and steps have **parameters** that can be set to change various
aspects of how they execute. To see what parameters are available for any given
pipeline or step, use the ``-h`` option on ``strun``. Some examples are:
::

   $ strun calwebb_detector1 -h
   $ strun jwst.dq_init.DQInitStep -h

To set a parameter, simply specify it on the command line. For example, to have
:ref:`calwebb_detector1 <calwebb_detector1>` save the calibrated ramp files, the
``strun`` command would be as follows:
::

   $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits --save_calibrated_ramp=true

To specify parameter values for an individual step when running a pipeline
use the syntax ``--steps.<step_name>.<parameter>=value``.
For example, to override the default selection of a dark current reference
file from CRDS when running a pipeline:
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.dark_current.override_dark='my_dark.fits'

If there is need to re-use a set of parameters often, parameters can be stored
in **parameter files**. See :ref:`parameter_files` for more information.

Pipeline/Step Parameters
========================



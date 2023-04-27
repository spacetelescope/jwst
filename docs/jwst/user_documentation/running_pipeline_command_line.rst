.. _run_from_strun:

=============================================================
Running the JWST pipeline: Command Line Interface (``strun``)
=============================================================

.. note::

   For seasoned users who are familiar with using ``collect_pipeline_cfgs`` and
   running pipelines by the default configuration (CFG) files, please note that
   this functionality has been deprecated. Please read :ref:`CFG Usage
   Deprecation Notice<cfg_usage_deprecation_notice>`.

Individual steps and pipelines (consisting of a series of steps) can be run
and configured from the command line using the ``strun`` command.
``strun`` is one of two options for running the pipeline. See
:ref:`here <run_from_python>` for an overview of the alternative Python
interface.


CRDS Environment Variables
--------------------------

The CRDS environment variables need to be defined *before* running a pipeline
or step with ``strun`` to allow the pipeline to access reference and parameter
files. See :ref:`crds` for more information.


Overview of Running the Pipeline with ``strun``
-----------------------------------------------

The first argument to ``strun`` must be one of either a pipeline name, Python
class of the step or pipeline to be run, or the name of a parameter file for the
desired step or pipeline (see :ref:`parameter_files`). The second argument to
``strun`` is the name of the input data file to be processed.

::

    $ strun <pipeline_name, class_name, or parameter_file> <input_file>


Pipeline classes also have a **pipeline name**, or **alias**, that can be used
instead of the full class specification. For example, ``jwst.pipeline.Detector1Pipeline``
has the alias ``calwebb_detector1`` and can be run as
::

  $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits

A full list of pipeline aliases can be found in :ref:`Pipeline Stages <pipelines>`.

.. _exit_status:

Exit Status
-----------
``strun`` produces the following exit status codes:

- 0: Successful completion of the step/pipeline
- 1: General error occurred
- 64: No science data found

The "No science data found" condition is returned by the ``assign_wcs`` step of
the ``calwebb_spec2`` pipeline when, after successfully determining the WCS
solution for a file, the WCS indicates that no science data will be found. This
condition most often occurs with NIRSpec's Multi-object Spectroscopy (MOS) mode:
There are certain optical and MSA configurations in which dispersion will not
cross one or the other of NIRSpec's detectors.

.. _configuring_pipeline_strun:

Configuring a Pipeline/Step with ``strun``
==========================================

By default, pipeline parameters and reference files are chosen by CRDS based on 
instrument, observing mode, date, etc. If set to the most current :ref:`crds_context`,
these represent the 'best' set of parameters and reference files for the pipeline
as determined by the JWST instrument teams.

A Pipeline/Step can be configured for custom processing. Pipeline-level and
step-level parameters can be changed, output file behavior can be set, references
files can be overridden, and pipeline steps can be skipped if desired. This
section will be a general overview on how to configure the pipeline when running with `strun`,
and the following sections will elaborate on each of these possible customizations
and demonstrate usage.

**When running command line with ``strun``, there are two ways two configure a Pipeline/Step.**

1. By passing in arguments to a pipeline/step on the command line
2. By using a :ref:`parameter file<parameter_files>` and passing this in as an argument on the command line


A combination of arguments and custom parameter files can be used
for configuration, but keep in mind the hierarchy of :ref:`parameter precedence<Parameter Precedence>`
to keep track of which value will get used if set in multiple locations.


.. _setting_parameters_strun:

Setting Step Parameters on a Pipeline or Individual Step
--------------------------------------------------------

All pipelines and steps have parameters that can be set to change various aspects
of how they execute (e.g switching on and off certain options in a step,
setting thresholds). By default, the values of these parameters are set in
the CRDS-chosen parameter file, but they can be overridden when running the
pipeline with ``strun``. As mentioned, this can either be done by passing in command line
arguments or by passing in a custom parameter file - both methods will be described in this
section.

**Using Command Line Arguments**

When running a pipeline, step-level parameters can be changed by passing in a command
line argument to that step. For example, to change the ``rejection_threshold`` parameter of
the jump detection step when running the full Detector1Pipeline:

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.jump.rejection_threshold=12.0


When running a standalone step, command line arguments do not need to be nested within
``steps``. For example, to change the parameter ``rejection_threshold`` for the jump detection
step when running the step individually:

::

    $ strun jump jw00017001001_01101_00001_nrca1_uncal.fits --rejection_threshold=12.0


**Using a Parameter File**

Alternatively, if using a :ref:`parameter file<parameter_files>`, edit the
file to add the following snippet (in this example, to a file named
'my_config_file.asdf' in the current working directory):

::

  steps:
  - class: jwst.jump.jump_step.JumpStep
    name: jump
    parameters:
      rejection_threshold : 12

And pass in the modified file to ``strun``:

::

	$ strun my_config_file.asdf jw00017001001_01101_00001_nrca1_uncal.fits

.. _override_ref_strun:

Overriding Reference Files
--------------------------
By default, when the pipeline or step is run, CRDS will determine the best set of 
reference files based on file metadata and the current CRDS mapping (also known
as 'context'). It is possible to override these files and use a custom reference file,
or one not chosen by CRDS.

**Using Command Line Arguments**

For any step that uses a calibration reference file you always have the
option to override the automatic selection of a reference file from CRDS and
specify your own file to use. Parameters for this are of the form
``--override_<ref_type>``, where ``ref_type`` is the name of the reference file
type, such as ``mask``, ``dark``, ``gain``, or ``linearity``. When in doubt as to
the correct name, just use the ``-h`` argument to ``strun`` to show you the list
of available override parameters.

To override the use of the default linearity reference file selection with a custom
file in the current working directory called `my_lin.fits`, for example,
you would do:
::

  $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.linearity.override_linearity='my_lin.fits'

Or, if running the step individually, to override the reference file:
::

  $ strun linearity jw00017001001_01101_00001_nrca1_uncal.fits
          --override_linearity='my_lin.fits'


**Using a Parameter File**

If  using a :ref:`parameter file<parameter_files>` for configuration, to override
a reference edit the file to add the following snippet (in this example, to a
file named 'my_config_file.asdf' in the current working directory):
::

  steps:
  - class: jwst.saturation.saturation_step.SaturationStep
    name: saturation
    parameters:
      override_saturation: '/path/to/new_saturation_ref_file.fits'


And pass in the modified file to ``strun``:

::

	$ strun my_config_file.asdf jw00017001001_01101_00001_nrca1_uncal.fits

To use an entire set of past reference files from a previous CRDS mapping,
see :ref:`here<crds_context>`.

.. _skip_step_strun:

Skipping a Pipeline Step
------------------------

.. note::

   Some steps in a pipeline expect certain previous steps to have been run
   beforehand, and therefore won't run if that expected previous correction
   has not been applied. Proceed with caution when skipping steps.

When running a pipeline with `strun`, one or several steps within that pipeline
can be skipped.

**Using Command Line Arguments**

Every step in a pipeline has a ``skip`` parameter that when set to true, will entirely
skip that step. For example, to skip the saturation step in the Detector1Pipeline:
::

	$ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
	    --steps.saturation.skip=True

**Using a Parameter File**

The equivalent to the above example can be done by adding the following snippet
to your parameter file (in this example, to a file named 'my_config_file.asdf'
in the current working directory):

::

	steps:
	- class: jwst.saturation.saturation_step.SaturationStep
	  parameters:
	    skip: true

And pass in the modified file to the ``config_file`` argument:

::

	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                 config_file='my_config_file.asdf')

.. _strun_outputs:

Controlling Output File Behavior with ``strun``
===============================================

By default, when running the pipeline with ``strun``, the final outputs of a pipeline
(or final outputs when running an individual step) will be written out to a file
in the current working directory. The base name of these final output files is
derived from the input file name, by default. Additionally, no intermediate step
results will be saved. This behavior can be modified to change output file names, 
locations, and specify that intermediate results from a step in a pipeline should
be written out to a file.

.. _strun_intermediate_outputs: 

Saving Intermediate Pipeline Results to a File
----------------------------------------------

The ``stpipe`` infrastructure automatically passes the output data model from
one step to the input of the next step, without saving any intermediate results
to disk.  If you want to save the results from individual steps, you have two options:

  - Specify ``save_results`` on an individual step within the pipeline.
    This option will save the results of the step, using a filename
    created by the step.

  - Specify a file name using ``output_file <basename>`` for an individual step.
    This option indicated that results should be saved, and to use the name specified.

For example, to save the result from the dark current step of ``Detector1Pipeline``
(using the :ref:`alias <pipelines>` name ``calwebb_detector1``):

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.dark_current.save_results=true

This will create the file ``jw00017001001_01101_00001_dark_current.fits`` in the 
current working directory.



Setting Output File Name
------------------------

As demonstrated in the :ref:`section above <strun_intermediate_outputs>`, the ``output_file``
parameter is used to specify the desired name for output files. When done at the
step-level as shown in those examples, the intermediate output files from steps
within a pipeline are saved with the specified name.

You can also specify a particular file name for saving the end result of
the entire pipeline using the ``--output_file`` parameter:

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --output_file='stage1_processed'

In this situation, using the default configuration, three files are created:

  - ``stage1_processed_trapsfilled.fits``
  - ``stage1_processed_rate.fits``
  - ``stage1_processed_rateints.fits``

When running a standalone step, setting ``--output_file`` at the top-level
will determine the name of the final output product for that step, overriding
the default based on input name:

::

    $ strun linearity jw00017001001_01101_00001_nrca1_uncal.fits
        --output_file='intermediate_linearity'


Similarly, to save the result from a step within a pipeline (for example,
the dark current step of ``calwebb_detector1``) with a different file name:

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.dark_current.output_file='intermediate_result'

A file, ``intermediate_result_dark_current.fits``, will then be created. Note
that the name of the step will be appended as the file name suffix


Setting Output File Directory
-----------------------------

To change the output directory of the final pipeline products from the default of the
current working directory, use the ``output_dir`` option.

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.dark_current.output_dir='calibrated'

When this is run, all three final output products of `Detector1Pipeline` will
be saved within the subdirectory ``calibrated``.


Setting ``output_dir`` at the step-level indicates that the step's result should
be saved (so, also setting ``save_results`` is redundant), and that the files
should be saved in the directory specified instead of the current working directory.
For example, to save the intermediate results of ``DarkCurrentStep`` when running
``Detector1Pipeline`` in a subdirectory ``/calibrated``:

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.dark_current.output_dir='calibrated'


Similarly, when `output_dir` is set on an individual step class, this will indicate
that the result from that step should be saved to the specified directory:

::

    $ strun dark_current jw00017001001_01101_00001_nrca1_uncal.fits --output_dir='calibrated'


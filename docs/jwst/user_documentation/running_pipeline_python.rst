.. _run_from_python:

===========================================
Running the JWST pipeline: Python Interface
===========================================

.. note:: 
	The use of the ``run`` method to run a pipeline or step is not
	recommended. By default, using the ``pipeline.run()`` method defaults to
	pipeline and step-level coded defaults, ignoring parameter files,
	unless explicitly overridden. Please see :ref:`python_run_vs_call` for more details.

The Python interface is one of two options for running the pipeline.
See :ref:`here <run_from_strun>` for an overview of the alternative command line
interface. 


Overview of Running the Pipeline in Python
==========================================

When using the Python interface to the JWST pipeline, each ``pipeline`` and
``step`` is available as a module that can be imported into your Python session,
configured (either directly as arguments/attributes or with a
:ref:`parameter file<parameter_files>`), and used to process input data. The
following section will describe the necessary steps to run a pipeline or step in
Python.

CRDS Environment Variables
--------------------------

The CRDS environment variables need to be defined *before* importing anything
from ``jwst`` or ``crds`` to allow access to reference and parameter files.
These environment variables can be set in the shell, or
in a Python session by using `os.environ`. See :ref:`python_crds_variables`
for more information.

.. _importing_from_python:

Importing and Running Pipelines and Steps in Python
---------------------------------------------------

All full pipeline stages can be imported by name from the `pipeline` module::

	from jwst.pipeline import Image3Pipeline
	from jwst.pipeline import Spec2Pipeline

Individual pipeline steps can be imported by name from their respective module
in ``jwst``::

	from jwst.saturation import SaturationStep
	from jwst.ramp_fitting import RampFitStep

Details of all the available pipeline modules and their names can be found at
:ref:`pipeline-modules`.

Once :ref:`imported <importing_from_python>`, you can execute a pipeline or a
step from within Python by using the .call() method of the class. The input can
be either a string path to a file on disk or an open ``DataModel`` object. Note
that the .run() class method is also available for use, but is discouraged and
should be used only with caution (see :ref:`here <python_run_vs_call>` for
more information).


**Example: Running a Pipeline or Step with Default Parameters and Reference Files**
::

	# running a full pipeline stage, input is path to file
	from jwst.pipeline import Detector1Pipeline
	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits')

	# running a single pipeline step, input is datamodel object
	from jwst.linearity import LinearityStep
	import stdatamodels.jwst.datamodels as dm
	input_model = dm.open('jw00001001001_01101_00001_mirimage_uncal.fits')
	result = LinearityStep.call(input_model)


In the examples above, the returned value ``result``, is a ``Datamodel``
containing the corrected data - no files are written out, by default.
See :ref:`python_outputs` for information on how to control the generation of
output files.

Additionally in both examples above, there are no arguments other than the input
data being passed in to the ``call`` method, so the appropriate parameter files
and reference files are chosen by CRDS based on the :ref:`current context <crds_context>`.
The :ref:`following section <configuring_pipeline_python>` will
show how to configure the pipeline to override these defaults.

.. _configuring_pipeline_python:

Configuring a Pipeline/Step in Python
=====================================

By default when using the ``.call()`` method to run a pipeline/step, pipeline/step
parameters and reference files are chosen by CRDS based on instrument,
observing mode, date, etc. If set to the most current :ref:`context <crds_context>`,
these represent the 'best' set of parameters and reference files for the dataset
passed in, as determined by the JWST instrument teams.

To override parameter and reference file defaults, a pipeline/step can be
configured for custom processing. Pipeline-level and step-level parameters can be
changed, output file behavior can be set, reference files can be overridden,
and pipeline steps can be skipped if desired. This section will be a general
overview on how to configure the pipeline when running in Python, and the
following sections will elaborate on each of these options.

**When running in Python, there are two ways two configure a Pipeline/Step.**

1. By passing in keyword arguments to a pipeline/step's ``call`` method
2. By using a :ref:`parameter file<parameter_files>`

A combination of both keyword arguments and custom parameter files can be used
for configuration, but keep in mind the hierarchy of
:ref:`parameter precedence <Parameter Precedence>` to keep track of which value
will get used if set in multiple locations.


**Example: Configuring a pipeline/step with keyword arguments**

::

	# configuring a pipeline and the steps within the pipeline with keyword arguments
	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                save_results=False,
	                                steps={'jump': {'rejection_threshold': 12.0, 'save_results':True}})
	# configuring a pipeline step with keyword arguments
	result = JumpStep.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                       save_results=True, 'rejection_threshold'=12.0)

Both examples above show how to configure the jump detection step with the same
settings - the ``rejection_threshold`` set to 12.0, and ``save_results`` set to True to indicate
the result from the step should be written to an output file.

The first example shows when the jump step is run inside a pipeline - because a
pipeline consists of many steps, parameters for a substep are specified within
the ``steps`` argument, a nested dictionary keyed by each substep and again by each
possible parameter for each substep. Pipeline-level arguments (in this case,
``save_results``) are passed in individually as keyword arguments. Note that in this
example, the 'save_results' argument within ``steps`` will override the
pipeline-level 'save_results' argument.

The second example shows the same configuration to the jump step, but this time
when the step is run standalone. Here, there is no ``steps`` dictionary argument
and all arguments can be passed to the step directly since it is now at the step level.

**Example: Configuring a pipeline/step with a parameter file**

To use a custom parameter file, set the ``config_file`` parameter:

::

	# passing a custom parameter file to a pipeline
	result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits",\
	                                config_file='calwebb_detector1.asdf')

Again, note the :ref:`parameter precedence<Parameter Precedence>` rules. If an
override parameter file passed in does not contain the full set of required
parameters for a step, the others will be obtained according to those rules and
may grab values from the CRDS-chosen parameter file as well. If a custom
parameter file is passed in to ``config_file`` AND an argument is passed directly
to the pipeline/step class, the value in the parameter file is overridden.

.. _setting_parameters_python:

Setting Step Parameters on a Pipeline or Individual Step
--------------------------------------------------------

All steps have parameters that can be set to change various aspects
of how they execute (e.g switching on and off certain options in a step,
setting thresholds). By default, the values of these parameters are set in
the CRDS-chosen parameter file (and if absent, defer to the coded defaults),
but they can be overridden if desired.

**As Arguments to a Pipeline / Step**

As discussed in :ref:`above<configuring_pipeline_python>`, when setting a
step-level parameter when that step is a substep of a pipeline, it must be passed
to the `steps` argument dictionary. For example, to change the ``rejection_threshold``
parameter of the jump detection step when running the full Detector1Pipeline:

::

	from jwst.pipeline import Detector1Pipeline
	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                 steps={'jump': {'rejection_threshold':12.0)}})

When running a single step, step-level parameters can be passed in directly as
keyword arguments. For example, to change the parameter
``rejection_threshold`` for the jump detection step when running the step individually:

::

	from jwst.jump import JumpStep
	result = JumpStep.call('jw00017001001_01101_00001_nrca1_uncal.fits', rejection_threshold=12.0)


**Using a Parameter File**

Alternatively, if using a :ref:`parameter file<parameter_files>`, edit the
file to add the following snippet (in this example, to a file named
`my_config_file.asdf` in the current working directory):

::

  steps:
  - class: jwst.jump.jump_step.JumpStep
    parameters:
      rejection_threshold : 12

And pass in the modified file to the ``config_file`` argument:

::

	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                 config_file='my_config_file.asdf')

Disabling all CRDS Step Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Retrieval of Step parameters from CRDS can be completely disabled by setting the
STPIPE_DISABLE_CRDS_STEPPARS environment variable to TRUE. This can be done in the shell, or
using the os.environ() command:

::

	os.environ["STPIPE_DISABLE_CRDS_STEPPARS"] = 'True'

.. _override_ref_python:

Overriding Reference Files
--------------------------

To override the reference file for a step selected by CRDS:

**As Arguments to a Pipeline / Step**

To override a reference file for a step within a pipeline, for example the ``saturation``
step in the Detector1Pipeline the ``override_saturation`` argument can be set in the
``saturation`` section of the ``steps`` argument.

::

	# To override a reference file of a step within a pipeline
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                 steps={"saturation" : {"override_saturation": '/path/to/new_saturation_ref_file.fits'}})

Multiple reference file overrides can be provided, for example:

::

	# To override a reference file for multiple steps within a pipeline
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	 								 steps={"saturation": {"override_saturation": '/path/to/new_saturation_ref_file.fits'},
	 								       {"jump" : {"override_jump": '/path/to/new_jump_ref_file.fits'}})



To override a reference file for a standalone step, "override\_<stepname>"
can be passed directly as a keyword argument to that step's `call` method: 

::

	# To override a reference file when running a standalone step
	 from jwst.linearity import SaturationStep
	 SaturationStep.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	 					 override_saturation='/path/to/new_saturation_ref_file.fits')


**Using a Parameter File**

If  using a :ref:`parameter file<parameter_files>` for configuration, to override
a reference edit the file to add the following snippet (in this example, to a file named
`my_config_file.asdf` in the current working directory):
::

  steps:
  - class: jwst.linearity.saturation_step.SaturationStep
    parameters:
      override_saturation: '/path/to/new_saturation_ref_file.fits'


And pass in the modified file to the ``config_file`` argument:

::

	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                 config_file='my_config_file.asdf')

To use an entire set of past reference files from a previous CRDS mapping,
see :ref:`here<crds_context>`.

.. _skip_step_python:

Skipping a Pipeline Step
------------------------

.. note::

   Some steps in a pipeline expect certain previous steps to have been run
   beforehand, and therefore won't run if that expected previous correction
   has not been applied. Proceed with caution when skipping steps.

When using the Python interface you wish to run a pipeline but skip one or some
of the steps contained in that pipeline, this can be done in two different ways:

**As Arguments to a Pipeline / Step**

Every step in a pipeline has a ``skip`` parameter that when set to true, will entirely
skip that step. For example, to skip the saturation step in the Detector1Pipeline:
::

	# To set a step parameter on a step within a pipeline
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', steps={"saturation": {"skip": True}})

**Using a Parameter File**

The equivalent to the above example can be done by adding the following snippet
to your parameter file (in this example, to a file named `my_config_file.asdf`
in the current working directory):

::

	steps:
	- class: jwst.linearity.linearity_step.LinearityStep
	  parameters:
	    skip: true

And pass in the modified file to the ``config_file`` argument:

::

	result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	                                 config_file='my_config_file.asdf')

.. _python_outputs:

Controlling Output File Behavior
================================

By default, when running in Python, all outputs are returned in-memory
(typically as a `Datamodel`) and no output files are written - even the final
result of a pipeline. To control this behavior, and other aspects of output file
generation like directory and file name, certain pipeline and step-level parameters
can be set. 

**Output file behavior can be modified with the ``save_results``, ``output_file``, and ``output_dir`` parameters**

Saving Final Pipeline Results
-----------------------------

The ``save_results`` parameter, when set at the pipeline-level, indicates
that the final pipeline output products should be saved to a file. The output
files will be in the current working directory, and be named based on the
input file name and the appropriate file suffix. Note that setting ``save_results`` 
at the pipeline-level will not save the results from each step, only the final
results from the full pipeline.

::

	# To save the final results from a pipeline to a file
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', save_results=True)


In this example, the following output files will be written in the current working directory:
	- ``jw00017001001_01101_00001_nrca1_trapsfilled.fits``
	- ``jw00017001001_01101_00001_nrca1_rate.fits``
	- ``jw00017001001_01101_00001_nrca1_rateints.fits``

**Changing Output File Name**

Setting ``output_file`` at the pipeline-level indicates that the pipeline's final result
should be saved (so, also setting ``save_results`` is redundant), and that a new file
base name should be used with the appropriate file suffix appended. For example,
to save the intermediate result from the saturation step when running
``Detector1Pipeline`` with a file name based on the string `detector_1_final` instead
of `jw00017001001_01101_00001_nrca1`:

::

	# saving the final results from running a pipeline with a custom output file basename
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', output_file='detector_1_final_result')

In this example, the following output files will be written in the current working directory

	- ``detector_1_final_result_trapsfilled.fits``
	- ``detector_1_final_result_rate.fits``
	- ``detector_1_final_result_rateints.fits``

**Changing Output File Directory**

When set at the pipeline level, the ``output_dir`` parameter will set where the final
pipeline output products are placed. The default is the current working directory.
For example, to save the results from Detector1Pipeline in a subdirectoy ``/calibrated``:

Setting ``output_dir`` at the pipeline-level indicates that the pipeline's final
results should be saved (so, also setting ``save_results`` is redundant), and that
the files should be saved in the directory specified instead of the current working
directory. For example, to save the intermediate results of ``Detector1Pipeline``
in a subdirectory ``/calibrated``:

::

	# to save the final result of a pipeline in a different specified output directory
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', output_dir='calibrated')


Saving Intermediate Step Results
--------------------------------

When the ``save_results`` parameter is set at the step-level (either within a pipeline,
or on a standalone step), it indicates that the result from that step should be
saved to a file. 

To save the intermediate output from a step within a pipeline:

::

	# To save the intermediate results of a step within a pipeline to a file
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	 								steps={"saturation": {"save_results": True}})


Similarly, when ``save_results`` is set on an individual step class, this will indicate
that the final result from that step should be saved.

::

	# To save the final results from SaturationStep when run standalone
	 from jwst.linearity import SaturationStep
	 SaturationStep.call('jw00017001001_01101_00001_nrca1_uncal.fits', save_results=True)


**Setting Output File Name** 

Setting ``output_file`` at the step-level indicates that the step's result should
be saved (so, also setting ``save_results`` is redundant), and that a new file
base name should be used with the appropriate file suffix appended. For example,
to save the intermediate result from the saturation step when running
``Detector1Pipeline`` with a file name based on the string `saturation_result` instead
of `jw00017001001_01101_00001_nrca1`:

::

	# To save the intermediate results of a step within a pipeline to a file with a custom name
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	 								steps={"saturation": {"output_file": 'saturation_result'})

Similarly, when `output_file` is set on an individual step class, this will indicate
that the result from that step should be saved to a file with that basename and the
appropriate suffix.

::

	# To save the final results from SaturationStep with a custom output file name when run standalone
	 from jwst.linearity import SaturationStep
	 SaturationStep.call('jw00017001001_01101_00001_nrca1_uncal.fits', output_file="saturation_result")

**Setting Output File Directory** 

Setting ``output_dir`` at the step-level indicates that the step's result should
be saved (so, also setting ``save_results`` is redundant), and that the files
should be saved in the directory specified instead of the current working directory.
For example, to save the intermediate results of ``DarkCurrentStep`` when running
``Detector1Pipeline`` in a subdirectory ``/calibrated``:

::

	# to save the intermediate step result in a different specified output directory
	 from jwst.pipeline import Detector1Pipeline
	 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
	 								 steps={'dark': {'output_dir': 'calibrated'}})


Similarly, when `output_dir` is set on an individual step class, this will indicate
that the result from that step should be saved to the specified directory:

::

	# to save the final result of a 
	 from jwst.pipeline import Detector1Pipeline
	 result = DarkCurrentStep.call('jw00017001001_01101_00001_nrca1_uncal.fits', output_dir='calibrated')


.. _python_run_vs_call:

Advanced use - `pipeline.run()` vs. `pipeline.call`
===================================================

Another option for running pipelines or steps is to use the `.run()` method
instead of the `.call()` method. **Using .run() is not recommended** and
considered advanced use, but it is an option to users.

The difference between ``.run()`` in ``.call()`` is in the retrieval and use
of parameters from CRDS parameter files. When the ``.call()`` method is invoked,
there is additional setup done to retrieve parameter and reference files and
reconcile those with any passed into the pipeline directly as an argument or in
a custom parameter file. When ``.call()`` is invoked, a new instance of the
pipeline/step class is created internally, and after parameters are determined the
``.run()`` method of that internal class is called. Because the actual processing
occurs on this new instance, attributes cannot be set directly on the original
pipeline/step class. They must be passed in as arguments to ``.call()`` or set
in the parameter file.

In contrast, when using the ``.run()`` method directly on a pipeline/step, the
additional logic to determine parameters and reference files is skipped. The pipeline
instance is being run as-is, and coded defaults for the pipeline and each intermediate step
will be used unless explicitly overridden individually. Because the instance created is
being run directly on the data, attributes can be set directly:

::

	from jwst.pipeline import Detector1Pipeline
	pipe = Detector1Pipeline()
	pipe.jump.rejection_threshold = 12
	pipe.ramp_fit.skip = True
	result = pipe.run('jw00017001001_01101_00001_nrca1_uncal.fits')

The ``pipe`` object created and the attributes set will persist and this object
can be reused within a Python session for processing data. Keep in mind that each
individual step parameter must be set when using this method, or else the coded
defaults will be used, which may be inappropriate for the dataset being processed.

See :ref:`call_examples` for more information.


.. _multiprocessing:

Multiprocessing
===============

Multiprocessing is supported to speed up certain computationally-intensive steps
in the pipeline, including the :ref:`jump detection <jump_step>`,
:ref:`ramp fitting <ramp_fitting_step>`, and
:ref:`WFSS contamination correction <wfss_contam_step>` steps. The examples below show how
multiprocessing can be enabled for these steps, as well as how to set up
multiprocessing to simultaneously run the entire pipeline on multiple observations.

Since the pipeline uses multiprocessing it is critical that any code using the pipeline adhere
to the guidelines described in the
`python multiprocessing documentation <https://docs.python.org/3/library/multiprocessing.html#multiprocessing-programming>`_.
The pipeline uses the `forkserver` start method internally and it is recommended that any
multiprocessing scripts that use the pipeline use the same start. As detailed in the
`python documentation <https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods>`_
this will require that code be "protected" with a ``if __name__ == '__main__':`` check as follows

::

    if __name__ = '__main__':
        [code used in multiprocessing]


There are a couple of scenarios to use multiprocessing with the pipeline:

1.  Multiprocessing within a pipeline step. At the moment, the steps that
support this are the :ref:`jump <jump_step>`,
:ref:`ramp_fitting <ramp_fitting_step>`,
and :ref:`wfss_contam <wfss_contam_step>` steps. To enable multiprocessing, the
optional parameter is `maximum_cores` for the ``jump``, ``ramp_fitting``, and
``wfss_contam`` steps. This parameter can be set to a numerical value given
as a string or it can be set to the words `quarter`, `half`, `all`,
or `none`, which is the default value.

The following example turns on a step's multiprocessing option. Notice only
one of the steps has multiprocessing turned on. We do not recommend
simultaneously enabling both steps to do multiprocessing, as this may likely
lead to running out of system memory.



::

    # SampleScript1

    import os, sys
    from jwst.pipeline import Detector1Pipeline

    uncal_file = 'jw0000_0000_uncal.fits'
    output_dir = '/my_project'

    def main():
        det1 = Detector1Pipeline()
        parameter_dict = {"ramp_fit": {"maximum_cores": 'all'}}
        det1.call(uncal_file, save_results=True, steps=parameter_dict, output_dir=output_dir)

    if __name__ = '__main__':
        sys.exit(main())


2. Calling the pipeline using multiprocessing. The following example uses this
option setting up a log file for each run of the pipeline and a text file with
the full traceback in case there is a crash. Notice that the ``import`` statement
of the pipeline is within the multiprocessing block that gets called by every
worker. This is to avoid a known memory leak.


::

    # SampleScript2

    import os, sys
    import traceback
    import configparser
    import multiprocessing
    from glob import glob

    def mk_stpipe_log_cfg(output_dir, log_name):
        """
        Create a configuration file with the name log_name, where
        the pipeline will write all output.
        Args:
            outpur_dir: str, path of the output directory
            log_name: str, name of the log to record screen output
        Returns:
            nothing
        """
        config = configparser.ConfigParser()
        config.add_section("*")
        config.set("*", "handler", "file:" + log_name)
        config.set("*", "level", "INFO")
        pipe_log_config = os.path.join(output_dir, "pipeline-log.cfg")
        config.write(open(pipe_log_config, "w"))

    def run_det1(uncal_file, output_dir):
        """
        Run the Detector1 pipeline on the given file.
        Args:
            uncal_file: str, name of uncalibrated file to run
            outpur_dir: str, path of the output directory
        Returns:
            nothing
        """
        log_name = os.path.basename(uncal_file).replace('.fits', '')
        mk_stpipe_log_cfg(output_dir, log_name+'.log')
        from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
        pipe_success = False
        try:
            det1 = Detector1Pipeline()
            det1.call(uncal_file, output_dir=output_dir, logcfg="pipeline-log.cfg", save_results=True)
            pipe_success = True
            print('\n * Pipeline finished for file: ', uncal_file, ' \n')
        except Exception:
            print('\n *** OH NO! The detector1 pipeline crashed! *** \n')
            pipe_crash_msg = traceback.print_exc()
        if not pipe_success:
            crashfile = open(log_name+'_pipecrash.txt', 'w')
            print('Printing file with full traceback')
            print(pipe_crash_msg, file=crashfile)

    def main():
        input_data_dir = '/my_project_dir'
        output_dir = input_data_dir

        # get the files to run
        files_to_run = glob(os.path.join(input_data_dir, '*_uncal.fits'))
        print('Will run the pipeline on {} files'.format(len(files_to_run)))

        # the output list should be the same length as the files to run
        outptd = [output_dir for _ in range(len(files_to_run))]

        # get the cores to use
        cores2use = int(os.cpu_count()/2)   # half of all available cores
        print('* Using ', cores2use, ' cores for multiprocessing.')

        # set the pool and run multiprocess
        with multiprocessing.Pool(cores2use) as pool:
            pool.starmap(run_det1, zip(files_to_run, outptd))

        print('\n * Finished multiprocessing! \n')

    if __name__ == '__main__':
        sys.exit(main())


.. warning::
    Although it is technically possible to call the pipeline with
    multiprocessing while also enabling this option in a step, we
    strongly recommend not to do this. This scenario would be the same as
    `SampleScript2` except with adding and calling the parameter dictionary
    `parameter_dict` in `SampleScript1`. However, Python will crash
    if both multiprocessing options are set to use all the cores or even
    less, because it is not permitted that a worker has children processes.
    We recommend not enabling step multiprocessing for parallel pipeline
    runs to avoid potentially running out of memory.

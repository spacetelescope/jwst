.. _call_examples:

Executing a pipeline or pipeline step via call()
================================================

The ``call`` method will create an instance and run a pipeline or pipeline step
in a single call.

::

 from jwst.pipeline import Detector1Pipeline
 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits')

 from jwst.linearity import LinearityStep
 result = LinearityStep.call('jw00001001001_01101_00001_mirimage_uncal.fits')


To set custom parameter values when using the ``call`` method, set the
parameters in the pipeline or parameter file and then supply the file using the
``config_file`` keyword: ::

 # Calling a pipeline
 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', config_file='calwebb_detector1.asdf')

 # Calling a step
 result = LinearityStep.call('jw00017001001_01101_00001_nrca1_uncal.fits', config_file='linearity.asdf')


When running a pipeline, parameter values can also be supplied in the call to ``call`` itself by using a nested dictionary of step and
parameter names:

::

 result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits", config_file='calwebb_detector1.asdf', steps={"jump":{"rejection_threshold": 200}})

When running a single step with ``call``, parameter values can be supplied more simply:

::

 result = JumpStep.call("jw00017001001_01101_00001_nrca1_uncal.fits", rejection_threshold=200)

Running steps and pipelines with ``call`` also allows for the specification of a logging
configuration file using the keyword ``logcfg``:

::

 result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits",
                                 config_file="calwebb_detector1.asdf",
                                 logcfg="my-logging-config.cfg")

Note that naming the logging configuration file "stpipe-log.cfg" will configure logging
without assignment of the ``logcfg`` keyword, as ``stpipe`` searches for this filename
in the local directory during execution. If the configuration should be used only when
specified, ensure your file is named something other than "stpipe-log.cfg"!

Where are the results?
----------------------

A fundamental difference between running steps and pipelines in Python as
opposed to from the command line using ``strun`` is whether files are created or
not. When using ``strun``, results are automatically saved to files because that
is the only way to access the results.

However, when running within a Python interpreter or script, the presumption is
that results will be used within the script. As such, results are not
automatically saved to files. It is left to the user to decide when to save.

If one wishes for results to be saved by a particular ``call``, use the
parameter ``save_results=True``::

 result = JumpStep.call("jw00017001001_01101_00001_nrca1_uncal.fits",
                        rejection_threshold=200, save_results=True)

If one wishes to specify a different file name, rather than a system-generated
one, set :ref:`output_file<intro_output_file>` and/or
:ref:`output_dir<intro_output_directory>`.

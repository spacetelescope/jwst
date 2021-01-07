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


To set custom parameter values when using the ``call`` method, set
the parameters in the pipeline or step configuration file and
then supply the file using the ``config_file`` keyword:
::

 # Calling a pipeline
 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', config_file='calwebb_detector1.cfg')

 # Calling a step
 result = LinearityStep.call('jw00017001001_01101_00001_nrca1_uncal.fits', config_file='linearity.cfg')


When running a pipeline, parameter values can also be supplied in the call to ``call`` itself by using a nested dictionary of step and
parameter names:

::

 result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits", config_file='calwebb_detector1.cfg', steps={"jump":{"rejection_threshold": 200}})

When running a single step with ``call``, parameter values can be supplied more simply:

::

 result = JumpStep.call("jw00017001001_01101_00001_nrca1_uncal.fits", rejection_threshold=200)

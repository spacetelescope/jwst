.. _run-examples

Calling a pipeline or pipeline step directly, or via run()
==========================================================

When calling a pipeline or step instance directly, or using the ``run`` method,
you can specify individual parameter values manually. In this case, configuration
files are not used (and will be ignored if provided).

::
 # Instantiate the class. Do not provide a configuration file.
 pipe = Detector1Pipeline()

 # Manually set any desired non-default parameter values
 pipe.refpix.skip = True
 pipe.jump.rejection_threshold = 5
 pipe.ramp_fit.override_gain = 'my_gain_file.fits'
 pipe.save_result = True
 pipe.output_dir = '/my/data/pipeline_outputs'

 # Run the pipeline
 result = pipe('jw00017001001_01101_00001_nrca1_uncal.fits')

 # Or, execute the pipeline using the run method
 result = pipe.run('jw00017001001_01101_00001_nrca1_uncal.fits')

 # Or to run a single step
 step = LinearityStep()
 step.override_linearity = 'my_linearity_coefficients.fits'
 step.save_results = True
 step.output_dir = '/my/data/linearized_data'

 # Execute by calling the instance directly
 result = step('jw00017001001_01101_00001_nrca1_superbias.fits')

 # Execute using the run method
 result = step.run('jw00017001001_01101_00001_nrca1_superbias.fits')

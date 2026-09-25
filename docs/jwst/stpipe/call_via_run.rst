.. _run_examples:

Executing a pipeline via run()
==============================

When calling a pipeline or step instance using the ``run`` method,
you can specify individual parameter values manually. In this case, parameter
files are not used. If you use ``run`` after instantiating with a parameter
file (as is done when using the :ref:`call <call_examples>` method), the
parameter file will be ignored::

    # Instantiate the class. Do not provide a parameter file.
    pipe = Detector1Pipeline()

    # Manually set any desired non-default parameter values
    pipe.refpix.skip = True
    pipe.jump.rejection_threshold = 5
    pipe.ramp_fit.override_gain = 'my_gain_file.fits'
    pipe.save_result = True
    pipe.output_dir = '/my/data/pipeline_outputs'

    # Execute the pipeline using the run method
    result = pipe.run('jw00017001001_01101_00001_nrca1_uncal.fits')

To run a single step::

    from jwst.jump import JumpStep

    # Instantiate the step
    step = JumpStep()

    # Set parameter values
    step.rejection_threshold = 5
    step.save_results = True
    step.output_dir = '/my/data/jump_data'

    # Execute using the run method
    result = step.run('jw00017001001_01101_00001_nrca1_linearity.fits')

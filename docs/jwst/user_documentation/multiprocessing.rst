.. _multiprocessing:

==========================================
Running the JWST pipeline: Multiprocessing
==========================================

Python multiprocessing could be used in the following ways,
which are *mutually exclusive*. You must not mix them together
because spawning multiple processes where each process itself
calls a step that also calls multiprocessing would result
in errors related to "daemon" or "recursion" and possibly
memory usage:

* :ref:`multiproc_within_pipeline_step`
* :ref:`multiproc_multiple-obs`

.. _multiproc_within_pipeline_step:

Multiprocessing within a pipeline step
======================================

This usage of multiprocessing is recommended when you want to speed up
the processing of a particular dataset running computationally-intensive steps:

* :ref:`jump <jump_step>` (jump detection)
* :ref:`ramp_fitting <ramp_fitting_step>`
* :ref:`wfss_contam <wfss_contam_step>` (WFSS contamination correction)

Unlike :ref:`multiproc_multiple-obs`, this usage is compatible with running
the pipeline within Jupyter Notebook/Lab.

To enable multiprocessing, the optional parameter is ``maximum_cores`` for
each of the step stated above. This parameter can be set to one of these options:

* a numerical value given as a string, e.g., ``'8'``
  (it is *not* recommended to set to ``'1'`` nor a number larger
  than the available number of cores)
* ``'quarter'``
* ``'half'``
* ``'all'`` (this is usually not recommended because it might
  tie up your CPU completely for the duration of the run)
* ``'none'`` (default)

The following example turns on a step's multiprocessing option.
Note that only one of the steps (``ramp_fit``) has multiprocessing
turned on.

.. note::
    For more details on how to adjust ``.call(...)`` inputs in the
    example below, please see :ref:`setting_parameters_python`.

::

    from jwst.pipeline import Detector1Pipeline

    uncal_file = 'jw0000_0000_uncal.fits'
    output_dir = '/path/to/my_project'
    parameter_dict = {
        "ramp_fit": {
            "maximum_cores": 'half'
        }
    }

    Detector1Pipeline.call(
        uncal_file,
        save_results=True,
        steps=parameter_dict,
        output_dir=output_dir
    )

Alternately, you can also run the equivalent call as above via
:ref:`strun <run_from_strun>`::

    strun calwebb_detector1 jw0000_0000_uncal.fits --steps.ramp_fit.save_results=true --steps.ramp_fit.maximum_cores=half --output_dir=/path/to/my_project

.. _multiproc_multiple-obs:

Multiprocessing on multiple observations
========================================

This usage of multiprocessing is to simultaneously run
the entire pipeline on multiple observations.
You must *not* use :ref:`multiproc_within_pipeline_step` if
you choose this option. It is recommended that you refer to the
:ref:`Python multiprocessing documentation <python:multiprocessing-programming>`
in order to follow its best practices; When in doubt, stick to
the pattern in the given example below.

The pipeline uses the ``spawn`` start method
(see :ref:`python:multiprocessing-start-methods`) internally and
it is recommended that any multiprocessing scripts that run
the pipeline use the same start.
As detailed in :ref:`python:multiprocessing-programming-spawn`,
this will require that code be "protected" with a
``if __name__ == '__main__':`` check as follows::

    if __name__ = '__main__':
        # code used in multiprocessing

Because the code has to be "protected" as explained above,
unlike :ref:`multiproc_within_pipeline_step`, you will not
be able to run this from within a Jupyter Notebook/Lab.

The following example runs the pipeline with multiprocessing via
a :py:meth:`multiprocessing.pool.Pool.starmap` method and
using :py:func:`zip` to pack the pipeline inputs.
The example also uses an option to set up a text file with the
full traceback for debugging, in case there is a crash.
Note that the ``import`` statement of the pipeline is within
the multiprocessing block that gets called by every worker (``run_det1``);
this is to avoid a known memory leak.

.. note::
    For more details on how to adjust ``.call(...)`` inputs in the
    example below, please see :ref:`setting_parameters_python`.

::

    # Save the code in a file named SampleScript2.py and then run it with
    #     python SampleScript2.py

    import os
    import sys
    import traceback
    import multiprocessing
    from glob import glob


    def run_det1(uncal_file, output_dir):
        """
        Run the Detector1 pipeline on the given file.

        Parameters
        ----------
        uncal_file : str
            Name of uncalibrated file to run.
        output_dir : str
            Path of the output directory.
        """
        # Local import of pipeline to avoid known memory leak
        from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

        log_name = os.path.basename(uncal_file).replace('.fits', '')
        pipe_success = False

        try:
            # Run the pipeline, turning off terminal logging messages
            Detector1Pipeline.call(
                uncal_file,
                output_dir=output_dir,
                save_results=True,
                configure_log=False
            )
            pipe_success = True
            print(f'\n * Pipeline finished for file: {uncal_file}\n')
        except Exception:
            print('\n *** OH NO! The detector1 pipeline crashed! *** \n')
            pipe_crash_msg = traceback.print_exc()
        if not pipe_success:
            with open(f'{log_name}_pipecrash.txt', 'w') as crashfile:
                print('Printing file with full traceback')
                print(pipe_crash_msg, file=crashfile)

    def main():
        input_data_dir = '/path/to/my_project_dir'
        output_dir = input_data_dir

        # get the files to run
        files_to_run = sorted(glob(os.path.join(input_data_dir, '*_uncal.fits')))
        n_files = len(files_to_run)
        print(f'Will run the pipeline on {n_files} files')

        # the output list should be the same length as the files to run
        outptd = [output_dir] * n_files

        # get the cores to use
        # (please adjust this according to your hardware, as there is no point
        # to take half of available cores if you only have 3 cores or less)
        n_cpu = os.cpu_count()
        cores2use = int(n_cpu / 2)  # half of all available cores
        print(f'* Using {cores2use}/{n_cpu} cores for multiprocessing.')

        # set the pool and run multiprocess
        with multiprocessing.Pool(cores2use) as pool:
            pool.starmap(run_det1, zip(files_to_run, outptd))

        print('\n * Finished multiprocessing! \n')

    if __name__ == '__main__':
        sys.exit(main())

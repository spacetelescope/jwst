.. _logging:

=======
Logging
=======

.. _cal_logs:

DataModel cal_logs
==================

Saved files that contain :ref:`jwst-data-models` will contain log messages
from the run that produced the file. These are stored in the ASDF portion
of the file and can be inspected by opening the file as a
:ref:`jwst-data-models`.

The ``cal_logs`` attribute contains log messages as lists of strings
organized by step or pipeline name. For example to see log messages from
:ref:`calwebb_detector1`:

::

        import stdatamodels.jwst.datamodels as dm
        model = dm.open("jw00001001001_01101_00001_mirimage_cal.fits")
        print(model.cal_logs.calwebb_detector1)

Files processed by a pipeline will contain all logs messages for that
run under the pipeline name (and not contain ``cal_logs`` for individual
steps that were part of the pipeline).

Log messages that contain sensitive information (user, hostname, paths,
IP addresses) are replaced with empty strings. Please see the console
logs for those messages.

Configuration
=============

By default, either the command line interface or the ``call`` function for
Step or Pipeline classes will log messages at the INFO level to the ``stderr``
stream.

The name of a file in which to save log information, as well as the desired
level of logging messages, can be specified in optional command line arguments
provided to ``strun``, or can be directly configured via the logging module
in Python code.

From the Command Line
---------------------

The available options for the command line are::

  --verbose, -v         Turn on all logging messages
  --log_level LOG_LEVEL
                        Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL). Ignored if 'verbose' is specified.
  --log_format LOG_FORMAT
                        A format string for the logger
  --log_file LOG_FILE   Full path to a file name to record log messages
  --log_stream LOG_STREAM
                        Log stream for terminal messages (stdout, stderr, or null).

For example::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --log_level=INFO --log_file=pipeline.log --log_stream=stdout

will log messages to the terminal at the INFO level in the ``stdout`` stream
and also record them to a file called "pipeline.log" in the current working directory.

To turn off all logging instead::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits --log_stream=null

CRDS messages may still display, since its logger is separately configured.

In Python Code
--------------

In a Python environment, more complex configuration is possible. For example,
to configure logging to print only the message at the INFO level and direct time-stamped
messages to a file at the DEBUG level::

    import logging
    import sys
    from jwst.pipeline import Detector1Pipeline

    # Set the log level to DEBUG to allow all messages
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)

    # Add a file handler for all messages, time-stamped
    handler = logging.FileHandler('pipeline.log')
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG)
    log.addHandler(handler)

    # Add a stream handler to just print messages at the INFO level
    handler = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    log.addHandler(handler)

    result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits")

To override the default logging configuration from Python code without directly
configuring the logger, turn it off in the ``call`` function with the ``configure_log`` option::

    Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits", configure_log=False)

This will print no log messages from the JWST pipeline code. As with the command line configuration,
CRDS messages may still display, since its logger is separately configured.

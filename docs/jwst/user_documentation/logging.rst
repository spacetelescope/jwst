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
:ref:`calwebb_detector1`::

    import stdatamodels.jwst.datamodels as dm
    model = dm.open("jw00001001001_01101_00001_mirimage_cal.fits")
    print(model.cal_logs.calwebb_detector1)

Files processed by a pipeline will contain all logs messages for that
run under the pipeline name (and will not contain separate ``cal_logs``
for individual steps that were part of the pipeline).

Log messages that contain sensitive information (user, hostname, paths,
IP addresses) are replaced with empty strings. Please see the console
logs for those messages.

Configuration
=============

By default, either the command line interface or the ``call`` method for
Step or Pipeline classes will log messages at the INFO level to the ``stderr``
stream.

The name of a file in which to save log information, as well as the desired
level of logging messages, can be specified in optional command line arguments
provided to ``strun``, or can be directly configured via the `logging` module
in Python code.

From the Command Line
---------------------

The available options for the command line are:

.. code-block:: text

  --verbose, -v         Turn on all logging messages
  --log-level LOG_LEVEL
                        Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL). Ignored if 'verbose' is specified.
  --log-file LOG_FILE   Full path to a file name to record log messages
  --log-stream LOG_STREAM
                        Log stream for terminal messages (stdout, stderr, or null).

For example:

.. code-block:: shell

    strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits --log-level=INFO --log-file=pipeline.log --log-stream=stdout

will log messages to the terminal at the INFO level in the ``stdout`` stream
and also record them to a file called "pipeline.log" in the current working directory.

To turn off all logging instead:

.. code-block:: shell

    strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits --log-stream=null

CRDS messages may still display, since its logger is separately configured.

In Python Code
--------------

Pipelines or steps may be run via the ``call`` method or the lower-level
``run`` method.  For more information on the difference, see :ref:`python_run_vs_call`.

Using ``call``
^^^^^^^^^^^^^^

In a Python environment, the default settings for the ``call`` method are the same as
for the command line, but more complex configuration is possible via the `logging` module
if the default configuration is disabled with the ``configure_log`` parameter.

For example, to configure the root logger to print only the log name and message at the INFO
level and direct time-stamped messages to a file at the DEBUG level::

    import logging
    import sys
    from jwst.pipeline import Detector1Pipeline

    # Make a stream handler to just print the log name and the
    # message at the INFO level
    stream_handler = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)

    # Make a file handler for all messages, time-stamped
    file_handler = logging.FileHandler('pipeline.log')
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    # Get the root logger and allow all messages through
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)

    # Attach the handlers
    log.addHandler(file_handler)
    log.addHandler(stream_handler)

    result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits", configure_log=False)


Using ``run``
^^^^^^^^^^^^^

Since it is a lower-level interface, the ``run`` method does not configure loggers
by default: no log messages will display when the ``run`` method is called unless the log
is directly configured.

To configure loggers for the ``run`` method, the above example for configuring the root logger
with the logging module will work exactly as it does for the ``call`` method.

For a minimum configuration that replicates the messages produced by the ``call`` method,
the root logger can be configured as follows::

    import logging
    log = logging.getLogger()
    log.setLevel('INFO')
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    log.addHandler(handler)

Then, for example::

    pipe = Detector1Pipeline()
    pipe.run("jw00017001001_01101_00001_nrca1_uncal.fits")

will produce similar log messages to the equivalent ``call`` method.

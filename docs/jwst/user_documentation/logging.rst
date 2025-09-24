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
Note that the log file is recreated every time the command is run, so subsequent
commands will overwrite earlier logs unless a new file name is provided.

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
for the command line. The following code, for example, will log messages at the INFO level to the
``stderr`` stream::

    from jwst.pipeline import Detector1Pipeline
    result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits")

The following code will keep the default logging to the ``stderr`` stream and also add a log file
called "pipeline.log" that will be appended to for every subsequent pipeline or step call::

    import logging
    from jwst.pipeline import Detector1Pipeline

    # Set up a file handler
    log = logging.getLogger()
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler = logging.FileHandler("pipeline.log")
    handler.setFormatter(formatter)
    log.addHandler(handler)

    result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits")

More complex configuration is possible via the `logging` module
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
    file_handler = logging.FileHandler("pipeline.log")
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


Migration guide for ``logcfg``
------------------------------

Prior to ``jwst`` version 1.20.0, the primary method for log configuration was
via a logging configuration file (``logcfg``).  Support for ``logcfg`` is now
deprecated and will be removed in a future build.

For Python code, non-default logging configuration should be implemented via
the `logging` module (see the examples above).  For the command line,
users can port logging configuration features from a config file to the new
command line configuration options as follows.

#. ``level``: The level at and above which logging messages will be
   displayed.  May be one of (from least important to most
   important): DEBUG, INFO, WARNING, ERROR or CRITICAL.

   **Via the command line, specify the log level with "--log-level".**

#. ``handler``: Defines where log messages are to be sent.  By
   default, they are sent to stderr.  However, one may also
   specify:

     - ``file:filename.log`` to send the log messages to the given
       file.

     - ``stdout`` to send log messages to stdout.

   Multiple handlers may be specified by putting the whole value in
   quotes and separating the entries with a comma.

   **Via the command line, specify a log file name with "--log-file".**
   **Specify the output stream with "--log-stream".**

These features, formerly supported by ``logcfg``, will no longer be available
via the command line:

#. ``break_level``: The level at and above which logging messages
   will cause an exception to be raised.  For instance, if you
   would rather stop execution at the first ERROR message (rather
   than continue), set ``break_level`` to ``ERROR``.

#. ``append:filename.log`` to append the log messages to the given file.

#. ``format``: Allows one to customize what each log message
   contains.

These advanced features may still be implemented in Python code via the
`logging` module.

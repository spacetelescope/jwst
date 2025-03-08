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

The name of a file in which to save log information, as well as the desired
level of logging messages, can be specified in an optional configuration file.
Two options exist - if the configuration file should be used for all instances
of the pipeline, the configuration file should be named "stpipe-log.cfg".
This file must be in the same directory in which you run the pipeline in order
for it to be used.

If instead the configuration should be active only when specified,
you should name it something other than "stpipe-log.cfg"; this filename should be
specified using either the ``--logcfg`` parameter to the command line ``strun`` or
using the ``logcfg`` keyword to a .call() execution of either a Step or Pipeline
instance.

If this file does not exist, the default logging mechanism is STDOUT,
with a level of INFO. An example of the contents of the stpipe-log.cfg file is:

::

    [*]
    handler = file:pipeline.log
    level = INFO

If there's no ``stpipe-log.cfg`` file in the working directory, which specifies
how to handle process log information, the default is to display log messages
to stdout.

For example:
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --logcfg=pipeline-log.cfg

Or in an interactive python environment:
::

    result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits",
                                    logcfg="pipeline-log.cfg")

and the file ``pipeline-log.cfg`` contains:
::

    [*]
    handler = file:pipeline.log
    level = INFO

In this example log information is written to a file called ``pipeline.log``.
The ``level`` argument in the log cfg file can be set to one of the standard
logging level designations of ``DEBUG``, ``INFO``, ``WARNING``, ``ERROR``, and
``CRITICAL``. Only messages at or above the specified level
will be displayed.

.. note::

   Setting up ``stpipe-log.cfg`` can lead to confusion, especially if it is
   forgotten about. If one has not run a pipeline in awhile, and then sees no
   logging information, most likely it is because ``stpipe-log.cfg`` is
   present. Consider using a different name and specifying it explicitly on the
   command line.

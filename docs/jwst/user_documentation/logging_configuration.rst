=====================
Logging Configuration
=====================

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

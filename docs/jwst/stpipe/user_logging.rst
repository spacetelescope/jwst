.. _user-logging:

=======
Logging
=======

Log messages are emitted from each Step at different levels of
importance.  The levels used are the standard ones for Python (from
least important to most important:

#. DEBUG
#. INFO
#. WARNING
#. ERROR
#. CRITICAL

By default, only messages of type WARNING or higher are displayed.
This can be controlled by providing a logging configuration file.

Logging configuration
=====================

A logging configuration file can be provided to customize what is
logged.

A logging configuration file is searched for in the following places.
The first one found is used *in its entirety* and all others are
ignored:

#. The file specified with the ``--logcfg`` option to the
   ``strun`` script.

#. The file specified with the ``logcfg`` keyword to a
   .call() execution of a Step or Pipeline.

#. A file called ``stpipe-log.cfg`` in the current working
   directory.

#. ``~/stpipe-log.cfg``

#. ``/etc/stpipe-log.cfg``

The logging configuration file is in the standard ini-file format.

Each section name is a Unix-style filename glob pattern used to match
a particular Step’s logger.  The settings in that section apply only
to that Steps that match that pattern.  For example, to have the
settings apply to all steps, create a section called ``[*]``.  To have
the settings apply only to a Step called ``MyStep``, create a section
called ``[*.MyStep]``.  To apply settings to all Steps that are
substeps of a step called ``MyStep``, call the section
``[*.MyStep.*]``.

In each section, the following may be configured:

#. ``level``: The level at and above which logging messages will be
   displayed.  May be one of (from least important to most
   important): DEBUG, INFO, WARNING, ERROR or CRITICAL.

#. ``break_level``: The level at and above which logging messages
   will cause an exception to be raised.  For instance, if you
   would rather stop execution at the first ERROR message (rather
   than continue), set ``break_level`` to ``ERROR``.

#. ``handler``: Defines where log messages are to be sent.  By
   default, they are sent to stderr.  However, one may also
   specify:

     - ``file:filename.log`` to send the log messages to the given
       file.

     - ``append:filename.log`` to append the log messages to the
       given file.  This is useful over ``file`` if multiple
       processes may need to write to the same log file.

     - ``stdout`` to send log messages to stdout.

   Multiple handlers may be specified by putting the whole value in
   quotes and separating the entries with a comma.

#. ``format``: Allows one to customize what each log message
   contains.  What this string may contain is described in the
   `logging module LogRecord Attributes
   <https://docs.python.org/3/library/logging.html#logrecord-attributes>`_
   section of the Python standard library.

Examples
========

The following configuration turns on all log messages and saves them
in the file myrun.log:

.. code-block:: ini

    [*]
    level = INFO
    handler = file:myrun.log

In a scenario where the user is debugging a particular Step, they may
want to hide all logging messages except for that Step, and stop when
hitting any warning for that Step:

.. code-block:: ini

    [*]
    level = CRITICAL

    [*.MyStep]
    level = INFO
    break_level = WARNING

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

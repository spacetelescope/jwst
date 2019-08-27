=====
Steps
=====

.. _configuring-a-step:

Configuring a Step
==================

This section describes how to instantiate a Step and set configuration
parameters on it.

Steps can be configured by either:

    - Writing a configuration file

    - Instantiating the Step directly from Python

Running a Step from a configuration file
========================================

A Step configuration file is in the well-known ini-file format.
stpipe uses the `ConfigObj
<https://configobj.readthedocs.io/en/latest/>`_ library to parse
them.

Every step configuration file must contain the ``name`` and ``class``
of the step, followed by parameters that are specific to the step
being run.

``name`` defines the name of the step.  This is distinct from the
class of the step, since the same class of Step may be configured in
different ways, and it is useful to be able to have a way of
distinguishing between them.  For example, when Steps are combined
into :ref:`stpipe-user-pipelines`, a Pipeline may use the same Step class
multiple times, each with different configuration parameters.

``class`` specifies the Python class to run.  It should be a
fully-qualified Python path to the class.  Step classes can ship with
``stpipe`` itself, they may be part of other Python packages, or they
exist in freestanding modules alongside the configuration file.  For
example, to use the ``SystemCall`` step included with ``stpipe``, set
``class`` to ``stpipe.subprocess.SystemCall``.  To use a class called
``Custom`` defined in a file ``mysteps.py`` in the same directory as
the configuration file, set ``class`` to ``mysteps.Custom``.

Below ``name`` and ``class`` in the configuration file are parameters
specific to the Step.  The set of accepted parameters is defined in
the Step’s spec member.  You can print out a Step’s configspec using
the ``stspec`` commandline utility.  For example, to print the
configspec for an imaginary step called `stpipe.cleanup`::

    $ stspec stpipe.cleanup
    # The threshold below which to apply cleanup
    threshold = float()

    # A scale factor
    scale = float()

    # The output file to save to
    output_file = output_file(default = None)

.. note::

    Configspec information can also be displayed from Python, just
    call ``print_configspec`` on any Step class.

.. doctest-skip::

  >>> from jwst.stpipe import cleanup
  >>> cleanup.print_configspec()
  >>> # The threshold below which to apply cleanup
  >>> threshold = float()
  >>> # A scale factor
  >>> scale = float()

Using this information, one can write a configuration file to use this
step.  For example, here is a configuration file (``do_cleanup.cfg``)
that runs the ``stpipe.cleanup`` step to clean up an image.

.. code-block:: ini

    name = "MyCleanup"
    class = "stpipe.cleanup"

    threshold = 42.0
    scale = 0.01


.. _strun:

Running a Step from the commandline
-----------------------------------
The ``strun`` command can be used to run Steps from the commandline.

The first argument may be either:

    - The path to a configuration file

    - A Python class

Additional configuration parameters may be passed on the commandline.
These parameters override any that are present in the configuration
file.  Any extra positional parameters on the commandline are passed
to the step's process method.  This will often be input filenames.

For example, to use an existing configuration file from above, but
override it so the threshold parameter is different::

    $ strun do_cleanup.cfg input.fits --threshold=86

To display a list of the parameters that are accepted for a given Step
class, pass the ``-h`` parameter, and the name of a Step class or
configuration file::

    $ strun -h do_cleanup.cfg
    usage: strun [--logcfg LOGCFG] cfg_file_or_class [-h] [--pre_hooks]
                 [--post_hooks] [--skip] [--scale] [--extname]

    optional arguments:
      -h, --help       show this help message and exit
      --logcfg LOGCFG  The logging configuration file to load
      --verbose, -v    Turn on all logging messages
      --debug          When an exception occurs, invoke the Python debugger, pdb
      --pre_hooks
      --post_hooks
      --skip           Skip this step
      --scale          A scale factor
      --threshold      The threshold below which to apply cleanup
      --output_file    File to save the output to

Every step has an `--output_file` parameter.  If one is not provided,
the output filename is determined based on the input file by appending
the name of the step.  For example, in this case, `foo.fits` is output
to `foo_cleanup.fits`.

Debugging
`````````

To output all logging output from the step, add the `--verbose` option
to the commandline.  (If more fine-grained control over logging is
required, see :ref:`user-logging`).

To start the Python debugger if the step itself raises an exception,
pass the `--debug` option to the commandline.

Running a Step in Python
------------------------

Running a step can also be done inside the Python interpreter and is as simple
as calling its `run()` or `call()` classmethods.

run()
`````

The `run()` classmethod will run a previously instantiated step class. This is
very useful if one wants to setup the step's attributes first, then run it::

    from jwst.flatfield import FlatFieldStep

    mystep = FlatFieldStep()
    mystep.override_sflat = ‘sflat.fits’
    output = mystep.run(input)

Using the `.run()` method is the same as calling the instance or class directly.
They are equivalent::

    output = mystep(input)

call()
``````

If one has all the configuration in a configuration file or can pass the
arguments directly to the step, one can use call(), which creates a new
instance of the class every time you use the `call()` method.  So::

    output = mystep.call(input)

makes a new instance of `FlatFieldStep` and then runs. Because it is a new
instance, it ignores any attributes of `mystep` that one may have set earlier,
such overriding the sflat.

The nice thing about call() is that it can take a configuration file, so::

    output = mystep.call(input, config_file=’my_flatfield.cfg’)

and it will take all the configuration from the config file.

Configuration parameters may be passed to the step by setting the `config_file`
kwarg in `call` (which takes a path to a configuration file) or as keyword
arguments.  Any remaining positional arguments are passed along to the step's
`process()` method::

    from jwst.stpipe import cleanup

    cleanup.call('image.fits', config_file='do_cleanup.cfg', threshold=42.0)

So use call() if you’re passing a config file or passing along args or kwargs.
Otherwise use run().



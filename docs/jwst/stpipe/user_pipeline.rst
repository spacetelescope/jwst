.. _stpipe-user-pipelines:

=========
Pipelines
=========

.. TODO: Rewrite using a real-world example

It is important to note that a Pipeline is also a Step, so everything
that applies to a Step in the :ref:`stpipe-user-steps` chapter also
applies to Pipelines.

Configuring a Pipeline
======================

This section describes how to set parameters on the individual steps
in a pipeline.  To change the order of steps in a pipeline, one must
write a Pipeline subclass in Python.  That is described in the
:ref:`devel-pipelines` section of the developer documentation.

Just as with Steps, Pipelines can by configured either by a
configuration file or directly from Python.

From a configuration file
-------------------------

A Pipeline configuration file follows the same format as a Step
configuration file: the ini-file format used by the `ConfigObj
<https://configobj.readthedocs.io/en/latest/>`_ library.

Here is an example pipeline configuration file for a `TestPipeline`
class:

.. code-block:: ini

    name = "TestPipeline"
    class = "stpipe.test.test_pipeline.TestPipeline"

    science_filename = "science.fits"
    flat_filename = "flat.fits"
    output_filename = "output.fits"

    [steps]
      [[flat_field]]
        config_file = "flat_field.cfg"
        threshold = 42.0

      [[combine]]
        skip = True

Just like a ``Step``, it must have ``name`` and ``class`` values.
Here the ``class`` must refer to a subclass of `stpipe.Pipeline`.

Following ``name`` and ``class`` is the ``[steps]`` section.  Under
this section is a subsection for each step in the pipeline.  To figure
out what configuration parameters are available, use the `stspec`
script (just as with a regular step):

.. code-block:: python

    > stspec stpipe.test.test_pipeline.TestPipeline
    science_filename = input_file()  # The input science filename
    flat_filename = input_file()     # The input flat filename
    skip = bool(default=False)   # Skip this step
    output_filename = output_file()  # The output filename
    [steps]
    [[combine]]
    config_file = string(default=None)
    skip = bool(default=False)   # Skip this step
    [[flat_field]]
    threshold = float(default=0.0)# The threshold below which to remove
    multiplier = float(default=1.0)# Multiply by this number
    skip = bool(default=False)   # Skip this step
    config_file = string(default=None)

For each Step’s section, the parameters for that step may either be
specified inline, or specified by referencing an external
configuration file just for that step.  For example, a pipeline
configuration file that contains:

.. code-block:: ini

    [steps]
      [[flat_field]]
        threshold = 42.0
        multiplier = 2.0

is equivalent to:

.. code-block:: ini

    [steps]
      [[flat_field]]
        config_file = myflatfield.cfg

with the file ``myflatfield.cfg`` in the same directory:

.. code-block:: ini

    threshold = 42.0
    multiplier = 2.0

If both a ``config_file`` and additional parameters are specified, the
``config_file`` is loaded, and then the local parameters override
them.

Any optional parameters for each Step may be omitted, in which case
defaults will be used.


From Python
-----------

A pipeline may be configured from Python by passing a nested
dictionary of parameters to the Pipeline’s constructor.  Each key is
the name of a step, and the value is another dictionary containing
parameters for that step.  For example, the following is the
equivalent of the configuration file above:

.. code-block:: python

    from stpipe.test.test_pipeline import TestPipeline

    steps = {
        'flat_field':   {'threshold': 42.0}
        }

    pipe = TestPipeline(
        "TestPipeline",
        config_file=__file__,
        science_filename="science.fits",
        flat_filename="flat.fits",
        output_filename="output.fits",
        steps=steps)

Running a Pipeline
==================

From the commandline
--------------------

The same ``strun`` script used to run Steps from the commandline can
also run Pipelines.

The only wrinkle is that any step parameters overridden from the
commandline use dot notation to specify the parameter name.  For
example, to override the ``threshold`` value on the ``flat_field``
step in the example pipeline above, one can do::

    > strun stpipe.test.test_pipeline.TestPipeline --steps.flat_field.threshold=48

From Python
-----------

Once the pipeline has been configured (as above), just call the
instance to run it.

    pipe()

Caching details
---------------

The results of a Step are cached using Python pickles.  This allows
virtually most of the standard Python data types to be cached.  In
addition, any FITS models that are the result of a step are saved as
standalone FITS files to make them more easily used by external tools.
The filenames are based on the name of the substep within the
pipeline.

Hooks
=====

Each Step in a pipeline can also have pre- and post-hooks associated.
Hooks themselves are Step instances, but there are some conveniences
provided to make them easier to specify in a configuration file.

Pre-hooks are run right before the Step.  The inputs to the pre-hook
are the same as the inputs to their parent Step.
Post-hooks are run right after the Step.  The inputs to the post-hook
are the return value(s) from the parent Step. The return values are
always passed as a list. If the return value from the parent Step is a
single item, a list of this single item is passed to the post hooks.
This allows the post hooks to modify the return results, if necessary.

Hooks are specified using the ``pre_hooks`` and ``post_hooks``
configuration parameter associated with each step.  More than one pre-
or post-hook may be assigned, and they are run in the order they are
given.  There can also be ``pre_hooks`` and ``post_hooks`` on the
Pipeline as a whole (since a Pipeline is also a Step).  Each of these
parameters is a list of strings, where each entry is one of:

   - An external commandline application.  The arguments can be
     accessed using {0}, {1} etc.  (See
     `stpipe.subproc.SystemCall`).

   - A dot-separated path to a Python Step class.

   - A dot-separated path to a Python function.

For example, here’s a ``post_hook`` that will display a FITS file in
the ``ds9`` FITS viewer the ``flat_field`` step has done flat field
correction on it:

.. code-block:: ini

    [steps]
      [[flat_field]]
        threshold = 42.0
        post_hooks = "ds9 {0}",

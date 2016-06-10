.. _devel-pipelines:

=========
Pipelines
=========

.. _writing-a-pipeline:

Writing a Pipeline
==================

There are two ways to go about writing a pipeline depending on how
much flexibility is required.

1. A linear pipeline defines a simple linear progression of steps
where each step has a single input and a single output flowing
directly into the next step.

2. A flexible pipeline allows the pipeline to be defined in Python
code and all of the tools that implies, such as loops, conditionals
and multiple inputs and/or outputs.

Linear pipeline
---------------

To create a linear pipeline, one inherits from the
`~stpipe.LinearPipeline` class and adds a special member
`pipeline_steps` to define the order of the steps::

    from jwst_lib.stpipe import LinearPipeline

    # Some locally-defined steps
    from . import FlatField, RampFitting

    class ExampleLinearPipeline(LinearPipeline):
        """
        This example linear pipeline has only two steps.
        """
        pipeline_steps = [
            ('flat_field', FlatField),
            ('ramp_fitting', RampFitting)
            ]

The `pipeline_steps` member is a list of tuples.  Each tuple is a pair
(*name*, *class*) where *name* is the name of the specific step, and
*class* is the step's class.  Both are required so the same step class
can be used multiple times in the pipeline.  The name is also used for
the section headings in the pipeline's configuration file.

Flexible pipeline
-----------------

The basics of writing a flexible Pipeline are just like
:ref:`writing-a-step`, but instead of inheriting from the
`~stpipe.Step` class, one inherits from the `~stpipe.Pipeline` class.

In addition, a Pipeline subclass defines what its Steps so that the
framework can configure parameters for the individual Steps.  This is
done with the ``step_defs`` member, which is a dictionary mapping step
names to step classes.  This dictionary defines what the Steps are,
but says nothing about their order or how data flows from one Step to
the next.  That is defined in Python code in the Pipeline’s
``process`` method. By the time the Pipeline’s ``process`` method is
called, the Steps in ``step_defs`` will be instantiated as member
variables.

For example, here is a pipeline with two steps: one that processes
each chip of a multi-chip FITS file, and another to combine the chips
into a single image::

    from jwst_lib.stpipe import Pipeline

    from jwst_lib.models import ImageModel

    # Some locally-defined steps
    from . import FlatField, Combine

    class ExamplePipeline(Pipeline):
        """
        This example pipeline demonstrates how to combine steps
        using Python code, in some way that it not necessarily
        a linear progression.
        """

        step_defs = {
            'flat_field': FlatField,
            'combine': Combine,
            }

        def process(self, input):
            with ImageModel(input) as science:

                flattened = self.flat_field(science, self.multiplier)

                combined = self.combine(flattened)

            return combined

        spec = """
        multiplier = float()     # A multiplier constant
        """

When writing the spec member for a Pipeline, only the parameters
that apply to the Pipeline as a whole need to be included.  The
parameters for each Step are automatically loaded in by the framework.

In the case of the above example, we define two new pipeline
configuration parameters for the flat field file and the output
filename.

The parameters for the individual substeps that make up the Pipeline
will be implicitly added by the framework.

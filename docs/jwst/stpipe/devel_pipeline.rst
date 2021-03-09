.. _devel-pipelines:

=========
Pipelines
=========

.. _writing-a-pipeline:

Writing a Pipeline
==================

The basics of writing a Pipeline are just like
:ref:`writing-a-step`, but instead of inheriting from the
`~stpipe.Step` class, one inherits from the `~stpipe.Pipeline` class.

In addition, a Pipeline subclass defines what its Steps are so that the
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

    from jwst.stpipe import Pipeline

    from jwst.datamodels import ImageModel

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

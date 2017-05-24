from ... import (Pipeline, Step)
from .... import datamodels
from ....datamodels import ImageModel


class AnotherDummyStep(Step):
    """
    This is a crazy do-nothing step that demonstrates configuration
    argument parsing.
    """

    spec = """
    par1 = float() # Control the frobulization

    par2 = string() # Reticulate the splines

    par3 = boolean(default=False) # Does it blend?

    [foo]

    """

    reference_file_types = ['flat_field']

    def process(self, a=0, b=0):
        self.log.info("Found a: {0}, b: {1}".format(a, b))
        self.log.info("par1: {0}".format(self.par1))
        self.log.info("par2: {0}".format(self.par2))
        self.log.info("par3: {0}".format(self.par3))

        return a + b


class OptionalRefTypeStep(Step):
    """
    This is a do-nothing step that demonstrates optionally omitting
    a reference file.
    """

    reference_file_types = ['to_be_ignored_ref_type']

    def process(self):
        ref_file = self.get_reference_file(datamodels.open(), 'to_be_ignored_ref_type')
        assert ref_file == ""


class StepWithModel(Step):
    """A step with a model"""

    spec = """
    """

    def process(self, *args):
        from ....datamodels import ImageModel

        model = ImageModel(args[0])

        return model


class SaveStep(Step):
    """
    Step with explicit save.
    """

    spec = """
    """

    def process(self, *args):
        model = ImageModel(args[0])

        self.log.info('Saving model as "processed"')
        self.save_model(model, 'processed')

        return model


class SavePipeline(Pipeline):
    """Save model in pipeline"""

    spec = """
    """

    step_defs = {
        'stepwithmodel': StepWithModel,
        'savestep': SaveStep
    }

    def process(self, *args):
        model = ImageModel(args[0])

        r = self.stepwithmodel(model)
        r = self.savestep(r)

        return r


class ProperPipeline(Pipeline):
    """Pipeline with proper output setupt"""

    spec = """
    """

    step_defs = {
        'stepwithmodel': StepWithModel,
        'another_stepwithmodel': StepWithModel,
    }

    def process(self, *args):

        self.output_basename = 'ppbase'
        self.suffix = 'pp'

        model = ImageModel(args[0])

        self.stepwithmodel.suffix = 'swm'
        r = self.stepwithmodel(model)
        self.another_stepwithmodel.suffix = 'aswm'
        r = self.another_stepwithmodel(r)

        return r

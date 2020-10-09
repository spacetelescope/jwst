from jwst.stpipe import (Pipeline, Step)
from jwst import datamodels
from jwst.datamodels import (
    ImageModel,
    ModelContainer,
)

class StepWithReference(Step):
    """Step that refers to a reference file"""

    reference_file_types = ['flat']

    def process(self, data):
        return data


class PipeWithReference(Pipeline):
    """Pipeline calling step with reference"""

    spec = """
    """

    step_defs = {'step_with_reference': StepWithReference}

    def process(self, data):
        result = self.step_with_reference(data)

        return result


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


class WithDefaultsStep(Step):
    """A step that contains defaults for each of its pars."""

    spec = """
    par1 = string(default='default par1 value')
    par2 = string(default='default par2 value')
    par3 = string(default='default par3 value')
    par4 = string(default='default par4 value')
    """

    def process(self, input):
        self.log.info("Parameters par1=%s, par2=%s, par3=%s, par4=%s",
                      self.par1, self.par2, self.par3, self.par3)

        return input


class MakeListStep(Step):
    """Make a list of all arguments and parameters"""

    spec = """
    par1 = float() # Control the frobulization
    par2 = string() # Reticulate the splines
    par3 = boolean(default=False) # Does it blend?
    """

    def process(self, a=None, b=None):
        self.log.info('Arguments a=%s b=%s', a, b)
        self.log.info('Parameters par1=%s, par2=%s, par3=%s',
                      self.par1, self.par2, self.par3)

        result = [
            item
            for item in [a, b, self.par1, self.par2, self.par3]
            if item is not None
        ]

        self.log.info('The list is %s', result)
        return result


class OptionalRefTypeStep(Step):
    """
    This is a do-nothing step that demonstrates optionally omitting
    a reference file.
    """

    reference_file_types = ['to_be_ignored_ref_type']

    def process(self):
        ref_file = self.get_reference_file(datamodels.open(), 'to_be_ignored_ref_type')
        assert ref_file == ""


class PostHookStep(Step):
    """A step to try out hooks"""

    spec = """
    """

    def process(self, *args):
        self.log.info('Received args: "{}"'.format(args))
        self.log.info('Self.parent = "{}"'.format(self.parent))

        args[0].post_hook_run = True
        self.parent.post_hook_run = True


class PostHookWithReturnStep(Step):
    """A step to try out hooks"""

    spec = """
    """

    def process(self, *args):
        self.log.info('Received args: "{}"'.format(args))
        self.log.info('Self.parent = "{}"'.format(self.parent))

        self.parent.post_hook_run = True
        return 'PostHookWithReturnStep executed'


class PreHookStep(Step):
    """A step to try out hooks"""

    spec = """
    """

    def process(self, *args):
        self.log.info('Received args: "{}"'.format(args))
        self.log.info('Self.parent = "{}"'.format(self.parent))

        args[0].pre_hook_run = True
        self.parent.pre_hook_run = True


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


class StepWithContainer(Step):
    """A step with a model"""

    spec = """
    """

    def process(self, *args):
        container = ModelContainer()
        model1 = ImageModel(args[0]).copy()
        model2 = ImageModel(args[0]).copy()
        model1.meta.filename = 'swc_model1.fits'
        model2.meta.filename = 'swc_model2.fits'
        container.append(model1)
        container.append(model2)

        return container


class StepWithModel(Step):
    """A step with a model"""

    spec = """
    """

    def process(self, *args):
        model = self.open_model(args[0])
        model = ImageModel(model)
        return model


class EmptyPipeline(Pipeline):
    """A pipeline that has no substeps"""

    spec = """
    par1 = string(default='Name the atomizer') # Control the frobulization
    """

    def process(self, *args):

        return args


class ProperPipeline(Pipeline):
    """Pipeline with proper output setup"""

    spec = """
    """

    step_defs = {
        'stepwithmodel': StepWithModel,
        'another_stepwithmodel': StepWithModel,
        'stepwithcontainer': StepWithContainer,
        'withdefaultsstep': WithDefaultsStep
    }

    def process(self, *args):

        self.suffix = 'pp'

        model = ImageModel(args[0])

        self.stepwithmodel.suffix = 'swm'
        r = self.stepwithmodel(model)
        self.another_stepwithmodel.suffix = 'aswm'
        r = self.another_stepwithmodel(r)
        self.stepwithcontainer.suffix = 'swc'
        r = self.stepwithcontainer(r)
        self.withdefaultsstep.suffix = 'wds'
        r = self.withdefaultsstep(r)

        return r


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


class MakeListPipeline(Pipeline):
    """A pipeline that calls MakeListStep"""

    spec = """
    par1 = string(default='Name the atomizer') # Control the frobulization
    """

    step_defs = {
        'make_list': MakeListStep,
    }

    def process(self, *args):

        return self.make_list(*args)

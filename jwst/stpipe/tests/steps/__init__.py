from jwst.stpipe import Pipeline, Step

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import ImageModel
from jwst.datamodels import ModelContainer


class StepWithReference(Step):
    """Step that refers to a reference file."""

    reference_file_types = ["flat"]

    def process(self, data):  # noqa: D102
        return data


class PipeWithReference(Pipeline):
    """Pipeline calling step with reference."""

    spec = """
    """

    step_defs = {"step_with_reference": StepWithReference}

    def process(self, data):  # noqa: D102
        result = self.step_with_reference.run(data)

        return result


class AnotherDummyStep(Step):
    """Do-nothing step that demonstrates configuration argument parsing."""

    spec = """
    par1 = float() # Control the frobulization

    par2 = string() # Reticulate the splines

    par3 = boolean(default=False) # Does it blend?

    [foo]

    """
    class_alias = "stpipe_dummy"

    reference_file_types = ["flat_field"]

    def process(self, a=0, b=0):  # noqa: D102
        self.log.info(f"Found a: {a}, b: {b}")
        self.log.info(f"par1: {self.par1}")
        self.log.info(f"par2: {self.par2}")
        self.log.info(f"par3: {self.par3}")

        return a + b


class WithDefaultsStep(Step):
    """A step that contains defaults for each of its pars."""

    spec = """
    par1 = string(default='default par1 value')
    par2 = string(default='default par2 value')
    par3 = string(default='default par3 value')
    par4 = string(default='default par4 value')
    """

    def process(self, input_data):  # noqa: D102
        self.log.info(
            "Parameters par1=%s, par2=%s, par3=%s, par4=%s",
            self.par1,
            self.par2,
            self.par3,
            self.par3,
        )

        return input_data


class MakeListStep(Step):
    """Make a list of all arguments and parameters."""

    spec = """
    par1 = float() # Control the frobulization
    par2 = string() # Reticulate the splines
    par3 = boolean(default=False) # Does it blend?
    """

    def process(self, a=None, b=None):  # noqa: D102
        self.log.info("Arguments a=%s b=%s", a, b)
        self.log.info("Parameters par1=%s, par2=%s, par3=%s", self.par1, self.par2, self.par3)

        result = [item for item in [a, b, self.par1, self.par2, self.par3] if item is not None]

        self.log.info("The list is %s", result)
        return result


class OptionalRefTypeStep(Step):
    """Do-nothing step that demonstrates optionally omitting a reference file."""

    reference_file_types = ["to_be_ignored_ref_type"]

    def process(self):  # noqa: D102
        ref_file = self.get_reference_file(datamodels.open(), "to_be_ignored_ref_type")
        assert ref_file == ""  # noqa: S101


class PostHookStep(Step):
    """A step to try out hooks."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        self.log.info(f'Received args: "{args}"')
        self.log.info(f'Self.parent = "{self.parent}"')

        args[0].post_hook_run = True
        self.parent.post_hook_run = True


class PostHookWithReturnStep(Step):
    """A step to try out hooks."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        self.log.info(f'Received args: "{args}"')
        self.log.info(f'Self.parent = "{self.parent}"')

        self.parent.post_hook_run = True
        return "PostHookWithReturnStep executed"


class PreHookStep(Step):
    """A step to try out hooks."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        self.log.info(f'Received args: "{args}"')
        self.log.info(f'Self.parent = "{self.parent}"')

        args[0].pre_hook_run = True
        self.parent.pre_hook_run = True


class SaveStep(Step):
    """Step with explicit save."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        model = ImageModel(args[0])

        self.log.info('Saving model as "processed"')
        self.save_model(model, "processed")

        return model


class StepWithContainer(Step):
    """A step with a model."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        container = []
        if isinstance(args[0], ModelContainer):
            model = args[0][0]
        else:
            model = args[0]
        model1 = ImageModel(model).copy()
        model2 = ImageModel(model).copy()
        model1.meta.filename = "swc_model1.fits"
        model2.meta.filename = "swc_model2.fits"
        container.append(model1)
        container.append(model2)

        return container


class StepWithModel(Step):
    """A step with a model."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        model = self.open_model(args[0])
        model = ImageModel(model)
        return model


class EmptyPipeline(Pipeline):
    """A pipeline that has no substeps."""

    spec = """
    par1 = string(default='Name the atomizer') # Control the frobulization
    """

    def process(self, *args):  # noqa: D102
        return args


class ProperPipeline(Pipeline):
    """Pipeline with proper output setup."""

    spec = """
    """

    step_defs = {
        "stepwithmodel": StepWithModel,
        "another_stepwithmodel": StepWithModel,
        "stepwithcontainer": StepWithContainer,
        "withdefaultsstep": WithDefaultsStep,
    }

    def process(self, *args):  # noqa: D102
        self.suffix = "pp"

        model = ImageModel(args[0])

        self.stepwithmodel.suffix = "swm"
        r = self.stepwithmodel.run(model)
        self.another_stepwithmodel.suffix = "aswm"
        r = self.another_stepwithmodel.run(r)
        self.stepwithcontainer.suffix = "swc"
        r = self.stepwithcontainer.run(r)
        self.withdefaultsstep.suffix = "wds"
        r = self.withdefaultsstep.run(r)

        return r


class SavePipeline(Pipeline):
    """Save model in pipeline."""

    spec = """
    """

    step_defs = {"stepwithmodel": StepWithModel, "savestep": SaveStep}

    def process(self, *args):  # noqa: D102
        model = ImageModel(args[0])

        r = self.stepwithmodel.run(model)
        r = self.savestep.run(r)

        return r


class MakeListPipeline(Pipeline):
    """A pipeline that calls MakeListStep."""

    spec = """
    par1 = string(default='Name the atomizer') # Control the frobulization
    """

    step_defs = {
        "make_list": MakeListStep,
    }

    def process(self, *args):  # noqa: D102
        return self.make_list.run(*args)


class CalLogsStep(Step):
    """Step for testing cal_logs."""

    class_alias = "cal_logs_step"

    def process(self, msg):  # noqa: D102
        from stdatamodels.jwst.datamodels import ImageModel

        self.log.info(msg)
        return ImageModel()


class CalLogsPipeline(Pipeline):
    """Pipeline for testing cal_logs."""

    class_alias = "cal_logs_pipeline"

    step_defs = {
        "a_step": CalLogsStep,
    }

    def process(self, msg):  # noqa: D102
        self.log.info(msg)
        return self.a_step.run(msg)

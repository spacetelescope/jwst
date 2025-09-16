import logging

from jwst.stpipe import Pipeline, Step

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import ImageModel
from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe.utilities import record_step_status


log = logging.getLogger("jwst.stpipe.tests.steps")


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
        log.info(f"Found a: {a}, b: {b}")
        log.info(f"par1: {self.par1}")
        log.info(f"par2: {self.par2}")
        log.info(f"par3: {self.par3}")

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
        log.info(
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
        log.info("Arguments a=%s b=%s", a, b)
        log.info("Parameters par1=%s, par2=%s, par3=%s", self.par1, self.par2, self.par3)

        result = [item for item in [a, b, self.par1, self.par2, self.par3] if item is not None]

        log.info("The list is %s", result)
        return result


class OptionalRefTypeStep(Step):
    """Do-nothing step that demonstrates optionally omitting a reference file."""

    reference_file_types = ["to_be_ignored_ref_type"]

    def process(self):  # noqa: D102
        ref_file = self.get_reference_file(datamodels.JwstDataModel(), "to_be_ignored_ref_type")
        assert ref_file == ""  # noqa: S101


class PostHookStep(Step):
    """A step to try out hooks."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        log.info(f'Received args: "{args}"')
        log.info(f'Self.parent = "{self.parent}"')

        args[0].post_hook_run = True
        self.parent.post_hook_run = True


class PostHookWithReturnStep(Step):
    """A step to try out hooks."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        log.info(f'Received args: "{args}"')
        log.info(f'Self.parent = "{self.parent}"')

        self.parent.post_hook_run = True
        return "PostHookWithReturnStep executed"


class PreHookStep(Step):
    """A step to try out hooks."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        log.info(f'Received args: "{args}"')
        log.info(f'Self.parent = "{self.parent}"')

        args[0].pre_hook_run = True
        self.parent.pre_hook_run = True


class SaveStep(Step):
    """Step with explicit save."""

    spec = """
    """

    def process(self, *args):  # noqa: D102
        model = ImageModel(args[0])

        log.info('Saving model as "processed"')
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

        log.info(msg)
        return ImageModel()


class CalLogsPipeline(Pipeline):
    """Pipeline for testing cal_logs."""

    class_alias = "cal_logs_pipeline"

    step_defs = {
        "a_step": CalLogsStep,
    }

    def process(self, msg):  # noqa: D102
        log.info(msg)
        return self.a_step.run(msg)


class PrepareOutputStep(Step):
    """Step to test the prepare_output method with defaults."""
    class_alias = "prepare_output"

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")

        result = self.prepare_output(input_data)
        record_step_status(result, "prepare_output", True)
        log.info(f"Output data is {type(result)}")

        return result


class PrepareOutputNoOpenStep(Step):
    """Step to test the prepare_output method, turning off open_models."""
    class_alias = "prepare_output_no_open"

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")

        result = self.prepare_output(input_data, open_models=False)
        log.info(f"Output data is {type(result)}")

        return result


class PrepareOutputForceCopyStep(Step):
    """Step to test the prepare_output method, forcing a copy."""
    class_alias = "prepare_output_force_copy"

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")

        result = self.prepare_output(input_data, make_copy=True)
        record_step_status(result, "prepare_output", True)
        log.info(f"Output data is {type(result)}")

        return result


class PrepareOutputForceCopyNoOpenStep(Step):
    """Step to test the prepare_output method, forcing a copy and not opening models."""
    class_alias = "prepare_output_force_copy_no_open"

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")
        result = self.prepare_output(input_data, make_copy=True, open_models=False)
        log.info(f"Output data is {type(result)}")

        return result


class PrepareOutputNoCopyStep(Step):
    """Step to test the prepare_output method, without making a copy."""
    class_alias = "prepare_output_no_copy"

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")

        result = self.prepare_output(input_data, make_copy=False)
        record_step_status(result, "prepare_output", True)
        log.info(f"Output data is {type(result)}")

        return result


class PrepareOutputPipeline(Pipeline):
    """Pipeline to test the prepare_output method."""
    class_alias = "prepare_output_pipeline"

    step_defs = {
        "step1": PrepareOutputNoOpenStep,
        "step2": PrepareOutputStep,
        "step3": PrepareOutputStep,
    }

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")

        result = self.prepare_output(input_data)
        log.info(f"Opened data is {type(result)}")

        result = self.step1.run(result)
        log.info(f"Intermediate data 1 is {type(result)}")

        result = self.step2.run(result)
        log.info(f"Intermediate data 2 is {type(result)}")

        result = self.step3.run(result)
        log.info(f"Output data is {type(result)}")

        return result


class PrepareOutputForceCopyPipeline(Pipeline):
    """Pipeline to test the prepare_output method with a step that forces a copy."""
    class_alias = "prepare_output_pipeline"

    step_defs = {
        "step1": PrepareOutputStep,
        "step2": PrepareOutputForceCopyStep,
    }

    def process(self, input_data):
        log.info(f"Input data is {type(input_data)}")

        result = self.prepare_output(input_data)
        log.info(f"Opened data is {type(result)}")

        result = self.step1.run(result)
        log.info(f"Intermediate data is {type(result)}")

        result = self.step2.run(result)
        log.info(f"Output data is {type(result)}")

        return result

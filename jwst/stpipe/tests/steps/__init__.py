from ... import Step
from .... import datamodels

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


class SaveStep(Step):
    """
    """

    spec = """
    """

    def process(self, *args):
        from ....datamodels import ImageModel

        model = ImageModel(args[0])

        self.save_model(model, 'processed')

        return model

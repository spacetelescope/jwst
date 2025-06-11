from jwst.stpipe import Step
import logging

log = logging.getLogger("FOO")


class DummyStep(Step):
    """Placeholder step."""

    spec = """
    foo = string()
    """

    def process(self, *args):  # noqa: D102
        from stdatamodels.jwst.datamodels import ImageModel

        log.info("Default logger")
        log.debug("Default logger")

        log.info(f"Foo: {self.foo}")
        log.debug("Debug!!!")

        return ImageModel(args[0])

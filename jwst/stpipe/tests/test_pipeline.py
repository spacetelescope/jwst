from __future__ import absolute_import, print_function

from os.path import dirname, join, abspath
import sys

import numpy as np
from numpy.testing import assert_allclose

from .. import Step, Pipeline, LinearPipeline
# TODO: Test system call steps


def library_function():
    import logging
    log = logging.getLogger()
    log.info("This is a library function log")


class FlatField(Step):
    """
    An example flat-fielding Step.
    """

    # Load the spec from a file

    def process(self, science, flat):
        from ... import datamodels

        self.log.info("Removing flat field")
        self.log.info("Threshold: {0}".format(self.threshold))
        library_function()

        output = models.ImageModel(data=science.data - flat.data)
        return output


class Combine(Step):
    """
    A Step that combines a list of images.
    """

    def process(self, images):
        from ... import datamodels

        combined = np.zeros((50, 50))
        for image in images:
            combined += image.data
        return models.ImageModel(data=combined)


class Display(Step):
    """
    A Step to display an image.
    """

    def process(self, image):
        pass


class MultiplyBy2(Step):
    """
    A Step that does the incredibly complex thing of multiplying by 2.
    """

    def process(self, image):
        from ... import datamodels

        with models.ImageModel(image) as dm:
            with models.ImageModel() as dm2:
                dm2.data = dm.data * 2
                return dm2


class TestPipeline(Pipeline):
    """
    A test pipeline.
    """

    step_defs = {
        'flat_field'  : FlatField,
        'combine'     : Combine,
        'display'     : Display
        }

    spec = """
    science_filename = input_file()  # The input science filename
    flat_filename = input_file(default=None)     # The input flat filename
    output_filename = output_file()  # The output filename
    """

    def process(self, *args):
        from ... import datamodels

        science = models.open(self.science_filename)
        if self.flat_filename is None:
            self.flat_filename = join(dirname(__file__), "data/flat.fits")
        flat = models.open(self.flat_filename)
        calibrated = []
        calibrated.append(self.flat_field(science, flat))
        combined = self.combine(calibrated)
        self.display(combined)
        dm = models.ImageModel(combined)
        dm.save(self.output_filename)
        return dm


def test_pipeline():
    pipeline_fn = join(dirname(__file__), 'steps', 'python_pipeline.cfg')
    pipe = Step.from_config_file(pipeline_fn)
    pipe.output_filename = "output.fits"

    assert pipe.flat_field.threshold == 42.0
    assert pipe.flat_field.multiplier == 2.0

    pipe.run()


def test_pipeline_python():
    steps = {
        'flat_field':   {'threshold': 42.0}
        }

    pipe = TestPipeline(
        "TestPipeline",
        config_file=__file__,
        steps=steps,
        science_filename=abspath(join(dirname(__file__), 'data', 'science.fits')),
        flat_filename=abspath(join(dirname(__file__), 'data', 'flat.fits')),
        output_filename="output.fits")

    assert pipe.flat_field.threshold == 42.0
    assert pipe.flat_field.multiplier == 1.0

    pipe.run()


class TestLinearPipeline(LinearPipeline):
    pipeline_steps = [
        ('multiply', MultiplyBy2),
        ('multiply2', MultiplyBy2),
        ('multiply3', MultiplyBy2)
        ]


def test_partial_pipeline():
    pipe = TestLinearPipeline() # output_filename="output.fits")

    pipe.end_step = 'multiply2'
    result = pipe.run(abspath(join(dirname(__file__), 'data', 'science.fits')))

    pipe.start_step = 'multiply3'
    pipe.end_step = None
    result = pipe.run(abspath(join(dirname(__file__), 'data', 'science.fits')))

    assert_allclose(np.sum(result.data), 9969.82514685, rtol=1e-4)


def test_pipeline_commandline():
    args = [
        abspath(join(dirname(__file__), 'steps', 'python_pipeline.cfg')),
        '--steps.flat_field.threshold=47'
        ]

    pipe = Step.from_cmdline(args)

    assert pipe.flat_field.threshold == 47.0
    assert pipe.flat_field.multiplier == 2.0

    pipe.run()


def test_pipeline_commandline_class():
    args = [
        'jwst.stpipe.tests.test_pipeline.TestPipeline',
        '--logcfg={0}'.format(
            abspath(join(dirname(__file__), 'steps', 'log.cfg'))),
        # The file_name parameters are *required*
        '--science_filename={0}'.format(
            abspath(join(dirname(__file__), 'data', 'science.fits'))),
        '--output_filename={0}'.format(
            'output.fits'),
        '--steps.flat_field.threshold=47'
        ]

    pipe = Step.from_cmdline(args)

    assert pipe.flat_field.threshold == 47.0
    assert pipe.flat_field.multiplier == 1.0

    pipe.run()


def test_pipeline_commandline_invalid_args():
    from io import StringIO

    args = [
        'jwst.stpipe.tests.test_pipeline.TestPipeline',
        # The file_name parameters are *required*, and one of them
        # is missing, so we should get a message to that effect
        # followed by the commandline usage message.
        '--flat_filename={0}'.format(
            abspath(join(dirname(__file__), 'data', 'flat.fits'))),
        '--steps.flat_field.threshold=47'
        ]

    __stdout = sys.stdout
    sys.stdout = buffer = StringIO()

    try:
        pipe = Step.from_cmdline(args)
    except ValueError:
        pass
    else:
        raise AssertionError("Wanted a ValueError")
    finally:
        sys.stdout = __stdout

    help = buffer.getvalue()
    print(help)

    # Make sure the comments made it into the help message
    assert "Multiply by this number" in help

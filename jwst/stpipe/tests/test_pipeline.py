from os.path import dirname, join, abspath
import sys

import numpy as np
import pytest

from jwst import datamodels
from jwst.stpipe import Step, Pipeline, crds_client

from .steps import PipeWithReference, StepWithReference


def library_function():
    import logging
    log = logging.getLogger()
    log.info("This is a library function log")


class FlatField(Step):
    """
    An example flat-fielding Step.
    """
    spec = """
        threshold = float(default=0.0)  # The threshold below which to remove
        multiplier = float(default=1.0) # Multiply by this number
    """

    # Load the spec from a file

    def process(self, science, flat):
        self.log.info("Removing flat field")
        self.log.info("Threshold: {0}".format(self.threshold))
        library_function()

        output = datamodels.ImageModel(data=science.data - flat.data)
        return output


class Combine(Step):
    """
    A Step that combines a list of images.
    """

    def process(self, images):
        combined = np.zeros((50, 50))
        for image in images:
            combined += image.data
        return datamodels.ImageModel(data=combined)


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
        with datamodels.ImageModel(image) as dm:
            dm2 = datamodels.ImageModel()
            dm2.data = dm.data * 2
            return dm2


class MyPipeline(Pipeline):
    """
    A test pipeline.
    """

    step_defs = {
        'flat_field': FlatField,
        'combine': Combine,
        'display': Display
        }

    spec = """
    science_filename = input_file()  # The input science filename
    flat_filename = input_file(default=None)     # The input flat filename
    output_filename = output_file()  # The output filename
    """

    def process(self, *args):
        science = datamodels.open(self.science_filename)
        if self.flat_filename is None:
            self.flat_filename = join(dirname(__file__), "data/flat.fits")
        flat = datamodels.open(self.flat_filename)
        calibrated = []
        calibrated.append(self.flat_field(science, flat))
        combined = self.combine(calibrated)
        self.display(combined)
        dm = datamodels.ImageModel(combined)
        dm.save(self.output_filename)
        science.close()
        flat.close()
        return dm


def test_pipeline(_jail):
    pipeline_fn = join(dirname(__file__), 'steps', 'python_pipeline.cfg')
    pipe = Step.from_config_file(pipeline_fn)
    pipe.output_filename = "output.fits"

    assert pipe.flat_field.threshold == 42.0
    assert pipe.flat_field.multiplier == 2.0

    pipe.run()


def test_pipeline_python(_jail):
    steps = {
        'flat_field': {'threshold': 42.0}
        }

    pipe = MyPipeline(
        "MyPipeline",
        config_file=__file__,
        steps=steps,
        science_filename=abspath(join(dirname(__file__), 'data', 'science.fits')),
        flat_filename=abspath(join(dirname(__file__), 'data', 'flat.fits')),
        output_filename="output.fits")

    assert pipe.flat_field.threshold == 42.0
    assert pipe.flat_field.multiplier == 1.0

    pipe.run()


def test_prefetch(_jail, monkeypatch):
    """Test prefetching"""

    # Setup mock to crds to flag if the call was made.
    class MockGetRef:
        called = False

        def mock(self, dataset_model, reference_file_types, observatory=None):
            if 'flat' in reference_file_types:
                self.called = True
            result = {
                reftype: 'N/A'
                for reftype in reference_file_types
            }
            return result
        __call__ = mock
    mock_get_ref = MockGetRef()
    monkeypatch.setattr(crds_client, 'get_multiple_reference_paths', mock_get_ref)

    # Create some data
    model = datamodels.ImageModel((19, 19))
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.detector = 'NRCA1'
    model.meta.instrument.filter = "F070W"
    model.meta.instrument.pupil = 'CLEAR'
    model.meta.observation.date = "2019-01-01"
    model.meta.observation.time = "00:00:00"

    # Run the pipeline with prefetch set.
    StepWithReference.prefetch_references = True
    PipeWithReference.call(model)
    assert mock_get_ref.called

    # Now run with prefetch unset.
    mock_get_ref.called = False
    StepWithReference.prefetch_references = False
    PipeWithReference.call(model)
    assert not mock_get_ref.called


def test_pipeline_commandline(_jail):
    args = [
        join(dirname(__file__), 'steps', 'python_pipeline.cfg'),
        '--steps.flat_field.threshold=47',
        ]

    pipe = Step.from_cmdline(args)

    assert pipe.flat_field.threshold == 47.0
    assert pipe.flat_field.multiplier == 2.0

    pipe.run()


def test_pipeline_commandline_class(_jail):
    args = [
        'jwst.stpipe.tests.test_pipeline.MyPipeline',
        f"--science_filename={join(dirname(__file__), 'data', 'science.fits')}",
        '--output_filename=output.fits',
        '--steps.flat_field.threshold=47'
        ]

    pipe = Step.from_cmdline(args)

    assert pipe.flat_field.threshold == 47.0
    assert pipe.flat_field.multiplier == 1.0

    pipe.run()


def test_pipeline_commandline_invalid_args():
    from io import StringIO

    args = [
        'jwst.stpipe.tests.test_pipeline.MyPipeline',
        # The file_name parameters are *required*, and one of them
        # is missing, so we should get a message to that effect
        # followed by the commandline usage message.
        '--flat_filename={0}'.format(
            abspath(join(dirname(__file__), 'data', 'flat.fits'))),
        '--steps.flat_field.threshold=47'
        ]

    sys.stdout = buffer = StringIO()

    with pytest.raises(ValueError):
        Step.from_cmdline(args)

    help = buffer.getvalue()
    assert "Multiply by this number" in help

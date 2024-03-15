import pytest
import os
from jwst.pipeline.calwebb_image2 import Image2Pipeline
from jwst.stpipe import Step
from jwst.datamodels import ImageModel


INPUT_FILE = "dummy_rate.fits"
OUTPUT_FILE = "custom_name.fits"

@pytest.fixture(scope='module')
def make_dummy_rate_file(tmp_cwd_module):
    '''
    Make and save a dummy rate file in the temporary working directory
    '''

    image = ImageModel((2048, 2048))
    image.data[:, :] = 1
    image.meta.instrument.name = 'NIRCAM'
    image.meta.instrument.filter = 'CLEAR'
    image.meta.exposure.type = 'NRC_IMAGE'
    image.meta.observation.date = '2019-02-27'
    image.meta.observation.time = '13:37:18.548'
    image.meta.date = '2019-02-27T13:37:18.548'
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1

    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]

    image.meta.instrument.channel = 'SHORT'
    image.meta.instrument.module = 'A'
    image.meta.instrument.detector = 'NRCA1'

    # bare minimum wcs info to get assign_wcs step to pass
    image.meta.wcsinfo.crpix1 = 693.5
    image.meta.wcsinfo.crpix2 = 512.5
    image.meta.wcsinfo.v2_ref = -453.37849
    image.meta.wcsinfo.v3_ref = -373.810549
    image.meta.wcsinfo.roll_ref = 272.3237653262276
    image.meta.wcsinfo.ra_ref = 80.54724018120017
    image.meta.wcsinfo.dec_ref = -69.5081101864959

    with image as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope='module', params=[OUTPUT_FILE])
def run_image2_pipeline(make_dummy_rate_file, request):
    '''
    Run pipeline, skipping most steps
    '''
    args = ["calwebb_image2", INPUT_FILE, 
            "--steps.flat_field.skip=true",
            "--steps.photom.skip=true",
            "--steps.resample.skip=true",
            f"--output_file={request.param}",]

    Step.from_cmdline(args)


def test_output_file_rename(run_image2_pipeline):
    '''
    Covers a bug where the output_file parameter was not being
    respected in calls to Image2Pipeline.
    '''
    assert os.path.exists(INPUT_FILE) #ensures tmp_cwd_module is working
    custom_stem = OUTPUT_FILE.split('.')[0]
    for extension in ['cal']:
        assert os.path.exists(f'{custom_stem}_{extension}.fits')
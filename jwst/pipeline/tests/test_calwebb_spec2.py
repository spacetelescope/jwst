import pytest
import os
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.stpipe import Step
from jwst.datamodels import IFUImageModel


INPUT_FILE = "dummy_rate.fits"
OUTPUT_FILE = "custom_name.fits"

@pytest.fixture(scope='module')
def make_dummy_rate_file(tmp_cwd_module):
    '''
    Make and save a dummy rate file in the temporary working directory
    Partially copied from test_background.py
    '''
    image = IFUImageModel((2048, 2048))
    image.data[:, :] = 1

    image.meta.instrument.name = 'NIRSPEC'
    image.meta.instrument.detector = 'NRS1'
    image.meta.instrument.filter = 'CLEAR'
    image.meta.instrument.grating = 'PRISM'
    image.meta.exposure.type = 'NRS_IFU'
    image.meta.observation.date = '2019-02-27'
    image.meta.observation.time = '13:37:18.548'
    image.meta.date = '2019-02-27T13:37:18.548'
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1
    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]
    image.meta.instrument.gwa_tilt = 37.0610

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
def run_spec2_pipeline(make_dummy_rate_file, request):
    '''
    Run pipeline, saving one intermediate step  
    and skipping most of the rest for runtime
    '''
    args = ["calwebb_spec2", INPUT_FILE, 
            "--steps.msa_flagging.skip=true",
            "--steps.nsclean.skip=true",
            "--steps.flat_field.skip=true",
            "--steps.pathloss.skip=true",
            "--steps.photom.skip=true",
            "--steps.pixel_replace.skip=true",
            "--steps.cube_build.save_results=true",
            "--steps.extract_1d.skip=true",
            f"--output_file={request.param}",]

    Step.from_cmdline(args)


def test_output_file_rename(run_spec2_pipeline):
    '''
    Covers a bug where the output_file parameter was not being
    respected in calls to Spec2Pipeline.
    _s3d file expected from cube_build save_results=true
    _cal file expected from pipeline finishing
    '''
    assert os.path.exists(INPUT_FILE)
    custom_stem = OUTPUT_FILE.split('.')[0]
    for extension in ['s3d', 'cal']:
        assert os.path.exists(f'{custom_stem}_{extension}.fits')


def test_filenotfounderror_raised(capsys):
    # Verify the failure is in the traceback message
    with pytest.raises(RuntimeError, match="FileNotFoundError"):
        Spec2Pipeline().run('file_does_not_exist.fits')

    # Verify the failure is printed to stderr
    captured = capsys.readouterr()
    assert 'FileNotFoundError' in captured.err
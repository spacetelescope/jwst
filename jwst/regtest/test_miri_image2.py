from glob import glob
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


@pytest.mark.bigdata
def test_foo(request, rtdata, fitsdiff_defaults, _jail):
    request.node.user_properties = [('output', rtdata.output)]
    rtdata.input_remote = "miri/image/jw00001001001_01101_00001_mirimage_rate.fits"
    input_data = rtdata.get_data(rtdata.input_remote)
    Image2Pipeline.call(input_data, save_results=True)

    rtdata.output = "jw00001001001_01101_00001_mirimage_cal.fits"
    rtdata.truth_local = rtdata.get_truth(rtdata.output)

    diff = FITSDiff(rtdata.output, rtdata.truth_local, **fitsdiff_defaults)
    assert diff.identical, diff.report


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, jail):
    """Run calwebb_image2 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.input_remote = "miri/image/jw00001001001_01101_00001_mirimage_rate.fits"
    input_data = rtdata.get_data(rtdata.input_remote)

    collect_pipeline_cfgs('config')
    config_file = os.path.join('config', 'calwebb_image2.cfg')
    Image2Pipeline.call(input_data, config_file=config_file)

    return rtdata


@pytest.mark.bigdata
def test_miri_image2_completion(run_pipeline):
    files = glob('*_cal.fits')
    files += glob('*_i2d.fits')
    # There should be 2 outputs
    assert len(files) == 2


@pytest.mark.bigdata
@pytest.mark.parametrize("output", [
    'jw00001001001_01101_00001_mirimage_cal.fits',
    'jw00001001001_01101_00001_mirimage_i2d.fits',],
    ids=['cal', 'i2d'],
)
def test_miri_image2(run_pipeline, request, fitsdiff_defaults, output):
    """
    Regression test of calwebb_image2 pipeline performed on MIRI data.
    """
    request.node.user_properties = [('output', output)]
    diff = FITSDiff(output, output, **fitsdiff_defaults)
    assert diff.identical, diff.report()

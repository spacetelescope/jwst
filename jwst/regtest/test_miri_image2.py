from glob import glob
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


@pytest.mark.bigdata
def test_foo(request, rtdata, fitsdiff_default_kwargs, _jail):
    rtdata.get_data("miri/image/jw00001001001_01101_00001_mirimage_rate.fits")

    Image2Pipeline.call(rtdata.input, save_results=True)
    rtdata.output = "jw00001001001_01101_00001_mirimage_cal.fits"
    # Boilerplate needed for every test that interacts with Artifactory
    request.node.user_properties = [('output', rtdata.output)]

    rtdata.get_truth("truth/test_foo/jw00001001001_01101_00001_mirimage_cal.fits")
    assert rtdata.output != rtdata.truth

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report
    assert 0


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, jail):
    """Run calwebb_image2 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw00001001001_01101_00001_mirimage_rate.fits")

    collect_pipeline_cfgs('config')
    config_file = os.path.join('config', 'calwebb_image2.cfg')
    Image2Pipeline.call(rtdata.input, config_file=config_file,
        save_results=True)

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
    ids=['cal', 'i2d'])
def test_miri_image2(run_pipeline, request, fitsdiff_default_kwargs, output):
    """
    Regression test of calwebb_image2 pipeline performed on MIRI data.
    """
    rtdata = run_pipeline
    rtdata.output = output
    # Boilerplate needed for every test that interacts with Artifactory
    request.node.user_properties = [('output', rtdata.output)]

    rtdata.get_truth(os.path.join("truth/test_foo", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
    assert 0

from glob import glob
import os

import pytest
from astropy.io import fits
from ci_watson.artifactory_helpers import get_bigdata

from jwst.pipeline import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from .helpers import fitsdiff_default_args


@pytest.fixture(scope="module")
def run_pipeline(request, artifactory_repos, jail):
    """Run calwebb_image2 pipeline on MIRI imaging data."""
    inputs_root, results_root = artifactory_repos
    envopt = request.config.getoption("env")
    input_data = get_bigdata(inputs_root, envopt, 'miri',
        'test_image2pipeline', 'jw00001001001_01101_00001_mirimage_rate.fits')

    collect_pipeline_cfgs('config')
    config_file = os.path.join('config', 'calwebb_image2.cfg')
    Image2Pipeline.call(input_data, config_file=config_file)

    # truth_data = get_bigdata(inputs_root, envopt, )

    return jail


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
def test_miri_image2(run_pipeline, request, output):
    """
    Regression test of calwebb_image2 pipeline performed on MIRI data.
    """
    request.node.user_properties = [('output', output)]
    diff = fits.diff.FITSDiff(output, output, **fitsdiff_default_args)
    assert diff.identical, diff.report()
    assert output == 'foo'

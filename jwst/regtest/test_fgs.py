"""Regression tests for FGS Guidestar in ID and FINEGUIDE modes"""
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope='module')
def run_pipeline(jail, rtdata_module):
    """Run the pipelines"""
    rtdata = rtdata_module
    rtdata.get_data('fgs/level1b/exptype_fgs_image_uncal.fits')

    collect_pipeline_cfgs('config')
    args = [
        'config/calwebb_guider.cfg',
        rtdata.input
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'exptype_fgs_image_cal.fits',
    ],
    ids=['cal']
)
def test_fgs(run_pipeline, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rtdata = run_pipeline
    rtdata.output = output

    rtdata.get_truth(os.path.join('fgs/truth/test_fgs', output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

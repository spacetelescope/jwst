from glob import glob
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline import Detector1Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
def test_miri_caldetector1_rate(_jail, rtdata, fitsdiff_default_kwargs):
    rtdata.get_data("miri/image/jw00001001001_01101_00001_MIRIMAGE_uncal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_detector1.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw00001001001_01101_00001_MIRIMAGE_rate.fits"

    rtdata.get_truth("truth/test_miri_caldetector1/jw00001001001_01101_00001_MIRIMAGE_rate.fits")
    assert rtdata.output != rtdata.truth

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, jail):
    """Run calwebb_detector1 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw00001001001_01101_00001_MIRIMAGE_uncal.fits")

    step = Detector1Pipeline()
    step.lastframe.save_results = True
    step.firstframe.save_results = True
    step.dq_init.save_results = True
    step.saturation.save_results = True
    step.rscd.save_results = True
    step.linearity.save_results = True
    step.dark_current.save_results = True
    step.refpix.save_results = True
    step.jump.save_results = True
    step.save_results = True
    step.save_calibrated_ramp = True

    step.run(rtdata.input)
    return rtdata

#    collect_pipeline_cfgs('config')
#    config_file = os.path.join('config', 'calwebb_detector1.cfg')
#    Detector1Pipeline.call(rtdata.input, config_file=config_file,
#        save_results=True, save_calibrated_ramp=True)


@pytest.mark.bigdata
def test_miri_caldetector1_completion(run_pipeline):
    files = glob('*_rate.fits')
    files += glob('*_rateints.fits')
    files += glob('*_dq_init.fits')
    files += glob('*_saturation.fits')
    files += glob('*_rscd.fits')
    files += glob('*_firstframe.fits')
    files += glob('*_lastframe.fits')
    files += glob('*_linearity.fits')
    files += glob('*_jump.fits')
    files += glob('*_dark_current.fits')
    files += glob('*_refpix.fits')
    files += glob('*_ramp.fits')
    # There should be 12 outputs
    assert len(files) == 12


@pytest.mark.bigdata
@pytest.mark.parametrize("output", [
    'jw80600012001_02101_00003_mirimage_rate.fits',
    'jw80600012001_02101_00003_mirimage_rateints.fits',
    'jw80600012001_02101_00003_mirimage_linearity.fits',
    'jw80600012001_02101_00003_mirimage_rscd.fits',
    'jw80600012001_02101_00003_mirimage_saturation.fits',
    'jw80600012001_02101_00003_mirimage_dark_current.fits',
    'jw80600012001_02101_00003_mirimage_refpixel.fits',
    'jw80600012001_02101_00003_mirimage_ramp.fits',],
    ids=['rate', 'rateints', 'linearity', 'rscd', 'saturation',
         'dark_current', 'ref_pixel', 'ramp'])
def test_miri_detector1(run_pipeline, request, fitsdiff_default_kwargs, output):
    """
    Regression test of calwebb_detector1 pipeline performed on MIRI data.
    """
    rtdata = run_pipeline
    rtdata.output = output
    rtdata.get_truth(os.path.join("truth/test_miri_detector1", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import jwst.datamodels as dm
from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

ROOT = 'nrs_verify_nrs1'


@pytest.fixture(scope='module')
def run_detector1(rtdata_module):
    """Run NRS_VERIFY through detector1"""
    rtdata = rtdata_module
    rtdata.get_data('nirspec/imaging/' + replace_suffix(ROOT, 'uncal') + '.fits')

    collect_pipeline_cfgs('config')

    args = [
        'config/calwebb_detector1.cfg', rtdata.input,
        '--steps.dq_init.save_results=True',
        '--steps.saturation.save_results=True',
        '--steps.superbias.save_results=True',
        '--steps.refpix.save_results=True',
        '--steps.linearity.save_results=True',
        '--steps.dark_current.save_results=True',
        '--steps.jump.save_results=True',
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope='module')
def run_image2(run_detector1, rtdata_module):
    """Run NRS_VERIFY through image2"""
    rtdata = rtdata_module
    rtdata.input = replace_suffix(ROOT, 'rate') + '.fits'

    args = [
        'jwst.pipeline.Image2Pipeline', rtdata.input,
        '--steps.assign_wcs.save_results=true',
        '--steps.flat_field.save_results=true',
        '--steps.photom.save_results=true',
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix', [
        'dq_init', 'saturation', 'superbias', 'refpix', 'linearity',
        'dark_current', 'jump', 'rate',
    ])
def test_verify_detector1(run_detector1, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test results of the detector1 and image2 processing"""
    rtdata = rtdata_module
    rtdata.input = replace_suffix(ROOT, 'uncal') + '.fits'
    output = replace_suffix(ROOT, suffix) + '.fits'
    rtdata.output = output
    rtdata.get_truth('truth/test_nirspec_verify/' + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix', [
        'assing_wcs', 'flat_field', 'cal',
    ])
def test_verify_image2(run_image2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test results of the detector1 and image2 processing"""
    rtdata = rtdata_module
    rtdata.input = replace_suffix(ROOT, 'uncal') + '.fits'
    output = replace_suffix(ROOT, suffix) + '.fits'
    rtdata.output = output
    rtdata.get_truth('truth/test_nirspec_verify/' + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

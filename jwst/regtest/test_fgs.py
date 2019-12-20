"""Regression tests for FGS Guidestar in ID and FINEGUIDE modes"""
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


def is_like_truth(rtdata, fitsdiff_default_kwargs, suffix, truth_path='truth/fgs/test_fgs'):
    """Compare step outputs with truth"""
    output = replace_suffix(
        os.path.splitext(os.path.basename(rtdata.input))[0], suffix
    ) + '.fits'
    rtdata.output = output

    rtdata.get_truth(os.path.join(truth_path, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture(scope='module')
def run_fgs_imaging_pipelines(jail, rtdata_module):
    """Run the pipelines"""
    rtdata = rtdata_module
    rtdata.get_data('fgs/level1b/exptype_fgs_image_uncal.fits')

    collect_pipeline_cfgs('config')
    args = [
        'config/calwebb_detector1.cfg',
        rtdata.input,
        '--steps.group_scale.save_results=true',
        '--steps.dq_init.save_results=true',
        '--steps.saturation.save_results=true',
        '--steps.superbias.save_results=true',
        '--steps.refpix.save_results=true',
        '--steps.firstframe.save_results=true',
        '--steps.lastframe.save_results=true',
        '--steps.linearity.save_results=true',
        '--steps.dark_current.save_results=true',
        '--steps.persistence.save_results=true',
        '--steps.jump.save_results=true',
        '--steps.ramp_fit.save_results=true',
        '--steps.gain_scale.save_results=true',
    ]
    Step.from_cmdline(args)

    return rtdata


file_roots = ['exptype_fgs_acq1', 'exptype_fgs_id_image', 'exptype_fgs_id_stack']
@pytest.fixture(scope='module', params=file_roots, ids=file_roots)
def run_guider_pipelines(jail, rtdata_module, request):
    """Run pipeline for guider data"""
    rtdata = rtdata_module
    rtdata.get_data('fgs/level1b/' + request.param + '_uncal.fits')

    collect_pipeline_cfgs('config')
    args = [
        'config/calwebb_guider.cfg',
        rtdata.input,
        '--steps.dq_init.save_results=true',
        '--steps.guider_cds.save_results=true',
        '--steps.flat_field.save_results=true',
    ]
    Step.from_cmdline(args)

    return rtdata

guider_suffixes = ['cal', 'dq_init', 'flat_field', 'guider_cds']
@pytest.mark.bigdata
@pytest.mark.parametrize('suffix', guider_suffixes, ids=guider_suffixes)
def test_fgs_guider(run_guider_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression for FGS Guider data"""
    is_like_truth(run_guider_pipelines, fitsdiff_default_kwargs, suffix)


imaging_suffixes = ['0_ramp_fit', '1_ramp_fit', 'dark_current', 'dq_init', 'gain_scale',
            'gain_scaleints', 'group_scale', 'jump', 'linearity', 'persistence',
            'rate', 'rateints', 'refpix', 'saturation', 'superbias', 'trapsfilled']
@pytest.mark.bigdata
@pytest.mark.parametrize('suffix', imaging_suffixes, ids=imaging_suffixes)
def test_fgs_imaging(run_fgs_imaging_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    is_like_truth(run_fgs_imaging_pipelines, fitsdiff_default_kwargs, suffix)

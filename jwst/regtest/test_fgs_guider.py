"""Regression tests for FGS Guidestar in ID and FINEGUIDE modes"""
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

from . import regtestdata as rt


file_roots = ['exptype_fgs_acq1', 'exptype_fgs_fineguide', 'exptype_fgs_id_image', 'exptype_fgs_id_stack']
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
    ]
    Step.from_cmdline(args)

    return rtdata

guider_suffixes = ['cal', 'dq_init', 'guider_cds']
@pytest.mark.bigdata
@pytest.mark.parametrize('suffix', guider_suffixes, ids=guider_suffixes)
def test_fgs_guider(run_guider_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression for FGS Guider data"""
    rt.is_like_truth(run_guider_pipelines, fitsdiff_default_kwargs, suffix,
                     'truth/fgs/test_fgs_guider', is_suffix=True)


@pytest.mark.bigdata
def test_fgs_toobig(rtdata, fitsdiff_default_kwargs, caplog, monkeypatch):
    """Test for the situation where the combined mosaic is too large"""

    # Set the environment to not allow the resultant too-large image.
    monkeypatch.setenv('DMODEL_ALLOWED_MEMORY', 1.0)

    rtdata.get_asn('fgs/image3/image3_asn.json')

    collect_pipeline_cfgs('config')
    args = ['config/calwebb_image3.cfg', rtdata.input]
    Step.from_cmdline(args)
    assert 'model cannot be instantiated' in caplog.text

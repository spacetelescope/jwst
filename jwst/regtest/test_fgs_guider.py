"""Regression tests for FGS Guidestar in ID and FINEGUIDE modes"""
import pytest

from jwst.resample.resample import OutputTooLargeError
from jwst.stpipe import Step

from jwst.regtest import regtestdata as rt


file_roots = ['exptype_fgs_acq1', 'exptype_fgs_fineguide', 'exptype_fgs_id_image', 'exptype_fgs_id_stack']


@pytest.fixture(scope='module', params=file_roots, ids=file_roots)
def run_guider_pipelines(jail, rtdata_module, request):
    """Run pipeline for guider data"""
    rtdata = rtdata_module
    rtdata.get_data('fgs/level1b/' + request.param + '_uncal.fits')

    args = [
        'calwebb_guider',
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
    monkeypatch.setenv('DMODEL_ALLOWED_MEMORY', "0.9")

    rtdata.get_asn('fgs/image3/image3_asn.json')

    args = ['calwebb_image3', rtdata.input]
    with pytest.raises(OutputTooLargeError):
        Step.from_cmdline(args)

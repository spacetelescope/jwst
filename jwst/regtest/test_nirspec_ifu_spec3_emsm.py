"""Regression tests for NIRSpec IFU"""
import pytest

from jwst.regtest import regtestdata as rt

# Define artifactory source and truth
INPUT_PATH = 'nirspec/ifu'
TRUTH_PATH = 'truth/test_nirspec_ifu'


@pytest.fixture(scope='module')
def run_spec3_multi_emsm(jail, rtdata_module):
    """Run Spec3Pipeline"""
    rtdata = rtdata_module

    step_params = {
        'input_path': 'nirspec/ifu/jw01249-o005_20230622t074431_spec3_00001_asn.json',
        'step': 'calwebb_spec3',
        'args': {
            '--steps.cube_build.save_results=true',
            '--steps.cube_build.weighting=emsm',
            '--steps.cube_build.output_file="nirspec_emsm"',
            '--steps.extract_1d.save_results=true',
        }
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'nirspec_emsm_g395h-f290lp_s3d.fits',
        'nirspec_emsm_g395h-f290lp_x1d.fits',
    ]
)
def test_spec3_multi_emsm(run_spec3_multi_emsm, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3_multi_emsm, fitsdiff_default_kwargs, output,
        truth_path=TRUTH_PATH,
        is_suffix=False
    )

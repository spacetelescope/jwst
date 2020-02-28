"""Regression tests for NIRSpec IFU"""
import pytest

from jwst.associations.asn_from_list import asn_from_list
from jwst.lib.suffix import replace_suffix

from . import regtestdata as rt

# Define artifactory source and truth
INPUT_PATH = 'nirspec/ifu'
TRUTH_PATH = 'truth/test_nirspec_ifu'


@pytest.fixture(scope='module')
def run_spec3(jail, run_spec2):
    """Run the Spec3Pipeline on the results from the Spec2Pipeline run"""
    rtdata = run_spec2

    # Create the level 3 `spec3` association from the output product
    # name of `run_spec2`
    product = run_spec2.asn['products'][0]['name']
    input_name = replace_suffix(product, 'cal') + '.fits'
    asn = asn_from_list([input_name], product_name=product)
    _, serialized = asn.dump()
    asn_name = product + '_spec3_asn.json'
    with open(asn_name, 'w') as asn_fh:
        asn_fh.write(serialized)

    # Set the input in the RegtestData instance to avoid download.
    rtdata.input = asn_name
    step_params = {
        'step': 'calwebb_spec3.cfg',
        'args': [
            '--steps.master_background.save_results=true',
            '--steps.mrs_imatch.save_results=true',
            '--steps.outlier_detection.save_results=true',
            '--steps.resample_spec.save_results=true',
            '--steps.cube_build.save_results=true',
            '--steps.extract_1d.save_results=true',
            '--steps.combine_1d.save_results=true',
        ]
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_multi(jail, rtdata_module):
    """Run Spec3Pipeline"""
    rtdata = rtdata_module

    step_params = {
        'input_path': 'nirspec/ifu/single_nrs1-nrs2_spec3_asn.json',
        'step': 'calwebb_spec3.cfg',
        'args': {
            '--steps.master_background.save_results=true',
            '--steps.mrs_imatch.save_results=true',
            '--steps.outlier_detection.save_results=true',
            '--steps.resample_spec.save_results=true',
            '--steps.cube_build.save_results=true',
            '--steps.extract_1d.save_results=true',
            '--steps.combine_1d.save_results=true',
        }
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.mark.skip(reason="single exposure processing is redundant with calspec2")
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw00626009002_02101_00001_nrs1_run_spec2_g395h-f290lp_s3d.fits',
        'jw00626009002_02101_00001_nrs1_run_spec2_g395h-f290lp_x1d.fits',
    ]
)
def test_spec3(run_spec3, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3, fitsdiff_default_kwargs, output,
        truth_path=TRUTH_PATH,
        is_suffix=False
    )


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw00626009002_02101_00001_nrs1_o009_crf.fits',
        'jw00626009002_02101_00001_nrs2_o009_crf.fits',
        'single_nrs1-nrs2_g395h-f290lp_s3d.fits',
        'single_nrs1-nrs2_g395h-f290lp_x1d.fits',
    ]
)
def test_spec3_multi(run_spec3_multi, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3_multi, fitsdiff_default_kwargs, output,
        truth_path=TRUTH_PATH,
        is_suffix=False
    )

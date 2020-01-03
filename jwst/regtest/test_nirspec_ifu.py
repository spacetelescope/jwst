"""Regression tests for NIRSpec"""
import pytest

from . import regtestdata as rt


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the pipelines"""
    step_params = {
        'input_path': 'nirspec/ifu/nrs_ifu_nrs1_rate.fits',
        'step': 'calwebb_spec2.cfg',
        'args': [
            '--steps.bkg_subtract.save_results=true',
            '--steps.assign_wcs.save_results=true',
            '--steps.imprint_subtract.save_results=true',
            '--steps.msa_flagging.save_results=true',
            '--steps.extract_2d.save_results=true',
            '--steps.flat_field.save_results=true',
            '--steps.srctype.save_results=true',
            '--steps.straylight.save_results=true',
            '--steps.fringe.save_results=true',
            '--steps.pathloss.save_results=true',
            '--steps.barshadow.save_results=true',
            '--steps.photom.save_results=true',
            '--steps.resample_spec.save_results=true',
            '--steps.cube_build.save_results=true',
            '--steps.extract_1d.save_results=true',
        ]
    }

    return rt.run_step_from_dict(rtdata_module, **step_params)


@pytest.fixture(scope='module')
def run_spec3(jail, rtdata_module):
    """Run Spec3Pipeline"""
    step_params = {
        'input_path': 'nirspec/ifu/nrs_ifu_spec3_asn.json',
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

    return rt.run_step_from_dict(rtdata_module, **step_params)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'cal', 'flat_field', 'msa_flagging', 'pathloss', 'photom', 's3d', 'srctype', 'x1d']
)
def test_spec2(run_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    rt.is_like_truth(run_spec2, fitsdiff_default_kwargs, suffix,
                     truth_path='truth/nirspec/test_nirspec_ifu')


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'nrs_ifu_g395h-f290lp_s3d.fits',
        'nrs_ifu_g395h-f290lp_x1d.fits',
        'nrs_ifu_nrs1_o028_crf.fits',
        'nrs_ifu_nrs2_o028_crf.fits',
    ]
)
def test_spec3(run_spec3, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3, fitsdiff_default_kwargs, output,
        truth_path='truth/nirspec/test_nirspec_ifu',
        is_suffix=False
    )

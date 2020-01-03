"""Regression tests for MIRI MRS modes"""
import pytest

from . import regtestdata as rt


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the pipelines"""
    step_params = {
        'input_path': 'miri/mrs/ifushort_ch12_rate.fits',
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


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'cal', 'flat_field', 'fringe', 'photom', 's3d', 'srctype', 'straylight', 'x1d']
)
def test_spec2(run_spec2, fitsdiff_default_kwargs, suffix):
    """Test ensuring the callwebb_spec2 is operating appropriately for MIRI MRS data"""
    rt.is_like_truth(run_spec2, fitsdiff_default_kwargs, suffix,
                     truth_path='truth/test_miri_mrs')

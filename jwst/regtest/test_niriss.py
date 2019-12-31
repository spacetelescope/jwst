"""Regression tests for NIRISS"""
import pytest

from . import regtestdata as rt


@pytest.fixture(scope='module')
def run_nis_wfss_spectral(jail, rtdata_module):
    """Run the pipelines"""
    step_params = {
        'input_path': 'niriss/level2a/nir_wfss_spec2_asn.json',
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
    ['assign_wcs', 'bsub', 'cal', 'extract_2d', 'flat_field', 'photom', 'srctype', 'x1d']
)
def test_nis_wfss_spectral(run_nis_wfss_spectral, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    fitsdiff_default_kwargs['ignore_keywords'].append('SCATFILE')
    rt.is_like_truth(run_nis_wfss_spectral, fitsdiff_default_kwargs, suffix)
